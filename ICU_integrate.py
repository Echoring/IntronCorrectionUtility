#!/usr/bin/env python3
"""
Intron Correction Utility step 3: Adjust gene record according to changed transcripts
"""

import sys
import argparse
import concurrent.futures
from collections import defaultdict
from typing import Dict, List, Tuple
import copy
import logging
import os
from dataclasses import dataclass

@dataclass
class gffrecord:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str]
    raw_line: str

class SeqProcessor:
    def __init__(self, seqid: str):
        self.seqid = seqid
        self.genes: dict[str, gffrecord] = {}  # {ID: gffrecord, ...}
        self.all_children: dict[str, gffrecord] = {}
        self.other: dict[str, gffrecord] = {} 
        self.features_by_parent: dict[str, dict[str, gffrecord]] = defaultdict(dict)  # {'parentID': {ID: gffrecord}, ...}}
                        
    def rebuild_parent_relationships(self):
        self.features_by_parent.clear()
        
        for featid, feat in self.all_children.items():
            parent_id = feat.attributes.get('Parent')
            if parent_id:
                for pid in parent_id.split(','):
                    pid = pid.strip()
                    if pid:
                        self.features_by_parent[pid][featid] = feat

    def update_feature_parent(self, update_children_id: list[str], old_parent: str, new_parent: str):
        """update given feature's parent in self.all_children"""
        for cid in update_children_id:
            feature = self.all_children[cid]

            current_parents = str(feature.attributes['Parent']).split(',')
            updated_parents = []
            
            for parent in current_parents:
                parent = parent.strip()
                if parent == old_parent:
                    updated_parents.append(new_parent)
                else:
                    updated_parents.append(parent)
            
            feature.attributes['Parent'] = ','.join(updated_parents)
            
    def get_output_lines(self, count: int, geneprefix: str) -> Tuple[List[str], int]:
        '''
        return str record list and start count for next processor. count=0 and geneprefix='' will skip rename
        will not rename non-gene record
        '''
        output_lines = []
        type_priority = {
            'mrna': 0,
            'transcript': 0,
            'cds': 1,
            'exon': 2,
            'intron': 3,
            'utr': 4,
            'three_prime_utr': 4,
            'five_prime_utr': 4,
            'start_codon': 5,
            'stop_codon': 5
        }
        
        toprecord = list(self.genes.items()) + list(self.other.items())
        for recordid, record in sorted(toprecord, key=lambda x: x[1].start):
            if record.type == 'gene':
                newgid = f'{geneprefix}{count:05d}' if count != 0 else ''
                output_lines.append(format_record(record, newgid, ''))
                count = count + 1 if count != 0 else 0
            else:
                newgid = ''
                output_lines.append(format_record(record, '', ''))
            
            # process all child
            stack = [(recordid, newgid, defaultdict(int))]
            while stack:
                current_original_id, current_new_id, counters = stack.pop()
                
                children = self.features_by_parent.get(current_original_id, {})
                if not children:
                    continue
                sorted_children = sorted(children.items(), key=lambda x: (x[1].start, type_priority.get(x[1].type.lower(), 999)))
                for child_id, child_feat in reversed(sorted_children):
                    counters[child_feat.type] += 1
                    new_child_id = f'{current_new_id}.{child_feat.type}{counters[child_feat.type]}' if current_new_id != '' else ''
                    output_lines.append(format_record(child_feat, new_child_id, current_new_id))
                    stack.append((child_id, new_child_id, defaultdict(int)))
        
        return output_lines, count
                    
    # ---------- Core ----------
    def adjust_gene_boundaries(self):
        """adjust gene boundary inconsistant with child features"""
        logger.debug(f'{self.seqid}: adjusting gene boundary inconsistant with child features...')
        
        genes_adjusted = 0
        
        for gene_id, gene in list(self.genes.items()):
            direct_features = self.features_by_parent.get(gene_id, {})
            
            if not direct_features:
                continue
            
            min_start = min(f.start for i, f in direct_features.items())
            max_end = max(f.end for i, f in direct_features.items())
            
            if gene.start == min_start and gene.end == max_end:
                continue
             
            # inconsistant gene   
            original_range = f"{gene.start}-{gene.end}"
            new_range = f"{min_start}-{max_end}"
            
            gene.start = min_start
            gene.end = max_end
            
            genes_adjusted += 1
            self.genes[gene_id] = gene
            logger.debug(f'{self.seqid}: {gene_id} {original_range} -> {new_range}')

        logger.debug(f"{self.seqid}: {genes_adjusted} genes were adjusted for boundary inconsistant with child features.")
    
    def process_strand_inconsistency(self):
        """Adjust gene strand inconsistant with child features, create a new gene for inconsistant child features"""
        logger.debug(f'{self.seqid}: adjusting gene strand inconsistant with child features...')
        
        genes_processed = 0
        
        for gene_id, gene in list(self.genes.items()):      
            direct_features = self.features_by_parent.get(gene_id, {})
            if not direct_features:
                continue
            
            direct_features_by_strand: defaultdict[str, dict[str, gffrecord]] = defaultdict(dict)
            for fid, feat in direct_features.items():
                direct_features_by_strand[feat.strand][fid] = feat
            
            gene_strand = gene.strand
            if len(direct_features_by_strand) == 1 and gene_strand in direct_features_by_strand:
                continue
            
            # inconsistant gene
            genes_processed += 1
            
            forward = len(direct_features_by_strand.get('+', {}))
            minus = len(direct_features_by_strand.get('-', {}))
            logger.debug(f'{self.seqid}: {gene_id} ({gene_strand}) has {forward} feature (+) and {minus} feature (-)')
            
            # Keep consistant feature
            same_strand_features = direct_features_by_strand.get(gene_strand, {})
            if same_strand_features:
                new_gene = copy.deepcopy(gene)
                new_start = min(f.start for i, f in same_strand_features.items())
                new_end = max(f.end for i, f in same_strand_features.items())
                new_gene.start = new_start
                new_gene.end = new_end
                logger.debug(f"{self.seqid}: {gene_id} ({gene_strand}) {gene.start}-{gene.end} -> {new_start}-{new_end}")
                self.genes[gene_id] = new_gene
            else:
                self.genes.pop(gene_id)
            
            # Create new gene for inconsistant feature
            different_strand = '+' if gene_strand == '-' else '-'
            different_strand_features = direct_features_by_strand[different_strand]
            
            new_gene = copy.deepcopy(gene)
            new_gene.strand = different_strand
            new_start = min(f.start for i, f in different_strand_features.items())
            new_end = max(f.end for i, f in different_strand_features.items())
            new_gene.start = new_start
            new_gene.end = new_end
            
            new_id = generate_unique_id([gene_id], "rc")
            new_gene.attributes['ID'] = new_id
            self.genes[new_id] = new_gene
            logger.debug(f" {self.seqid}: {gene_id} ({different_strand}) -> {new_id} {new_start}-{new_end}")
            
            # update child feature
            update_children_id = list(different_strand_features.keys())
            self.update_feature_parent(update_children_id, gene_id, new_id)
        
        logger.debug(f"{self.seqid}: {genes_processed} genes were adjusted for strand inconsistant with child features.")
    
    def process_gaps_without_features(self):
        """Split genes by discontinued intervals"""
        logger.debug(f"{self.seqid}: splitting genes by discontinued intervals...")
        
        genes_split = 0
        
        for gene_id, gene in list(self.genes.items()):
            direct_features = self.features_by_parent.get(gene_id, {})
            if not direct_features:
                continue
            
            direct_features = sorted(direct_features.items(), key=lambda x: x[1].start)
            
            covered_intervals = []
            current_start = direct_features[0][1].start
            current_end = direct_features[0][1].end
            for feat in direct_features[1:]:
                if feat[1].start <= current_end + 1:
                    current_end = max(current_end, feat[1].end)
                else:
                    covered_intervals.append((current_start, current_end))
                    # Start new interval
                    current_start = feat[1].start
                    current_end = feat[1].end
            # Append final
            covered_intervals.append((current_start, current_end))
            
            if len(covered_intervals) == 1:
                continue
            
            # Split gene
            genes_split += 1
            logger.debug(f'{self.seqid}: {gene_id} has {len(covered_intervals)} discontinued intervals')
            
            # group features by interval
            interval_features: defaultdict[int, list[tuple[str, gffrecord]]] = defaultdict(list)
            for feat in direct_features:
                for i, (interval_start, interval_end) in enumerate(covered_intervals):
                    if interval_start <= feat[1].start and feat[1].end <= interval_end:
                        interval_features[i].append(feat)
                        break
            
            # Create new genes for each interval
            for i, (start, end) in enumerate(covered_intervals):
                new_gene = copy.deepcopy(gene)
                new_gene.start = start
                new_gene.end = end
                
                new_id = generate_unique_id([gene_id], f"split{i+1}")
                new_gene.attributes['ID'] = new_id
                self.genes[new_id] = new_gene
                logger.debug(f"{self.seqid}: {gene_id} (interval {i+1}) -> {new_id} {start}-{end}")
                
                # update child feature
                features_in_interval = interval_features.get(i, [])
                updatefid = [fid for fid, feat in features_in_interval]
                self.update_feature_parent(updatefid, gene_id, new_id)

            self.genes.pop(gene_id)
        
        logger.debug(f"{self.seqid}: {genes_split} gene split for discontinued intervals.")
    
    def merge_overlapping_genes(self):
        """Merge overlapped genes on same strand"""
        logger.debug(f"{self.seqid}: merging overlapped genes on same strand...")
        genes_merged = 0
        merge_gene_list: list[list[tuple[str, gffrecord]]] = []
        
        # Group by strand
        genes_by_strand: defaultdict[str, list[tuple[str, gffrecord]]] = defaultdict(list)
        for gid, gene in self.genes.items():
            genes_by_strand[gene.strand].append((gid, gene))
        
        for strand, genes in genes_by_strand.items():
            genes.sort(key=lambda x: x[1].start)
            
            current_gene = genes[0]
            current_list = [current_gene]
            
            for gene in genes[1:]:
                # Overlapped
                if gene[1].start <= current_gene[1].end:
                    current_gene[1].end = max(current_gene[1].end, gene[1].end)
                    current_list.append(gene)
                else:
                    if len(current_list) > 1:
                        merge_gene_list.append(current_list)
                    # Start new interval
                    current_gene = gene
                    current_list = [current_gene]
            
            # final interval
            if len(current_list) > 1:
                merge_gene_list.append(current_list)
                
        # Create merged genes    
        for one_merged_gene in merge_gene_list:
            merged_ids = [x[0] for x in one_merged_gene]
            new_id = generate_unique_id(merged_ids, "merge")
            new_gene = copy.deepcopy(one_merged_gene[0][1])
            new_gene.start = one_merged_gene[0][1].start
            new_gene.end = one_merged_gene[-1][1].end
            new_gene.score = '.' 
            new_gene.attributes['ID'] = new_id
            # Merge other attribute except ID
            for attrs in [x[1].attributes for x in one_merged_gene[1:]]:
                for k, v in attrs.items():
                    if k != 'ID':
                        if k in new_gene.attributes:
                            new_gene.attributes[k] += f',{v}'
                        else:
                            new_gene.attributes[k] = v
                            
            merged_idstr = ', '.join(merged_ids)
            logger.debug(f"{self.seqid}: {merged_idstr} are overlapped")
            ranges = [f"{g[0]} {g[1].start}-{g[1].end}" for g in one_merged_gene]
            logger.debug(f"{self.seqid}: {' + '.join(ranges)} -> {new_id} {new_gene.start}-{new_gene.end}")
            self.genes[new_id] = new_gene
            # update children and remove previous gene
            for gid in merged_ids:
                feat = list(self.features_by_parent[gid].keys())
                self.update_feature_parent(feat, gid, new_id)
                self.genes.pop(gid)
                genes_merged += 1
        
        logger.debug(f"{self.seqid}: {genes_merged} overlaped genes are merged.")
    
    def merge_identical_mrna(self):
        """Merge identical mRNA in gene"""
        logger.debug(f"{self.seqid}: merging identical mRNA in gene...")
        
        mrna_merged = 0
        
        for gene_id, gene in self.genes.items():
            direct_features = self.features_by_parent.get(gene_id, {})
            if not direct_features:
                continue
            
            mrnas = {i: f for i, f in direct_features.items() if f.type.lower() in ['mrna', 'transcript']}
            
            if len(mrnas) <= 1:
                continue
            
            # Use record tuple (except source, score and attribute) as key to find duplicates
            mrna_groups: defaultdict[tuple, list[tuple[str, gffrecord]]] = defaultdict(list)
            for mrnaid, mrna in mrnas.items():
                feature_fingerprint = get_feature_fingerprint(mrna)
                
                child_features = self.features_by_parent.get(mrnaid, {})
                child_fingerprints = []
                for cid, child in child_features.items():
                    if child.type.lower() == 'cds':
                        child_fingerprint = get_feature_fingerprint(child)
                        child_fingerprints.append(child_fingerprint)
                child_fingerprints.sort()
                fingerprint = (feature_fingerprint, tuple(child_fingerprints))
                
                mrna_groups[fingerprint].append((mrnaid, mrna))
            
            for fingerprint, mrna_list in mrna_groups.items():
                if len(mrna_list) == 1:
                    continue
                
                # Duplicated, create new mRNA and remove old
                original_ids = [x[0] for x in mrna_list]
                merged_idstr = ', '.join(original_ids)
                logger.debug(f"{self.seqid}: {merged_idstr} are duplicated")
                
                new_id = generate_unique_id(original_ids, "eq")
                new_mrna = copy.deepcopy(mrna_list[0][1])
                new_mrna.score = '.' 
                new_mrna.attributes['ID'] = new_id
                # Merge other attribute except ID and parent
                for attrs in [x[1].attributes for x in mrna_list[1:]]:
                    for k, v in attrs.items():
                        if k != 'ID' and k != 'Parent':
                            if k in new_mrna.attributes:
                                new_mrna.attributes[k] += f',{v}'
                            else:
                                new_mrna.attributes[k] = v
                self.all_children[new_id] = new_mrna
                for mid in original_ids:
                    self.all_children.pop(mid)
                    mrna_merged += 1
                logger.debug(f"{self.seqid}: {merged_idstr} -> {new_id}")
                                        
                # Similarly process duplicated children
                childdictitem: List[tuple[str, gffrecord]] = []
                for mid in original_ids:
                    childdictitem.extend(self.features_by_parent[mid].items())

                childgroup: defaultdict[tuple, list[tuple[str, gffrecord]]] = defaultdict(list)
                for cid, child in childdictitem:
                    childfingerprint = get_feature_fingerprint(child)
                    childgroup[childfingerprint].append((cid, child))
                    
                for childfingerprint, child_list in childgroup.items():
                    # Duplicated child
                    if len(child_list) > 1:
                        original_cids = [x[0] for x in child_list]
                        merged_cidstr = ', '.join(original_cids)
                        logger.debug(f"{self.seqid}: {merged_cidstr} are duplicated")
                        new_cid = generate_unique_id(original_cids, "eq")
                        new_child = copy.deepcopy(child_list[0][1])
                        new_child.score = '.' 
                        new_child.attributes['ID'] = new_cid
                        new_child.attributes['Parent'] = new_id
                        # Merge other attribute except ID and parent
                        for attrs in [x[1].attributes for x in child_list[1:]]:
                            for k, v in attrs.items():
                                if k != 'ID' and k != 'Parent':
                                    if k in new_child.attributes:
                                        new_child.attributes[k] += f',{v}'
                                    else:
                                        new_child.attributes[k] = v
                        self.all_children[new_cid] = new_child
                        for cid in original_cids:
                            self.all_children.pop(cid)
                        logger.debug(f"{self.seqid}: {merged_cidstr} -> {new_cid}")
                    # Not duplicated child (non-CDS)
                    else:
                        self.update_feature_parent([child_list[0][0]], child_list[0][1].attributes['Parent'], new_id)
        
        logger.debug(f"{self.seqid}: {mrna_merged} identical mRNA merged.")
  
# ---------- Util ---------- 
def generate_unique_id(base_id: List[str], operation: str) -> str:
    """generate ID like <originID.operation.originID2...>"""
    ids_sorted = sorted(set(base_id))
    if len(ids_sorted) == 1:
        new_id = f"{ids_sorted[0]}.{operation}"
    else:
        new_id = f".{operation}.".join(ids_sorted)
    return new_id

def get_feature_fingerprint(feature: gffrecord) -> tuple:
    """return tuple fingerprint, ignore source, score and attributes"""
    fingerprint = (
        feature.seqid,
        feature.type.lower(),
        feature.start,
        feature.end,
        feature.strand,
        feature.phase
    )
    return fingerprint

def format_record(record: gffrecord, newid: str, newparent: str) -> str:
    '''if newid/newparent = '', skip rename'''
    if not record.attributes:
        attr_str = '.'
    else:
        attrs = record.attributes.copy()
        if newid:
            attrs['ID'] = newid
        if newparent:
            attrs['Parent'] = newparent
        attr_str = ';'.join([f"{k}={v}" for k, v in attrs.items()])
    recordstr = f"{record.seqid}\t{record.source}\t{record.type}\t{record.start}\t{record.end}\t{record.score}\t{record.strand}\t{record.phase}\t{attr_str};"
    return recordstr

# ---------- input ---------- 
def parse_gff3_by_Seq(input_file: str) -> Tuple[Dict[str, SeqProcessor], List[str], List[str]]:
    """
    return {'seqid': SeqProcessor}, ['seqid'], ['header lines']
    """
    logger.debug("Phasing annotation...")
    
    genes_data = defaultdict(dict)
    children_data = defaultdict(dict)
    other_data = defaultdict(dict)
    all_header_lines = []
    usedid = set()
    
    genecount = 0
    noidcount = 0
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                all_header_lines.append(line)
                continue
                
            fields = line.split('\t')
            if len(fields) < 9:
                all_header_lines.append(line)
                continue
            
            seqid = fields[0]
            attrs = {}
            attr_string = fields[8]
            
            for pair in attr_string.split(';'):
                pair = pair.strip()
                if not pair:
                    continue
                if '=' in pair:
                    key, value = pair.split('=', 1)
                    attrs[key.strip()] = value.strip()
                    
            record = gffrecord(
                seqid=seqid,
                source=fields[1],
                type=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                score=fields[5],
                strand=fields[6],
                phase=fields[7],
                attributes=attrs,
                raw_line=line 
            )
            
            recordid = attrs.get('ID')
            parentid = attrs.get('Parent')
            
            if not recordid:
                noidcount += 1
                recordid = f'NoID{noidcount:05d}'  
            if recordid not in usedid:
                usedid.add(recordid)
            else:
                logger.error(f'Duplicated record ID: {recordid}')
                sys.exit(1)    
                        
            if record.type.lower() == 'gene':
                genes_data[seqid][recordid] = record
                genecount += 1
            elif parentid:
                children_data[seqid][recordid] = record
            else:
                other_data[seqid][recordid] = record
    
    seqidlist = list(genes_data.keys())
    seqid2processor = {}
    for seqid in seqidlist:
        processor = SeqProcessor(seqid)
        processor.genes = genes_data[seqid]
        processor.all_children = children_data[seqid]
        processor.other = other_data[seqid]
        seqid2processor[seqid] = processor
    
    logger.debug(f"{len(seqidlist)} sequences found: {', '.join(seqidlist)}")
    logger.debug(f"{genecount} genes imported.")
    
    # Remove duplicated headers. Some gff3 use a lot of ### to split record, which will be considered as header.
    all_header_lines = list(dict.fromkeys(all_header_lines)) 
    
    return seqid2processor, seqidlist, all_header_lines

# ---------- Main process ---------- 
def process_Seq(seqid: str, processor: SeqProcessor) -> SeqProcessor:
    logger.debug(f"{seqid}: process started")
    
    processor.rebuild_parent_relationships()
    
    # Core process
    processor.adjust_gene_boundaries()
    
    processor.process_strand_inconsistency()
    processor.rebuild_parent_relationships()
    
    processor.process_gaps_without_features()
    processor.rebuild_parent_relationships()
    
    processor.merge_overlapping_genes()
    processor.rebuild_parent_relationships()
    
    processor.merge_identical_mrna()
    processor.rebuild_parent_relationships()
    
    logger.debug(f"{seqid}: process finished")
    return processor

# ---------- Main ----------
def integrate(inarg=None):
    '''
    Intron Correction Utility step 3: Adjust gene record according to changed transcripts
    return integrated.gff3
    '''
    parser = argparse.ArgumentParser(
        description='Intron Correction Utility step 3: Adjust gene record according to changed transcripts',
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    
    # Positional
    parser.add_argument('fixed_gff3', 
                        help='Intermediate GFF3 file which new mRNA not integrated into gene record')
    
    # Optional
    parser.add_argument('-o', '--out', type=str, default='ICU.integrated.gff3',
                       help='Generate output file named (default: ICU.integrated.gff3)')
    parser.add_argument('--rename', type=str, default='',
                       help='Rename all gene and its child record ID begin with prefix+number; if not set, join original ID with operation (default: not set)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Threads limit to use. If set as 0, will use as much as possible (default: 1); Max vaild thread equals sequence counts')
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Print DEBUG level log')
    
    args = parser.parse_args() if inarg is None else parser.parse_args(inarg)
    
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stderr),
        ]
    )
    global logger
    logger = logging.getLogger(__name__)
    
    logger.debug("="*60)
    logger.debug("ICU_integrate started")
    logger.debug("="*60)
    
    num_threads = os.cpu_count() if args.threads == 0 else args.threads
    
    logger.debug(f"Input: {args.fixed_gff3}")
    logger.debug(f"Options:")
    logger.debug(f"  - Output: {args.out}")
    logger.debug(f"  - Rename prefix: {args.rename}")
    logger.debug(f"  - Threads: {num_threads}")      
    logger.debug(f"  - debug log: {args.debug}")
    
    # Main process
    logger.debug("="*60)
    logger.debug(f"Main process started")
    logger.debug("="*60)
    try:
        # Split GFF3 record by sequence
        seqid2processor, seqidlist, headers = parse_gff3_by_Seq(args.fixed_gff3)
        
        if not seqidlist:
            logger.error("No sequence found.")
            sys.exit(1)
        
        # Multi-threading each sequence
        processors: dict[str, SeqProcessor] = {}

        actual_threads = min(num_threads, len(seqidlist))
        logger.debug(f"Actual threads used: {actual_threads}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=actual_threads) as executor:
            future_to_seq = {
                executor.submit(process_Seq, seqid, processor): seqid # Core
                for seqid, processor in seqid2processor.items()
            }
            
            for future in concurrent.futures.as_completed(future_to_seq):
                seqid = future_to_seq[future]
                try:
                    processor = future.result()
                    processors[seqid] = processor
                except Exception as e:
                    logger.exception(f"Unexcepted error occurred in processing {seqid}: {e}")
                    sys.exit(1)
                    

        # Write output
        logger.debug(f'Writing output file {args.out} ...')
        
        with open(args.out, 'w') as f:
            for line in headers:
                f.write(line + '\n')
            rcount = 1 if args.rename else 0
            for seq in seqidlist:
                if seq in processors:
                    processor = processors[seq]
                    output_lines, rcount = processor.get_output_lines(rcount, args.rename)
                    for line in output_lines:
                        f.write(line + '\n') 
        
        logger.debug("="*60)
        logger.debug("ICU_integrate completed.")
        logger.debug("="*60)
        logger.debug("Outputs:")
        logger.debug(f"  - Gene record integrated with new transcripts GFF3: {args.out}")
        return args.out
        
    except Exception as e:
        logger.exception(f"Unexcepted error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    integrate()