#!/usr/bin/env python3
"""
Intron Correction Utility step 2: Fix invalid intron in transcripts by trusted annotation
"""

import re
import sys
import argparse
import concurrent.futures
import logging
import os
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

@dataclass
class TranscriptInfo:
    seqid: str
    start: int
    end: int
    strand: str
    source: str
    score: str # convert float to str to avoid '.'

class GeneInfo:
    def __init__(self, gene_id: str, gene_record = None):
        self.gene_id = gene_id
        self.gene_record = gene_record  # original record line of gene
        self.transcripts = {}  # {transcript_id: [self_records, subfeature_record, ...]
        self.other_records = []  # other records
        
    def add_record(self, record: str, feature: str, record_id = None, parent = None):
        if feature == 'gene':
            self.gene_record = record
        elif feature == 'mRNA':
            if record_id not in self.transcripts:
                self.transcripts[record_id] = []
            self.transcripts[record_id].append(record)
        else:
            if parent and parent in self.transcripts: # subfeature of mRNA
                self.transcripts[parent].append(record)
            else: # independent record
                self.other_records.append(record)
    
    def get_all_records(self) -> List[str]:
        all_records = []
        if self.gene_record:
            all_records.append(self.gene_record)
        
        for transcript_records in self.transcripts.values():
            all_records.extend(transcript_records)
        
        all_records.extend(self.other_records)
        return all_records

# ---------- Parse File ----------
def readFastaAsDict(fastafile: str) -> Dict[str, str]:
    fastaDict = {}
    with open(fastafile, 'r') as f:
        allline = f.read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

def parse_portcullis_gtf(file_path: str) -> Tuple[Dict[str, TranscriptInfo], List[str]]:
    """
    Parse portcullis markup GTF, find invalid transcripts
    return: {original ID: TranscriptInfo}, [invalid original transcript ID]
    """
    invalid_transcripts = []
    transcript_info = {}
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            if feature == 'transcript':
                attrs = {}
                attr_str = fields[8]
                pattern = r'(\S+)\s+"([^"]+)"'
                matches = re.findall(pattern, attr_str)
                for key, value in matches:
                    attrs[key] = value
                
                if 'transcript_id' not in attrs:
                    logger.warning(f'Skip transcript for transcript_id attribute is missing in line: {line}')
                    continue
                tid = attrs['transcript_id']
                
                transcript_info[tid] = TranscriptInfo(
                    seqid=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    source=fields[1],
                    score=fields[5]
                )
                
                if 'introns' not in attrs:
                    logger.warning(f'Skip transcript for introns attribute is missing in line: {line}') 
                    continue
                
                if 'invalid' in attrs['introns']:
                    invalid_transcripts.append(tid)
                                      
    return transcript_info, invalid_transcripts


def parse_trusted_gtf(file_path: str) -> Tuple[Dict[str, TranscriptInfo], Dict[str, List[Tuple[int, int, str]]]]:
    """
    return {trusted ID: TranscriptInfo}, {trusted ID: [(start, end, strand), ...#each exon]}
    """
    transcript_info = {}
    transcript_exons = defaultdict(list)
    linefieldlist: List[List[str]] = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            linefieldlist.append(fields)
    
    # unsorted annotation resolve
    order = {'transcript': 1, 'mrna': 1, 'exon': 2, 'cds': 2}
    linefieldlist.sort(key=lambda x: order.get(x[2].lower(), 3))    
    
    for fields in linefieldlist:
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        
        attrs = {}
        attr_str = fields[8]
        pattern = r'(\S+)\s+"([^"]+)"'
        matches = re.findall(pattern, attr_str)
        for key, value in matches:
            attrs[key] = value

        if 'transcript_id' not in attrs:
            line = "\t".join(fields)
            logger.warning(f'Skip record for transcript_id attribute is missing in line: {line}')
            continue
        
        tid = attrs['transcript_id']
        
        # transcript
        if order.get(feature.lower()) == 1:
            transcript_info[tid] = TranscriptInfo(
                seqid=fields[0],
                start=start,
                end=end,
                strand=strand,
                source=fields[1],
                score='.'
            )
        # exon
        elif order.get(feature.lower()) == 2: 
            if tid not in transcript_info:
                line = "\t".join(fields)
                logger.warning(f'skip exon for parent transcript {tid} missing in line: {line}')
                continue
            
            transcript_exons[tid].append((start, end, strand))

    return transcript_info, transcript_exons


def parse_gff3_to_genes(file_path: str) -> Tuple[Dict[str, GeneInfo], list[str]]:
    """
    parse GFF3 and build GeneInfo
    return: {ID: GeneInfo}, [header]
    """
    genes: OrderedDict[str, GeneInfo] = OrderedDict()
    mRNA2genes: dict[str, GeneInfo] = {}
    header = []
    linelist: List[str] = []
    
    with open(file_path, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if not stripped_line:
                continue
            if stripped_line.startswith('#'):
                header.append(stripped_line)
                continue
            fields = stripped_line.split('\t')
            if len(fields) < 9:
                header.append(stripped_line)
                continue
            linelist.append(stripped_line)
            
    # unsorted annotation resolve
    order = {'gene': 1, 'mrna': 2, 'transcript': 2}
    linelist.sort(key=lambda x: order.get(x.split('\t')[2].lower(), 3))
        
    for line in linelist:
        fields = line.split('\t')
        feature = fields[2]
        
        attrs: dict[str, str] = {}
        for pair in fields[8].split(';'):
            if '=' in pair:
                key, value = pair.split('=', 1)
                attrs[key] = value
        
        record_id = attrs.get('ID')
        parent = attrs.get('Parent')
        # gene
        if order.get(feature.lower()) == 1:
            genes[record_id] = GeneInfo(record_id, line)
        # mRNA
        elif order.get(feature.lower()) == 2:
            relatedgene = genes.get(parent)
            if not relatedgene:
                logger.warning(f'skip mRNA for parent gene {parent} missing in line: {line}')
                continue
            relatedgene.add_record(line, 'mRNA', record_id, parent)
            mRNA2genes[record_id] = relatedgene
        # subfeature
        else:
            relatedgene = mRNA2genes.get(parent)
            if not relatedgene:
                logger.warning(f'skip feature for parent {parent} missing in line: {line}')
                continue
            relatedgene.add_record(line, feature, record_id, parent)
         
    # Remove duplicated headers. Some gff3 use a lot of ### to split record, which will be considered as header.
    header = list(dict.fromkeys(header))     
                
    return genes, header

# ---------- Core ----------
def find_longest_orf_in_transcript(transcript_exons: List[Tuple[int, int, str]],
                                   genome_seq: Dict[str, str],
                                   seqid: str,
                                   strand: str,
                                   min_orf_length: int = 300) -> Optional[Tuple[int, int, int, int, int, int, int]]:
    """
    return: orf_start, orf_end, start_codon_start, start_codon_end, stop_codon_start, stop_codon_end, orf_len
    """
    if not transcript_exons:
        return None
    if seqid not in genome_seq:
        return None
    
    # Build transcript sequence and genomic position index
    sorted_exons = sorted(transcript_exons, key=lambda x: x[0])
    seq_parts = []
    genomic_positions = []  # index=position on transcript (forward), value=position on genome
    
    if strand == '+':
        for exon_start, exon_end, exon_strand in sorted_exons:
            exon_seq = genome_seq[seqid][exon_start-1:exon_end]
            seq_parts.append(exon_seq)
            for pos in range(exon_start, exon_end + 1):
                genomic_positions.append(pos)
    else:
        for exon_start, exon_end, exon_strand in reversed(sorted_exons):
            exon_seq = genome_seq[seqid][exon_start-1:exon_end]
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            exon_seq_rc = ''.join(complement.get(base, base) for base in exon_seq[::-1])
            seq_parts.append(exon_seq_rc)
            for pos in range(exon_end, exon_start - 1, -1):
                genomic_positions.append(pos)
    
    transcript_seq = ''.join(seq_parts)
    
    # find all possible ORF and get longest
    longest_orf_len = 0
    best_start_pos = -1
    best_stop_pos = -1
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    for i in range(0, len(transcript_seq) - 2):
        codon = transcript_seq[i:i+3]
        if codon == 'ATG':
            for j in range(i + 3, len(transcript_seq) - 2, 3):
                next_codon = transcript_seq[j:j+3]
                if next_codon in stop_codons:
                    orf_len = j + 3 - i
                    if orf_len >= min_orf_length and orf_len > longest_orf_len:
                        longest_orf_len = orf_len
                        best_start_pos = i
                        best_stop_pos = j
                    break
    
    if longest_orf_len == 0 or longest_orf_len < min_orf_length:
        return None
    
    # convert transcript position to genomic position
    if strand == '+':
        start_codon_start = genomic_positions[best_start_pos]
        start_codon_end = genomic_positions[best_start_pos + 2]
        stop_codon_start = genomic_positions[best_stop_pos]
        stop_codon_end = genomic_positions[best_stop_pos + 2]
        orf_start = start_codon_start
        orf_end = stop_codon_end
    else:
        stop_codon_start = genomic_positions[best_stop_pos + 2] # 3'
        stop_codon_end = genomic_positions[best_stop_pos] 
        start_codon_start = genomic_positions[best_start_pos + 2]  
        start_codon_end = genomic_positions[best_start_pos] # 5'
        orf_start = stop_codon_start
        orf_end = start_codon_end
    
    return (orf_start, orf_end, 
            start_codon_start, start_codon_end,
            stop_codon_start, stop_codon_end,
            longest_orf_len)
    
def find_best_transcript_for_invalid(invalid_tinfo: TranscriptInfo,
                                     trusted_transcripts: Dict[str, TranscriptInfo],
                                     trusted_exons: Dict[str, List[Tuple[int, int, str]]],
                                     genome_seq: Dict[str, str],
                                     min_orf_length: int = 300) -> Tuple[Optional[str], Optional[Tuple], Optional[str]]:
    """
    return: 
      best transcript ID to fix, 
      (orf_start, orf_end, start_codon_start, start_codon_end, stop_codon_start, stop_codon_end, longest_orf_len),
      strand 
    
    if no transcript satisfied requirement, return None, None, None
    """
    seqid = invalid_tinfo.seqid
    
    best_tid = None
    best_orf = None
    best_length = 0
    best_strand = None
    
    candidates = []
    for st_tid, st_info in trusted_transcripts.items():
        if st_info.seqid != seqid:
            continue
        # has overlap
        if invalid_tinfo.end >= st_info.start and invalid_tinfo.start <= st_info.end:
            candidates.append((st_tid, st_info))
    
    if not candidates:
        return None, None, None
    
    # search every candidate and get best one
    for st_tid, st_info in candidates:
        exons = trusted_exons[st_tid]
        
        # Sometimes stringtie can generate strand '.' record. 
        # In this case, search both.
        search_strands = []
        if st_info.strand == '.':
            search_strands = ['+', '-']
        else:
            search_strands = [st_info.strand]
        
        for search_strand in search_strands:
            orf_result = find_longest_orf_in_transcript(
                exons, genome_seq, seqid, search_strand, min_orf_length
            )
            
            if orf_result:
                orf_start = orf_result[0]
                orf_end = orf_result[1]
                if invalid_tinfo.end >= orf_start and invalid_tinfo.start <= orf_end:
                    orf_len = orf_result[6]
                    if orf_len > best_length:
                        best_tid = st_tid
                        best_orf = orf_result
                        best_length = orf_len
                        best_strand = search_strand
        
    return best_tid, best_orf, best_strand

def create_corrected_transcript_records(invalid_tid: str,
                                        gene_id: str,
                                        trusted_tid: str,
                                        trusted_info: TranscriptInfo,
                                        trusted_exons: List[Tuple[int, int, str]],
                                        orf_info: Tuple,
                                        found_strand: str) -> List[str]:
    """
    return [mRNArecord, subrecord, ...]
    """
    records = []
    
    orf_start, orf_end, start_codon_start, start_codon_end, stop_codon_start, stop_codon_end, orf_len = orf_info
    
    # Create mRNA record
    mrna_record = f"{trusted_info.seqid}\t{trusted_info.source}\tmRNA\t{orf_start}\t{orf_end}\t{trusted_info.score}\t{found_strand}\t.\tID={invalid_tid};Parent={gene_id};fixby={trusted_tid};"
    records.append(mrna_record)
    
    # For NGS transcriptome, near-boundary reads can be rare and fail to be assembled. 
    # To avoid misleading boundary, do not generate UTR record.
    filtered_exons = [exon for exon in trusted_exons if not (exon[0] > orf_end or exon[1] < orf_start)]
    filtered_exons.sort(key=lambda x: x[0])
    
    # Create exon and intron record
    for i, (exon_start, exon_end, exon_strand) in enumerate(filtered_exons, 1):
        # exon
        exon_id = f"{invalid_tid}.exon{i}"
        exon_start = max(exon_start, orf_start)
        exon_end = min(exon_end, orf_end)
        exon_record = f"{trusted_info.seqid}\t{trusted_info.source}\texon\t{exon_start}\t{exon_end}\t.\t{found_strand}\t.\tID={exon_id};Parent={invalid_tid};"
        records.append(exon_record)
        
        # intron
        if i < len(filtered_exons):
            next_exon_start = filtered_exons[i][0]
            intron_start = exon_end + 1
            intron_end = next_exon_start - 1
            intron_id = f"{invalid_tid}.intron{i}"
            intron_record = f"{trusted_info.seqid}\t{trusted_info.source}\tintron\t{intron_start}\t{intron_end}\t.\t{found_strand}\t.\tID={intron_id};Parent={invalid_tid};"
            records.append(intron_record)
    
    # Create start and stop codon record
    sc_start, sc_end = (start_codon_start, start_codon_end) if start_codon_start <= start_codon_end else (start_codon_end, start_codon_start)
    start_codon_id = f"{invalid_tid}.start_codon1"
    start_codon_record = f"{trusted_info.seqid}\t{trusted_info.source}\tstart_codon\t{sc_start}\t{sc_end}\t.\t{found_strand}\t0\tID={start_codon_id};Parent={invalid_tid};"
    records.append(start_codon_record)

    tc_start, tc_end = (stop_codon_start, stop_codon_end) if stop_codon_start <= stop_codon_end else (stop_codon_end, stop_codon_start)
    stop_codon_id = f"{invalid_tid}.stop_codon1"
    stop_codon_record = f"{trusted_info.seqid}\t{trusted_info.source}\tstop_codon\t{tc_start}\t{tc_end}\t.\t{found_strand}\t0\tID={stop_codon_id};Parent={invalid_tid};"
    records.append(stop_codon_record)
    
    # Create CDS record
    # UTR is not included, so exon equals CDS.
    exon_cds_pairs = []
    for exon_start, exon_end, exon_strand in filtered_exons:
        cds_start = max(exon_start, orf_start)
        cds_end = min(exon_end, orf_end)
        if cds_start <= cds_end:
            exon_cds_pairs.append((cds_start, cds_end))
    
    # calculate CDS phase
    nextphase = 0
    if found_strand == '+':
        for i, (cds_start, cds_end) in enumerate(exon_cds_pairs, 1):
            cds_id = f"{invalid_tid}.CDS{i}"
            cds_record = f"{trusted_info.seqid}\t{trusted_info.source}\tCDS\t{cds_start}\t{cds_end}\t.\t{found_strand}\t{nextphase}\tID={cds_id};Parent={invalid_tid};"
            records.append(cds_record)
            nextphase = - (nextphase + cds_end - cds_start + 1) % 3
    else:
        for i, (cds_start, cds_end) in reversed(list(enumerate(exon_cds_pairs, 1))):
            cds_id = f"{invalid_tid}.CDS{i}"
            cds_record = f"{trusted_info.seqid}\t{trusted_info.source}\tCDS\t{cds_start}\t{cds_end}\t.\t{found_strand}\t{nextphase}\tID={cds_id};Parent={invalid_tid};"
            records.append(cds_record)
            nextphase = - (nextphase + cds_end - cds_start + 1) % 3
    
    return records

def process_single_invalid_transcript(invalid_tid: str,
                                      invalid_tinfo: TranscriptInfo,
                                      genes: Dict[str, GeneInfo],
                                      trusted_transcripts: Dict[str, TranscriptInfo],
                                      trusted_exons: Dict[str, List[Tuple[int, int, str]]],
                                      genome_seq: Dict[str, str],
                                      min_orf_length: int = 300) -> Tuple[Optional[str], Optional[List[str]], Optional[str], Optional[str]]:
    """
    process single invalid transcript
    return: (gene_id, [record, ...], fixed_by_trusted_tid, strand) 
    if no associated gene found，return (None, None, None, None)
    if no satisfied transcript can be used to fix, return (gene_id, None, None, None)
    """
    logger.debug(f'Process {invalid_tid}: started')
    
    # get associated gene
    gene_id = None
    for gid, gene in genes.items():
        if invalid_tid in gene.transcripts:
            gene_id = gid
            break
    
    if not gene_id:
        logger.debug(f'Process {invalid_tid}: no associated gene found')
        return None, None, None, None
    logger.debug(f'Process {invalid_tid}: associated gene is {gid}')
    
    # find a best trusted transcript
    best_tid, best_orf, found_strand = find_best_transcript_for_invalid(
        invalid_tinfo, 
        trusted_transcripts, 
        trusted_exons, 
        genome_seq,
        min_orf_length
    )
    
    if not best_tid or not best_orf:
        logger.debug(f'Process {invalid_tid}: no satisfied transcript can be used to fix')
        return gene_id, None, None, None
    logger.debug(f'Process {invalid_tid}: found replacement ORF from {best_tid}; range: {best_orf[0]}-{best_orf[1]}; CDS length: {best_orf[6]}); strand: {found_strand}')
    
    # generate new record
    trusted_info = trusted_transcripts[best_tid]
    trusted_exons_list = trusted_exons[best_tid]
    
    new_records = create_corrected_transcript_records(
        invalid_tid, gene_id, best_tid, 
        trusted_info, trusted_exons_list, best_orf, found_strand
    )
    
    logger.debug(f'Process {invalid_tid}: finished')
    return gene_id, new_records, best_tid, found_strand

# ---------- Main process ----------
def process_invalid_transcripts(portcullis_file: str, 
                                trusted_file: str, 
                                genome_file: str,
                                gff3_file: str,
                                output_file: str,
                                min_orf_length: int = 300,
                                num_threads: int = 1):
    """Main process"""
    logger.debug("="*60)
    logger.debug("Main process started")
    logger.debug("="*60)
    
    logger.debug("Importing genome sequence...")
    genome_seq = readFastaAsDict(genome_file)
    logger.debug(f"  - {len(genome_seq)} sequences imported.")
    
    logger.debug("Parsing portcullis markup GTF...")
    portcullis_transcripts, invalid_transcripts = parse_portcullis_gtf(portcullis_file)
    logger.debug(f"  - {len(invalid_transcripts)} invalid transcripts found.")
    
    logger.debug("Parsing trusted GTF...")
    trusted_transcripts, trusted_exons = parse_trusted_gtf(trusted_file)
    logger.debug(f"  - {len(trusted_transcripts)} trusted transcripts imported.")
    
    logger.debug("Parsing reference GFF3...")
    genes, header = parse_gff3_to_genes(gff3_file)
    logger.debug(f"  - {len(genes)} genes imported.")
    
    # Core
    logger.debug(f"Processing invalid transcripts...")
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_tid = {
            executor.submit(
                process_single_invalid_transcript, 
                tid, 
                portcullis_transcripts[tid], 
                genes, 
                trusted_transcripts, 
                trusted_exons, 
                genome_seq, 
                min_orf_length
                ): tid for tid in invalid_transcripts
            }    
        
        processed_count = 0
        skipped_count = 0
        
        for future in concurrent.futures.as_completed(future_to_tid):
            tid = future_to_tid[future]
            try:
                result = future.result()
                if result:
                    gene_id, new_records, used_tid, found_strand = result
                    if gene_id:
                        if new_records:
                            results.append((tid, gene_id, new_records, used_tid, found_strand))
                            processed_count += 1
                        else:
                            skipped_count += 1
                    else:
                        skipped_count += 1
            except Exception as e:
                skipped_count += 1
                logger.exception(f"Process {tid}: skip process for unexpected error: {e}")
    
    logger.debug(f"All subprocess Completed.")
    logger.debug(f"  - {processed_count} invalid transcripts fixed.")
    logger.debug(f"  - {skipped_count} invalid transcripts remained.")
    
    # Replace invalid records
    logger.debug(f"Replacing invalid records...")
    for invalid_tid, gene_id, new_records, used_tid, found_strand in results:
        if gene_id in genes and new_records:
            # Remove old record
            if invalid_tid in genes[gene_id].transcripts:
                del genes[gene_id].transcripts[invalid_tid]
            
            # Add new record
            for record in new_records:
                fields = record.split('\t')
                if len(fields) >= 9:
                    feature = fields[2]
                    
                    attrs = {}
                    for pair in fields[8].split(';'):
                        if '=' in pair:
                            key, value = pair.split('=', 1)
                            attrs[key] = value
                    
                    record_id = attrs.get('ID')
                    parent = attrs.get('Parent')
                    
                    if feature == 'mRNA' and record_id:
                        genes[gene_id].add_record(record, 'mRNA', record_id, parent)
                    else:
                        genes[gene_id].add_record(record, feature, record_id, parent)
    
    # Write output file
    logger.debug(f'Writing output file {output_file} ...')
    with open(output_file, 'w') as f_out:
        for line in header:
            f_out.write(line + '\n')
        for gene in genes.values():
            for record in gene.get_all_records():
                f_out.write(record + '\n')
    
    logger.debug("="*60)
    logger.debug("ICU_fix completed.")
    logger.debug("="*60)
    logger.debug("Outputs:")
    logger.debug(f"  - Invalid transcript fixed intermediate GFF3: {output_file}")

# ---------- Entry ----------
def fix(inarg=None):
    '''
    Intron Correction Utility step 2: Fix invalid intron in transcripts by trusted annotation
    return fixed.gff3 (intermediate, only mRNA are adjusted, haven't intergrated into gene record)
    '''
    parser = argparse.ArgumentParser(
        description='Intron Correction Utility step 2: Fix invalid intron in transcripts by trusted annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    # Positional
    parser.add_argument('portcullis_gtf', help='Portcullis-[junctools gtf markup] marked transcript annotation, GTF format')
    parser.add_argument('trusted_gtf', help='trusted transcript annotation, GTF format')
    parser.add_argument('genome_fasta', help='Reference genome, FASTA format')
    parser.add_argument('input_gff3', help='Genome annotation, GFF3 format')
    
    # Optional
    parser.add_argument('-o', '--out', type=str, default='ICU.intermediate.gff3',
                       help='Generate output file named (default: ICU.intermediate.gff3)')
    parser.add_argument('--min-orf-length', type=int, default=300,
                       help='Avoid using ORF on trusted transcript shorter than (default: 300) bp to replace invalid transcript')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Threads limit to use. If set as 0, will use as much as possible (default: 1)')
    parser.add_argument('--debug', default=False, action='store_true',
                       help='Print DEBUG level logs')
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
    logger.debug("ICU_fix started")
    logger.debug("="*60)
    
    # Check input
    logger.debug("Checking input...")
    input_files = [
        args.portcullis_gtf,
        args.trusted_gtf,
        args.genome_fasta,
        args.input_gff3
    ]
    for f in input_files:
        if not os.path.exists(f):
            logger.error(f"Cannot find file: {f}")
            sys.exit(1)
        logger.debug(f"File found: {f}")
        
    logger.debug(f"Input:")
    logger.debug(f"  - Portcullis markup GTF: {args.portcullis_gtf}")
    logger.debug(f"  - Trusted annotation GTF: {args.trusted_gtf}")
    logger.debug(f"  - Genome FASTA: {args.genome_fasta}")
    logger.debug(f"  - Genome annotation GFF3: {args.input_gff3}")
        
    num_threads = os.cpu_count() if args.threads == 0 else args.threads  

    logger.debug(f"Options:")
    logger.debug(f"  - Output file: {args.out}")
    logger.debug(f"  - Minimal ORF length: {args.min_orf_length}")
    logger.debug(f"  - Threads: {num_threads}")      
    logger.debug(f"  - debug log: {args.debug}")
    
    # Main process entry
    try:
        process_invalid_transcripts(
            args.portcullis_gtf,
            args.trusted_gtf,
            args.genome_fasta,
            args.input_gff3,
            args.out,
            args.min_orf_length,
            num_threads
        )
        return args.out

    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(1)
    

if __name__ == "__main__":
    fix()