#!/usr/bin/env python3
"""
Intron Correction Utility step 4: Run junctools again and generate statistical data
"""

import subprocess
import sys
import os
import argparse
import logging
import re

def check_dependency(tool_name, tool_path=None):
    """
    If path not provided, try $PATH.
    Return an executable. 
    """
    if tool_path:
        if os.path.exists(tool_path) and os.access(tool_path, os.X_OK):
            return tool_path
        else:
            raise FileNotFoundError(f"{tool_name} not found or not executable: {tool_path}")
    else:
        try:
            result = subprocess.run(['which', tool_name], 
                                   capture_output=True, text=True, check=False)
            if result.returncode == 0:
                return tool_name
            else:
                raise FileNotFoundError(f"{tool_name} not found in $PATH")
        except Exception as e:
            raise FileNotFoundError(f"Unexpected error when checking dependency {tool_name}: {e}")

def count_gene_mrna(gff3_path):
    """
    return gene_count, mrna_count
    """
    gene_count = 0
    mrna_count = 0
    
    with open(gff3_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            feature_type = parts[2].lower()
            
            if feature_type == 'gene':
                gene_count += 1
            elif feature_type in ['mrna', 'transcript']:
                mrna_count += 1
    return gene_count, mrna_count

def run_command(cmd, description, shell=False):
    logger.debug(f"subprocess {description} started, CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    try:
        if shell:
            process = subprocess.run(cmd, shell=True, check=True, 
                                   capture_output=True, text=True)
        else:
            process = subprocess.run(cmd, check=True, 
                                   capture_output=True, text=True)
        logger.debug(f"subprocess {description} completed")
        output = process.stdout
        return output
    except subprocess.CalledProcessError as e:
        logger.error(f"subprocess {description} failed (CalledProcessError), return code: {e.returncode}")
        logger.error(f"subprocess {description} error info: {e}")
        raise RuntimeError(f"subprocess {description} failed") from e
    except Exception as e:
        logger.exception(f"subprocess {description} failed (unexpected error), error info: {e}")
        raise RuntimeError(f"subprocess {description} failed") from e

def run_gffread(annotation_gff3, prefix, gffread_path, overwrite):
    """
    Run gffread, convert GFF3 annotation to GTF format
    return annotation.gtf
    """
    gtf_file = f"{prefix}.annotation.gtf"
    cmd = [
        gffread_path,
        '-T',
        annotation_gff3,
        '-o', gtf_file,
        '--force-exons' # if gff3 has only CDS record without exon, junctools gtf markup will fail.
    ]
    if not os.path.exists(gtf_file) or overwrite:
        run_command(cmd, "gffread")
    else:
        logger.warning(f"{gtf_file} existed. Try to reuse.")
    return gtf_file

def run_junctools(annotation_gtf, portcullis_out, prefix, junctools_path):
    """
    Run junctools gtf markup, 
    return output string
    """
    junctions_file = os.path.join(portcullis_out, 
                                  "3-filt", 
                                  "portcullis_filtered.pass.junctions.tab")
    
    if not os.path.exists(junctions_file):
        raise FileNotFoundError(f"portcullis_filtered.pass.junctions.tab not found in {portcullis_out}")
    
    output_gtf = f"{prefix}.markup.gtf"
    cmd = [
        junctools_path, 'gtf', 'markup',
        annotation_gtf,
        '-j', junctions_file,
        '-o', output_gtf
    ]
    
    output = run_command(cmd, "junctools gtf markup")
    return output

def extract_values_from_stdout(stdout_text):
    """
    extract statistical data from junctools gtf markup stdout
    return supported_junction, total_junction, invalid_multi_exonic_transcripts, total_multi_exonic_transcripts
    """
    junction_pattern = r'(\d+)\s*/\s*(\d+)\s+\(\d+\.?\d*\)\s+junctions supported'
    junction_match = re.search(junction_pattern, stdout_text)
    
    invalid_pattern = r'(\d+)\s*/\s*(\d+)\s+invalid multi-exonic'
    invalid_match = re.search(invalid_pattern, stdout_text)
    
    supported_junction = int(junction_match.group(1))
    total_junction = int(junction_match.group(2))

    invalid_multi_exonic_transcripts = int(invalid_match.group(1))
    total_multi_exonic_transcripts = int(invalid_match.group(2))
    
    return supported_junction, total_junction, invalid_multi_exonic_transcripts, total_multi_exonic_transcripts

def stat(inarg=None):
    '''
    Intron Correction Utility step 4: Run junctools again and generate statistical data
    return statfile
    '''
    parser = argparse.ArgumentParser(
        description='Intron Correction Utility step 4: Run junctools again and generate statistical data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    # Positional
    parser.add_argument('previous_gff3', help='Genome annotation (before process), GFF3 format')
    parser.add_argument('corrected_gff3', help='Genome annotation (after process), GFF3 format')
    parser.add_argument('portcullis_out', help='portcullis output folder')
    
    # Optional
    parser.add_argument('--prefix', type=str, default='ICU',
                       help='Generate output file named started with (default: ICU)')
    parser.add_argument('--debug', default=False, action='store_true',
                       help='Print DEBUG level logs')
    
    # Dependency, default $PATH
    parser.add_argument('--junctools', help='junctools executable path, default: $PATH')
    parser.add_argument('--gffread', help='gffread executable path, default: $PATH')
    
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
    logger.debug("ICU_stat started")
    logger.debug("="*60)
    
    logger.debug(f"Input:")
    logger.debug(f"  - Previous_gff3: {args.previous_gff3}")
    logger.debug(f"  - Corrected_gff3: {args.corrected_gff3}")
    logger.debug(f"  - portcullis_out: {args.portcullis_out}")   
    logger.debug(f"Options:")  
    logger.debug(f"  - Prefix: {args.prefix}")
    logger.debug(f"  - debug log: {args.debug}")

    # Check dependency
    logger.debug("Checking dependency...")
    try:
        junctools_path = check_dependency('junctools', args.junctools)
        gffread_path = check_dependency('gffread', args.gffread)
        
        logger.debug(f"  - junctools: {junctools_path}")
        logger.debug(f"  - gffread: {gffread_path}")
    except FileNotFoundError as e:
        logger.error(f"Fail to check dependency: {e}")
        sys.exit(1)
    
    # Check input
    logger.debug("Checking input...")
    required_files = [args.previous_gff3, args.corrected_gff3]
    for file_path in required_files:
        if not os.path.exists(file_path):
            logger.error(f"  - Cannot find file: {file_path}")
            sys.exit(1)
        logger.debug(f"  - File found: {file_path}")
    
    # Main process
    logger.debug("="*60)
    logger.debug(f"Main process started")
    logger.debug("="*60)
    
    try:
        previous_genecount, previous_mrnacount = count_gene_mrna(args.previous_gff3)
        previous_gtf = run_gffread(args.previous_gff3, args.prefix, gffread_path, False)  
        previous_stat = run_junctools(previous_gtf, args.portcullis_out, args.prefix, junctools_path)
        previous_supported_junction, previous_total_junction, previous_invalid_transcripts, previous_total_transcripts = extract_values_from_stdout(previous_stat)
        
        corrected_genecount, corrected_mrnacount = count_gene_mrna(args.corrected_gff3)
        corrected_gtf = run_gffread(args.corrected_gff3, args.prefix+'.corrected', gffread_path, True)  
        corrected_stat = run_junctools(corrected_gtf, args.portcullis_out, args.prefix+'.corrected', junctools_path)
        corrected_supported_junction, corrected_total_junction, corrected_invalid_transcripts, corrected_total_transcripts = extract_values_from_stdout(corrected_stat)
    
        statfile = f'{args.prefix}.stat'
        with open(statfile, 'w') as s:
            logger.debug(f"Gene counts: {previous_genecount} -> {corrected_genecount}")
            s.write(f'Gene counts: {previous_genecount} -> {corrected_genecount}\n')
            logger.debug(f'mRNA counts: {previous_mrnacount} -> {corrected_mrnacount}')
            s.write(f'mRNA counts: {previous_mrnacount} -> {corrected_mrnacount}\n')
            
            previous_supported_junction_ratio = round(previous_supported_junction / previous_total_junction * 100, 2)
            corrected_supported_junction_ratio = round(corrected_supported_junction / corrected_total_junction * 100, 2)
            logger.debug(f'Supported junctions: {previous_supported_junction} / {previous_total_junction} ({previous_supported_junction_ratio}%) ' +
                         f'-> {corrected_supported_junction} / {corrected_total_junction} ({corrected_supported_junction_ratio}%)')
            s.write(f'Supported junctions: {previous_supported_junction} / {previous_total_junction} ({previous_supported_junction_ratio}%) ' +
                    f'-> {corrected_supported_junction} / {corrected_total_junction} ({corrected_supported_junction_ratio}%)\n')
            
            previous_valid_transcripts = previous_total_transcripts - previous_invalid_transcripts
            previous_valid_transcripts_ratio = round(previous_valid_transcripts / previous_total_transcripts * 100, 2)
            corrected_valid_transcripts = corrected_total_transcripts - corrected_invalid_transcripts
            corrected_valid_transcripts_ratio = round(corrected_valid_transcripts / corrected_total_transcripts * 100, 2)
            
            logger.debug(f'Valid multi-exonic transcripts: {previous_valid_transcripts} / {previous_total_transcripts} ({previous_valid_transcripts_ratio}%) ' +
                         f'-> {corrected_valid_transcripts} / {corrected_total_transcripts} ({corrected_valid_transcripts_ratio}%)')
            s.write(f'Valid multi-exonic transcripts: {previous_valid_transcripts} / {previous_total_transcripts} ({previous_valid_transcripts_ratio}%) ' +
                    f'-> {corrected_valid_transcripts} / {corrected_total_transcripts} ({corrected_valid_transcripts_ratio}%)\n')    
    
    except Exception as e:
        logger.exception(f"Unexpected error oucurred in main process: {e}")
        sys.exit(1)
    
    logger.debug("="*60)
    logger.debug("ICU_stat completed.")
    logger.debug("="*60)
    logger.debug("Outputs:")
    logger.debug(f"  - statistical data: {statfile}")
    
    return statfile

if __name__ == "__main__":
    stat()