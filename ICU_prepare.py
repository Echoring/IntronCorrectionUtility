#!/usr/bin/env python3
"""
Intron Correction Utility step 1: Run portcullis and stringtie
"""

import subprocess
import sys
import os
import argparse
import logging

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

def run_command(cmd, description, shell=False):
    logger.debug(f"subprocess {description} started, CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    
    try:
        if shell:
            process = subprocess.Popen(
                cmd, 
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8',
                errors='replace',
                bufsize=1,
                universal_newlines=True
            )
        else:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8',
                errors='replace',
                bufsize=1,
                universal_newlines=True
            )
        
        output_lines = []
        for line in iter(process.stdout.readline, ''):
            if line:
                clean_line = line.rstrip('\n')
                logger.debug(f"[{description} output] {clean_line}")
                output_lines.append(clean_line)
        
        process.stdout.close()
        return_code = process.wait()
        
        if return_code == 0:
            logger.debug(f"subprocess {description} completed")
            return True
        else:
            output = "\n".join(output_lines)
            logger.error(f"subprocess {description} failed (non 0 return), return code: {return_code}")
            logger.error(f"subprocess {description} outputs:\n{output}")
            
            raise subprocess.CalledProcessError(return_code, cmd, output=output)
            
    except subprocess.CalledProcessError as e:
        logger.error(f"subprocess {description} failed (CalledProcessError), return code: {e.returncode}")
        logger.error(f"subprocess {description} error info: {e}")
        if hasattr(e, 'output') and e.output:
            logger.error(f"subprocess {description} outputs:\n{e.output}")
        raise RuntimeError(f"subprocess {description} failed") from e
        
    except Exception as e:
        logger.exception(f"subprocess {description} failed (unexpected error), error info: {e}")
        raise RuntimeError(f"subprocess {description} failed") from e

def run_stringtie(bam_file, min_junc, prefix, process, stringtie_path):
    """Return stringtie.gtf"""
    output_gtf = f"{prefix}_{os.path.basename(bam_file)}.stringtie.gtf"
    cmd = [
        stringtie_path,
        bam_file,
        '-o', output_gtf,
        '-p', str(process),
        '--conservative',
        '-l', f'{prefix}_STRG',
        '-j', str(min_junc)
    ]
    if not os.path.exists(output_gtf):
        run_command(cmd, "stringtie")
    else:
        logger.warning(f"{output_gtf} existed. Try to reuse.")
    
    return output_gtf

def run_stringtie_pipeline(bam_files, min_junc, prefix, process, stringtie_path):
    "Run stringtie for each bam input, then run stringtie merge, return merged.gtf"
    allgtf = []
    merged_gtf = f"{prefix}.merged.stringtie.gtf"

    if not os.path.exists(merged_gtf):
        for bam in bam_files:
            onegtf = run_stringtie(bam, min_junc, prefix, process, stringtie_path)
            allgtf.append(onegtf)
        cmd = [
            stringtie_path,
            '--merge',
            '-o', merged_gtf,
            '-l', f'{prefix}_MSTRG',
            '-f', '0.5',
            *allgtf
        ]
        run_command(cmd, "stringtie --merge")
    else:
        logger.warning(f"{merged_gtf} existed. Try to reuse.")  
    logger.debug("Stringtie pipeline completed.")
    return merged_gtf
    

def run_portcullis(genome_fasta, bam_files, min_junc, prefix, process, portcullis_path):
    """
    Run portcullis full
    return output folder
    """
    outfolder = f'{prefix}_portcullis_out'
    cmd = [
        portcullis_path, 'full',
        '-t', str(process),
        '-o', f'{prefix}_portcullis_out',
        '--min_cov', str(min_junc),
        genome_fasta
    ]
    cmd.extend(bam_files)
    if not os.path.exists(outfolder):
        run_command(cmd, "portcullis full")
    else:
        logger.warning(f"{outfolder} existed. Try to reuse.")
    return outfolder

def run_gffread(annotation_gff3, prefix, gffread_path):
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
    if not os.path.exists(gtf_file):
        run_command(cmd, "gffread")
    else:
        logger.warning(f"{gtf_file} existed. Try to reuse.")
    return gtf_file

def run_junctools(annotation_gtf, portcullis_out, prefix, junctools_path):
    """
    Run junctools gtf markup, 
    return markup.gtf
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
    if not os.path.exists(output_gtf):
        run_command(cmd, "junctools gtf markup")
    else:
        logger.warning(f"{output_gtf} existed. Try to reuse.")
    return output_gtf

def run_portcullis_pipeline(genome_fasta, annotation_gff3, bam_files, min_junc, prefix,
                           process, portcullis_path, gffread_path, junctools_path):
    """
    Run portcullis full, gffread, junctools gtf markup, 
    return portcullis_out(folder), annotation.gtf, markup.gtf
    """
    portcullis_out = run_portcullis(genome_fasta, bam_files, min_junc, prefix, process, portcullis_path)
    annotation_gtf = run_gffread(annotation_gff3, prefix, gffread_path)
    markup_gtf = run_junctools(annotation_gtf, portcullis_out, prefix, junctools_path)
    
    logger.debug("Portcullis pipeline completed.")
    return portcullis_out, annotation_gtf, markup_gtf

# ---------- Main ----------
def prepare(inarg=None):
    '''
    Intron Correction Utility step 1: Run portcullis and stringtie
    return stringtie.merged.gtf, invalid.markup.gtf, portcullis_out(folder)
    '''
    parser = argparse.ArgumentParser(
        description='Intron Correction Utility step 1: Run portcullis and stringtie',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    # Positional
    parser.add_argument('genome_fasta', help='Reference genome, FASTA format')
    parser.add_argument('input_gff3', help='Genome annotation, GFF3 format')
    parser.add_argument('bam_files', nargs='+', help='splice-alignment(s) against reference genome, BAM format')
    
    # Optional
    parser.add_argument('--min-junc', type=int, default=10,
                       help='Ignore splice-alignment that supported by lower than (default: 10) read pairs')
    parser.add_argument('--prefix', type=str, default='ICU',
                       help='Generate output file named started with (default: ICU)')
    parser.add_argument('-p', '--process', type=int, default=1,
                       help='Process limit to use for stringtie and portcullis. If set as 0, will use as much as possible (default: 1)')
    parser.add_argument('--debug', default=False, action='store_true',
                       help='Print DEBUG level logs')
    
    # Dependency, default $PATH
    parser.add_argument('--portcullis', help='portcullis executable path, default: $PATH')
    parser.add_argument('--stringtie', help='stringtie executable path, default: $PATH')
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
    logger.debug("ICU_prepare started")
    logger.debug("="*60)
    
    logger.debug(f"Input:")
    logger.debug(f"  - Genome: {args.genome_fasta}")
    logger.debug(f"  - Annotation: {args.input_gff3}")
    logger.debug(f"  - Splice-alignment: {', '.join(args.bam_files)}")   
    logger.debug(f"Options:")  
    logger.debug(f"  - minimal junction-support read pairs: {args.min_junc}")
    logger.debug(f"  - prefix: {args.prefix}")
    logger.debug(f"  - process: {args.process}")
    logger.debug(f"  - debug log: {args.debug}")
        
    # Check dependency
    logger.debug("Checking dependency...")
    try:
        portcullis_path = check_dependency('portcullis', args.portcullis)
        stringtie_path = check_dependency('stringtie', args.stringtie)
        junctools_path = check_dependency('junctools', args.junctools)
        gffread_path = check_dependency('gffread', args.gffread)
        
        logger.debug(f"  - portcullis: {portcullis_path}")
        logger.debug(f"  - stringtie: {stringtie_path}")
        logger.debug(f"  - junctools: {junctools_path}")
        logger.debug(f"  - gffread: {gffread_path}")
    except FileNotFoundError as e:
        logger.error(f"Fail to check dependency: {e}")
        sys.exit(1)

    # Check input
    logger.debug("Checking input...")
    required_files = [args.genome_fasta, args.input_gff3] + args.bam_files
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
        num_process = os.cpu_count() if args.process == 0 else args.process
        stringtie_gtf = run_stringtie_pipeline(args.bam_files, args.min_junc, args.prefix, num_process, stringtie_path)
        portcullis_out, annotation_gtf, markup_gtf = run_portcullis_pipeline(args.genome_fasta, args.input_gff3, args.bam_files, args.min_junc, args.prefix, num_process, portcullis_path, gffread_path, junctools_path)
            
    except Exception as e:
        logger.exception(f"Unexpected error oucurred in main process: {e}")
        sys.exit(1)
    
    logger.debug("="*60)
    logger.debug("ICU_prepare completed.")
    logger.debug("="*60)
    logger.debug("Major outputs:")
    logger.debug(f"  - stringtie merged GTF: {stringtie_gtf}")
    logger.debug(f"  - junctools markup GTF: {markup_gtf}")
    logger.debug(f"  - portcullis output DIR: {portcullis_out}")
    
    return stringtie_gtf, markup_gtf, portcullis_out
    
if __name__ == "__main__":
    prepare()