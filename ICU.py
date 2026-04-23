#!/usr/bin/env python3
"""
Intron Correction Utility
"""
import ICU_prepare
import ICU_fix
import ICU_integrate
import ICU_stat

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Intron Correction Utility v0.0.2',
        prog='ICU.py', 
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument('subcommand', choices=['prepare', 'fix', 'integrate', 'stat', 'pipe'], help=
'''Availble subcommand:
    prepare   run portcullis and stringtie
    fix       fix invalid intron in transcripts by trusted annotation
    integrate adjust gene record according to changed transcripts
    stat      run junctools again and generate statistical data
    pipe      run all four steps as pipeline
''')
    
    args, command = parser.parse_known_args()
    
    subcommandref = {'prepare': ICU_prepare.prepare, 'fix': ICU_fix.fix, 'integrate': ICU_integrate.integrate, 'stat': ICU_stat.stat}
    if args.subcommand in subcommandref:
        subcommandref[args.subcommand](command)
    elif args.subcommand == 'pipe':
        pipeparser = argparse.ArgumentParser(
            description='Intron Correction Utility Pipeline', 
            prog='ICU.py pipe',
            formatter_class=argparse.RawTextHelpFormatter
        )
        # Positional
        pipeparser.add_argument('genome_fasta', help='Reference genome, FASTA format')
        pipeparser.add_argument('input_gff3', help='Genome annotation, GFF3 format')
        pipeparser.add_argument('bam_files', nargs='+', help='splice-alignment(s) against reference genome, BAM format')
        
        # Optional
        # Force str to pass to sub argument parser
        pipeparser.add_argument('--min-junc', type=str, default='10',
                        help='Ignore splice-alignment that supported by lower than (default: 10) read pairs')
        pipeparser.add_argument('--min-orf-length', type=str, default='300',
                        help='Avoid using ORF on StringTie assembled transcript shorter than (default: 300) bp to replace invalid transcript')
        pipeparser.add_argument('--rename', type=str, default='',
                        help='Rename all gene and its child record ID begin with prefix+number; if not set, join original ID with operation (default: not set)')
        pipeparser.add_argument('--prefix', type=str, default='ICU',
                        help='Generate output file named started with (default: ICU)')
        pipeparser.add_argument('-p', '--process', type=str, default='1',
                        help='Process limit for stringtie and portcullis.  If set as 0, will use as much as possible (default: 1)')
        pipeparser.add_argument('--debug', default=False, action='store_true',
                        help='Print DEBUG level logs')
        
        # Dependency, default $PATH
        pipeparser.add_argument('--portcullis', help='portcullis executable path, default: $PATH')
        pipeparser.add_argument('--stringtie', help='stringtie executable path, default: $PATH')
        pipeparser.add_argument('--junctools', help='junctools executable path, default: $PATH')
        pipeparser.add_argument('--gffread', help='gffread executable path, default: $PATH')
        
        pipeargs = pipeparser.parse_args(command)
        
        stringtie_gtf, portcullis_gtf, portcullis_out = ICU_prepare.prepare(
            [pipeargs.genome_fasta, pipeargs.input_gff3] + pipeargs.bam_files + 
            ['--min-junc', pipeargs.min_junc] +
            ['--prefix', pipeargs.prefix] +
            ['-p', pipeargs.process] +
            (['--debug'] if pipeargs.debug else []) +
            ['--portcullis', pipeargs.portcullis] +
            ['--stringtie', pipeargs.stringtie] +
            ['--junctools', pipeargs.junctools] +
            ['--gffread', pipeargs.gffread]
        )
        intermediate_gff3 = ICU_fix.fix(
            [portcullis_gtf, stringtie_gtf, pipeargs.genome_fasta, pipeargs.input_gff3] +
            ['-o', f'{pipeargs.prefix}.intermediate.gff3'] +
            ['--min-orf-length', pipeargs.min_orf_length] +
            (['--debug'] if pipeargs.debug else [])
        )
        corrected_gff3 = ICU_integrate.integrate(
            [intermediate_gff3] +
            ['-o', f'{pipeargs.prefix}.integrated.gff3'] +
            ['--rename', pipeargs.rename] +
            (['--debug'] if pipeargs.debug else [])
        )
        ICU_stat.stat(
            [pipeargs.input_gff3, corrected_gff3, portcullis_out] +
            ['--prefix', pipeargs.prefix] +
            (['--debug'] if pipeargs.debug else []) +
            ['--junctools', pipeargs.junctools] +
            ['--gffread', pipeargs.gffread]
        )

if __name__ == "__main__":
    main()