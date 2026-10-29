## All pipeline functions

## ----- Libraries ----- ##
## General functionality
import os
import sys
import re
import subprocess
import logging
from pathlib import Path
import threading
import glob
import requests
import collections
import hashlib
import time
import warnings
import random
import gzip

## Multiprocessing
import multiprocessing
import functools

## General bioinformatics
import pysam
import polars as pl
import numpy as np
import pandas as pd

## Specialized functionality
## Note: this list does not include packages called in the command line
import editdistance
import parasail
import rapidfuzz
import itertools
from umi_tools import UMIClusterer

## Importing from the lightweight EPI2ME script
from ._epi2me import parse_vsearch, create_stranded_reads

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################
## ----- Basic utility ----- ##
###############################

def setup_workdir(
    work_dir: str | Path | None = None,
    clean_start: bool = False,
    log_file_name: str = 'log.log',
):
    r"""Generic function for creating working directory and logger.

    Parameters
    ----------
    work_dir
        String or Path specifying working directory. Defaults to `'~/WORK'`.
    clean_start:
        Boolean determining if log file should be clean. Otherwise, logger appends to the log file. 
    log_file_name:
        String specifying log file name. Defaults to `log.log`.
    """
    if work_dir is None:
        work_dir = f'{os.path.expanduser('~')}/WORK'
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    logger = logging.getLogger(__name__) 
    logging.basicConfig(
        filename = f'{work_dir}/{log_file_name}',
        filemode = 'a',
        level = logging.DEBUG, 
        format = '%(asctime)s - %(funcName)s - %(levelname)s - %(message)s')
    if clean_start:
        if os.path.exists(f'{work_dir}/{log_file_name}'):
            os.remove(f'{work_dir}/{log_file_name}')
    if not os.path.exists(f'{work_dir}/{log_file_name}'):
        Path(f'{work_dir}/{log_file_name}').touch()
    return logger, work_dir


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################################
## ----- Step 1: Genome handling ----- ##
#########################################

def get_genome(
    genome_dir: str | Path,
    genome_fasta: str | Path | None = None,
    genome_gtf: str | Path | None = None,
    threads: int | None = None,
    work_dir: str | Path | None = None,
):
    r"""Get genome information from user-specified folder.

    Parameters
    ----------
    genome_dir
        String or Path indicating where folder with genome information is.
    genome_fasta
        String or Path indicating which file is our genome fasta file. If not specified, inferred from `genome_dir`.
    genome_gtf
        String or Path indicating which file is our gene annotation file. If not specified, inferred from `genome_dir`.
    threads
        Integer indicating number of threads for `minimap2` to build the `.mmi` index. Defaults to `threading.active_count()`.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    """
    ## Creating working directory and start logger
    logger, work_dir = setup_workdir(work_dir)
    
    ## If files not specified, infer
    if genome_fasta is None or genome_gtf is None:
        if not os.path.exists(genome_dir):
            logger.error(f"{genome_dir} does not exist")
            sys.exit(1)
        genome_files = []
        for root, _, files in os.walk(genome_dir):
            for file_name in files:
                genome_files.append(os.path.join(root, file_name))
    if genome_fasta is None:
        genome_fasta_list = [x for x in genome_files if re.search('fa.gz$|fa$|fna$|fna.gz$', x)]
        genome_fasta = genome_fasta_list[0]
        if len(genome_fasta_list) > 1:
            logger.warning(f'More than one possible fasta file detected. Choosing first: {genome_fasta}')
    if genome_gtf is None:
        genome_gtf_list = [x for x in genome_files if re.search('gtf$|gtf.gz', x)]
        genome_gtf = genome_gtf_list[0]
        if len(genome_gtf_list) > 1:
            logger.warning(f'More than one possible gtf file detected. Choosing first: {genome_gtf}')
    if not os.path.exists(genome_fasta):
        logger.error(f'{genome_fasta} does not exist')
        sys.exit(1)
    if not os.path.exists(genome_gtf):
        logger.error(f'{genome_gtf} does not exist')
        sys.exit(1)

    ## Prep the files
    out_fasta = f"{work_dir}/GENOME.fa"
    out_faidx = f"{work_dir}/GENOME.fa.fai"
    out_chrsizes = f"{work_dir}/GENOME.chr.sizes"
    out_gtf = f"{work_dir}/GENOME.gtf"
    out_mmi = f"{work_dir}/GENOME.mmi"
    out_gtfbed = f"{work_dir}/GENOME.bed"
    out_genecsv = f"{work_dir}/genes.csv"
    logger.info(f'Moving relevant genome files to work directory')
    if threads is None:
        threads = threading.active_count()
        
    if not os.path.exists(out_gtf):
        logger.info(f'Copying gtf into work directory')
        if bool(re.search('.gz$', genome_gtf)):
            proc = subprocess.run(f'zcat {genome_gtf} > {out_gtf}', shell=True, capture_output=True)
        else:
            proc = subprocess.run(f'cat {genome_gtf} > {out_gtf}', shell=True, capture_output=True)
        ## And now provide chromosome-specific gtfs for later . 
        ## Annoyingly for Stringtie, this needs to only contain entries where transcript_id="" are removed
        ## See https://www.biostars.org/p/488062/ for how to do this with awk: '$3 != "gene" '
        ## However, it turns out that the issue is that T2T genome has other entries beyond genes that also have transcript_id=""
        ## Found a workaround with grep
        grep_string = 'transcript_id "";'
        proc1 = subprocess.Popen(['grep', '-v', grep_string, out_gtf], stdout = subprocess.PIPE)
        awk_string = '/^[^#]/ {print > (out_dir "/" $1 "chr.gtf") }'
        proc2 = subprocess.Popen(['awk', '-v', f'out_dir={os.path.dirname(out_gtf)}', awk_string], stdin=proc1.stdout)
        proc1.stdout.close()
        proc2.communicate()
            
    if not os.path.exists(out_fasta):
        logger.info(f'Copying reference fasta into work directory')
        if bool(re.search('.gz$', genome_fasta)):
            proc = subprocess.run(f'zcat {genome_fasta} > {out_fasta}', shell=True, capture_output=True)
        else:
            proc = subprocess.run(f'cat {genome_fasta} > {out_fasta}', shell=True, capture_output=True)
            
    if not os.path.exists(out_faidx):
        logger.info(f'Indexing reference fasta')
        proc = subprocess.run(f'samtools faidx {out_fasta} --output {out_faidx}', shell=True, capture_output=True)

    if not os.path.exists(out_chrsizes):
        logger.info(f'Getting chromosome sizes')
        proc = subprocess.run(f'cat {out_faidx} | cut -f1,2 > {out_chrsizes}', shell=True, capture_output=True)
        
    if not os.path.exists(out_mmi):
        logger.info(f'Building .mmi index from reference fasta')
        ## From https://github.com/epi2me-labs/wf-single-cell/blob/master/modules/local/common.nf
        ## Caps at 16GB per batch
        proc = subprocess.run(f'minimap2 -t {threads} -I 16G -d "{out_mmi}" "{out_fasta}"', shell=True, capture_output=True)

    if not os.path.exists(out_gtfbed):
        logger.info(f'Converting .gtf to bed file for minimap2')
        ## paftools.js comes with minimap2 installation
        proc = subprocess.run(f'paftools.js gff2bed -j {out_gtf} > {out_gtfbed}', shell=True, capture_output=True)

    if not os.path.exists(out_genecsv):
        logger.info(f'Consolidating gene information from bed file')
        gene_info = pl.read_csv(out_gtfbed, separator = '\t', has_header=False, 
                        new_columns=['chr', 'start', 'end', 'meta', 'score', 'strand'])['meta'].to_list()
        gene_info = pl.DataFrame(
            [{'id': x.split('|')[0], 'biotype': x.split('|')[1], 'gene': x.split('|')[2]} for x in gene_info]
        ).unique().sort('gene')
        gene_info.write_csv(out_genecsv, separator=',')
        
    return


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################
## ----- Step 2: Cellular barcode extraction ----- ##
#####################################################

def fuzz_assign(
    barcode: str,
    whitelist: list,
    max_barcode_fuzzdist: int = 2,
    mindelta_barcode_fuzzdist: int = 2,
):
    r"""Rapidfuzz-based assignment of barcodes based on whitelist.

    Realised I was using this loop a lot, so have made this its own function.

    Parameters
    ----------
    barcode
        String indicating barcode.
    whitelist
        List of strings, indicating universe of possible sequences.
    max_barcode_fuzzdist
        Integer indicating max tolerable edit distance between barcode and whitelist.
    mindelta_barcode_fuzzdist
        Integer indicating min difference in edit distances between best and second best whitelist match. Indicates how certain top match is true.
    """
    ## assigned_barcode only gets overwritten if successfully passes one of the criteria
    assigned_barcode = ""
    barcode_corrected = False

    if barcode in whitelist:
        ## Check if barcode is already in whitelist
        assigned_barcode = barcode
    else:
        ## If not, then do a fuzzy match
        fuzzy_match_results = rapidfuzz.process.extract(
            barcode, 
            whitelist,
            scorer = rapidfuzz.distance.Levenshtein.distance,
            score_cutoff = max_barcode_fuzzdist + mindelta_barcode_fuzzdist + 1
        )
        
        if len(fuzzy_match_results)==0:
            ## No confident barcode assigned, drop this read (return "")
            return assigned_barcode, barcode_corrected
            
        if len(fuzzy_match_results)==1:
            if fuzzy_match_results[0][1] > max_barcode_fuzzdist:
                ## Match exceeds max score, drop (return "")
                return assigned_barcode, barcode_corrected
            assigned_barcode = fuzzy_match_results[0][0]
            barcode_corrected = True
            
        if len(fuzzy_match_results) > 1:
            if fuzzy_match_results[0][1] > max_barcode_fuzzdist:
                ## Top match exceeds max score, drop (return "")
                return assigned_barcode, barcode_corrected
            diff_scores = fuzzy_match_results[1][1] - fuzzy_match_results[0][1]
            if diff_scores < mindelta_barcode_fuzzdist:
                ## Difference between 2 assignments too small, drop
                return assigned_barcode, barcode_corrected
            assigned_barcode = fuzzy_match_results[0][0]
            barcode_corrected = True
    
    return assigned_barcode, barcode_corrected

## ~~~ ##

def get_10x_whitelist(
    kit: str,
    whitelist_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
):
    r"""Function to download relevant 10X barcode whitelist from the EPI2ME github.

    Pulled this out of extract_10x_barcodes() to avoid clashes between workers when running multiprocessing.

    Parameters
    ----------
    kit
        String indicating what kit has been used.
    whitelist_dir
        String or Path indicating where whitelist directory is. Defaults to `'~/WORK'`.
        Alternatively, if a file is supplied, will assume that the file is the whitelist file itself.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    """
    logger, work_dir = setup_workdir(work_dir)
    
    ## Whitelist: see https://kb.10xgenomics.com/s/article/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist
    ## Note, also won't bother with Visium here
    kit_whitelists = {
        '3prime:v4': ['3M-3pgex-may-2023.txt.gz', '3M-3pgex-may-2023_TRU.txt.gz'],
        '3prime:v3': ['3M-febrary-2018.txt.gz', '3M-february-2018_TRU.txt.gz'],
        '3prime:v2': ['737k-august-2016.txt'],
        '3prime:v1': ['737k-april-2014_rc.txt'],
        
        '5prime:v3': ['3M-5pgex-jan-2023.txt.gz'],
        '5prime:v2': ['737k-august-2016.txt'],
        '5prime:v1': ['737k-august-2016.txt'],

        # 'multiome:v1': ['737k-arc-v1.txt.gz'],
        # 'atac:v1': ['737-cratac-v1.txt.gz']
    }
    allowed_kits = kit_whitelists.keys()
    if kit not in allowed_kits:
        logger.error(f'Unsupported kit specified ({kit}). We only accept: {', '.join(list(allowed_kits))}')
        sys.exit(1)
    whitelists = kit_whitelists[kit]
    if whitelist_dir is None:
        whitelist_dir = work_dir
    whitelist_ready = False
    chosen_whitelist = ''
    for whitelistx in whitelists:
        if whitelist_ready:
            continue
            
        ## If is file
        if not os.path.isdir(whitelist_dir):
            if not os.path.isfile(whitelist_dir):
                logger.error(f'Unable to determine if {whitelist_dir} is a directory or a file.')
                sys.exit(1)
            whitelist_ready = True
            chosen_whitelist = f'{whitelist_dir}'
            continue
            
        ## If is directory
        if os.path.exists(f'{whitelist_dir}/{whitelistx}'):
            whitelist_ready = True
            chosen_whitelist = f'{whitelist_dir}/{whitelistx}'
            continue
        repo_url = 'https://github.com/epi2me-labs/wf-single-cell/raw/refs/heads/master/data/'
        url = f'{repo_url}/{whitelistx}'
        wget_result = subprocess.run(f'wget --spider -q --method=HEAD {url}', shell=True, check=False, text=True, capture_output=True)
        if wget_result.returncode != 0:
            continue
        logger.info(f'Downloading whitelist ({whitelistx})...')
        proc = subprocess.run(f'wget -q -nc -P {whitelist_dir} {url}', shell=True, capture_output=True)
        whitelist_ready = True
        chosen_whitelist = f'{whitelist_dir}/{whitelistx}'
    if not whitelist_ready:
        logger.error(f'Unable to find whitelist')
        sys.exit(1)
    if not os.path.exists(chosen_whitelist):
        logger.error(f'Unable to find chosen_whitelist: {chosen_whitelist}')
        sys.exit(1)
    out_whitelist = f'{work_dir}/WHITELIST.txt'
    if not os.path.exists(out_whitelist):
        logger.info(f'Moving whitelist to work directory')
        if bool(re.search('.gz$', chosen_whitelist)):
            proc = subprocess.run(f'zcat {chosen_whitelist} > {out_whitelist}', shell=True, capture_output=True)
        else:
            proc = subprocess.run(f'cat {chosen_whitelist} > {out_whitelist}', shell=True, capture_output=True)
    return out_whitelist

## ~~~ ## 

def extract_10x_barcodes(
    ## Global params
    fastq_file: str | Path,
    kit: str = '3prime:v4',

    ## Vsearch params
    min_adapter_id: float = 0.7,
    threads: int | None = None,

    ## 10X kit params (needed for parasail)
    adaptor_suffix_len: int = 10,
    polyt_prefix_len: int = 10,
    barcode_len: int = 16,
    umi_len: int = 12,
    window_len: int = 100,

    ## Parasail params
    match_score: int = 5,
    mismatch_score: int = -1,
    vn_match_score: int = 1,
    tn_match_score: int = 1,
    gap_open_penalty: int = 2,
    gap_extend_penalty: int = 4,

    ## Post-parasail params
    max_adaptor_edit: int = 3,

    ## Rapidfuzz params
    rapidfuzz_error_correct: bool = False,
    max_barcode_fuzzdist: int = 2,
    mindelta_barcode_fuzzdist: int = 2,
    
    ## Directories
    whitelist_dir: str | Path | None = None,
    work_dir: str | Path | None = None,

    ## Clean up at the end?
    cleanup: bool = True
):
    r"""A single worker that assigns 10X barcodes from FASTQ files for each read. 

    Performs the functions of `adapter_scan_vsearch`, and `extract_barcode` of the original EPI2ME pipeline.
    
    Briefly: 
    1. Copy required information into work directory, 
    2. Select relevant metadata, e.g. adaptor sequences and whitelists
    3. Run `vsearch` to get reads with the adaptor sequences and trim as needed
    4. Use `parasail` to align an expected probe sequence to each read
    5. Check that barcode is in whitelist, and if not, attempt to find closest match with `rapidfuzz`

    Parameters
    ----------
    fastq_file
        String or Path indicating FASTQ file to parse.
    kit
        String indicating what kit has been used.
    min_adapter_id
        Float indicating minimum adaptor alignment for `vsearch`.
    threads
        Integer indicating number of threads to run `vsearch`.
    adaptor_suffix_len
        Integer indicating number of bases from end of the adaptor to include in `parasail` search.
    polyt_prefix_len
        Integer indicating number of bases from start of polyT sequence to include in `parasail` search.
    barcode_len
        Integer indicating number of bases of 10X barcode. Used in `parasail` search.
    umi_len
        Integer indicating number of bases of 10X UMI. Used in `parasail` search.
    window_len
        Integer indicating number of bases to scan for. Used in `parasail` search.
    match_score
        Integer indicating `parasail` score assigned to a match.
    mismatch_score
        Integer indicating `parasail` score assigned to a mismatch.
    vn_match_score
        Integer indicating `parasail` score assigned to a match between A/C/G and N.
    tn_match_score
        Integer indicating `parasail` score assigned to match between T and N.
    gap_open_penalty
        Integer indicating penalty score for opening new gaps (inferred by similar parameter name for `Clustal`). Increasing this makes gaps less frequent. Used by `parasail`.
    gap_extend_penalty
        Integer indicating penalty score for extending existing gaps by a base (inferred by similar parameter name for `Clustal`). Increasing this makes gaps shorter. Used by `parasail`.
    max_adaptor_edit
        Integer indicating max tolerable edit distance of adaptor sequence 1 in read after `parasail` alignment.
    rapidfuzz_error_correct
        Boolean indicating whether `rapidfuzz` error correction to be done. Defaults to False. 
        Recommended to only do this with a shortlisted whitelist, otherwise comparisons with >1 mio possible barcodes would be slow. 
    max_barcode_fuzzdist
        Integer indicating max tolerable edit distance between barcode and whitelist.
    mindelta_barcode_fuzzdist
        Integer indicating min difference in edit distances between best and second best whitelist match. Indicates how certain top match is true.
    whitelist_dir
        String or Path indicating where whitelist directory is. Defaults to `'~/WORK'`.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    ## Creating working directory and start logger
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []

    ## Prep params
    kit = kit.lower()

    ## Get relevant whitelist
    if whitelist_dir is None:
        out_whitelist = get_10x_whitelist(
            kit = kit,
            whitelist_dir = work_dir,
            work_dir = work_dir,
        )
    else:
        out_whitelist = get_10x_whitelist(
            kit = kit,
            whitelist_dir = whitelist_dir,
            work_dir = work_dir,
        )
    if not os.path.exists(f'{out_whitelist}'):
        logger.error(f'Unable to find whitelist file!')
        sys.exit(1)
    # files_to_delete.append(out_whitelist)

    ## Now get relevant adaptor
    ## Taken from sc_util
    kit_adaptors = {
        '3prime': {
            'adapter1': 'CTACACGACGCTCTTCCGATCT',
            'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
        },
        'multiome': {
            'adapter1': 'CTACACGACGCTCTTCCGATCT',
            'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
        },
        'visium': {
            'adapter1': 'CTACACGACGCTCTTCCGATCT',
            'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
        },
        '5prime': {
            'adapter1': 'CTACACGACGCTCTTCCGATCT',
            'adapter2': 'GTACTCTGCGTTGATACCACTGCTT'
        }
    }
    chosen_adaptor = [str(x) for x in kit_adaptors.keys() if str(x) in kit]
    if len(chosen_adaptor) != 1:
        logger.error(f'Cannot identify which adaptor to use. This should NEVER happen. Get in touch.')
        sys.exit(1)
    chosen_adaptor = chosen_adaptor[0]
    adaptor_seqs = kit_adaptors[chosen_adaptor]
    
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")

    ## Write a fasta with the relevant adaptors
    adaptor_fasta_info = {
        'adapter1_f': adaptor_seqs['adapter1'],
        'adapter1_r': adaptor_seqs['adapter1'][::-1].translate(complement_map),
        'adapter2_f': adaptor_seqs['adapter2'],
        'adapter2_r': adaptor_seqs['adapter2'][::-1].translate(complement_map),
    }
    out_adaptor = f'{work_dir}/ADAPTORS.fa'
    if not os.path.exists(out_adaptor):
        logger.info(f'Creating an ADAPTORS fasta')
        with open(out_adaptor, 'w') as fh:
            for name, seq in adaptor_fasta_info.items():
                fh.write(f">{name}\n{seq}\n")
    # files_to_delete.append(out_adaptor)

    ## Vsearch to find adaptors for each read
    ## Original uses .Popen instead of .run: changed to .run because need this to complete before moving on
    if threads is None:
        threads = threading.active_count()
    data_name = re.sub('.fastq$|.fq$|fastq.gz$|fq.gz$', '', Path(fastq_file).stem)
    if not os.path.exists(f'{work_dir}/{data_name}'):
        os.makedirs(f'{work_dir}/{data_name}', exist_ok=True)
    out_vsearch = f'{work_dir}/{data_name}/VSEARCH_ADAPTORS_{data_name}.tsv'
    if not os.path.exists(out_vsearch):
        logger.info(f'Finding adaptors with vsearch')
        ## Note the - after --usearch_global. This is likely for compatability with the pipe (don't delete).
        vsearch_cmd = f"seqkit fq2fa {fastq_file} | vsearch --usearch_global - \
            --db {out_adaptor} --minseqlength 20 --maxaccepts 5 --id {min_adapter_id} \
            --strand plus --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
            'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
            --userout {out_vsearch} --threads {threads}"
        logger.info(vsearch_cmd)
        proc = subprocess.Popen(vsearch_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.communicate()
    files_to_delete.append(out_vsearch)
    if os.path.getsize(out_vsearch)==0:
        logger.warning(f'Assuming no matches found...Skipping...')
        sys.exit()
    
    ## Now parse the Vsearch results (for explanation of column names, read the vsearch manual)
    ## Will cheat a little here and use the epi2me functions
    out_validreads = f'{work_dir}/{data_name}/VALID_READS_{data_name}.fastq'
    if not os.path.exists(out_validreads):
        logger.info(f'Parsing vsearch results')
        read_info = parse_vsearch(out_vsearch)
        logger.info(f'Creating fastq for {len(read_info)} reads')
        read_data = create_stranded_reads(fastq_file, read_info, chosen_adaptor, fl_only=False)
        with open(out_validreads, 'w') as f:
            for read in read_data:
                f.write(f'{read}')
    files_to_delete.append(out_validreads)
    
    ## Extract uncorrected barcodes
    ## Uses parasail, which needs a scoring matrix (mat) and a proposed search sequence (probe_seq)
    logger.info(f'Preparing for extracting raw barcodes')
    wl = pd.read_csv(out_whitelist, header=None).iloc[:,0].values
    mat = parasail.matrix_create('ACGTN', match_score, mismatch_score)
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        mat.pointer[0].matrix[i] = vn_match_score
    mat.pointer[0].matrix[22] = tn_match_score
    mat.pointer[0].matrix[27] = tn_match_score
    adp_seq_1 = adaptor_seqs['adapter1'][-adaptor_suffix_len:]
    ## Should not need to check again, so won't bother here
    if chosen_adaptor in ['5prime']:
        probe_seq = "{a1}{bc}{umi}{tso}".format(
            a1  = adp_seq_1,
            bc  = "N" * barcode_len,
            umi = "N" * umi_len,
            tso = "TTTCTTATATGGG"
        )
    else:
        probe_seq = "{a1}{bc}{umi}{pt}".format(
            a1  = adp_seq_1,
            bc  = "N" * barcode_len,
            umi = "N" * umi_len,
            pt  = "T" * polyt_prefix_len
        )

    ## Now the big step -- this will check to see if read can be assigned a barcode post-parasail
    logger.info(f'Evaluating individual reads and assigning barcodes')
    out_final = f'{work_dir}/{data_name}/TRIMMED_READS_{data_name}.fastq'
    out_tags  = f'{work_dir}/{data_name}/TRIMMED_READS_{data_name}_TAGS.tsv'
    with pysam.FastxFile(out_validreads, persist = False) as fastq_fh, open(out_final, 'w') as fq_out, open(out_tags, 'w') as tab_out:
        tab_out.write(
            '\t'.join([
                'read_id', 
                'assigned_barcode',
                'barcode_corrected',
                'original_barcode',
                'barcode_qscore_ascii',
                'barcode_qscore',
                'original_umi',
                'umi_qscore_ascii',
                'umi_qscore',
            ]) + '\n')
        for read in fastq_fh:
            logger.info(f'Processing {read.name}')
            if chosen_adaptor in ['5prime']:
                trim_right = False
                prefix_seq = read.sequence[:window_len]
                prefix_qv = read.quality[:window_len]
            else:
                trim_right = True
                prefix_seq = read.sequence[::-1].translate(complement_map)[:window_len]
                prefix_qv = read.quality[::-1].translate(complement_map)[:window_len]
    
            ## Align read
            p_alignment = parasail.sw_trace(s1=prefix_seq, s2=probe_seq, open=gap_open_penalty, extend=gap_extend_penalty, matrix=mat)
            ref_alignment = p_alignment.traceback.ref
            query_alignment = p_alignment.traceback.query
            start_of_barcode = ref_alignment.find('N')
    
            ## Parse results of parasail alignment
            adapter1_ed = len(adp_seq_1) # Levenshtein distance between adaptors
            barcode_no_ins = ""          # Barcode without insert
            umi_no_ins = ""              # UMI without insert
            bc_q_ascii = ""              # Barcode Q-score 
            umi_q_ascii = ""             # UMI Q-score
            if start_of_barcode > -1:
                adapter1_ed = editdistance.eval(query_alignment[0:start_of_barcode], adp_seq_1) 
                barcode_no_ins = query_alignment[start_of_barcode:(start_of_barcode + barcode_len)].replace('-', '') 
                umi_no_ins = query_alignment[(start_of_barcode + barcode_len):(start_of_barcode + barcode_len + umi_len)].replace('-', '')
                bc_start = prefix_seq.find(barcode_no_ins)
                bc_q_ascii = prefix_qv[bc_start:(bc_start + len(barcode_no_ins))]
                umi_start = prefix_seq.find(umi_no_ins)
                umi_q_ascii = prefix_qv[umi_start:(umi_start + len(umi_no_ins))]
    
            ## Decide whether read is worth keeping
            ## While original script has barcode correction take place AFTER mapping, figured, might as well do it here
            if (adapter1_ed <= max_adaptor_edit) & (len(barcode_no_ins) > 0) & (len(umi_no_ins) >0):
    
                ## Get numeric barcode score for determining confidence in barcode call
                numeric_barcode_score = min(list(map(lambda x: ord(x) - 33, bc_q_ascii)))
                numeric_umi_score = min(list(map(lambda x: ord(x) - 33, umi_q_ascii)))
    
                ## Change the ASCII Q-scores into something that can be saved in a .tsv
                umi_q_quoted = umi_q_ascii.replace('"', '""')
                barcode_q_quoted = bc_q_ascii.replace('"', '""')
                
                ## Check if barcode is in whitelist
                if rapidfuzz_error_correct:
                    assigned_barcode, barcode_corrected = fuzz_assign(
                        barcode = barcode_no_ins,
                        whitelist = wl,
                        max_barcode_fuzzdist = max_barcode_fuzzdist,
                        mindelta_barcode_fuzzdist = mindelta_barcode_fuzzdist,
                    )
                    if assigned_barcode == "":
                        continue
                else:
                    assigned_barcode = "NOT_ASSIGNED"
                    barcode_corrected = False
                        
                ## If passes barcode assignment, trim read
                ## WARNING: this overwrites the read information
                if trim_right:
                    trim_pos = p_alignment.end_query - 1
                    trim_pos -= polyt_prefix_len
                    new_sequence = read.sequence[:-trim_pos]
                    new_quality = read.quality[:-trim_pos]
                else:
                    trim_pos = p_alignment.end_query
                    new_sequence = read.sequence[trim_pos:]
                    new_quality = read.quality[trim_pos:]
                if new_sequence:
                    fq_out.write(f'@{read.name}\n{new_sequence}\n+\n{new_quality}\n')
                    tab_out.write('\t'.join([
                        f'"{read.name}"', 
                        f'"{assigned_barcode}"',
                        f'"{barcode_corrected}"',
                        f'"{barcode_no_ins}"',
                        f'"{barcode_q_quoted}"',
                        f'"{numeric_barcode_score}"',
                        f'"{umi_no_ins}"', 
                        f'"{umi_q_quoted}"',
                        f'"{numeric_umi_score}"',
                    ]) + '\n')

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
    return out_final


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################
## ----- Step 3: Read sorting ----- ##
######################################

def group_reads(
    ## Inputs
    trimmed_fastq_file: str | Path,
    trimmed_tags_tsv: str | Path | None = None,

    ## Groups to consider
    groups: dict = {
        'STICR': ['AGCAAACTGGGGCACAAGC', 'GAATTATCGAATT', 'GAATTGAGGAATT'],
        'CURIO': ['TCTTCAGCGTTCCCGAGA', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT']
    },
    default_group: str = 'RNA',
    all_matches_mode: bool = False,

    ## Vsearch params
    min_sequence_length: int = 8,
    min_adapter_id: float = 0.7,
    threads: int | None = None,

    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
):
    r"""Groups reads according whether an index is present.

    Meant to be run _after_ extracting 10X cellular barcodes. 
    Reads are grouped into separate FASTQs and .tsvs according to group names specified by `groups`.
    Reads will also be reoriented according to whether they match sequences in `groups`, or their reverse complement.

    If a read belongs to multiple groups or has ambiguous orientation, will be assigned to a new group "PROBLEM".
    If `all_matches_mode`, then "PROBLEM" will also include reads that match only some -- but not all -- sequences in a given group.
    
    Parameters
    ----------
    trimmed_fastq_file
        String or Path indicating the fastq file with 10X cellular barcodes and adaptor sequences trimmed off. 
        Should be a `TRIMMED_READS_{sample_name}.fq` file in `work_dir`.
    trimmed_tags_tsv
        String or Path indicating the tsv file with tag information for every read. 
        If not specified, infers from `trimmed_fastq_file`: i.e. matching `TRIMMED_READS_{sample_name}_TAGS.tsv` in `work_dir`.
    groups
        Dictionary of lists indicating what _forward_ sequence(s) are needed for a read to be assigned to any given group. 
        Name of the group as keys. Lists must be mutually exclusive. 
        "PROBLEM" cannot be a key in this dictionary -- this designation is used internally for problem reads.
    default_group
        String indicating name of group for all reads that do not get assigned to any of the groups in `groups`.
    all_matches_mode
        Boolean indicating whether -- if multiple sequences per group are provided -- to only assign it to a group if _all_ sequences in the list appear in the read.
        If False, then as long as _any_ of the sequences per group appears, then the read gets assigned to that group.
    min_sequence_length
        Integer of minimum length of reads to consider mapping against. Should be shorter than the shortest sequence in `groups`.
    min_adapter_id
        Float indicating minimum adaptor alignment for `vsearch`.
    threads
        Integer indicating number of threads to run `vsearch`. Defaults to `threading.active_count()`.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    ## Basic utilities
    logger, work_dir = setup_workdir(work_dir)
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")
    files_to_delete = []

    ## Check that file exists -- otherwise, don't bother
    ## We will delete this file at the end
    if not os.path.exists(trimmed_fastq_file):
        logger.error(f'Cannot find {trimmed_fastq_file}')
        sys.exit(1)
    files_to_delete.append(trimmed_fastq_file)

    ## Now make sure that fastq and matching tsv files exist
    if trimmed_tags_tsv is None:
        trimmed_tags_tsv = re.sub('.fq$|.fq.gz$|.fastq$|.fastq.gz$', '_TAGS.tsv', trimmed_fastq_file)
        if not os.path.exists(trimmed_tags_tsv):
            logger.error(f'Unable to find inferred .tsv file: {trimmed_tags_tsv}. Make sure that 10x cell barcodes have been already extracted!')
            sys.exit(1)
    files_to_delete.append(trimmed_tags_tsv)

    ## Okay, can now start with the algorithm
    ## First, check that the dictionary is correctly formated.
    all_seqs = []
    fasta_dict = {}
    for key, value in groups.items():
        if key == "PROBLEM":
            logger.error(f'groups cannot have "PROBLEM" as a key name -- this is used internally for problem reads')
            sys.exit(1)
        if not isinstance(value, list):
            if not isinstance(value, str):
                logger.error(f'groups needs to be a dictionary of lists')
                sys.exit(1)
            logger.warning(f'Got a string for {key}: converting to list of 1 element.')
            value = list(value)
            groups[key] = value
        all_seqs.extend(value)
        for validx in range(len(value)):
            val = value[validx]
            fasta_dict[f'{key}_{validx}_f'] = val
            fasta_dict[f'{key}_{validx}_r'] = val[::-1].translate(complement_map)
    if len(all_seqs) != len(set(all_seqs)):
        logger.error(f'Duplicated sequence(s) found: sequences must be unique to each group!')
        sys.exit(1)
        
    ## Generate a fasta of forwards and reverse sequences
    out_groupfasta = f'{work_dir}/GROUPS.fa'
    if not os.path.exists(out_groupfasta):
        logger.info(f'Making fasta of group sequences')
        with open(out_groupfasta, 'w') as fh:
            for name, seq in fasta_dict.items():
                fh.write(f">{name}\n{seq}\n")

    ## Files exist, so now VSEARCH: output to the same folder as the fastq
    out_dir = str(os.path.dirname(trimmed_fastq_file))
    data_name = re.sub(work_dir, '', out_dir)
    data_name = re.sub('[/]', '', data_name)
    out_vsearch_group = f'{out_dir}/VSEARCH_GROUPS_{data_name}.tsv'
    if threads is None:
        threads = threading.active_count()
    if not os.path.exists(out_vsearch_group):
        logger.info(f'Finding sequence hits for each group with vsearch')
        vsearch_cmd = f"seqkit fq2fa {trimmed_fastq_file} | vsearch --usearch_global - \
        --db {out_groupfasta} --minseqlength {min_sequence_length} --maxaccepts 5 --id {min_adapter_id} \
        --strand plus --wordlength 3 --minwordmatches 10 --userfields \
        'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
        --userout {out_vsearch_group} --threads {threads}"
        logger.info(vsearch_cmd)
        # proc = subprocess.run(vsearch_cmd, shell=True, capture_output=True)
        proc = subprocess.Popen(vsearch_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.communicate()
    files_to_delete.append(out_vsearch_group)

    ## Parse the vsearch results and update the .tsv
    logger.info(f'Sorting reads')
    cb_info = pl.read_csv(trimmed_tags_tsv, separator = '\t', has_header = True)
    vsearch_colnames = [
        "query", "target", "id", "alnlen", "mism", "opens",
        "qilo", "qihi", "qstrand", "tilo", "tihi",
        "ql", "tl"]
    schema = {
        'query': pl.Utf8,
        'target': pl.Categorical,
        'qilo': pl.UInt32,
        # 'qihi': pl.UInt32,
        # 'ql': pl.UInt32
    }
    columns_idxs = [vsearch_colnames.index(x) for x in schema.keys()]
    vs_info = pl.read_csv(
        source = out_vsearch_group,
        has_header = False,
        separator = '\t',
        schema_overrides = schema,
        columns = columns_idxs,
        new_columns = list(schema.keys())
    ).sort(["query", "qilo"])
    vs_info = vs_info.filter(pl.col('target')!='*')

    with pysam.FastxFile(trimmed_fastq_file, persist = False) as infile:
        for read in infile:
            ## Parse info
            logger.info(f'Processing {read.name}')
            read_info = cb_info.filter(pl.col('read_id') == read.name)
            read_info_str = '\t'.join([f'"{str(x)}"' for x in read_info.row(0)]) 
            sub = vs_info.filter(pl.col('query') == read.name)
    
            ## The default case: no hits, assume that it belongs to default_group
            if sub.shape[0]==0:
                group_membership = default_group
                group_outname = default_group
                out_fastq = f'{out_dir}/{group_membership}_{data_name}.fastq'
                out_tags = f'{out_dir}/{group_membership}_{data_name}_TAGS.tsv'
                if not os.path.exists(out_tags):
                    with open(out_tags, 'w') as tabfile:
                        tabfile.write('\t'.join(cb_info.columns) + '\tgroup' + '\n')
                with open(out_fastq, 'a') as fq_out:
                    fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
                with open(out_tags, 'a') as tag_out:
                    tag_out.write(read_info_str + f'\t"{group_outname}"\n')
                continue
    
            ## If there is a vsearch hit, now figure out what group it belongs to
            ## and the orientation of the read
            orient_find = [bool(re.search(f'_f$', x)) for x in sub['target'].to_list()]
            group_find = []
            for key in groups.keys():
                group_find.append([bool(re.search(f'^{key}', x)) for x in sub['target'].to_list()])
            group_find = np.array(group_find)
            orient_find = np.array(orient_find)
    
            ## Problem 1: Since sequences are mutually exclusive, read should never belong to multiple groups
            if sum(group_find.sum(1)>0) != 1:
                group_membership = 'PROBLEM'
                group_outname = 'MULTI_GROUPS'
                out_fastq = f'{out_dir}/{group_membership}_{data_name}.fastq'
                out_tags = f'{out_dir}/{group_membership}_{data_name}_TAGS.tsv'
                if not os.path.exists(out_tags):
                    with open(out_tags, 'w') as tabfile:
                        tabfile.write('\t'.join(cb_info.columns) + '\tgroup' + '\n')
                with open(out_fastq, 'a') as fq_out:
                    fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
                with open(out_tags, 'a') as tag_out:
                    tag_out.write(read_info_str + f'\t"{group_outname}"\n')
                continue
    
            ## Problem 2: Inconsistent read orientation
            if not (orient_find == orient_find[0]).all():
                group_membership = 'PROBLEM'
                group_outname = 'AMBI_STRAND'
                out_fastq = f'{out_dir}/{group_membership}_{data_name}.fastq'
                out_tags = f'{out_dir}/{group_membership}_{data_name}_TAGS.tsv'
                if not os.path.exists(out_tags):
                    with open(out_tags, 'w') as tabfile:
                        tabfile.write('\t'.join(cb_info.columns) + '\tgroup' + '\n')
                with open(out_fastq, 'a') as fq_out:
                    fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
                with open(out_tags, 'a') as tag_out:
                    tag_out.write(read_info_str + f'\t"{group_outname}"\n')
                continue
    
            ## Above problems dealt with, can now assign
            ## True orientation = forward, false = reverse
            assigned_group = np.array(list(groups))[group_find.sum(1)>0]
            assigned_group = assigned_group.item()
            orientation_is_forwards = (orient_find.sum()>0).item()
            
            ## Problem 3: Not all hits present -- only consider if all_matches_mode is true
            if all_matches_mode:
                if len(groups[assigned_group]) != len(sub['target']):
                    group_membership = 'PROBLEM'
                    group_outname = 'TOO_FEW_HITS'
                    out_fastq = f'{out_dir}/{group_membership}_{data_name}.fastq'
                    out_tags = f'{out_dir}/{group_membership}_{data_name}_TAGS.tsv'
                    if not os.path.exists(out_tags):
                        with open(out_tags, 'w') as tabfile:
                            tabfile.write('\t'.join(cb_info.columns) + '\tgroup' + '\n')
                    with open(out_fastq, 'a') as fq_out:
                        fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
                    with open(out_tags, 'a') as tag_out:
                        tag_out.write(read_info_str + f'\t"{group_outname}"\n')
                    continue
    
            ## Assume there are no more problems
            group_membership = assigned_group
            group_outname = assigned_group
            out_fastq = f'{out_dir}/{group_membership}_{data_name}.fastq'
            out_tags = f'{out_dir}/{group_membership}_{data_name}_TAGS.tsv'
            if not os.path.exists(out_tags):
                with open(out_tags, 'w') as tabfile:
                    tabfile.write('\t'.join(cb_info.columns) + '\tgroup' + '\n')
            with open(out_fastq, 'a') as fq_out:
                if orientation_is_forwards:
                    fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
                else:
                    ## If not forwards, then write reverse complement
                    fq_out.write(f'@{read.name}\n{read.sequence[::-1].translate(complement_map)}\n+\n{read.quality[::-1]}\n')
            with open(out_tags, 'a') as tag_out:
                tag_out.write(read_info_str + f'\t"{group_outname}"\n')
            
    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
        
    return out_dir


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###################################################################    
## ----- Step 4: STICR and CurioTrekker barcode extraction ----- ##
###################################################################

def process_curio(
    ## Inputs
    curio_fastq: str | Path,
    curio_whitelist: str | Path | list,

    ## Parasail params
    match_score: int = 5,
    mismatch_score: int = -1,
    vn_match_score: int = 1,
    tn_match_score: int = 1,
    gap_open_penalty: int = 2,
    gap_extend_penalty: int = 4,

    ## Post-parasail params
    max_linker_edit: int = 3,

    ## Rapidfuzz params
    max_barcode_fuzzdist: int = 2,
    mindelta_barcode_fuzzdist: int = 2,
    
    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
):
    r"""Extracts and assigns CurioTrekker barcodes.

    These are of the format (5' to 3'):

    ```
    [PCR HANDLE]                              [LINKER]                  [POLYA]
    ----------------------------------        ------------------        -------
    GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT - N8 - TCTTCAGCGTTCCCGAGA - N6 - N7VV30A
    ```
    
    We want the N8 and N6 barcodes, which when concatenated, returns the N14 barcodes provided in the whitelist.
    
    Parameters
    ----------
    curio_fastq
        String or Path indicating FASTQ file with CURIO reads.
    curio_whitelist
        String or Path (or list of these) indicating the location of `.txt` files with bead barcode coordinate information (whitelist).
        This is provided by the CurioTrekker portal (https://knowledgebase.curiobioscience.com/bioinformatics/tile-barcode/) 
        as a `.zip` upon keying in the relevant slide ID used. Comprises three columns: N14 barcode, x coordinate, y coordinate.
        If more than one provided, the union of all N14 barcodes are used as the whitelist.
    match_score
        Integer indicating `parasail` score assigned to a match.
    mismatch_score
        Integer indicating `parasail` score assigned to a mismatch.
    vn_match_score
        Integer indicating `parasail` score assigned to a match between A/C/G and N.
    tn_match_score
        Integer indicating `parasail` score assigned to match between T and N.
    gap_open_penalty
        Integer indicating penalty score for opening new gaps (inferred by similar parameter name for `Clustal`). Increasing this makes gaps less frequent. Used by `parasail`.
    gap_extend_penalty
        Integer indicating penalty score for extending existing gaps by a base (inferred by similar parameter name for `Clustal`). Increasing this makes gaps shorter. Used by `parasail`.
    max_linker_edit
        Integer indicating max tolerable edit distance of linker sequence in read after `parasail` alignment.
    max_barcode_fuzzdist
        Integer indicating max tolerable edit distance between barcode and whitelist.
    mindelta_barcode_fuzzdist
        Integer indicating min difference in edit distances between best and second best whitelist match. Indicates how certain top match is true.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []

    linker_seq = 'TCTTCAGCGTTCCCGAGA'
    probe_seq = "{pcr_handle}{bc_n8}{linker}{bc_n6}{polya}".format(
        pcr_handle  = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        bc_n8  = "N" * 8,
        linker = linker_seq,
        bc_n6 = "N" * 6,
        polya = ("N" * 7) + "VV" + ("A" * 30)
    )
    
    ## Start prepping for parasail alignment
    logger.info(f'Evaluating individual reads and assigning CURIO barcodes')
    if not isinstance(curio_whitelist, list):
        curio_whitelist = [curio_whitelist] #Assume a list of 1
    wl = []
    for wlfile in curio_whitelist:
        if not os.path.exists(wlfile):
            logger.error(f'Cannot find {wlfile}')
            sys.exit(1)
        wlx = pl.read_csv(
            wlfile, 
            has_header = False, 
            separator = '\t',
            new_columns = ['bc', 'x', 'y']
        )['bc'].to_list()
        wl.extend(wlx)

    ## Parasail scoring matrix
    mat = parasail.matrix_create('ACGTN', match_score, mismatch_score)
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        mat.pointer[0].matrix[i] = vn_match_score
    mat.pointer[0].matrix[22] = tn_match_score
    mat.pointer[0].matrix[27] = tn_match_score

    out_dir = str(os.path.dirname(curio_fastq))
    data_name = re.sub('^CURIO_', '', os.path.basename(curio_fastq).split('.')[0])
    out_tags = f'{out_dir}/CURIO_{data_name}_ASSIGNS.tsv'
    files_to_delete.append(curio_fastq)
    with pysam.FastxFile(curio_fastq, persist = False) as fastq_fh, open(out_tags, 'w') as tab_out:
        tab_out.write(
            '\t'.join([
                'read_id', 
                'assigned_barcode',
                'barcode_corrected',
                'original_barcode',
                'barcode_qscore_ascii',
                'barcode_qscore',
                # 'query_alignment' ## Not really necessary, let's keep the file small
            ]) + '\n')
        for read in fastq_fh:
            logger.info(f'Processing {read.name}')
            
            ## Align
            p_alignment = parasail.sw_trace(s1=read.sequence, s2=probe_seq, open=gap_open_penalty, extend=gap_extend_penalty, matrix=mat)
            ref_alignment = p_alignment.traceback.ref
            query_alignment = p_alignment.traceback.query
            start_of_barcode = ref_alignment.find('N')
            
            ## Parse results of parasail alignment
            linker_ed = len(linker_seq)
            barcode_no_ins = ""          # N8 + N6 barcode without insert
            bc_q_ascii = ""         # N8 + N6 Q-score 
            if start_of_barcode > -1:
                qc_start = read.sequence.index(query_alignment.replace('-', '')) 
                n8_no_ins = query_alignment[start_of_barcode:(start_of_barcode + 8)].replace('-', '')
                n8_q_ascii = read.quality[(qc_start + start_of_barcode):(qc_start + start_of_barcode + len(n8_no_ins))]
                linker_no_ins = query_alignment[(start_of_barcode + 8):(start_of_barcode + 8 + len(linker_seq))].replace('-', '')
                linker_ed = editdistance.distance(linker_no_ins, linker_seq)
                n6_no_ins = query_alignment[(start_of_barcode + 8 + len(linker_seq)):(start_of_barcode + 8 + len(linker_seq) + 6)].replace('-', '')
                n6_q_ascii = read.quality[(qc_start + start_of_barcode + len(n8_no_ins) + len(linker_no_ins)):(qc_start + start_of_barcode + len(n8_no_ins) + len(linker_no_ins) + len(n6_no_ins))]
                barcode_no_ins = f'{n8_no_ins}{n6_no_ins}'
                bc_q_ascii = f'{n8_q_ascii}{n6_q_ascii}'
    
            if (linker_ed <= max_linker_edit) & (len(barcode_no_ins) > 0):
                ## Get numeric barcode score for determining confidence in barcode call
                numeric_barcode_score = min(list(map(lambda x: ord(x) - 33, bc_q_ascii))) 
                
                ## Change the ASCII Q-scores into something that can be saved in a .tsv
                barcode_q_quoted = bc_q_ascii.replace('"', '""')
                    
                ## Check if barcode is in whitelist
                assigned_barcode, barcode_corrected = fuzz_assign(
                    barcode = barcode_no_ins,
                    whitelist = wl,
                    max_barcode_fuzzdist = max_barcode_fuzzdist,
                    mindelta_barcode_fuzzdist = mindelta_barcode_fuzzdist,
                )
                if assigned_barcode == "":
                    continue
    
                ## If no issues, then write to file:
                tab_out.write('\t'.join([
                    f'"{read.name}"', 
                    f'"{assigned_barcode}"',
                    f'"{barcode_corrected}"',
                    f'"{barcode_no_ins}"',
                    f'"{barcode_q_quoted}"',
                    f'"{numeric_barcode_score}"',
                    # f'"{query_alignment}"'
                ]) + '\n')

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
    return

## ~~~ ##

def process_sticr(
    ## Inputs
    sticr_fastq: str | Path,
    sticr_whitelist_dir: str | Path | None = None,
    indices_to_consider: str | list = ['1', '3', 'e'],
    sticr_whitelist: dict = {
        'Bit1': 'Bit1.F.with.Empty.fasta',
        'Bit2': 'Bit2.F.with.Empty.fasta',
        'Bit3': 'Bit3.F.with.Empty.fasta'
    },

    ## Parasail params
    match_score: int = 5,
    mismatch_score: int = -1,
    vn_match_score: int = 1,
    tn_match_score: int = 1,
    gap_open_penalty: int = 2,
    gap_extend_penalty: int = 4,

    ## Post-parasail params
    max_adaptor_edit: int = 3,

    ## Rapidfuzz params
    max_barcode_fuzzdist: int = 2,
    mindelta_barcode_fuzzdist: int = 2,
    
    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
):
    r"""Extracts and assigns STICR barcodes.
    
    These are of the format (5' to 3'):

    ```
    [ADAPTOR]             [INDEX 1]              [LK 1]         [LK 2]         [LK 3]   [bGH polyA]
    -------------------   --------------         ------         ------         ------   -----------
    AGCAAACTGGGGCACAAGC - GAATTATCGAATTC - N15 - GGATCC - N14 - GCTAGC - N14 - CTCGAG - ...........
    ```

    Barcode comprises of individual barcode segments or "bits", each with their own whitelist file. 

    If user is uncertain of index used, all indices ('1', '3', 'e') are mapped against as default. 
    The top performing match is then brought forwards for barcode assignment.
    
    Parameters
    ----------
    sticr_fastq
        String or Path indicating FASTQ file with STICR reads.
    sticr_whitelist_dir
        String or Path indicating directory with Bit1-3 files. If `None` (default), assumes absolute paths was provided in the `sticr_whitelist` dictionary.
    indices_to_consider
        String or list of strings indicating which indices to consider during mapping. Valid values are '1', '3', and/or 'e'. 
    sticr_whitelist
        Dictionary indicating filenames of Bit1, Bit2, and Bit3 whitelist files. 
    match_score
        Integer indicating `parasail` score assigned to a match.
    mismatch_score
        Integer indicating `parasail` score assigned to a mismatch.
    vn_match_score
        Integer indicating `parasail` score assigned to a match between A/C/G and N.
    tn_match_score
        Integer indicating `parasail` score assigned to match between T and N.
    gap_open_penalty
        Integer indicating penalty score for opening new gaps (inferred by similar parameter name for `Clustal`). Increasing this makes gaps less frequent. Used by `parasail`.
    gap_extend_penalty
        Integer indicating penalty score for extending existing gaps by a base (inferred by similar parameter name for `Clustal`). Increasing this makes gaps shorter. Used by `parasail`.
    max_adaptor_edit
        Integer indicating max tolerable edit distance of STICR adaptor sequence in read after `parasail` alignment.
    max_barcode_fuzzdist
        Integer indicating max tolerable edit distance between barcode and whitelist.
    mindelta_barcode_fuzzdist
        Integer indicating min difference in edit distances between best and second best whitelist match. Indicates how certain top match is true.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    ## Prelim checks and basic utility
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    if not isinstance(indices_to_consider, list):
        indices_to_consider = [indices_to_consider]
    indices_to_consider = [re.sub('index| ', '', str(x).lower()) for x in indices_to_consider]
    for x in indices_to_consider:
        if x not in ['3', '1', 'e']:
            logger.error(f'Invalid STICR index specified {x}: only accepts 1 or 3 or e')
            sys.exit(1)
    if sticr_whitelist_dir is None:
        sticr_whitelist_dir = ""
    if len(sticr_whitelist) != 3 or not np.array(list(sticr_whitelist.keys())) in np.array(['Bit1', 'Bit2', 'Bit3']):
        logger.error(f'Invalid sticr_whitelist provided: ensure this is a dictionary with entries for Bit1, Bit2, and Bit3 only')
        sys.exit(1)

    ## Parse thw whitelist
    wl_bit_dict = dict()
    for key, value in sticr_whitelist.items():
        full_wl_path = f'{sticr_whitelist_dir}/{value}'
        if not os.path.exists(full_wl_path):
            logger.error(f'Unable to find {full_wl_path}: please ensure absolute paths are provided!')
            sys.exit(1)
        fasta_seqs = []
        with open(full_wl_path, 'r') as wlfile:
            for line in wlfile:
                line = line.strip()
                if line.startswith('>'):
                    continue
                ## Annoyingly, the raw whitelists have a constant linker sequence upstream of the actual barcode -- trim
                if key == 'Bit1':
                    if len(line) > 15:
                        line = line[-15:]
                else:
                    if len(line) > 14:
                        line = line[-14:]
                fasta_seqs.append(line)
        wl_bit_dict[key] = fasta_seqs

    ## Build the probes for parasail
    adaptor_seq = 'AGCAAACTGGGGCACAAGC'
    index_1 = 'GAATTATCGAATTC'
    index_3 = 'GAATTGAGGAATTC'
    index_e = 'GAATTC'
    index_dict = {
        '1' : index_1,
        '3' : index_3,
        'e' : index_e
    }
    probe_dict = dict()
    for key, value in index_dict.items():
        if key not in indices_to_consider:
            continue
        probe_dict[key] = '{adaptor_seq}{index}{bc_1}{linker_1}{bc_2}{linker_2}{bc_3}{linker_3}'.format(
            adaptor_seq = adaptor_seq,
            index = value,
            bc_1 = "N" * 15,
            linker_1 = 'GGATCC',
            bc_2 = "N" * 14,
            linker_2 = 'GCTAGC',
            bc_3 = "N" * 14,
            linker_3 = 'CTCGAG'
        )

    ## Build the scoring matrix
    mat = parasail.matrix_create('ACGTN', match_score, mismatch_score)
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        mat.pointer[0].matrix[i] = vn_match_score
    mat.pointer[0].matrix[22] = tn_match_score
    mat.pointer[0].matrix[27] = tn_match_score

    ## Parasail will be a little complicated here: have to consider which index is best
    logger.info('Evaluating and extracting STICR barcodes')
    out_dir = str(os.path.dirname(sticr_fastq))
    data_name = re.sub('^STICR_', '', os.path.basename(sticr_fastq).split('.')[0])
    out_tags = f'{out_dir}/STICR_{data_name}_ASSIGNS.tsv'
    files_to_delete.append(sticr_fastq)
    with pysam.FastxFile(sticr_fastq, persist = True) as fastq_fh, open(out_tags, 'w') as tab_out:
        tab_out.write(
            '\t'.join([
                'read_id', 
                'assigned_index',
                'sticr_barcode',
                
                'assigned_bit1',
                'original_bit1',
                'bit1_corrected',
                'bit1_qscore_ascii',
                'bit1_qscore',

                'assigned_bit2',
                'original_bit2',
                'bit2_corrected',
                'bit2_qscore_ascii',
                'bit2_qscore',

                'assigned_bit3',
                'original_bit3',
                'bit3_corrected',
                'bit3_qscore_ascii',
                'bit3_qscore',
                
            ]) + '\n')
        for read in fastq_fh:
            logger.info(f'Processing {read.name}')
            
            ## Align
            p_aligns = dict()
            p_alignment = None
            p_align_best_score = 0
            p_best_index = ''
            for key, probe_seq in probe_dict.items():
                p_align_x = parasail.sw_trace(s1=read.sequence, s2=probe_seq, open=gap_open_penalty, extend=gap_extend_penalty, matrix=mat)
                p_aligns[key] = p_align_x # Not really necessary, but useful for downstream troubleshooting if needed
                if p_align_x.score > p_align_best_score:
                    p_align_best_score = p_align_x.score
                    p_alignment = p_align_x
                    p_best_index = key
            
            if p_best_index == "" or p_align_best_score == 0 or p_alignment is None:
                logger.warning(f'Problem with parasail aligning for {read.name}')
                continue
            
            ref_alignment = p_alignment.traceback.ref
            query_alignment = p_alignment.traceback.query
            start_of_barcode = ref_alignment.find('N')
            
            ## Parse results of parasail alignment
            adaptor_ed = len(adaptor_seq)
            bit1 = ""
            bit1_q_ascii = ""
            bit2 = ""
            bit2_q_ascii = ""
            bit3 = ""
            bit3_q_ascii = ""        
            if start_of_barcode > -1:
                ## Adaptor info
                aligned_adaptor = ref_alignment[0:len(adaptor_seq)].replace('-', '')
                adaptor_ed = editdistance.distance(aligned_adaptor, adaptor_seq)
                
                ## Bit 1 information
                bit1 = query_alignment[start_of_barcode:(start_of_barcode + 15)].replace('-', '')
                bit1_start = read.sequence.find(bit1)
                bit1_q_ascii = read.quality[(bit1_start):(bit1_start + len(bit1))]
            
                ## Bit 2 information
                ref_align_sub = ref_alignment[(start_of_barcode + 15):]
                quer_align_sub = query_alignment[(start_of_barcode + 15):]
                start_sub = ref_align_sub.find('N')
                if start_sub < 0:
                    continue #Unable to find bit 2
                bit2 = quer_align_sub[start_sub:(start_sub + 14)].replace('-', '')
                bit2_start = read.sequence[(bit1_start + len(bit1)):].find(bit2) + (bit1_start + len(bit1))
                bit2_q_ascii = read.quality[bit2_start:(bit2_start + len(bit2))]
            
                ## Bit 3 information
                ref_align_sub = ref_align_sub[(start_sub + 14):]
                quer_align_sub = quer_align_sub[(start_sub + 14):]
                start_sub = ref_align_sub.find('N')
                if start_sub < 0:
                    continue #Unable to find bit 3
                bit3 = quer_align_sub[start_sub:(start_sub + 14)].replace('-', '')
                bit3_start = read.sequence[(bit2_start + len(bit2)):].find(bit3) + (bit2_start + len(bit2))
                bit3_q_ascii = read.quality[bit3_start:(bit3_start + len(bit3))]
            
            if (adaptor_ed <= max_adaptor_edit) & (len(bit1) > 0) & (len(bit2) > 0) & (len(bit3) > 0):
                ## Get numeric barcode score for determining confidence in each bit's call
                bit1_q_score = min(list(map(lambda x: ord(x) - 33, bit1_q_ascii))) 
                bit2_q_score = min(list(map(lambda x: ord(x) - 33, bit2_q_ascii))) 
                bit3_q_score = min(list(map(lambda x: ord(x) - 33, bit3_q_ascii))) 
             
                ## Change the ASCII Q-scores into something that can be saved in a .tsv
                bit1_q_quoted = bit1_q_ascii.replace('"', '""')
                bit2_q_quoted = bit2_q_ascii.replace('"', '""')
                bit3_q_quoted = bit3_q_ascii.replace('"', '""')
            
                ## For each bit, check if they are in whitelist
                assigned_bit1, bit1_corrected = fuzz_assign(
                    barcode = bit1,
                    whitelist = wl_bit_dict['Bit1'],
                    max_barcode_fuzzdist = max_barcode_fuzzdist,
                    mindelta_barcode_fuzzdist = mindelta_barcode_fuzzdist,
                )
                if assigned_bit1 == "":
                    logger.info('No bit 1 found')
                    continue
            
                assigned_bit2, bit2_corrected = fuzz_assign(
                    barcode = bit2,
                    whitelist = wl_bit_dict['Bit2'],
                    max_barcode_fuzzdist = max_barcode_fuzzdist,
                    mindelta_barcode_fuzzdist = mindelta_barcode_fuzzdist,
                )
                if assigned_bit2 == "":
                    logger.info('No bit 2 found')
                    continue
            
                assigned_bit3, bit3_corrected = fuzz_assign(
                    barcode = bit3,
                    whitelist = wl_bit_dict['Bit3'],
                    max_barcode_fuzzdist = max_barcode_fuzzdist,
                    mindelta_barcode_fuzzdist = mindelta_barcode_fuzzdist,
                )
                if assigned_bit3 == "":
                    logger.info('No bit 3 found')
                    continue
            
                ## If no issues, then write to file
                final_barcode = f'{assigned_bit1}_{assigned_bit2}_{assigned_bit3}'
                tab_out.write('\t'.join([
                    f'"{read.name}"', 
                    f'"{p_best_index}"',
                    f'"{final_barcode}"',
                    
                    f'"{assigned_bit1}"',
                    f'"{bit1}"',
                    f'"{bit1_corrected}"',
                    f'"{bit1_q_quoted}"',
                    f'"{bit1_q_score}"',
            
                    f'"{assigned_bit2}"',
                    f'"{bit2}"',
                    f'"{bit2_corrected}"',
                    f'"{bit2_q_quoted}"',
                    f'"{bit2_q_score}"',
            
                    f'"{assigned_bit3}"',
                    f'"{bit3}"',
                    f'"{bit3_corrected}"',
                    f'"{bit3_q_quoted}"',
                    f'"{bit3_q_score}"',
                ]) + '\n')

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
        
    return


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################
## ----- Step 5: Mapping RNA reads ----- ##
###########################################

def map_rna(
    ## Inputs
    rna_fastq: str | Path,
    rna_tags: str | Path | None = None,
    
    ## Minimap2 params
    threads: int | None = None,
    
    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
):
    r"""Maps FASTQ reads to a genome using minimap2.

    Assumes that get_genome() has already been run and that GENOME files already exist in the work directory.

    Outputs `RNA_ALN_{sample_name}.bam` files with alignment information. These will be merged and deleted in subsequent steps (e.g. `feature_assign`).
    
    Parameters
    ----------
    rna_fastq
        String or Path indicating FASTQ file with RNA reads.
    rna_fastq
        String or Path indicating `_TAGS.tsv` file with matching tag information for reads in `rna_fastq`. Defaults to matching `RNA_{sample name}_TAGS.tsv` file in the same direcotry.
    threads
        Integer indicating number of threads to run `minimap2`. Defaults to `threading.active_count()`.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    ## Prelim checks and basic utility
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    if not os.path.exists(rna_fastq):
        logger.error(f'{rna_fastq} not found')
        sys.exit(1)
    files_to_delete.append(rna_fastq)
    out_dir = str(os.path.dirname(rna_fastq))
    data_name = re.sub('^RNA_', '', os.path.basename(rna_fastq).split('.')[0])
    files_to_delete = []
    genome_file_checks = ['GENOME.fa', 'GENOME.fa.fai', 'GENOME.gtf', 'GENOME.mmi']
    for file_check in genome_file_checks:
        if not os.path.exists(f'{work_dir}/{file_check}'):
            logger.error(f'{work_dir}/{file_check} not found: ensure get_genome() has already been run!')
            sys.exit(1)
    if rna_tags is None:
        rna_tags = f'{out_dir}/RNA_{data_name}_TAGS.tsv'
    if not os.path.exists(rna_tags):
        logger.error(f"Unable to find TAGS file: {rna_tags}")
        sys.exit(1)
    tag_df = pl.read_csv(rna_tags, has_header=True, separator='\t')

    ## Check if final output exists: tagged and aligned .bam file
    out_tagged_bam = f'{out_dir}/RNA_ALN_{data_name}.bam'
    if os.path.exists(out_tagged_bam):
        return

    ## The minimap2 command
    ref_mmi = f'{work_dir}/GENOME.mmi'
    ref_gtfbed = f'{work_dir}/GENOME.bed'
    out_bam = f'{out_dir}/TEMPORARY_{data_name}.bam' #Highly recommend to delete this after
    if threads is None:
        threads = threading.active_count()
    if not os.path.exists(out_bam):
        logger.info(f'Mapping with minimap2')
        ## Map with minimap2 --> convert sam to bam with samtools --> 
        ## Exclude reads with secondary alignments (flag 256) -->
        ## Sort and index --> Save the sorted .bam and its .csi index (akin to .bai)
        minimap_cmd = f'minimap2 -ax splice -uf --MD \
            -t {threads} -K 10M \
            --junc-bed {ref_gtfbed} \
            --cap-kalloc 100m \
            {ref_mmi} {rna_fastq} | \
            samtools view -uh --no-PG -F 256 - | \
            samtools sort -o {out_bam} -' #Won't write index yet
        logger.info(minimap_cmd)
        proc = subprocess.run(minimap_cmd, shell = True, capture_output=True)
    if not os.path.exists(out_bam):
        logger.error(f'Something has gone wrong with minimap2 alignment')
        sys.exit(1)
    files_to_delete.append(out_bam)

    ## When this is done, can now annotate with previously extracted barcode and UMI information
    ## Note we have not yet assigned UMIs -- these will just be raw UMIs
    with pysam.AlignmentFile(out_bam, 'rb') as infile:
        with pysam.AlignmentFile(out_tagged_bam, 'wb', header = infile.header) as outfile:
            for read in infile:
                tag_info_for_read = tag_df.filter(pl.col('read_id')==read.query_name).unique()
                if tag_info_for_read.shape[0] != 1:
                    logger.warning(f'Unable to tag read {read.query_name}')
                    continue
                read.set_tag("CB", tag_info_for_read['assigned_barcode'].item(), "Z")
                read.set_tag("CR", tag_info_for_read['original_barcode'].item(), "Z")
                read.set_tag("CY", tag_info_for_read['barcode_qscore_ascii'].item(), "Z")
                read.set_tag("UR", tag_info_for_read['original_umi'].item(), "Z")
                read.set_tag("UY", tag_info_for_read['umi_qscore_ascii'].item(), "Z")
                outfile.write(read)
    
    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
    return


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################
## ----- Step 6: Feature assignment ----- ##
############################################

def fa_stringtie(
    ## Inputs from wrapper
    bam_file: str | Path,
    work_dir: str | Path | None = None,

    ## kwargs user can define
    stringtie_opts: str = "-c 2",
    threads: int | None = None,
    min_mapping_quality: int = 30, 
    min_transcript_coverage: float = 0.4,
    min_read_coverage: float = 0.4,
    use_gffcompare: bool = True,
):
    r"""String-tie based approach for feature assignment.

    Splits input BAM file into individual chromosome subsets before aligning to `stringtie`-constructed transcript sequences.
    
    If `use_gffcompare = True`, then the `stringtie` transcript sequences are first checked according to their `class_code`
    (see https://ccb.jhu.edu/software/stringtie/gffcompare.shtml for details on what class codes mean).
    
    `minimap2` with the `map-ont` setting is then used to remap reads to each `stringtie` transcript.
    Best performing mapping (based on transcript and read coverage) is returned as the annotation for that read.

    Parameters
    ----------
    bam_file
        String or Path indicating location of BAM file to process.
    work_dir
        String or Path indicating where working directory is.  Defaults to `'~/WORK'`.
    stringtie_opts
        String indicating `stringtie` parameters to use. Defaults to `-c 2`, meaning a minimum coverage of 2 reads.
    threads
        Integer indicating number of CPUs used for other software: `stringtie`, `minimap2`, and `samtools`. Defaults to `threading.active_count()`.
    min_mapping_quality
        Integer indicating minimum mapping quality score for alignment to be considered for downstream assignment.
    min_transcript_coverage
        Float indicating minimum fraction (out of 1) of transcript sequence aligned after `minimap2`.
    min_read_coverage
        Float indicating minimum fraction (out of 1) of read sequence aligned after `minimap2`.
    use_gffcompare
        Boolean indicating if `gffcompare` should be used to QC `stringtie` sequences. 
        Defaults to True, but can be turned off if `gffcompare` is not installed on the machine. 
    """
    ## Generic checks -- should be used for every fa_ function.
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    if not os.path.exists(bam_file):
        logger.error(f'Cannot find {bam_file}')
        sys.exit(1)
    data_name = re.sub('.bam$', '', str(os.path.basename(bam_file)))
    out_dir = str(os.path.dirname(bam_file))

    ## Function-specific processing goes here -----------------

    ## Step 0: Get all chr.gtf splits (previously done in get_genome)
    chr_gtf_files = glob.glob(f'{work_dir}/*chr.gtf')
    genome_fa = f'{work_dir}/GENOME.fa'
    
    ## Step 1: Sort and write index for bam_file
    tmp_bam = f'{out_dir}/TEMP_{data_name}.bam'
    cmd = f'samtools sort --write-index -o {tmp_bam} {bam_file}'
    proc = subprocess.run(cmd, shell=True, capture_output=True)
    files_to_delete.append(tmp_bam)
    files_to_delete.append(f'{tmp_bam}.csi')

    ## Loop through each chromosome in merged bam and process
    ## Note I have hard coded the default stringtie_opts into the cmd itself
    if threads is None:
        threads = threading.active_count()

    transcript_assigns_out = f'{out_dir}/TRANSCRIPT_ASSIGNS_{data_name}.tsv'
    with open(transcript_assigns_out, 'w') as final_out_file:
        final_out_file.write(
            '\t'.join([
                'read_name',
                'ref_name',
                'map_score',
                'aln_score',
                'transcript_coverage',
                'read_coverage',
                'gene_id',
                'gene_name',
                'chr',
            ]) + '\n')
                            
    for i in range(len(chr_gtf_files)):

        chr_gtf_file = chr_gtf_files[i]
        chrx = re.sub('chr.gtf$', '', os.path.basename(chr_gtf_file))
        logger.info(f'Processing {chr_gtf_file} (chr: {chrx})')
        
        ## Step 2: Build a gff from stringtie. This contains information about exons etc.
        logger.info(f'Building stringtie reference')
        stringtie_out = f'{out_dir}/{data_name}_chr{chrx}.gff'
        cmd = f'samtools view -h {tmp_bam} {chrx} \
        | stringtie -L {stringtie_opts} -p {threads} -G {chr_gtf_file} \
        -l "{chrx}.stringtie" -o {stringtie_out} -'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        files_to_delete.append(stringtie_out)
        if os.stat(stringtie_out).st_size==0: ## No results
            continue

        ## Optional step (2.5): QC check of stringtie outputs
        ## Drop transcripts that do not have successful class_code values
        if use_gffcompare:
            logger.info(f'Running gffcompare to QC stringtie sequences')
            ## First, rename the existing stringtie to another file name
            raw_stringtie_out = f'{out_dir}/{data_name}_chr{chrx}_RAW.gff'
            if os.path.exists(stringtie_out) and not os.path.exists(raw_stringtie_out):
                cmd = f'mv {stringtie_out} {raw_stringtie_out}'
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
            files_to_delete.append(raw_stringtie_out)
            gff_prefix = f'{out_dir}/GFFCOMPARE_{data_name}_chr{chrx}'
            cmd = f'gffcompare -o {gff_prefix} -r {chr_gtf_file} {raw_stringtie_out}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
            files_to_delete.extend(glob.glob(f'{gff_prefix}*'))
            tmap_file = [x for x in glob.glob(f'{gff_prefix}*') if '.gff.tmap' in x][0]
            if not os.path.exists(tmap_file):
                continue
            if os.stat(tmap_file)==0:
                continue
            gff_info = pl.read_csv(tmap_file, separator='\t')
            n_lines_original = gff_info.shape[0]
            gff_info = gff_info.filter(~pl.col('class_code').is_in(['i', 'y', 'p', 's']))
            logger.info(f'After class code filtering, gffcompare results have changed from {n_lines_original} to {gff_info.shape[0]} rows')
            with open(stringtie_out, 'w') as stringtie_out_file, open(raw_stringtie_out, 'r') as stringtie_in_file:
                for line in stringtie_in_file:
                    if line.startswith('#'):
                        stringtie_out_file.write(line)
                        continue
                    field = 'transcript_id'
                    field_info = re.findall(rf'{field} "([^"]*)"', line)
                    if len(field_info)!=1:
                        continue
                    if field_info[0] in gff_info['qry_id']:
                        stringtie_out_file.write(line)
            if os.stat(stringtie_out).st_size==0: ## No results
                continue
            
        ## Step 3: Get transcriptome sequences from stringtie .gff
        logger.info(f'Getting transcriptome sequences')
        transcriptome_out = f'{out_dir}/{data_name}_chr{chrx}_transcriptome.fa'
        cmd = f'gffread -g {genome_fa} -w "{transcriptome_out}" "{stringtie_out}"'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        files_to_delete.append(transcriptome_out)
        if os.stat(transcriptome_out).st_size==0: ## No results
            continue
            
        ## Step 4: Minimap2 re-mapping
        ## Choice of minimap2 params come from https://github.com/epi2me-labs/wf-single-cell/blob/master/subworkflows/assign_features.nf
        logger.info(f'Minimap2 mapping of reads to transcriptome')
        new_out_bam = f'{out_dir}/{data_name}_chr{chrx}_tralign.bam'
        cmd = f'samtools view -h {tmp_bam} {chrx} \
        | samtools fastq | minimap2 -ax map-ont \
        --cap-kalloc 100m --cap-sw-mem 50m \
        --end-bonus 10 -p 0.9 -N 3 -t {threads} \
        {transcriptome_out} - | samtools view -h -@ {threads} -b -F 2052 - \
        | samtools sort -n -@ {threads} --no-PG - > {new_out_bam}'
        proc = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
        proc.communicate()
        files_to_delete.append(new_out_bam)
        if os.stat(new_out_bam).st_size==0: ## No results
            continue

        ## Step 5A: Parse the GFF file StringTie outputs -- will use this for cross-referencing to gene
        ## See https://www.ensembl.org/info/website/upload/gff.html for column names
        tr_info = pl.read_csv(stringtie_out, has_header = False, comment_prefix='#', separator='\t',
                              new_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']).unique()
        fields_of_interest = ['transcript_id', 'reference_id', 'ref_gene_id']
        for field in fields_of_interest:
            ## See https://gist.github.com/rf-santos/22f521c62ca2f85ac9582bf0d91e4054: but modify for polars
            tr_info = tr_info.with_columns(
                (tr_info['attribute'].map_elements(lambda x: re.findall(rf'{field} "([^"]*)"', x)[0] if rf'{field} "' in x else '')).alias(field)
            )
        ## Note: the following line will drop any transcript_id that is not matched to a reference_id (gene)
        ## Will assume that we only care about reads that were successfully mapped to a gene
        ## Alternatively, can create a pseudo-gene using coordinate info (but how to interpret?)
        tr_info = tr_info.select(fields_of_interest).unique().sort('transcript_id').filter(pl.col('reference_id')!="")

        ## Step 5B: Parse the new BAM file -- we want to know what gene got assigned to what transcript
        records = {}
        with pysam.AlignmentFile(new_out_bam, 'rb', check_sq = False) as bamfile: #Not sure if check_sq is needed, but kept here to better match EPI2ME
            ## Get reference lengths for all reference sequences
            ref_lengths = dict(zip(bamfile.header.references, bamfile.header.lengths))
            for align in bamfile.fetch(until_eof=True):
                if align.query_name not in records.keys():
                    records[align.query_name] = []
                try:
                    aln_score = align.get_tag('AS')
                except KeyError:
                    aln_score = 0
                if align.reference_name not in tr_info['transcript_id']:
                    continue
                if align.mapping_quality < min_mapping_quality:
                    continue
                tr_info_sub = tr_info.filter(pl.col('transcript_id')==align.reference_name)
                reference_id = tr_info_sub.unique()['reference_id'][0]
                ref_gene_id = tr_info_sub.unique()['ref_gene_id'][0]
                transcript_coverage = float(align.reference_length / ref_lengths[align.reference_name])
                read_coverage = float(align.query_alignment_length / align.infer_read_length())
                if transcript_coverage < min_transcript_coverage:
                    continue
                if read_coverage < min_read_coverage:
                    continue
                info_dict = {
                    'read_name': align.query_name,
                    'ref_name': align.reference_name,
                    'map_score': align.mapping_quality,
                    'aln_score': aln_score,
                    'transcript_coverage': transcript_coverage,
                    'read_coverage': read_coverage,
                    'gene_id': reference_id,
                    'gene_name': ref_gene_id,
                    'chr': chrx,
                }
                records[align.query_name].append(info_dict)
    
        ## Step 6: Disambiguate assignments
        ## If there is ambiguity (multiple transcripts), then pick the best coverage
        final_records = []
        for key, value in records.items():
            if len(value)==0:
                continue
            if len(value)==1:
                final_records.append(value[0])
                continue
                
            ## Ambiguity: pick the best assignment
            read_assign_info = pl.DataFrame(value).sort(['read_coverage', 'transcript_coverage'], descending=True)
            value = read_assign_info[0].to_dicts()
            final_records.append(value[0])

        ## Step 6: tabulate the identity of each read GN and TR
        ## This will be used for cross-referencing
        with open(transcript_assigns_out, 'a') as final_out_file:
            pl.DataFrame(final_records).write_csv(final_out_file, include_header=False, separator='\t')

        ## Optional step (7): If want, can merge and save these intermediate outputs
        ## I will just save the transcriptome gff files -- easy enough to get back the .fa from them if needed
        ## There are commented lines that will split every .gff
        merged_gff = f'{out_dir}/{data_name}_alltranscripts_stringtie.gff'
        cmd = f'cat {stringtie_out} >> {merged_gff}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)

    ## Function-specific processing ends here -----------------
    
    ## All of these functions need to return:
    ## 1. Files to delete (should cleanup be true)
    ## 2. Name of file that assigns an identity to every read
    return files_to_delete, transcript_assigns_out

## ~~~ ##

def feature_assign(
    ## Inputs
    bam_dir: str | Path | list,
    merge_bams: bool = True,
    merged_bam_name: str | None = None,
    method: str = 'stringtie',
    
    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
    **kwargs
):
    r"""Wrapper function for different algorithms for feature assignment.

    Briefly:
    1. Merge BAMs by sample. This is the MERGED.bam. 
    2. Call genes and transcripts per sample. Algorithm loops over each chromosome (for better memory conservation). 
    3. Return a MERGED_annotated.bam. Delete the unnecessary MERGED.bam if `cleanup = True`.

    The MERGED_annotated.bam will then be used for UMI deduplication in the next step.

    Currently only supports the "stringtie" protocol from EPI2ME.

    Parameters
    ----------
    bam_dir
        String or Path (or list of these) indicating location(s) of `RNA_ALN_*.bam` files. If multiple files found (e.g. in subdirectories),
        will merge bams together if `merge_bam = True`; otherwise, will run the feature assignment file by file.
    merge_bam
        Boolean indicating if bams are to be merged into one big merged `.bam` before feature calling.
    merged_bam_name
        String indicating name of the merged bam. If not specified, defaults to "MERGED_*.bam", where "*" is the hex of the current time. 
    method
        String indicating method for calling features per read. Accepts "stringtie". 
        TO DO: allow implementation of other approaches in future.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    **kwargs
        Additional keyword arguments passed to feature assignment method of choice.
    """
    ## Basic utility and checks
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    supported_methods = ['stringtie']
    if method.lower() not in supported_methods:
        logger.error(f'{method} is not supported')
        sys.exit(1)
    if not isinstance(bam_dir, list):
        bam_dir = [bam_dir]
    if merged_bam_name is None:
        merged_bam_name = f"MERGED_{hex(int(time.time()))}"

    ## Pop out kwargs users are never allowed to change
    kwargs.pop('bam_file', None)
    kwargs.pop('work_dir', None) #This actually can't happen, but put here to be safe

    ## For each bam_dir, find all the needed bam files
    bamf_list = []
    for bam_dirx in bam_dir:
        found_files = [str(x) for x in list(Path(bam_dirx).rglob(f'RNA_ALN_*.bam'))]
        bamf_list.extend(found_files)
        files_to_delete.extend(found_files)
    if len(bamf_list)==0:
        logger.warning(f'No RNA_ALN_*.bam files found in the specified bam_dir: assuming no successful mapping...stopping here...')
        return

    ## Specify one or more bam(s) to process
    bams_to_process = bamf_list.copy()
    if merge_bams:
        logger.info(f'Merging {len(bamf_list)} bams into {merged_bam_name}.bam in {work_dir}')
        merged_bam_out = f'{work_dir}/{merged_bam_name}.bam'

        ## Copy the first bam and call this the merged bam
        first_bam = bamf_list.pop(0)
        proc = subprocess.run(f'cp {first_bam} {merged_bam_out}', shell=True, capture_output=True) 

        ## Now that the merged file exists, merge with every remaining bam file with samtools
        ## If bamf_list is only 1 file big, end result should just be copying of the original bam into work_dir
        for bam_file in bamf_list:
            samtools_cmd = f"samtools merge -o {merged_bam_out} {merged_bam_out} {bam_file}"
            logger.info(samtools_cmd)
            proc = subprocess.run(samtools_cmd, shell = True, capture_output=True)

        ## Change bams_to_process to just contain the merged_bam_out
        bams_to_process = [merged_bam_out]

    ## Now the actual feature assignment
    final_bam_outs = []
    for bam_file in bams_to_process:
        data_name = re.sub('.bam$', '', str(os.path.basename(bam_file)))
        out_dir = str(os.path.dirname(bam_file))
        final_bam_out = f'{out_dir}/{data_name}_annotated.bam'

        if os.path.exists(final_bam_out):
            logger.warning(f'{final_bam_out} already exists...Will skip...')
            continue

        ## The feature assignment method chosen
        if method=='stringtie':
            files_to_delete_x, transcript_assigns_out = fa_stringtie(
                bam_file = bam_file,
                work_dir = work_dir,
                **kwargs,
            )
        ## Other methods go here

        ## Consolidate files to delete
        files_to_delete.extend(files_to_delete_x)

        ## Read the transcript_assigns_out file
        ## It is very important that the column names match up
        tr_assign_info = pl.read_csv(transcript_assigns_out, separator='\t')

        ## Now that we have the annotations, we will update the bam names
        with pysam.AlignmentFile(bam_file, 'rb') as infile:
            with pysam.AlignmentFile(final_bam_out, 'wb', header = infile.header) as outfile:
                for read in infile:
                    tag_info_for_read = tr_assign_info.filter(pl.col('read_name')==read.query_name).unique()
                    if tag_info_for_read.shape[0] != 1:
                        logger.warning(f'Unable to find transcript information for read {read.query_name}')
                        continue
                    read.set_tag("GI", tag_info_for_read['gene_id'].item(), "Z") #New field, gene ID
                    read.set_tag("GN", tag_info_for_read['gene_name'].item(), "Z")
                    read.set_tag("TR", tag_info_for_read['ref_name'].item(), "Z")
                    outfile.write(read)

        ## Once we have annotated file, don't need the original bam either
        ## We also no longer need the transcript assigns: but will keep them for QC for now
        files_to_delete.append(bam_file)
        files_to_delete.append(transcript_assigns_out)
        final_bam_outs.append(final_bam_out)
    

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
        
    return final_bam_outs


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################
## ----- Step 7: UMI deduplication ----- ##
###########################################

def get_adj_list_directional_lev(self, umis, counts, threshold=2):
    r"""Use Levenshtein distance for UMI clustering instead of Hamming.

    This function is to monkey-patch UMIClusterer._get_adj_list_directional.

    Copied from EPI2ME wf-single-cell.
    """
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if editdistance.eval(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)
    return adj_list


def umi_clusterer_call(self, umis, threshold):
    r"""To replace UMIClusterer.__call__.

    Use this method to monkey-patch the UMICluterer.__call__ in order to remove the
    nessesity for all UMIs to be the same length, allowing for deletions in the UMIs.

    Copied from EPI2ME wf-single-cell.
    """
    counts = umis
    umis = list(umis.keys())
    self.positions += 1
    number_of_umis = len(umis)
    self.total_umis_per_position += number_of_umis
    if number_of_umis > self.max_umis_per_position:
        self.max_umis_per_position = number_of_umis
    adj_list = self.get_adj_list(umis, counts, threshold)
    clusters = self.get_connected_components(umis, adj_list, counts)
    final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]
    return final_umis

# Monkey-patch the umi-tools clusterer with a modified method
UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
UMIClusterer.__call__ = umi_clusterer_call
## Let's also keep track of whether the clusterer is updated
UMIClusterer.is_monkey_patched = True

## ~~~ ##

def cluster_umi(
    umis: list | dict,
    umi_clusterer: type | None = None,
    cluster_threshold: float = 2.0,
    return_corrected_umis: bool = False,
    logger: logging.Logger | None = None
):
    r"""Function for clustering a list of UMIs.

    Returns a list of UMI groups. But if `return_corrected_umis = True` and the original input was a list of UMIs (instead of a dictionary of frequencies),
    then a corrected list of UMIs will be returned instead.
    
    Parameters
    ----------
    umis
        A list of UMIs or a dictionary of UMIs and their frequency. If a list provided, will first tabulate UMI frequencies in list with `collections.Counter()`.
    umi_clusterer
        A umi-tools "functor" that will determine distance between UMIs. Defaults to the monkey-patched UMIClusterer.
    cluster_threshold
        Float indicating edit distance for clustering.
    return_corrected_umis
        Boolean indicating if want to instead return corrected list of UMIs. 
    logger
        A logging.Logger to keep track of warning messages. If not provided, messages are returned in console instead. 
    """
    ## Basic value checks
    umi_correction_possible = False
    umi_list = []
    if not isinstance(umis, dict):
        if not isinstance(umis, list):
            if logger:
                logger.error(f"umis is neither a dict nor a list")
            else:
                raise TypeError(f"umis is neither a dict nor a list")
            sys.exit(1)
        umi_list = umis.copy()
        umis = collections.Counter(umis)
        umi_correction_possible = True
    if return_corrected_umis:
        if not umi_correction_possible:
            if logger:
                logger.error(f"Cannot correct UMIs if umis is not a list because no positional information")
            else:
                raise TypeError(f"Cannot correct UMIs if umis is not a list because no positional information")
            sys.exit(1)

    ## Specify the clusterer
    if umi_clusterer is None:
        umi_clusterer = UMIClusterer(cluster_method = "directional")
    try:
        is_monkey_patched = umi_clusterer.is_monkey_patched
    except AttributeError:
        if logger:
            logger.warning(f"Specified umi_clusterer does not have a is_monkey_patched attribute: assuming this is NOT the default Levenshtein clusterer")
        else:
            warnings.warn(f"Specified umi_clusterer does not have a is_monkey_patched attribute: assuming this is NOT the default Levenshtein clusterer")

    ## Actual clustering        
    clustered_umis = umi_clusterer(umis, threshold = cluster_threshold)
    if not return_corrected_umis:
        return clustered_umis

    ## If return_corrected_umis
    umi_df = pl.DataFrame({'umi_list': umi_list})
    for umi_group in clustered_umis:
        umi_df = umi_df.with_columns(
            pl.when(pl.col('umi_list').is_in(umi_group)).then(pl.lit(umi_group[0])).otherwise('umi_list').alias('umi_list')
        )
    return umi_df['umi_list'].to_list()

## ~~~ ##    

def deduplicate_umi(
    ## Inputs
    annotated_bam: str | Path | list,
    keep_all_reads: bool = True,

    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
):
    
    r"""Deduplicate UMIs.

    Briefly:
    1. For every set of RNA reads from the same cell AND mapping to the same transcript, pull out their raw UMI sequences.
    2. How far their UMI sequences are in edit distance (Levenshtein), and cluster those that are very close.
    3. Remove duplicate reads (pick the highest alignment score). 
    Won't actually delete them, but will tag them with a "DP" tag for "duplicated".

    Note that de-duplication will be only done for RNA reads here; STICR and CURIO reads are handled separately.
    
    Parameters
    ----------
    annotated_bam
        String, Path, or list of these indicating the absolute path(s) of annotated bam file(s) returned by the previous `feature_assign` step. 
        If a list provided, will loop through each bam file and output to said file's directory.
    keep_all_reads
        Boolean indicating if all reads should be kept. 
        If True, then reads to be filtered can be identified by the DP:1 tag.
        If False, filtered reads are already removed and only DP:0 tags exist.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    """
    ## Prelim checks and basic utility
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    if not isinstance(annotated_bam, list):
        annotated_bam = [annotated_bam]
    files_to_delete.extend(annotated_bam)

    ## Loop through bam files
    final_bam_outs = []
    for abam in annotated_bam:
        logger.info(f'Processing {abam}')
        data_name = re.sub('_annotated.bam$', '', str(os.path.basename(abam)))
        out_dir = str(os.path.dirname(abam))
        final_bam_out = f'{out_dir}/{data_name}_annotated_deduplicated.bam'
        final_geneexpr_out = f'{out_dir}/{data_name}_geneexpr.csv'

        if os.path.exists(final_bam_out) and os.path.exists(final_geneexpr_out):
            logger.warning(f'Both {final_bam_out} and {final_geneexpr_out} already exists...Will skip...')
            continue

        ## Parse annotated bam file
        awk_string = '{for (i=12; i<=NF; ++i) {if (substr($i, 1, 2) == TAG) {print substr($i, 6); next}}}'
        desired_tags = ['CB', 'TR', 'UR', 'UY']
        align_info = {}
        for desired_tag in desired_tags:
            cmd = f"samtools view {abam} | awk -v TAG='{desired_tag}' '{awk_string}' " 
            info = subprocess.run(cmd, shell = True, capture_output=True).stdout.decode('utf-8').split('\n')
            info.remove('')
            align_info[desired_tag] = info
            
        ## Now figure out which is reads to deduplicate
        gp = list(zip(align_info['CB'], align_info['TR']))
        gp_quant = collections.Counter(gp)
        ur = align_info['UR']
        uy = [min(list(map(lambda x: ord(x) - 33, s))) for s in align_info['UY']]
        umi_phred_all = pl.DataFrame({'idx': list(range(len(uy))), 'uy': uy})
        umi_df = pl.DataFrame({
            'umi_list': ur.copy(),
            'is_dup': False * len(ur),
            'ur': align_info['UR'],
            'cb': align_info['CB'],
            'tr': align_info['TR'],
        })
        for gpx in set(gp):
            if gp_quant[gpx]<=1:
                continue
            bool_vector = pl.Series(gp).is_in([gpx])
            original_umis = pl.Series(ur).filter(bool_vector).to_list()
            corrected_umis = cluster_umi(original_umis, return_corrected_umis=True, logger=logger)
            umi_phred = umi_phred_all.filter(bool_vector)
            best_idx = umi_phred.filter(pl.col('uy') == umi_phred['uy'].max())['idx'].to_list()[0]
            duplicated_idx = umi_phred.filter(umi_phred['idx']!=best_idx)['idx'].to_list()
            umi_df = umi_df.with_columns(
                pl.when(bool_vector).then(pl.lit(corrected_umis[0])).otherwise(pl.col('umi_list')).alias('umi_list'),
                pl.when(umi_phred_all['idx'].is_in(duplicated_idx)).then(True).otherwise(pl.col('is_dup')).alias('is_dup')
            )

        ## Now add duplicated tag if duplicated and update the UB
        cell_gene_counts = {}
        with pysam.AlignmentFile(abam, 'rb') as infile:
            with pysam.AlignmentFile(final_bam_out, 'wb', header = infile.header) as outfile:
                i = 0
                for read in infile:
                    logger.info(f'Processing {read.query_name}')
                    if i > umi_df.shape[0]:
                        logger.error(f'Number of reads in {abam} does not match rows in umi_df')
                        sys.exit(1)
                    tag_info_for_read = umi_df[i]
                    if read.get_tag('CB') != tag_info_for_read['cb'].item():
                        logger.error(f'CB in umi_df does not match alignment')
                        sys.exit(1)
                    if read.get_tag('TR') != tag_info_for_read['tr'].item():
                        logger.error(f'TR in umi_df does not match alignment')
                        sys.exit(1)
                    if read.get_tag('UR') != tag_info_for_read['ur'].item():
                        logger.error(f'UR in umi_df does not match alignment')
                        sys.exit(1)
                    read.set_tag("UB", tag_info_for_read['umi_list'].item(), "Z")
                    if tag_info_for_read['is_dup'].item()!=0:
                        if not keep_all_reads:
                            continue
                        read.set_tag("DP", 1, "i")
                    else:
                        read.set_tag("DP", 0, "i")
                        cell_gene_tuple = (read.get_tag('CB'), read.get_tag('GN'))
                        if cell_gene_tuple in cell_gene_counts.keys():
                            cell_gene_counts[cell_gene_tuple] += 1
                        else:
                            cell_gene_counts[cell_gene_tuple] = 1
                            
                    outfile.write(read)
                    i += 1

        logger.info(f'{umi_df['is_dup'].sum()} reads duplicated in {abam}')
        final_bam_outs.append(final_bam_out)

        ## Consolidate cell_gene_counts as a long dataframe
        with open(final_geneexpr_out, 'w') as outgenefile:
            outgenefile.write(
                ','.join(["cell", "gene", "numi"]) + '\n'
            )
            for (cell, gene), counts in cell_gene_counts.items():
                outgenefile.write(
                    ','.join([str(cell), str(gene), str(counts)]) + '\n'
                )
        
    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)
    return final_bam_outs


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################
## ----- Step 8: SNV calling ----- ##
#####################################

## Below are some internal functions for use inside sc_longshot()
## This is to facilitate parallel processing
def _sc_longshot_bulk(
    chrx: str,
    exon_split_out: str | None = None,
    genome_fa: str | None = None,
    cell_dir: str | None = None,
    data_name: str | None = None,
    work_dir: str | None = None,
):
    r"""Internal function for chromosome-wise SNV calling of bulk BAM with longshot."""
    logger, work_dir = setup_workdir(work_dir)
    bulk_vcf_out = f'{cell_dir}/{data_name}_chr{chrx}_bulkchr_candidate.vcf'
    if not os.path.exists(f'{bulk_vcf_out}.gz.csi'):
        if not os.path.exists(f'{bulk_vcf_out}.gz'):
            if not os.path.exists(bulk_vcf_out):
                ## Note the slightly higher criteria for calling from bulk
                cmd = f'longshot \
                    --bam {exon_split_out} \
                    --ref {genome_fa} \
                    --min_alt_count 2 \
                    --min_cov 2 \
                    --out {bulk_vcf_out} \
                    --region {chrx}'
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
            cmd = f'bgzip {bulk_vcf_out}' #Needs tabix
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
        cmd = f'bcftools index {bulk_vcf_out}.gz'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
    return bulk_vcf_out

##
def _sc_longshot_scfirst(
    cell_bam: str,
    genome_fa: str | None = None,
    cell_dir: str | None = None,
    data_name: str | None = None,
    work_dir: str | None = None,
):
    r"""Internal function for first-pass cell-wise SNV calling with longshot."""
    logger, work_dir = setup_workdir(work_dir)

    logger.info(f'Processing {cell_bam}')
    cell_name = re.sub('_cell.bam$', '', os.path.basename(cell_bam))
    
    ## First, index the file, which longshot needs
    ## -M probably not needed here since we're doing one file at a time
    ## But will keep to be safe
    cell_index = f'{cell_bam}.bai'
    # files_to_delete.append(cell_index)
    if not os.path.exists(cell_index):
        cmd = f'samtools index -M {cell_bam}' 
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
    ## Now that index file exists, we can run longshot
    sc_vcf_out = f'{cell_dir}/{cell_name}_cell_candidate.vcf'
    # files_to_delete.append(f'{sc_vcf_out}.gz')
    # files_to_delete.append(f'{sc_vcf_out}.gz.csi')
    if not os.path.exists(f'{sc_vcf_out}.gz.csi'):
        if not os.path.exists(f'{sc_vcf_out}.gz'):
            if not os.path.exists(sc_vcf_out):
                cmd = f'longshot \
                    --bam {cell_bam} \
                    --ref {genome_fa} \
                    --min_alt_count 1 \
                    --min_cov 1 \
                    --sample_id {cell_name} \
                    --out {sc_vcf_out}'
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
            ## For downstream merging, need to bgzip and index (aka tabix in the original script)
            cmd = f'bgzip {sc_vcf_out}' #Needs tabix
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
        cmd = f'bcftools index {sc_vcf_out}.gz'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
    return sc_vcf_out

## 
def _sc_longshot_processmerged(
    chrx: str,
    sample_split: bool = True,
    tmp1: str | None = None,
    tmp_file_prefix: str = 'temp02_merged',
    cell_dir: str | None = None,
    data_name: str | None = None,
    work_dir: str | None = None,
):
    r"""Internal function for chromosome-wise processing of VCFs post-merging."""
    logger, work_dir = setup_workdir(work_dir)
    
    tmp2 = f'{cell_dir}/{tmp_file_prefix}_chr{chrx}.vcf.gz'
    if not os.path.exists(tmp2):
        if sample_split:
            cmd = f'bcftools norm --regions {chrx} --atomize --multiallelics - {tmp1} -Ou | \
            bcftools +fill-tags -- -t all | bcftools view -s SAMPLE -Oz -o {tmp2}'
        else:
            cmd = f'bcftools norm --regions {chrx} --atomize --multiallelics - {tmp1} -Ou | \
            bcftools +fill-tags -- -t all | bgzip > {tmp2}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    if not os.path.exists(tmp2):
        logger.error(f'Did not generate tmp2 file: {proc.stderr}')
        sys.exit(1)
    if not os.path.exists(f'{tmp2}.csi'):
        cmd = f'bcftools index {tmp2}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    return tmp2

##
def _sc_longshot_scsecond(
    cell_bam: str,
    scsecond_full: bool = False,
    vcf_candidates_out: str | None = None,
    chr_sizes_file: str | None = None,
    max_target_depth: int | None = None,
    genome_fa: str | None = None,
    cell_dir: str | None = None,
    threads: int | None = None, 
    data_name: str | None = None,
    work_dir: str | None = None,
    seed: int = 12345,
):
    r"""Internal function for second-pass cell-wise SNV calling with longshot using candidates."""
    logger, work_dir = setup_workdir(work_dir)
    
    logger.info(f'Processing {cell_bam}')
    cell_name = re.sub('_cell.bam$', '', os.path.basename(cell_bam))
    cell_clipped_bam = f'{cell_dir}/{cell_name}_cellclipped.bam'
    cell_final_out = f'{cell_dir}/{cell_name}_finalgenotypes.vcf.gz'

    if threads is None:
        threads = threading.active_count()

    if os.path.exists(cell_final_out):
        if not os.path.exists(f'{cell_final_out}.csi'):

            ## Check that file is actually gzip
            try:
                with gzip.open(cell_final_out, 'rb') as f:
                    f.read(1)
            except gzip.BadGzipFile:
                logger.warning(f'Not actually a gzipped file: correcting...')
                cell_final_out_unzipped = f'{cell_dir}/{cell_name}_finalgenotypes.vcf'
                cmd = f'mv {cell_final_out} {cell_final_out_unzipped} && bgzip {cell_final_out_unzipped}'
                print(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
                if not os.path.exists(cell_final_out):
                    logger.error(f'Something has gone wrong: {cell_final_out} no longer exists')
                    sys.exit(1)

            ## Check again after fixing
            try:
                with gzip.open(cell_final_out, 'rb') as f:
                    f.read(1)
            except gzip.BadGzipFile:
                logger.error(f'Unable to solve gzip problem')

            ## If no issues, then try bcftools index
            proc = subprocess.run(
                f'bcftools index {cell_final_out}',
                shell=True, capture_output=True)
            
            if not os.path.exists(f'{cell_final_out}.csi'):
                 logger.error(f'{cell_final_out} exists but cannot be indexed')
                 sys.exist(1)
            
        logger.warning(f'{cell_final_out} already exists...Skipping...')
        return cell_final_out
    
    if not os.path.exists(cell_clipped_bam):
        ## For better memory usage, just only consider VCFs that overlap a given cell's BAM 
        tmp_cell_candidate_out = f'{cell_dir}/{cell_name}_cellshortlist.vcf.gz'
        cmd = f'bedtools intersect -header -u -a {vcf_candidates_out} -b {cell_bam} \
        | bgzip > {tmp_cell_candidate_out}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        n_vcfs = subprocess.run(f'zcat {tmp_cell_candidate_out} | grep -c -v "^#"', shell=True, capture_output=True)
        n_vcfs = int(n_vcfs.stdout.decode('utf-8'))
        logger.info(f'{n_vcfs} candidate VCFs for {cell_name}')
        if n_vcfs == 0:
            logger.warning(f"No VCFs found...Cleaning up and moving to next BAM...")
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
            return
    
        ## Index the file
        if not os.path.exists(f'{tmp_cell_candidate_out}.csi'):
            cmd = f'bcftools index {tmp_cell_candidate_out}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
    
        ## Depth correction of BAM so that max depth is set to 200x
        ## First, make a BED of locations for depth clipping
        ## There is some dilation by 100 bp in both directions (-b).
        tmp_bed_regions = f'{cell_dir}/{cell_name}_tmpcoverage.bed'
        awk_string = r'BEGIN {OFS="\t"} {print $1, $2, $2}'
        cmd = f"zgrep -v '^#' {tmp_cell_candidate_out} | awk '{awk_string}' | uniq \
        | bedtools slop -i - -g {chr_sizes_file} -b 100 \
        | bedtools merge > {tmp_bed_regions}"
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        if not os.path.exists(tmp_bed_regions):
            logger.warning(f"No {tmp_bed_regions} found...Assuming no regions to handle...Cleaning up and moving to next BAM...")
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
            return
        if os.path.getsize(tmp_bed_regions)==0:
            logger.warning(f"{tmp_bed_regions} has no data...Assuming no regions to handle...Cleaning up and moving to next BAM...")
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
            return
    
        ## Assuming the bed file creation was successful, now figure out the depth per region with mosdepth
        ## The tmp_bed_regions file is no longer needed after this step, so just delete
        tmp_mosdepth_prefix = f'{cell_dir}/{cell_name}_md'
        cmd = f'mosdepth --no-per-base --by {tmp_bed_regions} --threads {threads} --fast-mode {tmp_mosdepth_prefix} {cell_bam}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        proc = subprocess.run(f'rm {tmp_bed_regions}', shell=True, capture_output=True)
        mos_file = f'{tmp_mosdepth_prefix}.regions.bed.gz'
        if not os.path.exists(mos_file):
            logger.warning(f"No {mos_file} generated...Assuming no depth information available...Cleaning up and moving to next BAM...")
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
            for mdfile in glob.glob(f'{tmp_mosdepth_prefix}*'):
                if os.path.exists(mdfile):
                    os.remove(mdfile)
            return
    
        ## Have a look at the bed.gz. If there exist reads > 110% of target coverage, then clip
        ## In the EPI2ME, they used csvtk for CSV handling -- but not necessary (we have pl.)
        try:
            coverage_df = pl.read_csv(mos_file, has_header=False, separator='\t', new_columns=['chrom', 'start', 'end', 'depth'])
        except Exception as e:
            logger.warning(f'Coverage depth csv read file error: {e}')
            logger.warning(f'Assuming insufficient coverage for this cell for SNV calling...Skipping {cell_name} entirely...')
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
            proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
            return

        ## If no exception, then csv is populated and SNVs can be called
        if not os.path.exists(cell_clipped_bam):
            if any(coverage_df['depth'] > round(max_target_depth * 1.1)):
                
                logger.info(f'Performing depth clipping...')
        
                ## Start by getting high coverage regions as a bed
                high_cov_bed = f'{tmp_mosdepth_prefix}.high_coverage.bed'
                coverage_df.filter(
                    coverage_df['depth'] > round(max_target_depth * 1.1)
                ).write_csv(
                    high_cov_bed, include_header=False, separator='\t'
                )
        
                ## Get read_ids in high coverage regions
                ## If read spans multiple regions (it shouldn't), we take the first unique value
                clip_bed = f'{tmp_mosdepth_prefix}.to_clip.bed'
                cmd = f"bedtools intersect -a {cell_bam} -b {high_cov_bed} -wb -bed | cut -f 4,16 | sort -u -k1,1 > {clip_bed}"
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
                clip_bed_df = pl.read_csv(clip_bed, has_header = False, separator='\t', new_columns=['read_id', 'depth'])
        
                ## See https://github.com/epi2me-labs/wf-single-cell/blob/master/bin/workflow_glue/clip_depth.py for depth clipping
                ## Essentially they just randomly generate a number between 0 and 1, and if it it exceeds target_depth / window_depth, drop
                random.seed(seed)
                with pysam.AlignmentFile(cell_bam, 'rb', check_sq=False) as inbam:
                    with pysam.AlignmentFile(cell_clipped_bam, 'wb', template=inbam) as outbam:
                        for aln in inbam.fetch(until_eof=True):
                            if aln.query_name in clip_bed_df['read_id']:
                                window_depth = clip_bed_df.filter(pl.col('read_id')==aln.query_name)['depth'].to_list()
                                if len(window_depth) > 1:
                                    logger.warning(f'Multiple entries for {aln.query_name} in {clip_bed}: defaulting to first entry')
                                window_depth = window_depth[0]
                                if random.random() > max_target_depth / window_depth:
                                    continue # Discard the read by probability of target_depth / window_depth
                            outbam.write(aln)
            else:
                logger.info(f'No depth clipping necessary')
                proc = subprocess.run(f'cp {cell_bam} {cell_clipped_bam}', shell=True, capture_output=True)
        
        ## Delete all the unneeded mosdepth files
        for mdfile in glob.glob(f'{tmp_mosdepth_prefix}*'):
            if os.path.exists(mdfile):
                os.remove(mdfile)
    
    ## This is the final genotyping
    ## Note that the command outputs second_call_out, 
    ## but we'll actually be taking the results inside variant_debug_dir
    ## So we can delete the second_call_out almost immediately
    logger.info(f'Running second pass of VCF calling on {cell_clipped_bam}...')  
    if not os.path.exists(f'{cell_clipped_bam}.bai'):
        ## Index does not exist -- so create
        cmd = f'samtools index {cell_clipped_bam}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        # files_to_delete.append(cell_clipped_bam)
    if scsecond_full:
        logger.warning(f'Running full scsecond approach to get whole genome -- high memory needed and slow')
        second_call_out = f'{cell_dir}/{cell_name}_cellsecondpass.vcf'
        variant_debug_dir = f'{cell_dir}/{cell_name}_LONGSHOT_DEBUG/'
        cmd = f'longshot \
        --potential_variants {vcf_candidates_out} \
        --bam {cell_clipped_bam} \
        --ref {genome_fa} \
        --min_alt_count 1 \
        --min_cov 1 \
        --sample_id {cell_name} \
        --force_overwrite \
        --out {second_call_out} \
        --variant_debug_dir {variant_debug_dir}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    
        ## Now go into variant_debug_dir and get the file we want
        logger.info(f'Creating and indexing {cell_final_out} from {variant_debug_dir}/4.0.final_genotypes.vcf')
        if not os.path.exists(f'{variant_debug_dir}/4.0.final_genotypes.vcf'):
            proc = subprocess.run(
                f'mv {variant_debug_dir}/1.0.potential_SNVs.vcf {variant_debug_dir}/4.0.final_genotypes.vcf', 
                shell=True, capture_output=True)
        proc = subprocess.run(
            f'bgzip {variant_debug_dir}/4.0.final_genotypes.vcf', 
            shell=True, capture_output=True)
        proc = subprocess.run(
            f'cp {variant_debug_dir}/4.0.final_genotypes.vcf.gz {cell_final_out}', 
            shell=True, capture_output=True)
        
        proc = subprocess.run(f'rm {second_call_out}', shell=True, capture_output=True)
        proc = subprocess.run(f'rm -rf {variant_debug_dir}', shell=True, capture_output=True)
    else:
        ## Faster to not consider variants debug dir
        cell_final_out_unzipped = f'{cell_dir}/{cell_name}_finalgenotypes.vcf'
        cmd = f'longshot \
        --potential_variants {vcf_candidates_out} \
        --bam {cell_clipped_bam} \
        --ref {genome_fa} \
        --min_alt_count 1 \
        --min_cov 1 \
        --sample_id {cell_name} \
        --force_overwrite \
        --out {cell_final_out_unzipped}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        cmd = f'bgzip {cell_final_out_unzipped}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
    proc = subprocess.run(
        f'bcftools index {cell_final_out}', 
        shell=True, capture_output=True)
    # files_to_delete.append(cell_final_out)
    # files_to_delete.append(f'{cell_final_out}.csi')
    
    ## Mandatory cleanup
    proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
    proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
    proc = subprocess.run(f'rm {tmp_bed_regions}', shell=True, capture_output=True)

    return cell_final_out

## ~~~ ##

def sc_longshot(
    ## General inputs
    dbam: str | Path,
    work_dir: str | Path | None = None,
    n_workers_snv: int | None = None,
    n_vcfs_per_merge: int = 20000,

    ## Specific inputs
    scsecond_full: bool = False,
    samtools_split: bool = False,
    file_cap: int | None = None,
    max_target_depth: int = 200,
    threads: int | None = None,
    seed: int = 12345
):
    r"""Longshot protocol for calling SNVs.

    Briefly:
    1. Remove the duplicated reads from the BAM.
    2. Exon split the BAM by exon using `gatk` (which requires the `reference.fa` to be procesed
    by `gatk CreateSequenceDictionary`). 
    3. Split the exon-splitted BAM into cellwise BAMs. EPI2ME defaults to samtools, but `pysam` work just as well (and can report which read being processed).  
    4. For each cell, run the SNV calling protocol. May also need to call SNVs from the bulk.
    5. Consolidate SNV calls into an easy-to-read format.

    Parameters
    ---------- 
    dbam
        A single string or Path indicating the bulk BAM file with "DP" tags. Will deduplicate internally.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    n_workers_snv
        Integer indicating number of workers per `multiprocessing.Pool`. If not specified, defaults to `os.cpu_count()`.
    n_vcfs_per_merge
        Integer indicating number of VCF files processed per loop during `bcftools merge`. This is done to avoid overloading memory.
    scsecond_full
        Boolean indicating whether to use the slow method to get VCFs per cell or not. Slow method returns full genome. Fast method (default) only returns SNVs that are not 0/0 per cell.
        Full set of SNVs can be obtained from the _allcandidates.vcf.gz file. 
    samtools_split
        Boolean indicating if should use `samtools` to split the file by cellular barcodes (CB) tags. If no, then will loop through the bam file line by line (with `pysam`), and sort.
        `samtools` approach works ok when n cells / n files is small. Recommend to do this line by line otherwise.
    file_cap
        Integer indicating maximum number of open files allowed. Should be set to more than the number of cells, since every cell will yield a single bam file. 
        If unspecified, defaults to the current `ulimit -n`.  
    target_depth
        Integer indicating max number of deduplicated reads per region to target. 
    threads
        Integer indicating number of threads for `samtools` etc.
    seed
        Integer indicating random seed.
    """
    ## Basic utility
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []

    ## There is a problem with polars working across multiprocessing pools.
    ## Note that this is a non-issue for Python 3.14, which sets default multiprocessing behaviour to "spawn".
    mp = multiprocessing.get_context("spawn")
    if n_workers_snv is None:
        n_workers_snv = os.cpu_count()
    n_workers = n_workers_snv
    
    ## Prep information
    data_name = re.sub('_annotated_deduplicated.bam$', '', str(os.path.basename(dbam)))
    out_dir = str(os.path.dirname(dbam))
    cell_dir = f'{out_dir}/{data_name}_CELLSPLITBAMS/'
    os.makedirs(cell_dir, exist_ok=True)
    if threads is None:
        threads = threading.active_count() # High RAM requirement -- keep threads as 1
        # threads = round(n_workers * 0.8)
    
    ## Check we have needed extra files
    chr_sizes_file = f'{work_dir}/GENOME.chr.sizes'
    if not os.path.exists(chr_sizes_file):
        logger.error(f'Cannot find {chr_sizes_file}: file that contains chromosome lengths')
        sys.exit(1)
    
    ## Step 0: Use GATK to generate GENOME.dict
    ## This is really only necessary if GATK is used, so will keep this here instead of get_genome()
    genome_fa = f'{work_dir}/GENOME.fa'
    genome_dict = f'{work_dir}/GENOME.dict'
    if not os.path.exists(genome_fa):
        logger.error(f"Cannot find {genome_fa}")
        sys.exit(1)
    if not os.path.exists(genome_dict):
        logger.info(f'Creating sequence dictionary for {genome_fa}')
        cmd = f'gatk CreateSequenceDictionary --REFERENCE {genome_fa}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    
    working_bam = dbam
    exon_split_out = f'{cell_dir}/{data_name}_exonsplit.bam'
    if not os.path.exists(exon_split_out):
    
        ## Step 1: Deduplicate.
        ## This is only needed if DP:1 tag exists in the file.
        ## Check that the DP tag exists:
        check_dp_tag = subprocess.run(f'samtools view {dbam} | head -n 100 | grep -q "DP:"', shell=True, capture_output=True, text=True).stdout
        if len(check_dp_tag) > 0:
            if subprocess.run(f'samtools view -d DP:1 {dbam}', shell=True, capture_output=True).stdout and (len(check_dp_tag)>0):
                logger.info(f'DP:1 tag found in file -- assuming deduplication is necessary')
                working_bam = f'{cell_dir}/{data_name}_WORKING_TEMP.bam'
                cmd = f'samtools view -h -d DP:0 {dbam} > {working_bam}'
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
            else:
                logger.info(f'No DP:1 tags found -- assuming no deduplication needed')
        else:
            logger.info(f'No DP tags found -- assuming no deduplication needed')
            
        ## Step 2: Exon split with GATK.
        ## Much faster if done on the whole BAM rather than individual BAMs.
        ## If tmp dir not specified, the function might quietly fail
        logger.info(f'Exon splitting of {working_bam}')
        tmp_dir = f'{cell_dir}TMP/'
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir, exist_ok=True)
        cmd = f'gatk SplitNCigarReads -R {genome_fa} \
        --input {working_bam} \
        --output {exon_split_out} \
        --tmp-dir {tmp_dir} \
        --java-options "-Xmx8g" \
        --max-reads-in-memory 100000 \
        --create-output-bam-index true \
        --verbosity ERROR'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    
        if working_bam != dbam:
            logger.info(f'Deleting {working_bam}')
            cmd = f"rm {working_bam}"
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
    
    files_to_delete.append(exon_split_out)
    files_to_delete.append(f'{exon_split_out}.bai')

    ## Check that exon_split_out is of correct size
    if os.path.getsize(exon_split_out)==0:
        logger.error(f'Writing {exon_split_out} has failed: a 0 byte object is returned instead. Likely a memory issue: ensure sufficient space on disk.')
        sys.exit(1)

    ## If file-cap is none, specify to ulimit -n
    if file_cap is None:
        actual_limit = int(subprocess.run(f'ulimit -n', shell=True, capture_output=True, text=True).stdout)
        logger.info(f'Setting file_cap to current max limit of {actual_limit} files allowed open')
        file_cap = actual_limit

    ## Step 3: Generate cell-wise BAMs
    file_check = f'{cell_dir}/cellwise_bams_generated.txt'
    if not os.path.exists(file_check):
        logger.info(f'Generating cell-wise BAMs in {cell_dir}')
    
        if samtools_split:
            logger.info(f'Splitting BAM using samtools')
            cell_bam_format = f'{cell_dir}/%!_cell.bam'
            if file_cap > 20000:
                logger.warning(f'High file_cap ({file_cap}): consider line by line instead')
            cmd = f'ulimit -n {file_cap}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
            actual_limit = int(subprocess.run(f'ulimit -n', shell=True, capture_output=True, text=True).stdout)
            logger.info(f'Current max limit of {actual_limit} files allowed to be open')
            if actual_limit < file_cap:
                logger.error(f'Max limit of files is less than specified file_cap of {file_cap}')
                sys.exit(1)
        
            cmd = f'samtools split --max-split {file_cap} -d CB --threads {threads} -f "{cell_bam_format}" --no-PG {exon_split_out}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
            
        else:
            logger.info(f'Splitting BAM line-by-line using pysam')
            with pysam.AlignmentFile(exon_split_out, 'rb' ) as fin:
                fouts_dict = {}
                for read in fin:
                    logger.info(f'Processing {read.query_name}')
                    try:
                        tag_value = read.get_tag('CB') #Cell barcode
                    except KeyError:
                        continue #If tag missing
                    if tag_value not in fouts_dict.keys():
                        output_file = f'{cell_dir}/{tag_value}_cell.bam'
                        fouts_dict[tag_value] = pysam.AlignmentFile(output_file, "wb", template=fin)
                    fouts_dict[tag_value].write(read)
                ## Once all reads done, close all files
                for fout in fouts_dict.values():
                    fout.close()
                    
        Path(file_check).touch()
    else:
        logger.info(f'{file_check} found: skipping cell-wise bam generation')

    cell_bams = glob.glob(f'{cell_dir}/*_cell.bam')
    cell_bams.sort()
    files_to_delete.extend(cell_bams)
    files_to_delete.append(file_check)
    
    ## Step 4A: Call bulk VCF with longshot
    ## Split by chromosome for better memory handling
    logger.info(f'First pass bulk VCF calling')
    chromosomes = [re.sub('chr.gtf$', '', os.path.basename(x)) for x in glob.glob(f'{work_dir}/*chr.gtf')]
    chromosomes.sort()
    # for chrx in chromosomes:
    #     bulk_vcf_out = f'{cell_dir}/{data_name}_chr{chrx}_bulkchr_candidate.vcf'
    #     files_to_delete.append(bulk_vcf_out)
    #     if not os.path.exists(bulk_vcf_out):
    #         ## Note the slightly higher criteria for calling from bulk
    #         cmd = f'longshot \
    #             --bam {exon_split_out} \
    #             --ref {genome_fa} \
    #             --min_alt_count 2 \
    #             --min_cov 2 \
    #             --out {bulk_vcf_out} \
    #             --region {chrx}'
    #         logger.info(cmd)
    #         proc = subprocess.run(cmd, shell=True, capture_output=True)
    pp_sc_longshot_bulk = functools.partial(
        _sc_longshot_bulk,
        exon_split_out = exon_split_out,
        genome_fa = genome_fa,
        cell_dir = cell_dir,
        data_name = data_name,
        work_dir = work_dir,
    )
    with mp.Pool(round(n_workers * 0.5)) as pool:
        bulk_vcf_outs = pool.map(pp_sc_longshot_bulk, chromosomes)
    files_to_delete.extend(bulk_vcf_outs)
    
    ## Step 4B: For each cell BAM, generate VCF
    logger.info(f'First pass single cell VCF calling')
    # for cell_bam in cell_bams:
    #     logger.info(f'Processing {cell_bam}')
    #     cell_name = re.sub('_cell.bam$', '', os.path.basename(cell_bam))
        
    #     ## First, index the file, which longshot needs
    #     ## -M probably not needed here since we're doing one file at a time
    #     ## But will keep to be safe
    #     cell_index = f'{cell_bam}.bai'
    #     files_to_delete.append(cell_index)
    #     if not os.path.exists(cell_index):
    #         cmd = f'samtools index -M {cell_bam}' 
    #         logger.info(cmd)
    #         proc = subprocess.run(cmd, shell=True, capture_output=True)
            
    #     ## Now that index file exists, we can run longshot
    #     sc_vcf_out = f'{cell_dir}/{cell_name}_cell_candidate.vcf'
    #     files_to_delete.append(f'{sc_vcf_out}.gz')
    #     files_to_delete.append(f'{sc_vcf_out}.gz.csi')
    #     if not os.path.exists(sc_vcf_out):
    #         cmd = f'longshot \
    #             --bam {cell_bam} \
    #             --ref {genome_fa} \
    #             --min_alt_count 1 \
    #             --min_cov 1 \
    #             --sample_id {cell_name} \
    #             --out {sc_vcf_out}'
    #         logger.info(cmd)
    #         proc = subprocess.run(cmd, shell=True, capture_output=True)
    #         ## For downstream merging, need to bgzip and index (aka tabix in the original script)
    #         cmd = f'bgzip {sc_vcf_out}' #Needs tabix
    #         logger.info(cmd)
    #         proc = subprocess.run(cmd, shell=True, capture_output=True)
    #         cmd = f'bcftools index {sc_vcf_out}.gz'
    #         logger.info(cmd)
    #         proc = subprocess.run(cmd, shell=True, capture_output=True)
    pp_sc_longshot_scfirst = functools.partial(
        _sc_longshot_scfirst,
        genome_fa = genome_fa,
        cell_dir = cell_dir,
        data_name = data_name,
        work_dir = work_dir,
    )
    with mp.Pool( round(n_workers * 0.5) ) as pool: #Reduce number of workers a little to avoid overloading memory
        sc_vcf_outs = pool.map(pp_sc_longshot_scfirst, cell_bams)
    files_to_delete.extend([f'{x}.gz' for x in sc_vcf_outs])
    files_to_delete.extend([f'{x}.gz.csi' for x in sc_vcf_outs])
    files_to_delete.extend([f'{x}.bai' for x in cell_bams])            
    
    ## Step 5: Merge together all the .vcfs generated so far into a "MERGED.vcf"
    ## Merge all the bulks first (assuming they are chromsome-specific)
    bulk_vcfs = glob.glob(f'{cell_dir}/*_bulkchr_candidate.vcf.gz')
    sc_vcfs = glob.glob(f'{cell_dir}/*_cell_candidate.vcf.gz')

    bulk_filelist_txt = f'{cell_dir}BULKS_FILELIST.txt'
    with open(bulk_filelist_txt, 'w') as f:
        for line in bulk_vcfs:
            f.write(line + '\n')

    files_to_delete.append(f'{bulk_filelist_txt}')

    ## Step 5A: Merge the bulks
    merged_vcf_bulk = f'{cell_dir}/{data_name}_allbulks.vcf.gz'
    if not os.path.exists(merged_vcf_bulk):
        logger.info(f'Merging first pass bulk VCF calls')
        cmd = f'bcftools concat --file-list {bulk_filelist_txt} --threads {threads} \
            --output-type z -o {merged_vcf_bulk}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        ## Again, for merging, need the index
        cmd = f'bcftools index {merged_vcf_bulk}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
    # files_to_delete.append(merged_vcf_bulk) #Might be good to keep the bulk
    files_to_delete.append(f'{merged_vcf_bulk}.csi')
    
    ## Step 5B: Merge all.
    ## Here, we are atomizing to split MNVs into consecutive SNVs, and split multiallelics into biallelic records
    ## fill-tags fills in missing tags
    vcf_candidates_out = f'{cell_dir}/{data_name}_allcandidates.vcf.gz'
    if not os.path.exists(vcf_candidates_out):
        logger.info(f'Merging first pass bulk and single cell VCF calls')
        ## Too much memory used
        # cmd = f'bcftools merge --file-list {all_filelist_txt} \
        # | bcftools norm --atomize --multiallelics - \
        # | bcftools +fill-tags - -- -t all \
        # | bcftools view -s SAMPLE - --threads {threads} --output-type z -o {vcf_candidates_out}'

        already_done_tracking = f'{cell_dir}/MERGED_DONE.txt'
        tmp1 = f'{cell_dir}/temp01_merged.vcf.gz'
        already_done_list = []
        tmp1new = f'{cell_dir}/temp01new_merged.vcf.gz'
        if not os.path.exists(already_done_tracking):
            cmd = f'bcftools concat --file-list {bulk_filelist_txt} -Oz -o {tmp1} && bcftools index {tmp1}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
        else:
            with open(already_done_tracking, 'r') as f:
                already_done_list = [line.rstrip('\n') for line in f]
        scvcfs_to_process = list(set(sc_vcfs))
        if len(already_done_list) > 0:
            scvcfs_to_process = list(set(sc_vcfs) - set(already_done_list))
        scvcfs_to_process.sort()
        tmp_filelist_txt = f'{cell_dir}TEMP_FILELIST.txt'
        while len(scvcfs_to_process) > 0:
            ## This is to control amount of memory being used -- we only process X cells at a time
            logger.info(f'Files left to process: {len(scvcfs_to_process)}')
            sub_scvcfs = scvcfs_to_process[:n_vcfs_per_merge]
            del scvcfs_to_process[:n_vcfs_per_merge]
            with open(tmp_filelist_txt, 'w') as f:
                f.write(tmp1 + '\n')
                for line in sub_scvcfs:
                    f.write(line + '\n')
            cmd = f'bcftools merge --threads {n_workers} --file-list {tmp_filelist_txt} -Oz -o {tmp1new} && mv {tmp1new} {tmp1} \
            && rm {tmp1}.csi && bcftools index {tmp1} && rm {tmp_filelist_txt}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)
            with open(already_done_tracking, 'a') as ft:
                for line in sub_scvcfs:
                    ft.write(line + '\n')
        if not os.path.exists(f'{tmp1}.csi'):
            cmd = f'bcftools index {tmp1}'
            proc = subprocess.run(cmd, shell=True, capture_output=True)

        ## The next steps can be done on each chromosome separately (for speed)
        # tmp2 = f'{cell_dir}/temp02_merged.vcf.gz'
        # if not os.path.exists(tmp2):
        #     cmd = f'bcftools norm --threads {n_workers} --atomize --multiallelics - {tmp1} -Oz -o {tmp2}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
        # tmp3 = f'{cell_dir}/temp03_merged.vcf.gz'
        # if not os.path.exists(tmp2):
        #     cmd = f'bcftools +fill-tags {tmp2} -- -t all -Oz -o {tmp3} --threads {n_workers}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
        # cmd = f'bcftools view -s SAMPLE {tmp3} --threads {n_workers} -Oz -o {vcf_candidates_out}'
        # logger.info(cmd)
        # proc = subprocess.run(cmd, shell=True, capture_output=True)

        pp1_sc_longshot_processmerged = functools.partial(
            _sc_longshot_processmerged,
            sample_split = True,
            tmp1 = tmp1,
            tmp_file_prefix = 'temp02a_merged',
            cell_dir = cell_dir,
            data_name = data_name,
            work_dir = work_dir,
        )
        with mp.Pool(n_workers) as pool:
            tmp2_vcf_outs = pool.map(pp1_sc_longshot_processmerged, chromosomes)

        tmpbulk_filelist_txt = f'{cell_dir}TMPMERGEDCHR_1_FILELIST.txt'
        with open(tmpbulk_filelist_txt, 'w') as f:
            for line in tmp2_vcf_outs:
                f.write(line + '\n')

        cmd = f'bcftools concat --file-list {tmpbulk_filelist_txt} --threads {threads} \
            --output-type z -o {vcf_candidates_out}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)

        if os.path.exists(vcf_candidates_out):
            os.remove(tmp1)
            for tmp2fx in tmp2_vcf_outs:
                if os.path.exists(tmp2fx):
                    os.remove(tmp2fx)
                if os.path.exists(f'{tmp2fx}.csi'):
                    os.remove(f'{tmp2fx}.csi')
            os.remove(already_done_tracking)
            os.remove(tmpbulk_filelist_txt)
        
    if not os.path.exists(f'{vcf_candidates_out}.csi'):
        cmd = f'bcftools index {vcf_candidates_out}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)
    # files_to_delete.append(vcf_candidates_out)
    files_to_delete.append(f'{vcf_candidates_out}.csi')
    
    ## Step 6: Now with all the candidates, attempt to re-call SNPs cellwise
    final_vcf_out = f'{cell_dir}/{data_name}_finalsnv_cells.vcf.gz'
    if not os.path.exists(final_vcf_out):
        logger.info(f'Second pass single cell VCF calls')
        # for cell_bam in cell_bams:
        #     logger.info(f'Processing {cell_bam}')
        #     cell_name = re.sub('_cell.bam$', '', os.path.basename(cell_bam))
    
        #     ## For better memory usage, just only consider VCFs that overlap a given cell's BAM 
        #     tmp_cell_candidate_out = f'{cell_dir}/{cell_name}_cellshortlist.vcf.gz'
        #     cmd = f'bedtools intersect -header -u -a {vcf_candidates_out} -b {cell_bam} \
        #     | bgzip > {tmp_cell_candidate_out}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
        #     n_vcfs = subprocess.run(f'zcat {tmp_cell_candidate_out} | grep -c -v "^#"', shell=True, capture_output=True)
        #     n_vcfs = int(n_vcfs.stdout.decode('utf-8'))
        #     logger.info(f'{n_vcfs} candidate VCFs for {cell_name}')
        #     if n_vcfs == 0:
        #         logger.warning(f"No VCFs found...Cleaning up and moving to next BAM...")
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
        #         continue
    
        #     ## Index the file
        #     cmd = f'bcftools index {tmp_cell_candidate_out}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
    
        #     ## Depth correction of BAM so that max depth is set to 200x
        #     ## First, make a BED of locations for depth clipping
        #     ## There is some dilation by 100 bp in both directions (-b).
        #     tmp_bed_regions = f'{cell_dir}/{cell_name}_tmpcoverage.bed'
        #     awk_string = r'BEGIN {OFS="\t"} {print \$1, \$2, \$2}'
        #     cmd = f"zgrep -v '^#' {tmp_cell_candidate_out} | awk '{awk_string}' | uniq \
        #     | bedtools slop -i - -g {chr_sizes_file} -b 100 \
        #     | bedtools merge > {tmp_bed_regions}"
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
        #     if not os.path.exists(tmp_bed_regions):
        #         logger.warning(f"No {tmp_bed_regions} found...Assuming no regions to handle...Cleaning up and moving to next BAM...")
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
        #         continue
        #     if os.path.getsize(tmp_bed_regions)==0:
        #         logger.warning(f"{tmp_bed_regions} has no data...Assuming no regions to handle...Cleaning up and moving to next BAM...")
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
        #         continue
    
        #     ## Assuming the bed file creation was successful, now figure out the depth per region with mosdepth
        #     ## The tmp_bed_regions file is no longer needed after this step, so just delete
        #     tmp_mosdepth_prefix = f'{cell_dir}/{cell_name}_md'
        #     cmd = f'mosdepth --no-per-base --by {tmp_bed_regions} --threads {threads} --fast-mode {tmp_mosdepth_prefix} {cell_bam}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
        #     proc = subprocess.run(f'rm {tmp_bed_regions}', shell=True, capture_output=True)
        #     mos_file = f'{tmp_mosdepth_prefix}.regions.bed.gz'
        #     if not os.path.exists(mos_file):
        #         logger.warning(f"No {mos_file} generated...Assuming no depth information available...Cleaning up and moving to next BAM...")
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
        #         proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
        #         for mdfile in glob.glob(f'{tmp_mosdepth_prefix}*'):
        #             if os.path.exists(mdfile):
        #                 os.remove(mdfile)
        #         continue
    
        #     ## Have a look at the bed.gz. If there exist reads > 110% of target coverage, then clip
        #     ## In the EPI2ME, they used csvtk for CSV handling -- but not necessary (we have pl.)
        #     coverage_df = pl.read_csv(mos_file, has_header=False, separator='\t', new_columns=['chrom', 'start', 'end', 'depth'])
        #     cell_clipped_bam = cell_bam
        #     if any(coverage_df['depth'] > round(max_target_depth * 1.1)):
                
        #         logger.info(f'Performing depth clipping...')
    
        #         ## Start by getting high coverage regions as a bed
        #         high_cov_bed = f'{tmp_mosdepth_prefix}.high_coverage.bed'
        #         coverage_df.filter(
        #             coverage_df['depth'] > round(max_target_depth * 1.1)
        #         ).write_csv(
        #             high_cov_bed, include_header=False, separator='\t'
        #         )
    
        #         ## Get read_ids in high coverage regions
        #         ## If read spans multiple regions (it shouldn't), we take the first unique value
        #         clip_bed = f'{tmp_mosdepth_prefix}.to_clip.bed'
        #         cmd = f"bedtools intersect -a {cell_bam} -b {high_cov_bed} -wb -bed | cut -f 4,16 | sort -u -k1,1 > {clip_bed}"
        #         clip_bed_df = pl.read_csv(clip_bed, has_header = False, separator='\t', new_columns=['read_id', 'depth'])
    
        #         ## See https://github.com/epi2me-labs/wf-single-cell/blob/master/bin/workflow_glue/clip_depth.py for depth clipping
        #         ## Essentially they just randomly generate a number between 0 and 1, and if it it exceeds target_depth / window_depth, drop
        #         random.seed(seed)
        #         cell_clipped_bam = f'{cell_dir}/{cell_name}_cellclipped.bam'
        #         with pysam.AlignmentFile(cell_bam, 'rb', check_sq=False) as inbam:
        #             with pysam.AlignmentFile(cell_clipped_bam, 'wb', template=inbam) as outbam:
        #                 for aln in inbam.fetch(until_eof=True):
        #                     if aln.query_name in clip_bed_df['read_id']:
        #                         window_depth = clip_bed_df.filter(pl.col('read_id')==aln.query_name)['depth'].to_list()
        #                         if len(window_depth) > 1:
        #                             logger.warning(f'Multiple entries for {aln.query_name} in {clip_bed}: defaulting to first entry')
        #                         window_depth = window_depth[0]
        #                         if random.random() > max_target_depth / window_depth:
        #                             continue # Discard the read by probability of target_depth / window_depth
        #                     outbam.write(aln)
        #     else:
        #         logger.info(f'No depth clipping necessary')
            
        #     ## Delete all the unneeded mosdepth files
        #     for mdfile in glob.glob(f'{tmp_mosdepth_prefix}*'):
        #         if os.path.exists(mdfile):
        #             os.remove(dfile_path)
            
        #     ## This is the final genotyping
        #     ## Note that the command outputs second_call_out, 
        #     ## but we'll actually be taking the results inside variant_debug_dir
        #     ## So we can delete the second_call_out almost immediately
        #     logger.info(f'Running second pass of VCF calling on {cell_clipped_bam}...')  
        #     if cell_clipped_bam != cell_bam:
        #         ## Index does not exist -- so create
        #         cmd = f'samtools index {cell_clipped_bam}'
        #         logger.info(cmd)
        #         proc = subprocess.run(cmd, shell=True, capture_output=True)
        #         files_to_delete.append(cell_clipped_bam)
        #     second_call_out = f'{cell_dir}/{cell_name}_cellsecondpass.vcf'
        #     variant_debug_dir = f'{cell_dir}/{cell_name}_LONGSHOT_DEBUG/'
        #     cmd = f'longshot \
        #     --potential_variants {vcf_candidates_out}\
        #     --bam {cell_clipped_bam} \
        #     --ref {genome_fa} \
        #     --min_alt_count 1 \
        #     --min_cov 1 \
        #     --sample_id {cell_name} \
        #     --force_overwrite \
        #     --out {second_call_out}'
        #     logger.info(cmd)
        #     proc = subprocess.run(cmd, shell=True, capture_output=True)
    
        #     ## Now go into variant_debug_dir and get the file we want
        #     cell_final_out = f'{cell_dir}/{cell_name}_finalgenotypes.vcf.gz'
        #     logger.info(f'Creating and indexing {cell_final_out} from {variant_debug_dir}/4.0.final_genotypes.vcf')
        #     if not os.path.exists(f'{variant_debug_dir}/4.0.final_genotypes.vcf'):
        #         proc = subprocess.run(
        #             f'mv {variant_debug_dir}/1.0.potential_SNVs.vcf {variant_debug_dir}/4.0.final_genotypes.vcf', 
        #             shell=True, capture_output=True)
        #     proc = subprocess.run(
        #         f'bgzip {variant_debug_dir}/4.0.final_genotypes.vcf', 
        #         shell=True, capture_output=True)
        #     proc = subprocess.run(
        #         f'cp {variant_debug_dir}/4.0.final_genotypes.vcf.gz {cell_final_out}', 
        #         shell=True, capture_output=True)
        #     proc = subprocess.run(
        #         f'bcftools index {cell_final_out}', 
        #         shell=True, capture_output=True)
        #     files_to_delete.append(cell_final_out)
        #     files_to_delete.append(f'{cell_final_out}.csi')
            
        #     ## Mandatory cleanup
        #     proc = subprocess.run(f'rm {tmp_cell_candidate_out}', shell=True, capture_output=True)
        #     proc = subprocess.run(f'rm {tmp_cell_candidate_out}.csi', shell=True, capture_output=True)
        #     proc = subprocess.run(f'rm {tmp_bed_regions}', shell=True, capture_output=True)
        #     proc = subprocess.run(f'rm {second_call_out}', shell=True, capture_output=True)
        #     proc = subprocess.run(f'rm -rf {variant_debug_dir}', shell=True, capture_output=True)
        pp_sc_longshot_scsecond = functools.partial(
            _sc_longshot_scsecond,
            vcf_candidates_out = vcf_candidates_out,
            scsecond_full = scsecond_full,
            chr_sizes_file = chr_sizes_file,
            max_target_depth = max_target_depth,
            genome_fa = genome_fa,
            cell_dir = cell_dir,
            threads = threads,
            data_name = data_name,
            work_dir = work_dir,
            seed = seed,
        )
        with mp.Pool( round(n_workers * 0.25) ) as pool: #Reduce number of workers to avoid overloading memory
            cell_final_outs = pool.map(pp_sc_longshot_scsecond, cell_bams)
        # files_to_delete.extend(cell_final_outs)
        files_to_delete.extend([f'{x}.csi' for x in cell_final_outs])
        
        ## Once all files have been processed, consolidate
        # cell_final_outs = glob.glob(f'{cell_dir}/*_finalgenotypes.vcf.gz')
        if cell_final_outs is None:
            logger.error(f'No single cell calls at all')
            sys.exit(1)
        cell_final_outs.sort()
        logger.info(f'{len(cell_final_outs)} cell-wise final VCF calls found...')

        if not os.path.exists(final_vcf_out):
            already_done_tracking = f'{cell_dir}/MERGED_DONE.txt'
            tmp1 = f'{cell_dir}/temp01b_merged.vcf.gz'
            already_done_list = []
            tmp1new = f'{cell_dir}/temp01bnew_merged.vcf.gz'
            if os.path.exists(already_done_tracking):
                with open(already_done_tracking, 'r') as f:
                    already_done_list = [line.rstrip('\n') for line in f]
            scvcfs_to_process = list(set(cell_final_outs))
            if len(already_done_list) > 0:
                scvcfs_to_process = list(set(cell_final_outs) - set(already_done_list))
            scvcfs_to_process.sort()
            tmp_filelist_txt = f'{cell_dir}TEMP_FILELIST.txt'
            while len(scvcfs_to_process) > 0:
                logger.info(f'Files left to process: {len(scvcfs_to_process)}')
                sub_scvcfs = scvcfs_to_process[:n_vcfs_per_merge]
                del scvcfs_to_process[:n_vcfs_per_merge]
                with open(tmp_filelist_txt, 'w') as f:
                    if os.path.exists(tmp1):
                        f.write(tmp1 + '\n')
                    for line in sub_scvcfs:
                        f.write(line + '\n')
                cmd = f'bcftools merge --threads {n_workers} --file-list {tmp_filelist_txt} -Oz -o {tmp1new} && mv {tmp1new} {tmp1} \
                && rm {tmp1}.csi && bcftools index {tmp1} && rm {tmp_filelist_txt}'
                logger.info(cmd)
                proc = subprocess.run(cmd, shell=True, capture_output=True)
                with open(already_done_tracking, 'a') as ft:
                    for line in sub_scvcfs:
                        ft.write(line + '\n')
            if not os.path.exists(f'{tmp1}.csi'):
                cmd = f'bcftools index {tmp1}'
                proc = subprocess.run(cmd, shell=True, capture_output=True)
            # tmp2 = f'{cell_dir}/temp02_merged.vcf.gz'
            # if not os.path.exists(tmp2):
            #     cmd = f'bcftools norm --threads {n_workers} --atomize --multiallelics - {tmp1} -Oz -o {tmp2}'
            #     logger.info(cmd)
            #     proc = subprocess.run(cmd, shell=True, capture_output=True)
            # cmd = f'bcftools +fill-tags {tmp2} -- -t all -Oz -o {final_vcf_out} --threads {n_workers}'
            # logger.info(cmd)
            # proc = subprocess.run(cmd, shell=True, capture_output=True)
            # os.remove(tmp1)
            # os.remove(tmp2)
            # os.remove(already_done_tracking)
            pp2_sc_longshot_processmerged = functools.partial(
                _sc_longshot_processmerged,
                sample_split = False,
                tmp1 = tmp1,
                tmp_file_prefix = 'temp02b_merged',
                cell_dir = cell_dir,
                data_name = data_name,
                work_dir = work_dir,
            )
            with mp.Pool(n_workers) as pool:
                tmp2_vcf_outs = pool.map(pp2_sc_longshot_processmerged, chromosomes)
    
            tmpbulk_filelist_txt = f'{cell_dir}TMPMERGEDCHR_2_FILELIST.txt'
            with open(tmpbulk_filelist_txt, 'w') as f:
                for line in tmp2_vcf_outs:
                    f.write(line + '\n')
            cmd = f'bcftools concat --file-list {tmpbulk_filelist_txt} --threads {threads} \
                --output-type z -o {final_vcf_out}'
            logger.info(cmd)
            proc = subprocess.run(cmd, shell=True, capture_output=True)

            if os.path.exists(final_vcf_out):
                os.remove(tmp1)
                for tmp2fx in tmp2_vcf_outs:
                    if os.path.exists(tmp2fx):
                        os.remove(tmp2fx)
                    if os.path.exists(f'{tmp2fx}.csi'):
                        os.remove(f'{tmp2fx}.csi')
                os.remove(already_done_tracking)
                os.remove(tmpbulk_filelist_txt)
        
    return files_to_delete, final_vcf_out

## ~~~ ##

def snv_call(
    ## Inputs
    dedup_bam: str | Path | list,
    method: str = 'longshot',

    ## General utlity
    work_dir: str | Path | None = None,
    cleanup: bool = True,
    **kwargs
):
    r"""Wrapper function for SNV calling.
    
    Currently only supports the "longshot" protocol from EPI2ME.

    Parameters
    ----------
    dedup_bam
        String, Path, or list of these indicating the location of BAM files with the "DP" tag for duplicated reads. 
        If a list provided, will loop through each bam file and output to said file's directory.
    method
        String indicating method for SNV calling. Currently only supports "longshot".
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    **kwargs
        Additional keyword arguments passed to SNV calling method of choice.
    """
    ## Prelim checks and basic utility
    logger, work_dir = setup_workdir(work_dir)
    files_to_delete = []
    if not isinstance(dedup_bam, list):
        dedup_bam = [dedup_bam]
    supported_methods = ['longshot']
    if method.lower() not in supported_methods:
        logger.error(
            f"{method} is not a supported method. \
            Currently acceptable values are: {', '.join(supported_methods)}"
        )
        sys.exit(1)

    ## Pop out kwargs users are never allowed to change
    kwargs.pop('dbam', None)
    kwargs.pop('work_dir', None) #This actually can't happen, but put here to be safe

    ## Loop through all the BAMs to process
    final_vcfs = []
    for dbam in dedup_bam:
        if method.lower() == 'longshot':
            files_to_delete_x, final_vcf = sc_longshot(dbam=dbam, work_dir=work_dir, **kwargs)
            files_to_delete.extend(files_to_delete_x)
        final_vcfs.append(final_vcf)

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                if not os.path.isdir(dfile_path):
                    os.remove(dfile_path)
    return final_vcfs


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################################
## ----- For integrating with EPI2ME ----- ##
#############################################

def get_missing_reads(
    fastq_file: str | Path,
    reads_to_exclude: list | set,
    work_dir: str | Path | None = None,
):
    r"""Get reads that are missing from a list.

    Returns a FASTQ file 

    Parameters
    ----------
    fastq_files
        A string or Path indicating the FASTQ file to examine.
    reads_to_exclude
        A list or set indicating which reads to ignore in the FASTQ.
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`. THis will also be where the final results will be stored.
    """
    ## Basic utility and checks
    logger, work_dir = setup_workdir(work_dir)
    if not isinstance(reads_to_exclude, set):
        reads_to_exclude = set(reads_to_exclude)
    if not os.path.exists(fastq_file):
        logger.error(f'{fastq_file} not found')
        sys.exit(1)
    file_name = re.sub('.fastq|.gz|.fq', '', os.path.basename(os.path.normpath(fastq_file)))
    out_file = f'{work_dir}/{file_name}_filtered.fastq'
    logger.info(f'Writing missing reads to {out_file}')
    
    ## Loop through the FASTQ
    
    ## Assume this is faster
    p = subprocess.Popen(['seqkit', 'fx2tab', fastq_file], text=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        read_id = line.split('\t', 1)[0].split(' ', 1)[0]
        if read_id in reads_to_exclude:
            continue
        read_info = line.strip('\n').split('\t')
        if len(read_info)!=3:
            logger.error(f'Incorrect parsing of {read_id}: expecting 3 fields to be returned')
            sys.exit(1)
        if len(read_info[1])!=len(read_info[2]):
            logger.error(f'Incorrect parsing of {read_id}: expecting sequence to be as long as QC info')
        read_out = f'@{read_info[0]}\n{read_info[1]}\n+\n{read_info[2]}\n'
        with open(out_file, 'a') as fout:
            fout.write(read_out)
    p.wait()

    ## Alternative: pysam approach
    # with pysam.FastxFile(fastq_file, persist = False) as fastq_fh, open(out_file, 'w') as fq_out:
    #     for read in fastq_fh:
    #         if read.name in reads_to_exclude:
    #             continue
    #         print(f'Keeping {read.name}')
    #         fq_out.write(f'@{read.name}\n{read.sequence}\n+\n{read.quality}\n')
    return out_file

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
