## Actual pipeline here
import multiprocessing
import functools
import shutil

from .pipeline_functions import *

def run_after_epi2me(
    # Output of EPI2ME folder location
    epi2me_output_dir: str | Path,
    
    # Global params: user-specified
    fastq_files: str | Path | list,
    genome_dir: str | Path,
    curio_whitelist: str | Path,
    sticr_whitelist_dir: str | Path,

    ## Global params: defaults given
    cellbc_whitelist_dir: str | Path | None = None,
    sample_name: str = "po3t",
    n_workers: int | None = None,

    ## Step 1: Getting genome information 
    ## Move genome_dir up into Global params
    genome_fasta: str | Path | None = None,
    genome_gtf: str | Path | None = None,
    
    ## Step 2: extracting 10x barcodes
    whitelist_dir: str | Path | None = None,
    kit: str = '3prime:v4',

    ## Step 3: Sorting reads into groups
    groups: dict = {
        'STICR': ['AGCAAACTGGGGCACAAGC', 'GAATTATCGAATT', 'GAATTGAGGAATT'],
        'CURIO': ['TCTTCAGCGTTCCCGAGA', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT']
    },
    default_group: str = 'RNA',
    all_matches_mode: bool = False,

    ## Step 4A: CURIO parsing
    ## Move curio_whitelist up into Global params
    ## No other params needed

    ## Step 4B: STICR parsing
    ## Move sticr_whitelist_dir up into Global params
    indices_to_consider: str | list = ['1', '3', 'e'],
    sticr_whitelist: dict = {
        'Bit1': 'Bit1.F.with.Empty.fasta',
        'Bit2': 'Bit2.F.with.Empty.fasta',
        'Bit3': 'Bit3.F.with.Empty.fasta'
    },

    ## Step 8: SNV calling
    snv_calling_method: str = 'longshot',

    ## General utility
    work_dir: str | Path | None = None,
    cleanup: bool = True,
    **kwargs,
):
    r"""Pipeline for processing FASTQ files after running the EPI2ME pipeline.

    The EPI2ME wf-single-cell is a fast method for mapping long read FASTQ files. After running, the outputs should be in a subdirectory in the `/output/` folder:
    e.g. `/output/sample_X/`. We will use the `.tagged.bam` for variant calling, and `.read_summary.tsv` for cross-referencing to sequences in `fastq_files`.
    Thus, to avoid ambiguity, users should specify the `sample_X` subdirectory and matching FASTQ files in `fastq_files`.
    
    Parameters
    ----------
    fastq_files
        A string, Path, or list of these indicating FASTQs belonging to the same biological sample.
    genome_dir
        String or Path indicating where folder with genome information is.
    curio_whitelist
        String or Path (or list of these) indicating the location of `.txt` files with bead barcode coordinate information (whitelist).
        This is provided by the CurioTrekker portal (https://knowledgebase.curiobioscience.com/bioinformatics/tile-barcode/) 
        as a `.zip` upon keying in the relevant slide ID used. Comprises three columns: N14 barcode, x coordinate, y coordinate.
        If more than one provided, the union of all N14 barcodes are used as the whitelist.
    sticr_whitelist_dir
        String or Path indicating directory with Bit1-3 files. If `None` (default), assumes absolute paths was provided in the `sticr_whitelist` dictionary.

    cellbc_whitelist_dir
        String of Path indicating directory with 10X barcodes whitelist. If unspecified (default) will download relevant file from online. 
    sample_name
        String indicating name of sample these datasets belong to. This will be used in prefixes of merged files.
    n_workers
        Integer indicating number of workers per `multiprocessing.Pool`. If not specified, defaults to `os.cpu_count()`.
        
    genome_fasta
        String or Path indicating which file is our genome fasta file. If not specified, inferred from `genome_dir`. Used in `get_genome()`.
    genome_gtf
        String or Path indicating which file is our gene annotation file. If not specified, inferred from `genome_dir`. Used in `get_genome()`.
        
    whitelist_dir
        String or Path indicating where 10X cellular barcode whitelist directory is. If unspecified, downloads the inferred whitelist into `work_dir`.
        Thus, default `whitelist_dir` is `'~/WORK/'`. Alternatively, if a file is supplied, will assume that the file is the whitelist file itself. 
        Used in `extract_10x_barcodes()`.  
    kit
        String indicating what 10X kit was used. Used in `extract_10x_barcodes()`.
        
    groups
        Dictionary of lists indicating what _forward_ sequence(s) are needed for a read to be assigned to any given group. 
        Name of the group as keys. Lists must be mutually exclusive. 
        "PROBLEM" cannot be a key in this dictionary -- this designation is used internally for problem reads.
        Used in `group_reads()`. If "STICR" and/or "CURIO" are missing from this dictionary, will not run the matching step 4.
    default_group
        String indicating name of group for all reads that do not get assigned to any of the groups in `groups`. Used in `group_reads()`.
    all_matches_mode
        Boolean indicating whether -- if multiple sequences per group are provided -- to only assign it to a group if _all_ sequences in the list appear in the read.
        If False, then as long as _any_ of the sequences per group appears, then the read gets assigned to that group. Used in `group_reads()`.

    indices_to_consider
        String or list of strings indicating which indices to consider during mapping. Valid values are '1', '3', and/or 'e'. Used in `process_sticr()`. 
    sticr_whitelist
        Dictionary indicating filenames of Bit1, Bit2, and Bit3 whitelist files. Used in `process_sticr()`. 

    snv_calling_method
        String indicating method for SNV calling. Currently only supports "longshot". Used in `snv_call()`.
        
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`. THis will also be where the final results will be stored.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    **kwargs
        Other keyword arguments passed to pipeline functions. 
    """
    ## Basic utility and checks
    logger, work_dir = setup_workdir(work_dir)
    if not os.path.exists(epi2me_output_dir):
        logger.error(f'{epi2me_output_dir} not found')
        sys.exit(1)
    epi2me_output_dir = os.path.normpath(epi2me_output_dir)
    sample_name = os.path.basename(os.path.normpath(epi2me_output_dir))
    logger.info(f'Inferred sample name to be {sample_name}')
    tagged_bam_file = f'{epi2me_output_dir}/{sample_name}.tagged.bam'
    tagged_bai_file = f'{tagged_bam_file}.bai'
    if not os.path.exists(tagged_bam_file):
        logger.error(f'{tagged_bam_file} does not exist: ensure you hae specified the correct sample subdirectory')
        sys.exit(1)
    if not isinstance(fastq_files, list):
        fastq_files = [fastq_files]
    for fastq_file in fastq_files:
        any_files_not_found = False
        if not os.path.exists(fastq_file):
            logger.error(f'{fastq_file} not found: ensure absolute paths are specified and spelling is correct')
            any_files_not_found = True
        if any_files_not_found:
            sys.exit(1)
    if n_workers is None:
        n_workers = os.cpu_count()

    ## Get annotated reads
    file_check = f'{work_dir}/step02_extract_10x_barcodes.txt'
    if not os.path.exists(file_check):
        read_summary_tsv_file = f'{epi2me_output_dir}/{sample_name}.read_summary.tsv'
        if not os.path.exists(read_summary_tsv_file):
            logger.warning(f'{read_summary_tsv_file} does not exist...Defaulting to getting read IDs from the BAM instead...')
            annotated_reads = subprocess.run(f'samtools view {tagged_bam_file} | cut -f 1', capture_output=True, shell=True).stdout.decode().split('\n')
        else:
            annotated_reads = pl.read_csv(f'{read_summary_tsv_file}', columns=['read_id'], separator='\t')['read_id'].to_list()
        annotated_reads = set(annotated_reads)
        while '' in annotated_reads:
            annotated_reads.remove('')
        logger.info(f'{len(annotated_reads)} annotated sub reads found')
        pattern = re.compile(r'_\d+$')
        parent_reads = [pattern.sub('', x) for x in annotated_reads]
        parent_reads = set(parent_reads)
        logger.info(f'{len(parent_reads)} annotated parent reads found')
 
    ## There is a problem with polars working across multiprocessing pools.
    ## Note that this is a non-issue for Python 3.14, which sets default multiprocessing behaviour to "spawn".
    mp = multiprocessing.get_context("spawn")

    #####################################
    ## ----- Core Pipeline START ----- ##
    #####################################
    
    ## Step 1: get genome, only needs to be done once
    file_check = f'{work_dir}/step01_get_genome.txt'
    if not os.path.exists(file_check):
        logger.info(f'Prepping genome information')
        get_genome(genome_dir=genome_dir, genome_fasta=genome_fasta, genome_gtf=genome_gtf, work_dir=work_dir, **kwargs)
        Path(file_check).touch()

    ## Step 1.5: Process fastq files so that they do not have processed reads
    file_check = f'{work_dir}/step01A_get_missing_reads.txt'
    if not os.path.exists(file_check):
        logger.info(f'Filtering out reads already processed')
        p_get_missing_reads = functools.partial(
            get_missing_reads,
            reads_to_exclude = parent_reads,
            work_dir = work_dir,
        )
        with mp.Pool(n_workers) as pool:
            chunk_files = pool.map(p_get_missing_reads, fastq_files)
        logger.info(f'{len(chunk_files)} filtered fastqs generated.')
        with open(file_check, 'w') as fclog:
            for line in chunk_files:
                fclog.write(line + '\n')
    else:
        with open(file_check, 'r') as fclog:
            chunk_files = [line.strip() for line in fclog]
            
    ## Step 2: Barcode extraction
    file_check = f'{work_dir}/step02_extract_10x_barcodes.txt'
    if not os.path.exists(file_check):
        if cellbc_whitelist_dir is None:
            cellbc_whitelist_dir = f'{work_dir}'
        out_whitelist = get_10x_whitelist(
            kit = kit,
            whitelist_dir = cellbc_whitelist_dir,
            work_dir = work_dir
        )
        logger.info(f'Obtained 10X whitelist: {out_whitelist}')
        
        p_extract_10x_barcodes = functools.partial(
            extract_10x_barcodes,
            kit = kit,
            rapidfuzz_error_correct = True,
            work_dir = work_dir,
            cleanup = cleanup,
            whitelist_dir = cellbc_whitelist_dir,
            **kwargs
        )
        with mp.Pool(n_workers) as pool:
            trimmed_fastqs = pool.map(p_extract_10x_barcodes, chunk_files)
        logger.info(f'{len(trimmed_fastqs)} trimmed fastqs generated.')
        with open(file_check, 'w') as fclog:
            for line in trimmed_fastqs:
                fclog.write(line + '\n')
    else:
        with open(file_check, 'r') as fclog:
            trimmed_fastqs = [line.strip() for line in fclog]

    if len(chunk_files) >0:
            for cfx in chunk_files:
                if os.path.exists(cfx):
                    os.remove(cfx)
    
    ## Step 3: Sort reads into groups
    file_check = f'{work_dir}/step03_group_reads.txt'
    if not os.path.exists(file_check):
        p_group_reads = functools.partial(
            group_reads,
            work_dir = work_dir,
            cleanup = cleanup,
            **kwargs
        )
        with mp.Pool(n_workers) as pool:
            chunked_dirs = pool.map(p_group_reads, trimmed_fastqs)
        logger.info(f'{len(chunked_dirs)} directories with sorted reads.')
        with open(file_check, 'w') as fclog:
            for line in chunked_dirs:
                fclog.write(line + '\n')
    else:
        with open(file_check, 'r') as fclog:
            chunked_dirs = [line.strip() for line in fclog]

    ## Step 4A: CURIO processing
    file_check = f'{work_dir}/step04A_process_curio.txt'
    if not os.path.exists(file_check):
        curio_files = []
        do_curio = True
        if "CURIO" not in groups.keys():
            logger.warning(f"No CURIO in groups.keys(): Assuming do not have to parse CURIO barcodes")
            do_curio = False
        else:
            for chunkx in chunked_dirs:
                fx = glob.glob(f'{chunkx}/CURIO_*.fastq')
                if len(fx)==1:
                    curio_files.extend(fx)
        if len(curio_files)==0:
            logger.warning(f"No CURIO_*.fastq files found: Assuming do not have to parse CURIO barcodes")
            do_curio = False
        if do_curio:
            p_process_curio = functools.partial(
                process_curio,
                curio_whitelist = curio_whitelist,
                work_dir = work_dir,
                cleanup = False,
                **kwargs
            )
            with mp.Pool(n_workers) as pool:
                pool.map(p_process_curio, curio_files) #No output
            logger.info(f'{len(curio_files)} CURIO FASTQs processed.')
        Path(file_check).touch()
            

    ## Step 4B: STICR processing
    file_check = f'{work_dir}/step04B_process_sticr.txt'
    if not os.path.exists(file_check):
        sticr_files = []
        do_sticr = True
        if "STICR" not in groups.keys():
            logger.warning(f"No STICR in groups.keys(): Assuming do not have to parse STICR barcodes")
            do_sticr = False
        else:
            for chunkx in chunked_dirs:
                fx = glob.glob(f'{chunkx}/STICR_*.fastq')
                if len(fx)==1:
                    sticr_files.extend(fx)
        if len(sticr_files)==0:
            logger.warning(f"No STICR_*.fastq files found: Assuming do not have to parse STICR barcodes")
            do_sticr = False
        if do_sticr:
            p_process_sticr = functools.partial(
                process_sticr,
                sticr_whitelist_dir = sticr_whitelist_dir,
                indices_to_consider = indices_to_consider,
                sticr_whitelist = sticr_whitelist,
                work_dir = work_dir,
                cleanup = False,
                **kwargs
            )
            with mp.Pool(n_workers) as pool:
                pool.map(p_process_sticr, sticr_files) #No output
            logger.info(f'{len(sticr_files)} STICR FASTQs processed.')
        Path(file_check).touch()

    ## Consolidate barcode parsing -- moved into the pipeline because SNV calling takes a while
    logger.info(f'Consolidating barcode calls')
    final_output_dir = f'{work_dir}/OUTPUT/'
    if not os.path.exists(final_output_dir):
        os.makedirs(final_output_dir, exist_ok=True)
    parsed_barcode_files = subprocess.run(['find', work_dir, '-type', 'f', '-name', "*_ASSIGNS.tsv"], capture_output=True, text=True).stdout.split('\n')
    parsed_barcode_files.remove('')
    group_outs = dict.fromkeys(groups)
    for assign_filex in parsed_barcode_files:
        for group_namey in list(group_outs.keys()):
            if group_namey in assign_filex:
                assign_tagx = re.sub('_ASSIGNS.tsv$', '_TAGS.tsv', assign_filex)
                if not os.path.exists(assign_tagx):
                    sys.exit(1)
                assign_df = pl.read_csv(assign_filex, separator='\t')
                if assign_df.shape[0]==0:
                    continue
                assign_tag = pl.read_csv(assign_tagx, separator='\t')
                assign_tag = assign_tag.filter(
                    pl.col('read_id').is_in(assign_df['read_id'].implode())
                )
                assign_final = assign_df.join(assign_tag, on='read_id', how='inner', suffix='_cellbc')
                if group_outs[group_namey] is None:
                    group_outs[group_namey] = [assign_final]
                else:
                    group_outs[group_namey].append(assign_final)
    for key,value in group_outs.items():
        group_outx = ''
        if value is None:
            logger.warning(f'No assigned reads for {key}')
            continue
        group_outx = pl.concat(value)    
        group_outx.write_csv(f'{final_output_dir}{sample_name}_{key}.csv')

    ## Jump ahead to SNV calling            
    ## Step 8: SNV calling (renamed to step 05)
    file_check = f'{work_dir}/step05_snv_call.txt'
    if not os.path.exists(file_check):
        ## Do some prep for file names
        dedup_bam_files = f'{work_dir}/{sample_name}_annotated_deduplicated.bam'
        dedup_bai_files = f'{dedup_bam_files}.bai'
        if not os.path.exists(dedup_bam_files):
            logger.info(f'Copying {tagged_bam_file} to {dedup_bam_files}')
            subprocess.run(['cp', tagged_bam_file, dedup_bam_files], check=True)
        if os.path.exists(tagged_bai_file):
            if not os.path.exists(dedup_bai_files):
                logger.info(f'Copying {tagged_bai_file} to {dedup_bai_files}')
                subprocess.run(['cp', tagged_bai_file, dedup_bai_files], check=True)
        else:
            logger.warning(f'No matching .bai file found for {tagged_bam_file}')
        dedup_bam_files = [dedup_bam_files]

        ## Now do the SNV calling -- note that this duplicates the ~40GB BAM file
        ## Ensure sufficient memory available
        logger.info(f'{len(dedup_bam_files)} BAM files to be used for SNV calling')
        
        ## Do some file cap checks before running
        ## First, check whether an actual limit exists
        actual_limit = int(subprocess.run(f'ulimit -n', shell=True, capture_output=True, text=True).stdout)
        logger.info(f'Current max limit of {actual_limit} files allowed to be open')

        ## Next check if there is a file_cap already user-specified, and check that it is <= actual_limit
        file_cap = kwargs.pop('file_cap', None)
        if file_cap is None:
            file_cap = actual_limit
        if file_cap > actual_limit:
            logger.error(f'Proposed file cap {file_cap} exceeds actual ulimit: recommend to edit /etc/security/limits.conf with sudo permissions to raise limits')
            sys.exit(1)

        ## Finally check that ulimit is changeable
        cmd = f'ulimit -n {file_cap}'
        logger.info(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True, check=True)
        final_limit = int(subprocess.run(f'ulimit -n', shell=True, capture_output=True, text=True).stdout)
        logger.info(f'Current max limit of {final_limit} files allowed to be open')
        
        vcf_files = snv_call(
            dedup_bam = dedup_bam_files,
            method = snv_calling_method,
            work_dir = work_dir,
            cleanup = cleanup,
            file_cap = file_cap,
            **kwargs
        )
        logger.info(f'{len(vcf_files)} VCF files created')
        Path(file_check).touch()

        logger.info(f'Deleting the {len(dedup_bam_files)} copied .bam files')
        for dedup_bam_filex in dedup_bam_files:
            if os.path.exists(dedup_bam_filex):
                os.remove(dedup_bam_filex)

    ## Step 9: Consolidate VCFs into OUTPUT folder
    file_check = f'{work_dir}/step06_snv_consolidated.txt'
    if not os.path.exists(file_check):
        logger.info(f'Copying the most important VCFs into OUTPUT/')
        vcfs_to_transfer = glob.glob( f'{work_dir}/*/*finalsnv_cells.vcf.gz')
        vcfs_to_transfer.extend( glob.glob( f'{work_dir}/*/*allcandidates.vcf.gz') )
        vcfs_to_transfer.extend( glob.glob( f'{work_dir}/*/*allbulks.vcf.gz') )
        for src in vcfs_to_transfer:
            try:
                shutil.copy2(src, final_output_dir)
                logger.info(f'Copied {src}')
            except shutil.SameFileError:
                logger.warning(f'{src} is unchanged...skipping...')
            except Exception as e:
                logger.error(f'Unexpected error: {e}')
                sys.exit(1)
        
        vcfs_to_transfer = glob.glob( f'{work_dir}/*/*finalgenotypes.vcf.gz')
        if not os.path.exists(f'{final_output_dir}/CELL_SNVS/'):
            os.makedirs(f'{final_output_dir}/CELL_SNVS/', exist_ok=True)
        for src in vcfs_to_transfer:
            try:
                shutil.copy2(src, f'{final_output_dir}/CELL_SNVS/')
                logger.info(f'Copied {src}')
            except shutil.SameFileError:
                logger.warning(f'{src} is unchanged...skipping...')
            except Exception as e:
                logger.error(f'Unexpected error: {e}')
                sys.exit(1)
        Path(file_check).touch()

    ###################################
    ## ----- Core Pipeline END ----- ##
    ###################################

    return

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Deprecated function
def run_pipeline_old(
    ## Global params: user-specified
    fastq_files: str | Path | list,
    genome_dir: str | Path,
    curio_whitelist: str | Path,
    sticr_whitelist_dir: str | Path,

    ## Global params: defaults given
    cellbc_whitelist_dir: str | Path | None = None,
    sample_name: str = "po3t",
    chunk_size: int = 100 * 1000,
    n_workers: int | None = None,

    ## Step 1: Getting genome information 
    ## Move genome_dir up into Global params
    genome_fasta: str | Path | None = None,
    genome_gtf: str | Path | None = None,
    
    ## Step 2: extracting 10x barcodes
    whitelist_dir: str | Path | None = None,
    kit: str = '3prime:v4',

    ## Step 3: Sorting reads into groups
    groups: dict = {
        'STICR': ['AGCAAACTGGGGCACAAGC', 'GAATTATCGAATT', 'GAATTGAGGAATT'],
        'CURIO': ['TCTTCAGCGTTCCCGAGA', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT']
    },
    default_group: str = 'RNA',
    all_matches_mode: bool = False,

    ## Step 4A: CURIO parsing
    ## Move curio_whitelist up into Global params
    ## No other params needed

    ## Step 4B: STICR parsing
    ## Move sticr_whitelist_dir up into Global params
    indices_to_consider: str | list = ['1', '3', 'e'],
    sticr_whitelist: dict = {
        'Bit1': 'Bit1.F.with.Empty.fasta',
        'Bit2': 'Bit2.F.with.Empty.fasta',
        'Bit3': 'Bit3.F.with.Empty.fasta'
    },

    ## Step 5: RNA mapping
    ## No other params needed

    ## Step 6: Feature assignment
    merge_bams: bool = True,
    feature_assign_method: str = 'stringtie',

    ## Step 7: UMI de-duplication
    keep_all_reads: bool = True,

    ## Step 8: SNV calling
    snv_calling_method: str = 'longshot',

    ## General utility
    work_dir: str | Path | None = None,
    cleanup: bool = True,
    **kwargs,
):
    r"""Deprecated multi-processing pipeline for processing FASTQ reads. 

    WARNING: SLOW. This is likely because `` Recommend to instead run EPI2ME's pipeline first and then run the `run_pipeline` function instead. 

    Chunks each FASTQ file in `fastq_files`, and processes them with our pipeline. Merging of RNA reads from all chunks for all FASTQ files specified then performed, 
    and used as input for `map_rna()`. UMI de-duplication, SNV calling etc will thus be done on this merged BAM -- it is crucial all data comes from the same biological sample.

    At the end, the most important outputs are:
    * `geneexpr.csv`: CSV of gene expression values (number of UMIs) per gene per cell. In long format: users will need to convert this into a matrix of genes x cells for most scRNA-seq workflows.
    * `sticr.csv`: CSV of information for all _reads_ assigned a STICR barcode. To get UMIs, filter according to the "duplicated" column. 
    * `curio.csv`: CSV of information for all _reads_ assigned a CURIO barcode. To get UMIs, filter according to the "duplicated" column. 
    * `annotated_deduplicated.bam`: BAM of final alignment information, e.g. "CB" (cellular barcode), "GN" (gene assigned), "TR" (transcript assigned), "DP" (duplicated UMI if 1, unique if 0).
    * `finalsnv.vcf.gz`: VCF of single-cell SNV calls.
    
    Other useful outputs:
    * `alltranscripts_stringtie.gff`: GFF of transcript (TR) coordinates with respect to the user-specified genome. Useful for cross-referencing to `annotated_deduplicated.bam`. 
    * `genes.csv`: CSV of all possible gene names in user-specified genome. Useful for cross-referencing to `geneexpr.csv`.
    * `log.log`: Pipeline's log.
    
    Does not offer a lot of flexibility: for users wishing to run their own pipeline parameters beyond what `**kwargs` can specify,
    copy this function and tweak as needed.

    Function primarily works within `work_dir` (defaults to `~/WORK`). If looping through samples, recommend to give unique `work_dir` values for each sample.
    
    Parameters
    ----------
    fastq_files
        A string, Path, or list of these indicating FASTQs belonging to the same biological sample.
    genome_dir
        String or Path indicating where folder with genome information is.
    curio_whitelist
        String or Path (or list of these) indicating the location of `.txt` files with bead barcode coordinate information (whitelist).
        This is provided by the CurioTrekker portal (https://knowledgebase.curiobioscience.com/bioinformatics/tile-barcode/) 
        as a `.zip` upon keying in the relevant slide ID used. Comprises three columns: N14 barcode, x coordinate, y coordinate.
        If more than one provided, the union of all N14 barcodes are used as the whitelist.
    sticr_whitelist_dir
        String or Path indicating directory with Bit1-3 files. If `None` (default), assumes absolute paths was provided in the `sticr_whitelist` dictionary.

    cellbc_whitelist_dir
        String of Path indicating directory with 10X barcodes whitelist. If unspecified (default) will download relevant file from online. 
    sample_name
        String indicating name of sample these datasets belong to. This will be used in prefixes of merged files.
    chunk_size
        Integer indicating number of reads per chunk.
    n_workers
        Integer indicating number of workers per `multiprocessing.Pool`. If not specified, defaults to `os.cpu_count()`.
        
    genome_fasta
        String or Path indicating which file is our genome fasta file. If not specified, inferred from `genome_dir`. Used in `get_genome()`.
    genome_gtf
        String or Path indicating which file is our gene annotation file. If not specified, inferred from `genome_dir`. Used in `get_genome()`.
        
    whitelist_dir
        String or Path indicating where 10X cellular barcode whitelist directory is. If unspecified, downloads the inferred whitelist into `work_dir`.
        Thus, default `whitelist_dir` is `'~/WORK/'`. Used in `extract_10x_barcodes()`.  
    kit
        String indicating what 10X kit was used. Used in `extract_10x_barcodes()`.
        
    groups
        Dictionary of lists indicating what _forward_ sequence(s) are needed for a read to be assigned to any given group. 
        Name of the group as keys. Lists must be mutually exclusive. 
        "PROBLEM" cannot be a key in this dictionary -- this designation is used internally for problem reads.
        Used in `group_reads()`. If "STICR" and/or "CURIO" are missing from this dictionary, will not run the matching step 4.
    default_group
        String indicating name of group for all reads that do not get assigned to any of the groups in `groups`. Used in `group_reads()`.
    all_matches_mode
        Boolean indicating whether -- if multiple sequences per group are provided -- to only assign it to a group if _all_ sequences in the list appear in the read.
        If False, then as long as _any_ of the sequences per group appears, then the read gets assigned to that group. Used in `group_reads()`.

    indices_to_consider
        String or list of strings indicating which indices to consider during mapping. Valid values are '1', '3', and/or 'e'. Used in `process_sticr()`. 
    sticr_whitelist
        Dictionary indicating filenames of Bit1, Bit2, and Bit3 whitelist files. Used in `process_sticr()`. 

    merge_bams
        Boolean indicating if bams are to be merged into one big merged `.bam` before feature calling. Used in `feature_assign()`.
    feature_assign_method
        String indicating method for calling features per read. Currently only supports "stringtie". Used in `feature_assign()`.

    keep_all_reads
        Boolean indicating if all reads should be kept. 
        If True, then reads to be filtered can be identified by the DP:1 tag.
        If False, filtered reads are already removed and only DP:0 tags exist.
        Used in `deduplicate_umi()`.

    snv_calling_method
        String indicating method for SNV calling. Currently only supports "longshot". Used in `snv_call()`.
        
    work_dir
        String or Path indicating where working directory is. Defaults to `'~/WORK'`. THis will also be where the final results will be stored.
    cleanup
        Boolean indicating if want to delete unnecessary files at the end.
    **kwargs
        Other keyword arguments passed to pipeline functions. 
    """
    ## Basic utility and checks
    logger, work_dir = setup_workdir(work_dir)
    logger.warning(f'This run_pipeline_slow() is a deprecated function: pipeline being run is very slow. Recommend to run EPI2ME wf-single-cell first then use the new run_pipeline().')
    if not isinstance(fastq_files, list):
        fastq_files = [fastq_files]
    for fastq_file in fastq_files:
        any_files_not_found = False
        if not os.path.exists(fastq_file):
            logger.error(f'{fastq_file} not found: ensure absolute paths are specified and spelling is correct')
            any_files_not_found = True
        if any_files_not_found:
            sys.exit(1)
    if n_workers is None:
        n_workers = os.cpu_count()
    

    ## Pop out invalid keywords
    invalid_kwargs = ['method', 'merged_bam_name']
    for invalid_kwarg in invalid_kwargs:
        popped_out = kwargs.pop(invalid_kwarg, None)
        if popped_out:
            logger.warning(f'Invalid keyword detected: removing {invalid_kwarg} from **kwargs')

    ## There is a problem with polars working across multiprocessing pools.
    ## Note that this is a non-issue for Python 3.14, which sets default multiprocessing behaviour to "spawn".
    mp = multiprocessing.get_context("spawn")

    
    #####################################
    ## ----- Core Pipeline START ----- ##
    #####################################
    
    ## Step 1: get genome, only needs to be done once
    logger.info(f'Prepping genome information')
    file_check = f'{work_dir}/step01_get_genome.txt'
    if not os.path.exists(file_check):
        get_genome(genome_dir=genome_dir, genome_fasta=genome_fasta, genome_gtf=genome_gtf, work_dir=work_dir, **kwargs)
        Path(file_check).touch()
    

    ## Loop through each fastq file
    for fastq_file in fastq_files:

        data_name = re.sub('.fastq$|.fq$|fastq.gz$|fq.gz$', '', Path(fastq_file).stem)

        ## For each fastq file, chunk so that each file has max chunk_size number of reads
        ## Each read is 4 lines in a fastq file, so split size  = chunk_size * 4
        chunk_files = glob.glob(f'{work_dir}/{data_name}_chunk_*.fastq')
        if len(chunk_files) == 0 and not os.path.exists(f'{work_dir}/step02_extract_10x_barcodes_{data_name}.txt'):
            n_lines = int(chunk_size) * 4
            logger.info(f'Chunking {fastq_file} into {n_lines}-line chunks')
            chunk_prefix = f'{data_name}_chunk_'
            if re.search('.gz$', fastq_file):
                proc = subprocess.run(
                f'zcat {fastq_file} | split -l {n_lines} - {work_dir}/{chunk_prefix} --additional-suffix .fastq', 
                shell=True, capture_output=True)
            else:
                proc = subprocess.run(
                f'split -l {n_lines} {fastq_file} {work_dir}/{chunk_prefix} --additional-suffix .fastq', 
                shell=True, capture_output=True)
            chunk_files = glob.glob(f'{work_dir}/{chunk_prefix}*.fastq')

        ## After chunking, we will process each chunk separately until it comes time to merge RNA
        ## To allow users to specify new kwargs, we will have to create a bunch of partial functions to update our default arguments
        ## That way, we only have to care about the first positional argument
        ## Since these base functions are already defined at the top level of module, should not be a problem to create partials _inside_ this function
        logger.info(f'Parallel processing of {len(chunk_files)} chunks for {fastq_file}')

        ## Step 2: Barcode extraction
        ## For some very odd reason, this tends to stall. But when restarting, there's no issue.
        file_check = f'{work_dir}/step02A_extract_10x_barcodes_{data_name}.txt'
        if not os.path.exists(file_check):
            if cellbc_whitelist_dir is None:
                cellbc_whitelist_dir = f'{work_dir}'
            out_whitelist = get_10x_whitelist(
                kit = kit,
                whitelist_dir = cellbc_whitelist_dir,
                work_dir = work_dir
            )
            logger.info(f'Obtained 10X whitelist: {out_whitelist}')
            
            p_extract_10x_barcodes = functools.partial(
                extract_10x_barcodes,
                kit = kit,
                rapidfuzz_error_correct = False,
                work_dir = work_dir,
                cleanup = cleanup,
                whitelist_dir = cellbc_whitelist_dir,
                **kwargs
            )
            with mp.Pool(n_workers) as pool:
                trimmed_fastqs = pool.map(p_extract_10x_barcodes, chunk_files)
            logger.info(f'{len(trimmed_fastqs)} trimmed fastqs generated.')
            with open(file_check, 'w') as fclog:
                for line in trimmed_fastqs:
                    fclog.write(line + '\n')
        else:
            with open(file_check, 'r') as fclog:
                trimmed_fastqs = [line.strip() for line in fclog]
        
        ## Delete the chunk files after processing
        if len(chunk_files) >0:
            for cfx in chunk_files:
                if os.path.exists(cfx):
                    os.remove(cfx)
                
        ## Step 3: Sort reads into groups
        file_check = f'{work_dir}/step03_group_reads_{data_name}.txt'
        if not os.path.exists(file_check):
            p_group_reads = functools.partial(
                group_reads,
                work_dir = work_dir,
                cleanup = cleanup,
                rapidfuzz_error_correct = True,
                **kwargs
            )
            with mp.Pool(n_workers) as pool:
                chunked_dirs = pool.map(p_group_reads, trimmed_fastqs)
            logger.info(f'{len(chunked_dirs)} directories with sorted reads.')
            with open(file_check, 'w') as fclog:
                for line in chunked_dirs:
                    fclog.write(line + '\n')
        else:
            with open(file_check, 'r') as fclog:
                chunked_dirs = [line.strip() for line in fclog]

        ## Step 4A: CURIO processing
        file_check = f'{work_dir}/step04a_process_curio_{data_name}.txt'
        if not os.path.exists(file_check):
            curio_files = []
            do_curio = True
            if "CURIO" not in groups.keys():
                logger.warning(f"No CURIO in groups.keys(): Assuming do not have to parse CURIO barcodes")
                do_curio = False
            else:
                for chunkx in chunked_dirs:
                    fx = glob.glob(f'{chunkx}/CURIO_*.fastq')
                    if len(fx)==1:
                        curio_files.extend(fx)
            if len(curio_files)==0:
                logger.warning(f"No CURIO_*.fastq files found: Assuming do not have to parse CURIO barcodes")
                do_curio = False
            if do_curio:
                p_process_curio = functools.partial(
                    process_curio,
                    curio_whitelist = curio_whitelist,
                    work_dir = work_dir,
                    cleanup = False,
                    **kwargs
                )
                with mp.Pool(n_workers) as pool:
                    pool.map(p_process_curio, curio_files) #No output
                logger.info(f'{len(curio_files)} CURIO FASTQs processed.')
            Path(file_check).touch()
                

        ## Step 4B: STICR processing
        file_check = f'{work_dir}/step04b_process_sticr_{data_name}.txt'
        if not os.path.exists(file_check):
            sticr_files = []
            do_sticr = True
            if "STICR" not in groups.keys():
                logger.warning(f"No STICR in groups.keys(): Assuming do not have to parse STICR barcodes")
                do_sticr = False
            else:
                for chunkx in chunked_dirs:
                    fx = glob.glob(f'{chunkx}/STICR_*.fastq')
                    if len(fx)==1:
                        sticr_files.extend(fx)
            if len(sticr_files)==0:
                logger.warning(f"No STICR_*.fastq files found: Assuming do not have to parse STICR barcodes")
                do_sticr = False
            if do_sticr:
                p_process_sticr = functools.partial(
                    process_sticr,
                    sticr_whitelist_dir = sticr_whitelist_dir,
                    indices_to_consider = indices_to_consider,
                    sticr_whitelist = sticr_whitelist,
                    work_dir = work_dir,
                    cleanup = False,
                    **kwargs
                )
                with mp.Pool(n_workers) as pool:
                    pool.map(p_process_sticr, sticr_files) #No output
                logger.info(f'{len(sticr_files)} STICR FASTQs processed.')
            Path(file_check).touch()

        ## Step 5: Align RNA reads
        file_check = f'{work_dir}/step05_map_rna_{data_name}.txt'
        if not os.path.exists(file_check):
            rna_files = []
            if "RNA" != default_group.upper():
                logger.warning(f"Assuming {default_group} is our RNA group -- aligning with respect to these FASTQ files")
            for chunkx in chunked_dirs:
                fx = glob.glob(f'{chunkx}/{default_group}_*.fastq')
                if len(fx)==1:
                    rna_files.extend(fx)
            if len(rna_files)==0:
                logger.error(f"No {default_group}_*.fastq files found: something has gone wrong")
                sys.exit(1) ## Fatal error, cannot do the other steps

            ## Running into problems with getting this to work with multiprocessing
            # p_map_rna = functools.partial(
            #     map_rna,
            #     work_dir = work_dir,
            #     cleanup = cleanup,
            #     **kwargs
            # )
            # with mp.Pool(n_workers) as pool:
            #     pool.map(p_map_rna, rna_files)

            ## So instead, let's go back to looping through files
            for rna_filex in rna_files:
                logger.info(f'Processing {rna_filex}')
                map_rna(
                    rna_fastq = rna_filex,
                    work_dir = work_dir,
                    cleanup = cleanup,
                    # threads = os.cpu_count() - 1, ## Does not speed up (the memory caps probably stopping speed up).
                    **kwargs
                )
            logger.info(f'{len(rna_files)} RNA FASTQs processed.')
            Path(file_check).touch()

    ## This concludes the multiprocessing portion
    ## Next step is (usually) to merge all these bams together

    ## Step 6: Feature assign
    ## NB it may make sense to actually do multiprocesssing for every chromosome, but tricky, will ignore for now
    file_check = f'{work_dir}/step06_feature_assign.txt'
    if not os.path.exists(file_check):
        bams_to_process = glob.glob(f'{work_dir}/**/RNA_ALN_*.bam', recursive=True)
        bams_to_process = [os.path.dirname(x) for x in bams_to_process] ## Stupid of me, but next step takes in directories instead of files
        if merge_bams:
            logger.info(f'{len(bams_to_process)} BAM files to be merged into {sample_name}.')
        annotated_bam_files = feature_assign(
            bam_dir = bams_to_process,
            merge_bams = merge_bams,
            merged_bam_name = sample_name,
            method = feature_assign_method,
            work_dir = work_dir,
            cleanup = cleanup,
            **kwargs
        )
        with open(file_check, 'w') as fclog:
            for line in annotated_bam_files:
                fclog.write(line + '\n')
    else:
        with open(file_check, 'r') as fclog:
            annotated_bam_files = [line.strip() for line in fclog]

    ## Step 7: UMI deduplication
    file_check = f'{work_dir}/step07_deduplicate_umi.txt'
    if not os.path.exists(file_check):
        logger.info(f'{len(annotated_bam_files)} BAM files to undergo UMI deduplication')
        dedup_bam_files = deduplicate_umi(
            annotated_bam = annotated_bam_files,
            keep_all_reads = keep_all_reads,
            work_dir = work_dir,
            cleanup = cleanup
        )
        with open(file_check, 'w') as fclog:
            for line in dedup_bam_files:
                fclog.write(line + '\n')
    else:
        with open(file_check, 'r') as fclog:
            dedup_bam_files = [line.strip() for line in fclog]

    ## Step 8: SNV calling
    file_check = f'{work_dir}/step08_snv_call.txt'
    if not os.path.exists(file_check):
        logger.info(f'{len(dedup_bam_files)} BAM files to be used for SNV calling')
        vcf_files = snv_call(
            dedup_bam = dedup_bam_files,
            method = snv_calling_method,
            work_dir = work_dir,
            cleanup = cleanup,
            **kwargs
        )
        logger.info(f'{len(vcf_files)} VCF files created')
        Path(file_check).touch()

    ###################################
    ## ----- Core Pipeline END ----- ##
    ###################################

    logger.info(f'Consolidating barcode calls')
    final_output_dir = f'{work_dir}/OUTPUT/'
    if not os.path.exists(final_output_dir):
        os.makedirs(final_output_dir, exist_ok=True)
    parsed_barcode_files = subprocess.run(['find', work_dir, '-type', 'f', '-name', "*_ASSIGNS.tsv"], capture_output=True, text=True).stdout.split('\n')
    parsed_barcode_files.remove('')
    group_outs = dict.fromkeys(groups)
    for assign_filex in parsed_barcode_files:
        for group_namey in list(group_outs.keys()):
            if group_namey in assign_filex:
                assign_tagx = re.sub('_ASSIGNS.tsv$', '_TAGS.tsv', assign_filex)
                if not os.path.exists(assign_tagx):
                    sys.exit(1)
                assign_df = pl.read_csv(assign_filex, separator='\t')
                if assign_df.shape[0]==0:
                    continue
                assign_tag = pl.read_csv(assign_tagx, separator='\t')
                assign_tag = assign_tag.filter(
                    pl.col('read_id').is_in(assign_df['read_id'].implode())
                )
                assign_final = assign_df.join(assign_tag, on='read_id', how='inner', suffix='_cellbc')
                if group_outs[group_namey] is None:
                    group_outs[group_namey] = [assign_final]
                else:
                    group_outs[group_namey].append(assign_final)
    for key,value in group_outs.items():
        group_outx = ''
        if value is None:
            logger.warning(f'No assigned reads for {key}')
            continue
        group_outx = pl.concat(value)    
        group_outx.write_csv(f'{final_output_dir}{sample_name}_{key}.csv')

    ## QoL: Delete all the GENOME files
    files_to_delete = glob.glob(f'{work_dir}/GENOME*')
    files_to_delete.extend(glob.glob(f'{work_dir}/*chr.gtf'))

    ## QoL: Merge together all the STICR and CURIO
    ## TO DO: FIGURE THIS OUT

    if cleanup:
        for dfile_path in files_to_delete:
            if os.path.exists(dfile_path):
                os.remove(dfile_path)

    logger.info(f'Results are in {work_dir}')
    return work_dir