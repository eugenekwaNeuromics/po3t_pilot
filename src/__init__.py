package_name = "PO3T"
package_description = "Pipeline for processing reads from Oxford Nanopore Technologies (ONT) to obtain single cell (10X), spatial (CurioTrekker), and lineage tracing (STICR) data"
__version__ = "0.0.0"

## May not be a good idea to import EVERYTHING from .pipeline_functions
# from .pipeline_functions import *

## Let's instead just load the necessary functions

## Step 1: Set-up and genome prep
from .pipeline_functions import setup_workdir, get_genome

## Step 2: 10X cellular barcodes
from .pipeline_functions import fuzz_assign, get_10x_whitelist, extract_10x_barcodes

## Step 3: Sorting reads into groups
from .pipeline_functions import group_reads

## Step 4: Parsing other barcodes
from .pipeline_functions import process_curio, process_sticr

## Step 5: RNA read mapping
from .pipeline_functions import map_rna

## Step 6: Feature assignment
# from .pipeline_functions import fa_stringtie
from .pipeline_functions import feature_assign

## Step 7: UMI deduplication
from .pipeline_functions import cluster_umi, deduplicate_umi

## Step 8: SNV calling
# from .pipeline_functions import sc_longshot
from .pipeline_functions import snv_call

## The full pipeline
from .run_pipeline import run_after_epi2me