## Script version: 251204

## To use:
## nohup python -u post_epi2me.py &> nohup.out &

## PARAMS -----
package_dir 	= '/location_of_po3t_directory'
data_dir 	= '/location_of_your_data'
genome_dir 	= '/location_of_genome'
samples_to_drop = [] #If need to drop a sample because was already processed/QC failed
curio_whitelist = {
    'sample_name': [f'/location/XXX_BeadBarcodes.txt'], 
}
sticr_whitelist_dir = f'/location_of_sticr_directory'
epi2me_parent_dir = '/location_of_EPI2ME_OUTPUTS/'

## Modules -----
import os
import re
import sys
import subprocess

## Setting up environment ----
name = 'po3t'
if not name in sys.modules:
    sys.path.insert(0, f'{package_dir}')
import po3t

data_files = subprocess.run(f'find {data_dir} -type f -name "*.fastq.gz"', shell=True, capture_output=True, text=True).stdout.split('\n')
data_files.remove('')
data_dict = {}
for data_file in data_files:
    sample_name = re.sub(f'{data_dir}|/$', '', os.path.dirname(data_file))
    if sample_name not in data_dict.keys():
        data_dict[sample_name] = []
    data_dict[sample_name].append(data_file)

## If desired, drop samples here
for sx in samples_to_drop:
    print(f'Dropping {sx}')
    del data_dict[sx]

print(f'{len(data_dict)} samples found...')

## Running pipeline -----
if __name__ == "__main__":
    result_dirs = []
    for sample_name,fastq_files in data_dict.items():
        desired_results_folder = f'{home_dir}/RESULTS_{sample_name}'
        
        print(f'Running pipeline on sample {sample_name}')
        
        epi2me_sample_dir = f'{epi2me_parent_dir}{sample_name}/'
        print(f'EPI2ME folder: {epi2me_sample_dir}')
        if not os.path.exists(epi2me_sample_dir):
            print(f'{epi2me_sample_dir} not found!')
            sys.exit(1)

        ## Would be faster if we only care about the most relevant cellular barcodes instead of all 7mio
        epi2me_cellbc_whitelist = f'{epi2me_sample_dir}{sample_name}.whitelist.tsv'
        print(f'EPI2ME 10X whitelist: {epi2me_cellbc_whitelist}')
        if not os.path.exists(epi2me_cellbc_whitelist):
            print(f'{epi2me_cellbc_whitelist} not found!')
            sys.exit(1)

        results_dir = po3t.run_after_epi2me(
            epi2me_output_dir = epi2me_sample_dir,
            fastq_files = fastq_files,
            genome_dir = genome_dir,
            curio_whitelist = curio_whitelist[sample_name],
            sticr_whitelist_dir = sticr_whitelist_dir,
            cellbc_whitelist_dir = epi2me_cellbc_whitelist,
            sample_name = sample_name,
            work_dir = desired_results_folder,
        )
        result_dirs.append(results_dir)
    print(f'{len(result_dirs)} result directories created.:')
    for result_dirx in result_dirs:
        print(result_dirx)
    print(f'Done!')