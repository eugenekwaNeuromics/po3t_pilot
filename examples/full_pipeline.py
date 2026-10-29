## Script version: 251111

## To use:
## nohup python -u full_pipeline.py &> nohup.out &

## PARAMS -----
package_dir 	= '/location_of_po3t_directory'
data_dir 	= '/location_of_your_data'
genome_dir 	= '/location_of_genome'
out_dir 	= '/location_of_output_dirs'
curio_whitelist = {
    'sample_name': [f'/location/XXX_BeadBarcodes.txt'], 
}
sticr_whitelist_dir = f'/location_of_sticr_directory'

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
samples_to_drop = []
for sx in samples_to_drop:
    print(f'Dropping {sx}')
    del data_dict[sx]

print(f'{len(data_dict)} samples found...')

## Running pipeline -----
## Day 55 (851) is 013 slide, while Day 80 (852) is 014 slide


if __name__ == "__main__":
    result_dirs = []
    for sample_name,fastq_files in data_dict.items():
        desired_results_folder = f'{out_dir}/RESULTS_{sample_name}'
        print(f'Running pipeline on sample {sample_name}')
        results_dir = po3t.run_pipeline(
            fastq_files = fastq_files,
            genome_dir = genome_dir,
            curio_whitelist = curio_whitelist[sample_name],
            sticr_whitelist_dir = sticr_whitelist_dir,
            sample_name = sample_name,
            work_dir = desired_results_folder,
            chunk_size = 20 * 1000,
        )
        result_dirs.append(result_dir)

    print(f'{len(result_dirs)} result directories created.')
    print(f'Done!')