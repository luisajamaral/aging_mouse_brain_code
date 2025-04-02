import csv
import subprocess

# Function to run the command for a single row and save error messages
def run_command(sample, bam_file, threads, te_annotation, prefix, output_dir, error_file):
    bam_file_with_suffix = bam_file + "/outs/gex_possorted_bam.bam"
    cmd = f"python /mnt/silencer2/home/lamaral/software/SoloTE/SoloTE_pipeline.py --threads {threads} --bam {bam_file_with_suffix} --teannotation {te_annotation} --outputprefix {prefix} --outputdir {output_dir}"
    with open(error_file, 'a') as error_output:
        subprocess.call(cmd, shell=True, stderr=error_output)

# Path to your CSV file
csv_file = '/mnt/silencer2/home/lamaral/ps-renlab2/projects/combined_all/tracksheet_with_paths_oct1723.csv'
error_log_file = 'error_log.txt'  # Define the name of the error log file

# Read CSV and process samples sequentially, saving error messages to the log file
with open(csv_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        run_command(row['Sample_name'], row['Path-luisa'], '15', '/mnt/silencer2/home/lamaral/software/SoloTE/mm10_rmsk.bed', row['Sample_name'], '/mnt/silencer2/home/lamaral/ps-renlab2/projects/fmouse_multiome/SoloTE_out', error_log_file)
