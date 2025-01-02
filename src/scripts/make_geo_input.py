import subprocess
import os
import re
import shutil

def main():
	# Define directories
	input_file = "files/meta_mapping.csv"
	geo_dir = "GEO_submission"
	data_dir = "results/fastqs"
	results_dir = "results"
	antibody_reference = "files/antibodies.csv"
	all_samples = get_samples(input_file)
	move_fastq_files(geo_dir, data_dir, all_samples)
	copy_data_files(geo_dir, results_dir, all_samples)
	copy_antibody_reference(geo_dir, antibody_reference)

def get_samples(input_file):
	"""
	Pulls out the old and new sample names for mapping
	"""
	sample_dict = {}
	with open(input_file, "r") as in_file:
		next(in_file)
		for line in in_file:
			# Make a dictionary with the name as the key and new id as the value
			row, sid, name = line.strip().split(",")
			sample_dict[name] = sid

	return(sample_dict)


def move_fastq_files(geo_dir, data_dir, all_samples):
	"""
	Finds the fastq files using the results/fastqs directory because these have all been renamed 
	to have consistent naming
	"""
	all_fastq = os.listdir(data_dir)

	# Define the pattern for the current sample names
	pattern = r"JH[0-9]*-[0-9]*_[A-Z]*"
	for file_name in all_fastq:
		sample_name = re.search(pattern, file_name).group()

		# Replace the original sample name with the new id from our dictionary made earlier
		new_name = re.sub(pattern, all_samples[sample_name], file_name)
		full_file_path = os.path.join(os.path.abspath(data_dir), file_name)
		new_file_path = os.path.join(geo_dir, new_name)

		# Check if the destination link already exists
		if os.path.islink(new_file_path):
			os.remove(new_file_path)
		
		# Create a new soft link in the appropriate directory with the appropriate names
		subprocess.run(["ln", "-s", full_file_path, new_file_path], check=True)

def copy_data_files(geo_dir, results_dir, all_samples):
	"""
	Finds and copies all of the 10x processed files
	"""

	for sample in all_samples.keys():
		# Now find the required files
		output_dir = os.path.join(results_dir, sample, "outs/per_sample_outs", sample)

		count_dir = os.path.join(output_dir, "count", "sample_filtered_feature_bc_matrix")

		for file_name in os.listdir(count_dir):
			rsync_command(file_name = file_name, source_dir = count_dir, target_dir = geo_dir, sample = all_samples[sample])

		# Copy all the vdj files
		vdj_dir = os.path.join(output_dir, "vdj_b")
		vdj_files = ["filtered_contig_annotations.csv", "concat_ref.bam", "concat_ref.bam.bai", "airr_rearrangement.tsv"]
		vdj_files = ["filtered_contig_annotations.csv", "airr_rearrangement.tsv"]

		for file_name in vdj_files:
			rsync_command(file_name = file_name, source_dir = vdj_dir, target_dir = geo_dir, sample = all_samples[sample])

def rsync_command(file_name, source_dir, target_dir, sample):
	source_file = os.path.join(os.path.abspath(source_dir), file_name)
	new_file = os.path.join(target_dir, sample + "_" + file_name)

	subprocess.run(["ln", "-s", source_file, new_file], check=True)

def copy_antibody_reference(geo_dir, antibody_reference):
	subprocess.run(["rsync", "-avz", antibody_reference, geo_dir])

if __name__ == "__main__":
	main()