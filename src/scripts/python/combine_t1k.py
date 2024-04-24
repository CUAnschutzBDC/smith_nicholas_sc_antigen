import pandas as pd
import argparse
import os
from collections import defaultdict

def main():
	opts = parse_args()
	sample_metadata = pd.read_excel(opts.sample_metadata, sheet_name = "New_metadata")
	t1k_dict = read_t1k(opts.samples, opts.results_dir)
	sample_dict = make_sample_dict(sample_metadata)
	make_output(opts.out_file, t1k_dict, sample_dict)


def read_t1k(samples, results_dir):
	"""
	Reads through the output of T1K and keeps only the HLA information we care about
	"""
	t1k_dict = defaultdict(dict)
	for sample in samples:
		read_file = os.path.join(results_dir, sample, "t1k", sample + "_genotype.tsv")
		with open(read_file, "r") as in_file:
			full_dict = {}
			for line in in_file:
				line = line.strip().split("\t")
				if line[0] in ["HLA-DRB1", "HLA-DQB1"]:
					hla_dict = {}
					locus, count, allele_1, abundance_1, quality_1, allele_2, abundance_2, quality_2 = line
					if int(quality_1) > 1:
						hla_dict[allele_1] = abundance_1
					if int(quality_2) > 1:
						hla_dict[allele_2] = abundance_2
					full_dict[locus] = hla_dict
			t1k_dict[sample] = full_dict

	return(t1k_dict)

def make_sample_dict(df):
	"""
	Pulls the sample name and HLA type from the sample info
	"""
	sample_dict = {}
	for index, row in df.iterrows():
		sample_dict[row['Sample Name']] =  row['HLA type']

	return(sample_dict)


def make_output(output_file, t1k_dict, sample_dict):
	with open(output_file, "w") as out_file:
		out_file.write("sample\thla_type\tallele_1_DR\tallele_1_DR_count\tallele_2_DR\tallele_2_DR_count\t" +
			"allele_1_DQ\tallele_1_DQ_count\tallele_2_DQ\tallele_2_DQ_count\n")
		for sample in t1k_dict:
			if sample == "JH313-15_JB":
				orig_sample = "JH313-15_JBM"
			else:
				orig_sample = sample
			out_file.write("{}\t{}".format(sample, sample_dict[orig_sample]))
			for hla in ["HLA-DRB1", "HLA-DQB1"]:
				for allele in t1k_dict[sample][hla]:
					out_file.write("\t{}\t{}".format(allele, t1k_dict[sample][hla][allele]))
				if(len(t1k_dict[sample][hla])) == 1:
					out_file.write("\t\t")
			out_file.write("\n")


#########
# Setup #
#########

def parse_args():
	# House keeping to read in arguments from the command line
	parser = argparse.ArgumentParser()

	parser.add_argument("-m", "--sample_metadata",
		help = "Path to the sample metadata.",
		action = "store",
		metavar = "\b")

	parser.add_argument("-r", "--results_dir",
		help = "Path to the results directory.",
		action = "store",
		metavar = "\b")

	parser.add_argument("-o", "--out_file",
		help = "Path to the output file.",
		action = "store",
		metavar = "\b")

	# 	default = "results")
	parser.add_argument("-s", "--samples",
		help = "List of samples to process.",
		action = "store",
		nargs = "+")

	args = parser.parse_args()

	return args


if __name__ == "__main__":
	main()