from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import re
from collections import defaultdict

def main():
	opts = parse_args()

	process_file(opts.file, opts.output)

#########
# Setup #
#########

def parse_args():
	# House keeping to read in arguments from the command line
	parser = argparse.ArgumentParser()

	parser.add_argument("-f", "--file",
		help = "Path to the sequence file.",
		action = "store",
		metavar = "\b")

	parser.add_argument("-o", "--output",
		help = "Path to output fasta",
		action = "store",
		metavar = "\b")

	args = parser.parse_args()

	return args

def process_file(file_path, out_file):
	locus_mapping = {"IGH": "VTSS", "IGK": "VEIK", "IGL": "LTVL"}

	counting_dict = defaultdict(int)
	all_records = []
	with open(file_path, "r") as in_file:
		header = next(in_file)
		header = header.strip().split("\t")
		for line in in_file:
			line = line.strip().split("\t")
			full_line = dict(zip(header, line))
			clone_id = full_line["clone_id"]
			locus = full_line["locus"]
			save_name = clone_id + "_" + locus
			full_sequence = full_line["sequence"]
			sequence_start = full_line["v_sequence_start"]
			sequence_end = full_line["j_sequence_end"]
			start_seq = full_line["fwr1"]
			seq_end = full_line["fwr4"]
			test_sequence = full_sequence[int(sequence_start) -1 :int(sequence_end) - 1]
			start_seq = start_seq.split("...")[0]
			start_position = test_sequence.find(start_seq)
			end_position = test_sequence.find(seq_end[:-1])
			end_position = end_position + len(seq_end)
			difference = end_position - len(test_sequence)
			description = "vdj_seq"
			counting_dict[locus] += 1
			if start_position != 0:
				description += "_incomplete_match"
				count_dots = start_seq.count(".")
				test_sequence = "N" * count_dots + test_sequence
				new_test = full_sequence[int(sequence_start) - 1 - count_dots:int(sequence_end) -1]
				print(new_test)
				print(test_sequence)
				print()
			if difference -1 != 0:
				description += "_incomplete_match"

			if len(test_sequence) % 3 != 0:
				difference = len(test_sequence) % 3
				if difference == 1:
					test_sequence += "N" * 2
				elif difference == 2:
					test_sequence += "N" *1
				description += "_incomplete_match"
			dna_seq = Seq(test_sequence)

			# Translate the DNA sequence to amino acids
			amino_acid_seq = dna_seq.translate()
			last_four = amino_acid_seq[-4:]
			expected = locus_mapping[locus]
			description += "_end_amino_acid_" + str(last_four)

			# if last_four != expected:
			# 	#print(test_sequence)
			# 	print(last_four)
			# 	print(expected)
			# print()

			record = SeqIO.SeqRecord(id = save_name, seq = dna_seq, description = description)
			all_records.append(record)

	SeqIO.write(all_records, out_file, "fasta")
	print(counting_dict)


if __name__ == "__main__":
	main()