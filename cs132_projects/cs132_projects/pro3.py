#!/usr/bin/python3

# Xavier Rolon
# Open Reading Frames (Program Option 8)
# Date Created: 11/28/23
#-----------------------------------------------------------------------------

# import sys module to access command line arguments
# import re module for regular expression search
import sys
import re

# function to create a reverse complement string of a DNA sequence
# adenine(a) - thymine(t), cytosine(c) - guanine(g)
def reverse_complement(dna):
	# dictionary defined where key values are bases:complement bases
	complement_dna = {'a': 't', 't': 'a', 'c': 'g', 'g':'c'}
	# reverse the order of bases in the dna sequence
	reversed_dna = reversed(dna)
	# replace each base with its complement using a for loop
	complemented_dna = []
	for base in reversed_dna:
		complemented_dna.append(complement_dna[base])
	# join the complemented bases together into a string
	reverse_complement_dna = "".join(complemented_dna)
	return reverse_complement_dna	# reverse complement of original DNA

# function to find all of the open reading frames in a given DNA sequence
# called for each possible reading frames (-3 to 3)
# note: for negativee reading frames, reverse complement of the sequence used
def find_reading_frames(dna, reading_frame):
	# start codon 'atg' and ends with stop codon 'taa', 'tag', 'tga'
	start_codon = 'atg'
	stop_codons = ['taa', 'tag', 'tga']

	# empty list to store reading frames
	reading_frames = []
	# loop over DNA sequence starting at given frame and stepping 3 bases at a time
	for i in range(reading_frame - 1, len(dna), 3):
		if dna[i:i + 3] == start_codon:	# if start codon, enter another loop
			# inner loop steps 3 bases at a time
			for j in range(i, len(dna), 3):
				# if stop codon found, the frame from start to stop added to open reading frame list
				if dna[j:j + 3] in stop_codons:
					reading_frames.append((reading_frame, dna[i:j + 3]))
					break # exit loop after adding frame
	# after all positions checks, return the list of reading frames
	# each frame represented as a tuple: frame is first element and the open reading frame is the second element
	return reading_frames

# function to translate a DNA sequence into an amino acid sequence
def translation(dna):
	# declare amino acid dictionary
	amino_acids = {
    "Phe": ["TTC", "TTT"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Met": ["ATG"],
    "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Tyr": ["TAT", "TAC"],
    "His": ["CAC", "CAT"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAT", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAT", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["TGT", "TGC", "TGC"],
    "Trp": ["TGG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "***": ["TAA", "TAG", "TGA"],
}
	protein = ''
	# go through dna sequence and set the codons in steps of 3
	for i in range(0, len(dna), 3):
		codon = dna[i:i + 3].upper()	# upper case since sequence is lower and dictionary is upper
		# inner loop to check for corresponding amino acid
		for amino_acid, codons in amino_acids.items():
			if codon in codons:	# if codon is the same as the codons in dictionary
				# set the empty protein string to corresponding amino acid
				protein += amino_acid + " "
				break	# exit for loop
	# return protein string containing amino acid sequence
	return protein.strip()


# main function
def main():
	# check if user provided valid pathname
	if len(sys.argv) < 2:
		print("Enter only one dna file to be read. Usage: openreadingframe.py <dnafile>")
		sys.exit()

	# get pathname of DNA file from argument
	filename = sys.argv[1]

	# open the file and read the sequence
	try:
		with open(filename, 'r') as file:
			# remove whitespace and newline characters
			dna = file.read().strip()

			# check for FASTA files, starting with '>' in header line
			if dna.startswith('>'):
				dna = dna.split('\n', 1)[1].replace('\n', '') # remove the '>' header line

			# validate dna sequence checking only a,c,g, and t characters
			if re.search('[^ACGTacgt]', dna):
				print("Invalid DNA sequence. Must only contain characters 'a', 'c', 'g', and 't'")
				sys.exit()
			# lower case sequence for processing
			dna = dna.lower()
	except FileNotFoundError:
		print("File not found", file = sys.stderr)
		sys.exit()


	# initialize list of frames to represent open reading frames
	frames = [1, 2, 3, -1, -2, -3]
	# create new variable if reading frame is found at a position
	found = False

	# loop through list of frames and call the function to find the reading frames of the dna sequence
	# for negative frames, the reverse complement function is called
	for frame in frames:
		if frame > 0:
			reading_frames = find_reading_frames(dna, frame)
		else:
			reading_frames = find_reading_frames(reverse_complement(dna), frame)

		# loop through the reading_frames
		for reading_frame in reading_frames:
			found = True	# frame found in read
			codons = []	# list of the codons
			# loop through each reading frame item in steps of 3
			for i in range(0, len(reading_frame[1]), 3):
				codon = reading_frame[1][i:i + 3]
				codons.append(codon)	# add each codon to the list of codons
			codon_sequence = ' '.join(codons).upper()	# add codons onto one sequence

			# print sequence and the reading frame
			print('-' * 75)
			print(codon_sequence)
			print(f"Reading Frame: {reading_frame[0]}")

			# translating the codon sequence
			protein = translation(reading_frame[1])
			print(f"Amino acid sequence: {protein} \n")
	if not found:
		print("No open reading frame found at position.")

# call main function
if __name__ == "__main__":
	main()
