#!/usr/bin/python3

# Xavier Rolon
# Genomic Variations: Search for Motifs in DNA (Option 5)
# Date Created: 11/24/23
# ---------------------------------------------------------------------------

# import sys module to access command line arguments
# import os module to check whether a path is an exisiting file (used for query)
import sys
import os

# function to read dna files and return the sequence as a string
def read_dnafile(filename):
	with open(filename, 'r') as file:
		dna_sequence = ''	# empty string to hold sequence
		for line in file:	# loop over each line
			if not line.startswith('>'):	# if line does not start with '>', it is part of the sequence
				dna_sequence += line.strip()	# remove newline characters and add the line to the sequence
		if not dna_sequence:
			print(f"Warning: The file {file} might be corrupted.")
		return dna_sequence.upper()

# function to compare two sequences and return their differences
def compare_sequences(reference_sequence, sequence):
	minimum_length = min(len(reference_sequence), len(sequence))	# find length of shorter sequence

	# create list of tuples: each tuple contains the index and differing bases
	differences = [(i, reference_sequence[i], sequence[i]) for i in range(minimum_length) if reference_sequence[i] != sequence[i]]
	return differences

# function to find all occurences of a query in a sequence
def find_query(sequence, query):
	# list of starting indices for all occurences of the query
	location = [i for i in range(len(sequence)) if sequence.startswith(query, i)]
	return location

# main function that runs on execution
def main():
	# accept at least two pathnames of dna files, could be FASTA files
	if len(sys.argv) < 3:
        	print("Enter at least two file arguments. Usage: ./gvmotifs.py <file1> <file2> ...")
        	sys.exit()

	# reference file is first argument read
	reference_file = sys.argv[1]

	# call function read_dnafile to read the reference file sequence
	reference_sequence = read_dnafile(reference_file)

	# check if the reference sequence has only A,C,G, and T and no other letter characters, otherwise, print an error and exit
	valid_bases = {"A", "C", "G", "T"}
	for base in reference_sequence:
		if base not in valid_bases:
			print(f"Error: ", reference_file, " contains invalid characters.")
			sys.exit()

# PART 1:
# print out the differences for two dna files, and where they are located
# sequences are of different lenghts, comparison up to length of shorter sequence
# presentation of differences and their locations should be done in organized manner
# most likely require string formatting

	# reference file name
	print(f"Part 1:\nReference: {reference_file}\n")
	print("Comparing DNA files:")

	# loop over the remaining arguments except query argument
	for file in sys.argv[2:]:
		# check for text files and FASTA files (to keep track of query argument
		if file.endswith('.fasta') or file.endswith('.txt'):
			print(f"    {file}")
			# read the file to be checked with the reference sequence
			sequence = read_dnafile(file)

			# if the file has invalid characters, error
			for base in sequence:
				if base not in valid_bases:
					print(f"Error: ", file, " contains invalid characters.")
					sys.exit()

			# call function to compare sequences between the reference, and current sequence of current file
			differences = compare_sequences(reference_sequence, sequence)
			print("\nLocations of differences: ")

			# loop through function to locate all differences and print in a readable format
			for difference in differences:
				print(f"    Position: {difference[0] + 1}, Reference File: {difference[1]}, Compared File: {difference[2]}")
			print("\n")
# PART 2:
# if optional query is given, print out all locations of query within all files
	# last argument is the query if it was given
	if len(sys.argv) > 3:
		query = sys.argv[-1]
		# check if  query argument is a file
		if os.path.isfile(query) or set(query.upper()).issubset({'A', 'C', 'G', 'T'}):
		# if True, read the query from the file
			if os.path.isfile(query):
				query = read_dnafile(query)
			# otherwise treat the query as a a string
			else:
				query = query.upper()


			# query to be given is at most 1k-bp (1000 base pairs long)
			if len(query) > 1000:
				print(f"Warning: The query must be at most 1000 base pairs long. Skipping comparison.")
			# if check fails, skip to next iteration and print the query
			else:
				print(f"Part 2:\nQuery: ", query)

			# find locations of query in reference sequence first
				query_locations = find_query(reference_sequence, query)
				print("Motif found at the following positions in", reference_file, ":")
				for location in query_locations:
					print(f"    {location + 1}")

				# loop over remaining arguments to find query in other DNA files
				for file in sys.argv[2:-1]:
					if file.endswith('.fasta') or file.endswith('.txt'):
					# read the sequences from dna files
						sequence = read_dnafile(file)
					# find the locations of the query
						query_locations = find_query(sequence, query)
						print(f"Motif found at the following positions in ", file,":")
						for location in query_locations:
							print(f"    {location + 1}")
				print("\n")
# run the main function when the script is executed
if __name__ == "__main__":
	main()
