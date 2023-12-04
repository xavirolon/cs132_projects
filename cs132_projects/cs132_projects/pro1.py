#!/usr/bin/python3

# Xavier Rolon
# DNA Translation (Program Option 1)
# Date Created: 10/25/23
# Program accepts pathname of a dna file as the only argument,
# reporting back the codon sequence and corresponding amino acid
# sequence. Must also report additional information.

#--------------------------------------------------------------------------------------------------------

# Import sys module to access command line arguments
import sys

# Define a dictionary of amino acids and their corresponding codons
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

# Reverse the dictionary for easy lookup
codon_to_amino_acid = {codon: aa for aa, codons in amino_acids.items() for codon in codons}

# Check if user provided valid pathname
if len(sys.argv) < 2:  # only 2 arguments, python file and dna file it will read
    print("There is no file present, Usage: python_dna.py <filename>")
    sys.exit()

# Get the pathname of the DNA file from the command line argument
filename = sys.argv[1]

# Try to open file and read a DNA sequence
try:
    with open(filename, "r") as file:
        # Read the whole file as a string and then close it
        dna = file.read().replace("\n", "").replace(" ", "").upper()
except FileNotFoundError:
    print("File not found", file=sys.stderr)
    sys.exit()

# Validate DNA sequence
valid_nucleotides = {"A", "T", "C", "G"}
for char in dna:
    # Check if the nucleotides in list
    if char not in valid_nucleotides:
        print("Invalid DNA sequence. Must only contain letters 'a', 'c', 'g', and 't'")
        sys.exit()

# Store codons
codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]

# Print codons
print("Codon Sequence:", ' '.join(codons))

# Create empty list for valid sequences in file
sequences = []
current_sequence = []

# Loop over codons in DNA
for codon in codons:
    # Check if codon is a start codon
    if codon == "ATG":
        # Begin a new sequence with the start codon
        current_sequence = [codon]
    # Check if codon is a stop codon
    elif codon in ["TAA", "TAG", "TGA"]:
        # Check if current sequence not empty
        if current_sequence:
            # Add codon to current sequence
            current_sequence.append(codon)
            # Append current sequence to valid sequences list
            sequences.append(current_sequence)
            # Set current sequence to an empty list
            current_sequence = []
    # Check if current sequence not empty
    elif current_sequence:
        current_sequence.append(codon)

# Print valid codon sequences
print("Valid codon sequences:")
for seq in sequences:
    print(' '.join(seq).lower())

# Create empty string for amino acid sequences
amino_acid_sequences = []

# Loop through valid sequences
for seq in sequences:
    # Empty string for amino acid sequence
    amino_acid_sequence = []
    # Loop over codons in codon_sequence
    for codon in seq:
        # If codon in list of codons
        if codon in codon_to_amino_acid:
            # Add amino acid to amino acid sequence
            amino_acid_sequence.append(codon_to_amino_acid[codon])
    amino_acid_sequences.append(amino_acid_sequence)

# Print amino acid sequences
for i, aa_seq in enumerate(amino_acid_sequences):
    print(f"\nAmino Acid Sequence {i+1}:")
    print(' '.join(aa_seq))

    # Count amino acids
    aa_counts = {}
    for aa in aa_seq:
        if aa not in aa_counts:
            aa_counts[aa] = 0
        aa_counts[aa] += 1
    del aa_counts["***"]  # Do not count stop codon as an amino acid

    # Print counts of each amino acid
    print("\nCounts of each amino acid:")
    for aa, count in sorted(aa_counts.items(), key=lambda x: (-x[1], x[0])):
        print(f"{aa}: {count}")

    # Print total number of different amino acids
    print(f"\nTotal number of different amino acids: {len(aa_counts)}")
