from Bio import SeqIO  # Importing SeqIO from BioPython for parsing sequences.
import subprocess

# This function converts a nucleotide to its binary representation.
# It could be optimized by using a dictionary instead of multiple if/elif statements.
def nucleotide_to_binary(nucleotide):
    nucleotide = nucleotide.upper()  # Convert nucleotide to upper case to remove redundancy
    if nucleotide == "A":
        return "1,0,0,0"
    elif nucleotide == "C":
        return "0,1,0,0"
    elif nucleotide == "G":
        return "0,0,1,0"
    elif nucleotide == "T":
        return "0,0,0,1"
    else:
        return "0,0,0,0"

# This function converts a multifasta file to a binary sequence file.
def fasta2bin(input_file, dir_out):
    outf0 = str(input_file).split("/")[-1]
    nameFile = str(outf0).split("_")[0]
    output_file = str(outf0).split("_")[0] + "_bin.txt"


    with open(input_file, 'r') as fasta_file, open(output_file, 'w') as output:
        for record in SeqIO.parse(fasta_file, 'fasta'):  # Iterate over each record in the fasta file.
            header = record.id
            sequence = record.seq
            # Convert each nucleotide in the sequence to its binary representation and join them by spaces.
            binary_sequence = ' '.join([nucleotide_to_binary(nt) for nt in sequence])
            output.write(header + " " + nameFile + " " + binary_sequence + "\n")  # Write the header and binary sequence to the output file.
    mv_command = ["mv", output_file, dir_out]
    subprocess.run(mv_command)
