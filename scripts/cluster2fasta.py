import os
import subprocess
from Bio import SeqIO


def txt2fasta(input_file, headers_file, mapped_output_file):
    # Read the list of headers
    with open(headers_file, "r") as f:
        headers = [line.strip() for line in f]

    # Create a set of headers for more efficient searching
    header_set = set(headers)

    # Open the output file
    with open(mapped_output_file, "w") as mapped_out_file:
        # Iterate over the sequences in the multifasta file
        for record in SeqIO.parse(input_file, "fasta"):
            # Check if the sequence's header is in the list
            if record.id in header_set:
                # Write the sequence to the mapped sequences file
                SeqIO.write(record, mapped_out_file, "fasta")

# os.chdir("../example_data")
# dirc = os.getcwd()

# txtCluster = [f for f in os.listdir(dirc) if f.endswith('-.txt') and f.startswith("test")]
# print(txtCluster)

# for headers_file in txtCluster:
#     mapped_output_file = headers_file.split(".")[0] + ".fa"
#     txt2fasta("test.fasta", headers_file, mapped_output_file)