import subprocess
import os

def mafftRaw(input_file):
    # Define the MAFFT command as a list of strings
    mafft_command = ["mafft", "--reorder", "--thread", "30", input_file]
    outf= str(input_file).split("_")[0] + "_mafftRaw.fasta"
    # Use subprocess to run the command and redirect the output to a file
    with open(outf, 'w') as outfile:
        subprocess.run(mafft_command, stdout=outfile)

def mafftFinal(list_files):
    for input_file in list_files:
        # Define the MAFFT command as a list of strings
        mafft_command = ["mafft", "--reorder", "--thread", "30", input_file]
        outf= str(input_file).split(".")[0] + "_mafftFinal.fasta"
        # Use subprocess to run the command and redirect the output to a file
        with open(outf, 'w') as outfile:
            subprocess.run(mafft_command, stdout=outfile)
