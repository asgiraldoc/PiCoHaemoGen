import subprocess

def mafftRaw(input_file, dir_out):
    # Define the MAFFT command as a list of strings
    mafft_command = ["mafft", "--reorder", "--thread", "30", input_file]
    outf0 = str(input_file).split("/")[-1]
    outf= str(outf0).split("_")[0] + "_mafftRaw.fasta"
    # Use subprocess to run the command and redirect the output to a file
    with open(outf, 'w') as outfile:
        subprocess.run(mafft_command, stdout=outfile)
    mv_command = ["mv", outf, dir_out]
    subprocess.run(mv_command)

def mafftFinal(list_files):
    for input_file in list_files:
        # Define the MAFFT command as a list of strings
        mafft_command = ["mafft", "--reorder", "--thread", "30", input_file]
        outf= str(input_file).split(".")[0] + "_mafftFinal.fasta"
        # Use subprocess to run the command and redirect the output to a file
        with open(outf, 'w') as outfile:
            subprocess.run(mafft_command, stdout=outfile)

def mafftHap(list_files):
    for input_file in list_files:
        # Define the MAFFT command as a list of strings
        mafft_command = ["mafft", "--thread", "30", input_file]
        outf= str(input_file).split("-")[0] + "-_aliHap.fasta"
        # Use subprocess to run the command and redirect the output to a file
        with open(outf, 'w') as outfile:
            subprocess.run(mafft_command, stdout=outfile)