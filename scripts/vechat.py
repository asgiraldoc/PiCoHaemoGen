import subprocess

def vechatRun(list_files):
    for input_file in list_files:
        inp = str(input_file).split("_")[0]
        vechat_command = ["vechat", "--platform", "pb", "--scrub", 
                        "-o", inp + "_corrected.fasta", ## output file
                        "-t", "31", 
                        input_file] ## input file
        subprocess.run(vechat_command)