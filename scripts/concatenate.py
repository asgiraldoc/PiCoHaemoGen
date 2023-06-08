import os
import subprocess

def concat(filesC0, filesC1):
    for i in zip(filesC0, filesC1):
        outf = str(i[0]).split("-")[0] + "-_RawHap.fasta"
        with open(outf, "w") as f:
            cat_command = ["cat", i[0], i[1]]
            subprocess.run(cat_command, stdout=f)

