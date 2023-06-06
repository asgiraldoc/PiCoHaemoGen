import sys
import os
import argparse

sys.path.append('scripts/')

## import functions
from readsRedirection import primerDetection
from mafft import mafftRaw, mafftFinal
from fasta2binary import fasta2bin
from VAEmethod import VAE_model
from DBScan import extract_cluster_labels_dbscan
from cluster2fasta import txt2fasta
from consensus import polishing
from rmTmp import remove_temp_files


def main():
    ## Arguments 
    parser = argparse.ArgumentParser(description='PiCoHaemoGen: PacBio-integrated Computational Pipeline for Haemosporida Mitochondrial Genome Analysis')
    ### main arguments
    parser.add_argument('-rR', '--rawReads', type=str, help="rawReads from PacBio sequencing", required=True)
    ### readsRedirection arguments
    parser.add_argument('-pF', '--primerF', type=str, help="primerF 5'-3', default=GATTCTCTCCACACTTCAATTCGTACTTC", default='GATTCTCTCCACACTTCAATTCGTACTTC')
    parser.add_argument('-pR', '--primerR', type=str, help="primerR 5'-3, default=GAAGTACGAATTGAAGTGTGGAGAGAATC", default='GAAGTACGAATTGAAGTGTGGAGAGAATC')
    parser.add_argument('-rF', '--RemoveFiles', type=str, help="Removing temporal files, default= yes", default='yes')

    args = parser.parse_args()

    ## files
    rawReads = args.rawReads
    nameSample = str(rawReads).split(".")[0]
    ## run readsRedirection
    primerF = args.primerF
    primerR = args.primerR
    primerDetection(rawReads, primerF, primerR)

    rmFiles = args.RemoveFiles
    print(rawReads, "##########")
    ## run raw aligment
    mafftRawOut= str(rawReads).split(".")[0] + "_nolong.fasta"
    mafftRaw(mafftRawOut)

    ## run convert DNA seq to Binary format
    fasta2binOut= str(rawReads).split(".")[0] + "_mafftRaw.fasta"
    fasta2bin(fasta2binOut)

    ## run VAE program and clustering
    VAErunOut = str(rawReads).split(".")[0] + "_bin.txt"
    mu, VAErunData = VAE_model(VAErunOut)
    nameCluster = str(rawReads).split(".")[0]
    extract_cluster_labels_dbscan (mu, VAErunData, nameCluster)

    txtCluster = [f for f in os.listdir() if f.endswith('-.txt') and f.startswith(nameCluster)]
    for headers_file in txtCluster:
        mapped_output_file = headers_file.split(".")[0] + ".fa"
        txt2fasta(mafftRawOut, headers_file, mapped_output_file)

    ## run final aligment
    filesM = [f for f in os.listdir() if f.endswith('.fa') and f.startswith(nameSample)]

    mafftFinal(filesM)

    ## polishing output
    filesFinal = [f for f in os.listdir() if f.endswith('_mafftFinal.fasta') and f.startswith(nameSample)]
    polishing(filesFinal)

    # removing temporal files
    remove_temp_files(rmFiles)

if __name__ == '__main__':
    main()