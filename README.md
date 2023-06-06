# PiCoHaemoGen: PacBio-integrated Computational Pipeline for Haemosporida Mitochondrial Genome Analysis 🖥️ 🦟 🧬 🦠


PiCoHaemoGen is a computational pipeline designed to analyze PacBio sequenced Haemosporida mitochondrial genomes. It executes a sequence of operations that involves correcting sequence direction, generating a consensus sequence for each genetic cluster discovered by a Variational Autoencoder (VAE), a machine learning model for DNA sequence clustering.

How to run it?

`python PiCoHaemoGen.py -rR folder/data`

More options:
`
optional arguments:
  -h, --help            show this help message and exit
  -rR RAWREADS, --rawReads RAWREADS
                        folder with rawReads from PacBio sequencing
  -pF PRIMERF, --primerF PRIMERF
                        primerF 5'-3', default=GATTCTCTCCACACTTCAATTCGTACTTC
  -pR PRIMERR, --primerR PRIMERR
                        primerR 5'-3, default=GAAGTACGAATTGAAGTGTGGAGAGAATC
  -rF REMOVEFILES, --RemoveFiles REMOVEFILES
                        Removing temporal files, default=no
`
Requirements: 
mafft (MAFFT - a multiple sequence alignment program) and some python modules (Biopython, tensorflow, numpy, etc).


                        +----------------------+
                        | Primer Detection and |
                        |   Read Redirection   |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        |   Execution of Raw   |
                        |   Alignment (MAFFT)  |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        | Conversion of DNA to |
                        |     Binary Format    |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        |  VAE Model Execution |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        | Cluster Extraction   |
                        |   (DBSCAN)           |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        |Conversion of Clusters|
                        |   to FASTA Format    |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        | Execution of Final   |
                        |   Alignment (MAFFT)  |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        | Polishing of Output  |
                        |   Sequences          |
                        +-----------+----------+
                                    |
                        +-----------v----------+
                        | Removal of Temporary |
                        |   Files              |
                        +----------------------+

