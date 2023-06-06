# PiCoHaemoGen: PacBio-integrated Computational Pipeline for Haemosporida Mitochondrial Genome Analysis


PiCoHaemoGen is a computational pipeline designed to analyze PacBio sequenced Haemosporida mitochondrial genomes. It executes a sequence of operations that involves correcting sequence direction, generating a consensus sequence for each genetic cluster discovered by a Variational Autoencoder (VAE), a machine learning model for DNA sequence clustering.

How to run it? 
python PiCoHaemoGen.py -rR folder/data

Requirements: 
mafft (MAFFT - a multiple sequence alignment program)
Some python module (Biopython, tensorflow, numpy, etc)


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

