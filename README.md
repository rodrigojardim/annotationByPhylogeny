# annotationByPhylogeny

This tool performs the functional annotation of nucleotide or amino acid sequences based on their position in a phylogenetic tree. The algorithm looks up to two levels above the sequence's position, giving it a new annotation based on its ancestors.

Requirements:
1. Python > 3.9
2. Biopython > 1.79

How to run:

python annotByPhylo.py -t <newick file> -f <original fasta file> -a ../data/headers_itol_annotation.txt -o saidaPhylo.fasta -u ../data/Unknown.txt
