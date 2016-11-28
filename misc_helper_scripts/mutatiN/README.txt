README last updated 2015 Dec 17.

mutatiN.py accepts a sequential nucleotide fasta alignment (-i <file name>) and pseudorandomly changes nucleotides in each sequence to "N". The proportion of mutated positions in each sequence is controlled by the frequency parameter (-f <floating point>). 

For reproducibility, set the psuedorandom seed (defaults to 1; -s <integer>).

Users can also optionally protect existing gaps in the alignment (prevent them from being mutated) with -g T (defaults to F).

For full parameter information, use ./mutatiN.py --help.

The file 0.fa is provided as a small test input dataset.

The output is a fasta file (default name mutated_output) in the working directory.
