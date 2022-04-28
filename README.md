# Optimized-sgRNA-Design

![](datafiles/CRISPR-Logo.jpg)

This is a set of software tools for analyzing the on-target and off-target effects of sgRNA spacer inputs

This program is still currently under development. Intended functionality is that this program can take in a metagenome consisting of a large number of organisms' genetic data, along with a collection of 20 base pair (no PAM) or 23 base pair (with PAM) CRISPR Cas9 sgRNA sequences, and to output those sgRNA sequences ranked in terms of optimality. This program has the following required dependencies: Python 3.8, pysam, bowtie.

-- INSTALLATION/ENVIRONMENT CREATION

-- TOOL USAGE INFORMATION

-- TOOL USAGE EXAMPLES

-- OUTPUT EXPLANATION

This program can be run in multi target input mode or single target input mode. In multi target input mode, the program can take a .TSV or .CSV file path representing the file from which the DNA string target sequences will be read; these strings can be 20 base pairs (just the targets), 23 base pairs (the targets + NGG PAMs), and in RNA or DNA format (though currently, no generic 'N' nucleotide functionality). The program also takes in a .FASTA file path representing the metaGenome with which the target sequences are to be found and analyzed, and (optionally) a .CSV or .TSV file path representing the file to which the target analysis data should be written (**if the save file already exists, its contents will not be appended to**). Each target analysis row will take the form:\n[Spacer (20 bp RNA), Target Guide Sequence (35 bp DNA), Off-Target Sequences (35 bp DNA), On-Target Score, Off-Target Scores, Heuristic]