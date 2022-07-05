# Optimized-sgRNA-Design

![](CRISPR-Logo.jpg)

INTRODUCTION
=====================
The Optimized sgRNA Design repository is a set of command line software tools for analyzing the on-target and off-target effects of sgRNA spacer inputs developed from the Pickering Lab at UC Berkeley.

This program is still currently under development. Intended functionality is that this program can take in a metagenome consisting of a large number of organisms' genetic data, along with a collection of 20 base pair (no PAM) or 23 base pair (with PAM) CRISPR Cas9 sgRNA sequences, and to output those sgRNA sequences ranked in terms of optimality. This program has the following required dependencies: Python 3.8, pysam, biopython, bowtie.

INSTALLATION/ENVIRONMENT CREATION
=====================
This tool will only work on Linux and MacOS due to dependency on bowtie and pysam, neither of which are available on Windows. The steps to install the necessary dependencies for this program into a conda environment are as follows:
1. Install Python 3 (see http://www.python.org)

Linux: https://www.python.org/downloads/source/
MacOS: https://www.python.org/downloads/macos/
2.  Install conda (rcommend installing Miniconda rather than Anaconda due to size bloat, but both work since they are both just distributions of Conda; see https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Linux: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

MacOS: https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html

After installing the miniconda installer for your operating system using one of the links above, open your terminal, navigate to where you downloaded the file, and run one of the following commands. For Linux:
```bash Miniconda3-latest-Linux-x86_64.sh```

For MacOS:
```bash Miniconda3-latest-MacOSX-x86_64.sh```

3. Create a conda environment for the program

In the command line, navigate to the folder you want to place the environment. Then run the command:
```conda create --name sgRNA```

4. Install biopython into the sgRNA environment (see https://biopython.org/)

In the command line, navigate to the same folder and run the command:
```conda install biopython -c conda-forge -n sgRNA```

5. Install pysam into the sgRNA environment (see https://pysam.readthedocs.io/en/latest/)

In the command line, navigate to the same folder and run the command:
```conda install pysam -c bioconda -n sgRNA```

6. Install bowtie into the sgRNA environment

Finally, after ensuring Conda is installed (either Miniconda or Anaconda), open your terminal and run:
```conda install bowtie -c bioconda -n sgRNA```

You are done! If all above packages are succesfully installed, you should now be able to run the tool whenever you have the sgRNA conda environment activated.

TOOL USAGE INFORMATION
=====================
This program can only run with the required dependencies, which must be manually activated through the conda sgRNA environment whenever the terminal is opened. This can be done by navigating to the same folder and running the command:
```conda activate sgRNA```
While the environment is active, each line of your terminal prompt should start with the text "(sgRNA)"; this means the tool can be run!

This program takes in (1) A single sgRNA sequence or multiple sgRNA sequences, 2. A genome/metagenome file, and 3. An output file, and returns compiled data on the possible sites each sgRNA sequence could bind to in the given genome/metagenome, classifying perfect matches as "on-target sites" and imperfect matches as "off-target sites." The program can be run in single-target mode (STI, for a single sgRNA input) or multi-target mode (MTI, for multiple sgRNA inputs). Further details on each input are below:
1. The sgRNA sequences are always 20 base pairs (the spacer sequence) or 23 base pairs (the spacer sequence + PAM) and can be either DNA or RNA. In single-target mode (STI), the sgRNA sequence is written as a string to the command line. In multi-target mode (MTI), the sgRNA sequences are given in the form of a compiled .CSV or .TSV file, where the program reads every cell of the given file, evaluates whether it is an sgRNA sequence, and analyzes every sequence found in the file, meaning it does not matter which rows/columns the sequences are placed in nor what other information is included in other cells, so long as only 1 sgRNA sequence is written in each cell.
2. The genome/metagenome can be any size and contain any number of individual sequences; this is what the sgRNA sequence(s) will be aligned against. The genome/metagenome is always given as a .FASTA file.
3. The output file is optional, but if given, should be a .CSV or .TSV file path. If the file already exists, it will be overwritten. The output will have 2 sections for each sgRNA sequence; the first, the on-target section, will have columns taking the form [On-Target + PAM (23 bp RNA), On-Target Sequence (35  bp DNA), On-Target Heuristic, On-Target Score], which lists each on-target for a given sgRNA sequence and gives a score evaluating how well the sgRNA would bind to it in isolation and a heuristic evaluating the same thing when considering the number of off-targets. The second section, the off-target section, will have columns taking the form [On-Target + PAM (23 bp RNA), Off-Target Sequence (35 bp DNA), Off-Target Score, Off-Target Count], which lists each off-target for a given sgRNA sequence and gives a score evaluating how well the sgRNA would bind to it in isolation as well as a counter for the number of times that off-target appears in the genome.


TOOL USAGE EXAMPLES
=====================
This program can be run with either of the below input schemes:

Single target input: The program takes in a string representing an RNA spacer or DNA target and a .FASTA file path and saves the target analysis data to the saveFile path. Example:
```python SpacerSequenceHandler.PY STI TGACTGACTGACTGACTGAC metaGenome.FASTA saveFile.CSV```

Multi target input: The program takes in a .CSV or .TSV file path containing spacer/target sequences and a .FASTA file path and saves a list of target analysis data to the saveFile path. Example:
```python SpacerSequenceHandler.PY MTI targetSequences.CSV metaGenome.FASTA saveFile.CSV```

Help: Print out a general or more detailed help message. Examples:
```python SpacerSequenceHandler.py help```
```python SpacerSequenceHandler.py help single-target-input```
```python SpacerSequenceHandler.py help multi-target-input```
