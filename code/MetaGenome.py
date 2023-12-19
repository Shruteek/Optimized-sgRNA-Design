from Bio import SeqIO
from GenomeTools import *
import os
import pysam


class MetaGenome:

    def __init__(self, metaGenomePath, name="Generic MetaGenome"):
        """Initialization method that takes in the local .FASTA filename of a metaGenome (and, optionally, a name for that
            metagenome) and instantiates."""
        self.__OriginalPath = os.path.abspath(metaGenomePath)
        self.__Name = name
        self.__Sequences = []
        if not isValidFasta(metaGenomePath):
            print("The given metaGenome file is not a valid .fasta string file name:" + metaGenomePath)

    def getName(self):
        """Getter method that returns the name of the MetaGenome."""
        return self.__Name

    def setName(self, name):
        """Takes in a string name and sets the Sequence name accordingly."""
        if isinstance(name, str) and len(name) > 0:
            self.__Name = name

    def size(self):
        """Method that returns the size of the Sequence array."""
        return len(self.__Sequences)

    def getSequence(self, key):
        """Returns a list containing either the string sequence at the given key index or the string sequences with the
            given string key sequence in __Sequences, ensuring the validity of the key input as either an index or a string
            sequence."""
        if isinstance(key, str):
            if not (isValidDNA(key) or isValidRNA(key)):
                # print("Subsequence entered is neither valid RNA nor valid DNA.")
                return []
            else:
                matchingSequences = []
                for sequence in self.__Sequences:
                    if sequence.hasSubsequence(key):
                        matchingSequences.append(sequence.getSequence())
                return matchingSequences
        elif isinstance(key, int):
            if key < 0 or key >= len(self.__Sequences):
                # print("Invalid given index of " + key + "for __Sequences size of " + len(self.__Sequences))
                return []
            else:
                return [self.__Sequences[key].getSequence()]
        else:
            # print("Invalid sequence key: must be int (index) or str (subsequence).")
            return []

    def getOriginalPath(self):
        """Returns the original path of the metagenome file as a string."""
        return self.__OriginalPath

    def bowtieFindTargetsFromSpacer(self, spacerSequence):
        """Method that takes an input 20 bp String RNA spacerSequence and uses bowtie to run alignment analysis on it
        with respect to the metagenome. Specified flags include: -a (include all alignments), -v 5 (limit to maximum
        of 5 mismatches), -n 1 (limit to maximum of 1 mismatch in seed sequence), -l 10 (assume seed sequence is 10
        base pairs), -c <SpacerSequence> (take direct sequence input rather than file), --np 0 (no penalty for N
        nucleotides in either the aligning or the reference sequence).

        The method returns the list 'foundTargets', which is populated by 2-element lists where each first element is
        a 35 base-pair DNA target sequence (either on-target or off-target) and the second element is its respective
        location in the reference metagenome."""
        foundTargets = []
        if not (len(spacerSequence) == 20 and (isValidRNA(spacerSequence) or isValidDNA(spacerSequence))):
            return foundTargets
        else:
            sequence_to_align = convertToDNA(spacerSequence)
            indexName = os.path.splitext(os.path.basename(self.__OriginalPath))[0]
            projectPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
            outputPath = os.path.join(projectPath, "Outputs")
            genedataPath = os.path.join(projectPath, "Genedata")
            bowtieIndexPath = os.path.join(genedataPath, indexName)
            outputFilePath = os.path.join(outputPath, indexName + "_" + sequence_to_align + ".sam")
            if (not os.path.exists(bowtieIndexPath + ".rev.2.bt2")) and (not os.path.exists(bowtieIndexPath + ".rev.2.bt2l")):
                print("Index " + indexName + " does not exist. Building...")
                os.system("bowtie2-build " + self.__OriginalPath + " " + bowtieIndexPath)
                print()
            bowtie2Options = " -a --mp 1,1 -L 5 -N 0 --np 0 --score-min L,-3,0 "
            if os.path.exists(bowtieIndexPath + ".rev.2.bt2") or os.path.exists(bowtieIndexPath + ".rev.2.bt2l"):
                print("Aligning spacer sequence " + sequence_to_align + " to Bowtie2 index at " + bowtieIndexPath)
                os.system("bowtie2 " + bowtie2Options
                          + " -x " + bowtieIndexPath
                          + " -c " + sequence_to_align
                          + " -S " + outputFilePath)
            else:
                print("Failed to find and to build index.")
            # Now that we have aligned the FASTA file into a SAM file, let's check each alignment.
            fastaFile = pysam.FastaFile(self.__OriginalPath)
            alignmentFile = pysam.AlignmentFile(outputFilePath)
            alignments = 0
            for alignedSegment in alignmentFile:
                # We check if the aligned segment lines up with all 20 base pairs of the spacer.
                if alignedSegment.is_mapped and alignedSegment.cigarstring == "20M":
                    # If so, we fetch the full reference sequence that contains the target.
                    referenceSequence = fastaFile.fetch(reference=alignedSegment.reference_name)
                    for alignedBlock in alignedSegment.get_blocks():
                        # The below conditional checks if there are 6 bps around the aligned sequence - we need that
                        # much spacing for our on-target analysis.
                        if alignedBlock[0] >= 6 and alignedBlock[1] <= len(referenceSequence) - 9:
                            # If so, let's grab the target sequence, the 3 bp PAM, 6 bp upstream, and 6 bp downstream.
                            surroundingTargetSequence = referenceSequence[alignedBlock[0] - 6:alignedBlock[1] + 9]
                            # We check if the target sequence has a PAM ('NGG') after the spacer, and if there are at
                            # most 5 mismatches
                            if surroundingTargetSequence[27:29] == "GG" \
                                    and sum(n1 != n2 for n1, n2 in zip(surroundingTargetSequence[6:26],sequence_to_align)) <= 5:
                                targetSequenceInfo = [surroundingTargetSequence, alignedSegment.reference_name]
                                foundTargets.append(targetSequenceInfo)
                                print("Aligned spacer " + sequence_to_align + " versus target " + surroundingTargetSequence)
                                alignments += 1
                            else:
                                # Otherwise, we check for the reverse complement of the target sequence:
                                surroundingTargetSequence =  referenceSequence[alignedBlock[0] - 9:alignedBlock[1] + 6]
                                surroundingTargetSequence = reverseComplementaryDNA(surroundingTargetSequence)
                                if surroundingTargetSequence[27:29] == "GG" \
                                        and sum(n1 != n2 for n1, n2 in zip(surroundingTargetSequence[6:26],sequence_to_align)) <= 5:
                                    targetSequenceInfo = [surroundingTargetSequence, alignedSegment.reference_name]
                                    foundTargets.append(targetSequenceInfo)
                                    print("Aligned spacer " + sequence_to_align + " versus target " + surroundingTargetSequence)
                                    alignments += 1
            fastaFile.close()
            alignmentFile.close()
        print("Identified " + alignments + " alignments for spacer sequence: " + convertToDNA(spacerSequence))
        return foundTargets


