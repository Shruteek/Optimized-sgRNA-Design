from Bio import SeqIO
from GenomeTools import *
from Sequence import Sequence
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
        if not (len(spacerSequence) == 20 and isValidRNA(spacerSequence)):
            return foundTargets
        else:
            sequence_to_align = spacerSequence + "NGG"
            indexName = os.path.splitext(os.path.basename(self.__OriginalPath))[0]
            projectPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
            outputPath = os.path.join(projectPath, "Outputs")
            if (not os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwt"))) and (not os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwtl"))):
                print("Index " + indexName + " does not exist. Building...")
                os.system("bowtie-build " + self.__OriginalPath + " " + os.path.join(outputPath,
                                                                         indexName) + " > /dev/null")
            print("Aligning spacer sequence: " + convertToDNA(sequence_to_align))
            if os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwt")):
                os.system("bowtie -a -v 3 --np 0 " + os.path.join(outputPath,
                                                           indexName) + " -c " + convertToDNA(sequence_to_align) + " -S " + os.path.join(
                    outputPath, indexName + convertToDNA(sequence_to_align) + ".sam"))
            elif os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwtl")):
                os.system("bowtie -a -v 3 --np 0 --large-index " + os.path.join(outputPath,
                                                                         indexName) + " -c " + convertToDNA(sequence_to_align) + " -S " + os.path.join(
                    outputPath, indexName + convertToDNA(sequence_to_align) + ".sam"))
            else:
                print("Failed to find and to build index.")
            fastaFile = pysam.FastaFile(self.__OriginalPath)
            alignmentFile = pysam.AlignmentFile(os.path.join(outputPath, indexName + convertToDNA(sequence_to_align) + ".sam"))
            alignments = 0
            for alignedSegment in alignmentFile:
                # We check if the aligned segment matches all 23 base pairs (spacer + PAM) with its reference sequence.
                if alignedSegment.is_mapped and alignedSegment.cigarstring == "23M":
                    referenceSequence = fastaFile.fetch(reference=alignedSegment.reference_name)
                    for alignedBlock in alignedSegment.get_blocks():
                        # The below conditional checks if there are 6 bps around the aligned sequence - we need that
                        # much spacing for our on-target analysis.
                        if alignedBlock[0] >= 6 and alignedBlock[1] <= len(referenceSequence) - 6:
                            alignedReferenceSequence = referenceSequence[alignedBlock[0]:alignedBlock[1]]
                            if alignedReferenceSequence[-3:] == convertToDNA(sequence_to_align)[-3:] \
                                    and correlateSequences(alignedReferenceSequence, convertToDNA(sequence_to_align)) >= 18:
                                fullTargetSequence = referenceSequence[alignedBlock[0] - 6:alignedBlock[1] + 6]
                                foundTargets.append(fullTargetSequence)
                                print("Aligned target " + fullTargetSequence + " versus spacer " + alignedSegment.get_forward_sequence())
                                alignments += 1
                            elif reverseComplementaryDNA(alignedReferenceSequence)[-3:] == convertToDNA(sequence_to_align)[-3:] \
                                    and correlateSequences(reverseComplementaryDNA(alignedReferenceSequence), convertToDNA(sequence_to_align)) >= 18:
                                fullTargetSequence = reverseComplementaryDNA(referenceSequence[alignedBlock[0] - 6:alignedBlock[1] + 6])
                                targetSequenceInfo = [fullTargetSequence, alignedSegment.reference_name]
                                foundTargets.append(targetSequenceInfo)
                                print("Aligned target " + fullTargetSequence + " versus spacer " + alignedSegment.get_forward_sequence())
                                alignments += 1
            fastaFile.close()
            alignmentFile.close()
        print("Identified " + alignments + " alignments for spacer sequence: " + convertToDNA(spacerSequence))
        return foundTargets

    def findTargetsFromSpacer(self, spacerSequence):
        """Method that takes an input String RNA spacerSequence and uses the findTargetsFromSpacer method of each
        Sequence in __Sequences, compiles them, and returns."""
        nestedOffTargets = []
        if not (len(spacerSequence) == 20 and isValidRNA(spacerSequence)):
            nestedOffTargets = []
        else:
            for sequence in self.__Sequences:
                offTargets = sequence.findTargetsFromSpacer(spacerSequence)
                for offTarget in offTargets:
                    nestedOffTargets.append(offTarget)
        return nestedOffTargets

    # def addSequences(self, sequencesFile):
    #     """Method that takes in the filename of a sequence file and adds sequences directly from it to __Sequences."""
    #     count = 1
    #     with open(sequencesFile) as fileHandler:
    #         for record in SeqIO.parse(fileHandler, "fasta"):
    #             self.__Sequences.append(Sequence(record, "Sequence " + str(count)))
    #             count = count + 1

