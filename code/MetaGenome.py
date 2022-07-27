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
        else:
            self.addSequences(metaGenomePath)

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
        """Method that takes an input 23 bp String RNA spacerSequence and uses bowtie to run alignment analysis on it with
        respect to the metagenome. Specified flags include: -a (include all alignments), -v 5 (limit to maximum of 5
        mismatches), -n 1 (limit to maximum of 1 mismatch in seed sequence), -l 10 (assume seed sequence is 10 base pairs), -c <SpacerSequence> (take direct sequence input rather than file), -B"""
        foundTargets = []
        if not (len(spacerSequence) == 23 and isValidRNA(spacerSequence)):
            return foundTargets
        else:
            indexName = os.path.splitext(os.path.basename(self.__OriginalPath))[0]
            projectPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
            outputPath = os.path.join(projectPath, "Outputs")
            if (not os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwt"))) and (not os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwtl"))):
                print("Index " + indexName + " does not exist. Building...")
                os.system("bowtie-build " + self.__OriginalPath + " " + os.path.join(outputPath,
                                                                         indexName) + " > /dev/null")

            if os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwt")):
                os.system("bowtie -a -v 3 " + os.path.join(outputPath,
                                                           indexName) + " -c " + convertToDNA(spacerSequence) + " -S " + os.path.join(
                    outputPath, indexName + convertToDNA(spacerSequence) + ".sam"))
            elif os.path.exists(os.path.join(outputPath, indexName + ".rev.2.ebwtl")):
                os.system("bowtie -a -v 3 --large-index " + os.path.join(outputPath,
                                                                         indexName) + " -c " + convertToDNA(spacerSequence) + " -S " + os.path.join(
                    outputPath, indexName + convertToDNA(spacerSequence) + ".sam"))
            else:
                print("Failed to find and to build index.")

            print("Aligning spacer sequence: " + convertToDNA(spacerSequence))
            fastaFile = pysam.FastaFile(self.__OriginalPath)
            alignmentFile = pysam.AlignmentFile(os.path.join(outputPath, indexName + convertToDNA(spacerSequence) + ".sam"))
            for alignedSegment in alignmentFile.head(10000):
                if alignedSegment.is_mapped:
                    print("Cigarstring: " + alignedSegment.cigarstring)
                    print("Aligned " + alignedSegment.get_forward_sequence() + " versus " + alignedSegment.get_reference_sequence())
                if alignedSegment.is_mapped and alignedSegment.cigarstring == "23M":
                    referenceSequence = fastaFile.fetch(reference=alignedSegment.reference_name)
                    for alignedBlock in alignedSegment.get_blocks():
                        if alignedBlock[0] >= 6 and alignedBlock[1] <= len(referenceSequence) - 6 and referenceSequence[alignedBlock[1] - 3:alignedBlock[1]] == convertToDNA(spacerSequence)[-3:]:
                            fullTargetSequence = referenceSequence[alignedBlock[0] - 6:alignedBlock[1] + 6]
                            foundTargets.append(fullTargetSequence)
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

    def addSequences(self, sequencesFile):
        """Method that takes in the filename of a sequence file and adds sequences directly from it to __Sequences."""
        count = 1
        with open(sequencesFile) as fileHandler:
            for record in SeqIO.parse(fileHandler, "fasta"):
                self.__Sequences.append(Sequence(record, "Sequence " + str(count)))
                count = count + 1

