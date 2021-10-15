from Bio import SeqIO
import SpacerSequence


class MetaGenome:
    __Name = ""
    __Sequences = []

    """Initialization function that takes in the local filename of a metaGenome (and, optionally, a name for that 
    metagenome) and instantiates."""
    def __init__(self, metaGenome, name="Generic MetaGenome"):
        self.__Name = name
        if self.isValidFasta(metaGenome):
            print("The given metaGenome file is not a valid .fasta string file name.")
        else:
            self.addSequences(metaGenome)

    """Getter method that returns the name of the MetaGenome."""
    def getName(self):
        return self.__Name

    """Method that takes an input String spacerSequence and uses the findOffTargets method of each Sequence in 
    __Sequences, compiles them, and returns."""
    def findOffTargets(self, spacerSequence):
        offTargets = []
        for sequence in self.__Sequences:
            if not offTargets:
                offTargets = sequence.findOffTargets(spacerSequence)
            else:
                offTargets = [offTargets, sequence.findOffTargets(spacerSequence)]
        return offTargets

    """Method that takes in the local filename of a sequence file and adds sequences directly from it to __Sequences."""
    def addSequences(self, sequencesFile):
        for record in SeqIO.parse(sequencesFile, "fasta"):
            if isinstance(record, SeqIO.SeqRecord):
                print("Successful SeqRecord opened!")
            self.__Sequences.append(record)

    """Returns the sequence at the given index in __Sequences, ensuring the validity of the index."""
    def getSequenceAtIndex(self, index):
        if index < 0 or index >= len(self.__Sequences):
            print("Invalid given index of " + index + "for __Sequences size of " + len(self.__Sequences))
            return None
        else:
            return self.__Sequences[index]

    """Returns the first sequence containing the given string in __Sequences, ensuring the validity of the sequence."""
    def getSequenceWithSubsequence(self, subsequence):
        if not SpacerSequence.isValidDNA(subsequence) and not SpacerSequence.isValidRNA(subsequence):
            print("Sequence entered is neither valid RNA nor valid DNA.")
            return None
        else:
            for sequence in self.__Sequences:
                if not sequence.__Record.seq.find(subsequence) == -1:
                    return sequence
            print("Sequence not find in __Sequences.")
            return None

    """Takes in a filename and returns whether the given variable is a valid fasta filename"""
    def isValidFasta(self, sequencesFile):
        return isinstance(sequencesFile,str) and ".fasta" in sequencesFile and \
               str.lower(sequencesFile[-6:]) == ".fasta" and len(sequencesFile) > 6


