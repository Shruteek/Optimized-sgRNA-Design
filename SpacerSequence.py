from MetaGenome import MetaGenome
from Sequence import Sequence

class SpacerSequence:
    """A class where each instance represents a unique 35-base pair sequence on an sgRNA. The class contains instance
    variables and methods that analyze the on-target and off-target efficiency and effects of the sequence."""
    HsuMatrix = [[1.6533, 0.9030, 1.5977, 0.9235, 0.8070, 0.9632, 1.0163, 0.2658, 0.7119, 1.5211, 0.6942, 1.0434, 0.5255, 0.8981, 0.7164, 0.8399, 0.5376, 0.2821, 0.6898],
                 [1.5142, 1.1597, 1.6582, 0.9924, 0.0247, 0.5522, 1.8687, 0.7737, 0.9270, 0.7292, 0.4842, 0.4824, 0.7060, 1.0221, 0.0181, 0.3496, 0.1811, 0.1362, 0.2700],
                 [1.3234, 1.4157, 1.2967, 1.2060, 0.9524, 0.2304, 1.0163, 0.8100, 1.1559, 0.7075, 1.5791, 0.3490, 0.0899, 0.0497, 0.0045, 0.2267, 0.2153, 0.5250, 0.4965],
                 [1.5366, 1.2771, 1.2689, 1.2197, 1.0645, 0.7791, 1.2445, 0.9885, 1.5319, 0.7184, 1.7381, 0.4166, 0.1285, 0.0720, 0.0549, 0.2261, 0.3119, 0.1343, 0.0601],
                 [1.7347, 0.8215, 1.1579, 1.1816, 0.7380, 0.9004, 0.8368, 0.2997, 0.6210, 1.1400, 0.3561, 0.6192, 0.1799, 0.2665, 0.2793, 0.2613, 0.1152, 0.1680, 0.3372],
                 [1.0186, 0.9649, 2.2504, 1.6222, 0.2405, 0.7561, 1.0651, 0.1102, 1.4293, 0.3533, 0.6178, 0.7269, 0.1165, 0.0367, 0.5013, 0.4147, 0.1786, 0.5315, 0.1664],
                 [1.6719, 0.9688, 1.0732, 1.0869, 0.6475, 1.0142, 0.8635, 0.3059, 0.4487, 0.9046, 0.4327, 0.5576, 0.1379, 0.0722, 0.3279, 0.2420, 0.0433, 0.1351, 0.4403],
                 [1.1662, 0.4544, 2.7867, 1.0461, 0.6036, 0.8132, 0.7875, 0.6882, 1.3655, 0.1240, 0.1953, 0.2497, 0.0132, 0.0227, 0.0478, 0.3682, 0.3175, 0.5621, 0.4588],
                 [1.1916, 1.0954, 2.8709, 1.1310, 0.5160, 0.6439, 1.0322, 0.5356, 1.2868, 0.0780, 0.2592, 0.2675, 0.0469, 0.0252, 0.0052, 0.0218, 0.1718, 0.6970, 0.2720],
                 [1.4786, 1.0820, 1.2952, 0.7450, 0.9763, 0.4912, 0.9272, 0.6022, 1.0375, 0.3047, 0.8210, 0.0303, 0.0365, 0.0592, 0.0253, 0.1553, 0.1006, 0.2175, 0.0275],
                 [0.0400, 0.9954, 1.6466, 1.3410, 0.0102, 0.5428, 2.3401, 0.4367, 0.2143, 0.3405, 0.2640, 0.0935, 0.0462, 0.0688, 0.0165, 0.3659, 0.0546, 0.0857, 0.2543],
                 [0.0345, 1.0478, 1.0507, 1.4075, 0.0540, 0.6396, 2.0810, 0.4585, 0.1555, 0.1369, 0.1026, 0.0417, 0.0105, 0.0458, 0.0099, 0.2114, 0.0552, 0.0253, 0.0596]]
    MismatchToHsuIndexDict = {'G->T': 0,'A->C': 1,'G->G': 2,'T->G': 3,'T->T': 4,'C->A': 5,'C->T': 6,'G->A': 7,'A->A': 8,'A->G': 9,'T->C': 10,'C->C': 11}
    NucleotideFeaturesDict = {'AA19': -0.0974, 'TT18': -0.0944, 'TT13': -0.0862, 'CT26': -0.0843, 'GC25': -0.0735, 'T21': -0.0687, 'TG23': -0.0664, 'AG23': -0.0543,
                              'G30': -0.0463, 'A4': -0.0422, 'AG34': -0.0419, 'GA34': -0.0378}
    """Initialization method that takes in a core sequence and an associated metaGenome class file and instantiates,
    checking to ensure the coreSequence is in the format 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp 
    downstream, as well as that genome is a MetaGenome."""
    def __init__(self, coreSequence, genome):
        if not self.isValidCoreSequence(coreSequence):
            self.__guideSequence = ""
        else:
            self.__guideSequence = coreSequence
        if not isinstance(genome, MetaGenome):
            print("Given genome is not a valid MetaGenome object.")
            self.__offTargetSequences = []
        else:
            self.__offTargetSequences = genome.findOffTargets()
        self.__oonTargetScore = self.__calcOnTargetScore(self.complementaryDNA(self.__guideSequence))
        self.__offTargetScores = []
        for candidate in self.__offTargetSequences:
            self.__offTargetScores.append(self.__calcOffTargetScore(candidate))
        self.__heuristic = self.__calcHeuristic()

    """Method that returns a boolean representing whether the input sequence is a valid DNA sequence."""
    def isValidDNA(self, DNASequence):
        if not isinstance(DNASequence, str):
            return False
        validDNABasePairs = "ACTGN"
        for nucleotide in DNASequence:
            if not validDNABasePairs.__contains__(nucleotide):
                return False
        return True

    """Method that returns a boolean representing whether the input sequence is a valid RNA sequence."""
    def isValidRNA(self, RNASequence):
        if not isinstance(RNASequence, str):
            return False
        validRNABasePairs = "ACUG"
        for nucleotide in RNASequence:
            if not validRNABasePairs.__contains__(nucleotide):
                return False
        return True

    """Method that returns whether a given sequence is a valid string representing a 35 base-pair RNA sequence in the 
    format 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream"""
    def isValidCoreSequence(self, coreSequence):
        if not isinstance(coreSequence, str):
            print("Given sequence is not a valid string: " + coreSequence)
            return False
        elif not self.isValidRNA(coreSequence):
            print("Given sequence is not valid RNA: " + coreSequence)
        elif not (len(coreSequence) == 35 and coreSequence[27:28] == "GG"):
            print("Given RNA sequence is not valid RNA: " + coreSequence)
            return False
        else:
            return True



    """Method that returns the complementary RNA strand sequence to the given input DNA sequence."""
    def complementaryRNA(self, sequence):
        complement = ""
        for nucleotide in sequence:
            if nucleotide == "A":
                complement = complement + "U"
            elif nucleotide == "C":
                complement = complement + "G"
            elif nucleotide == "T":
                complement = complement + "A"
            elif nucleotide == "G":
                complement = complement + "C"
        return complement

    """Method that returns the complementary DNA strand sequence to the given input DNA or RNA sequence."""
    def complementaryDNA(self, sequence):
        complement = ""
        for nucleotide in sequence:
            if nucleotide == "A":
                complement = complement + "T"
            elif nucleotide == "C":
                complement = complement + "G"
            elif nucleotide == "T":
                complement = complement + "A"
            elif nucleotide == "U":
                complement = complement + "A"
            elif nucleotide == "G":
                complement = complement + "C"
        return complement

    """Method that returns the calculated on-target score of the given target sequence, assuming a perfectly 
    complementary RNA guide sequence on the sgRNA. """
    def __calcOnTargetScore(self, targetSequence):
        if len(targetSequence) == 0:
            return 0
        return 1.0

    """Method that returns the calculated off-target score of the given sequence with respect to to the instance RNA 
    guide sequence. """
    def __calcOffTargetScore(self, targetSequence):
        return 0.0

    """Method that takes in an on-target score and a list of off-target scores, and calculates a general heuristic to 
    measure the effectiveness of the core sequence. """
    def __calcHeuristic(self):
        heuristic = self.__onTargetScore
        for offTargetScore in self.__offTargetScores:
            heuristic = heuristic - offTargetScore
        return heuristic

    """Getter method that returns the calculated on target score."""
    def getOnTargetScore(self):
        return self.__onTargetScore

    """Getter method that returns the calculated off target scores."""
    def getOffTargetScore(self):
        return self.__offTargetScores

    """Getter method that returns the calculated heuristic of the coreSequence."""
    def getHeuristic(self):
        return self.__heuristic