import math
from GenomeTools import *
from MetaGenome import MetaGenome
import time
import os


class TargetSequence:
    """A class where each instance represents a single target sequence (either an off-target or an on-target),
    including its target type (on vs. off), its surrounding 35 bp sequence, how many mismatches it has, its
    location in the metagenome, and its evaluation scores with respect to its SpacerSequence instance."""

    def __init__(self, spacer_sequence_object, full_sequence, metagenome_location):
        self.SpacerSequenceObject = spacer_sequence_object
        self.surrounding_sequence = full_sequence
        self.location = metagenome_location

        self.spacer_to_align = self.SpacerSequenceObject.__spacerSequence
        self.target_sequence = self.surrounding_sequence[6:26]
        if self.target_sequence == convertToDNA(self.spacer_to_align):
            self.sequence_type = "On-Target"
        else:
            self.sequence_type = "Off-Target"
        self.mismatches = self.countMismatches()
        self.off_target_score = None
        self.on_target_score = None

    def countMismatches(self):
        spacer_DNA = convertToDNA(self.spacer_to_align)
        mismatches = sum(1 for i in range(len(spacer_DNA)) if spacer_DNA[i] != self.target_sequence[i])
        return mismatches

    def calculateScores(self):
        self.on_target_score = self.SpacerSequenceObject.calcOnTargetScore(self.surrounding_sequence)
        if self.sequence_type == "On-Target":
            self.off_target_score = None
        else:
            if len(self.SpacerSequenceObject.__onTargetScores) == 0:
                self.off_target_score = self.SpacerSequenceObject.calcOffTargetEstimate(self.surrounding_sequence)
            else:
                self.off_target_score = self.SpacerSequenceObject.calcOffTargetScore(self.surrounding_sequence,
                                                                                     max(self.SpacerSequenceObject.__onTargetScores))


class SpacerSequence:
    """A class where each instance represents a unique 20-base pair sequence on an sgRNA. The class contains instance
    variables and methods that analyze the on-target and off-target efficiency and effects of the sequence.
    TERMS KEY:
    Guide sequence: 35 base-pair 3'-5' RNA sequence with 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream
    Target guide sequence: 35 base-pair 5'-3' DNA sequence with 6 bp upstream, 20 bp target sequence, 3 bp PAM-match,
        and 6 bp downstream
    Spacer sequence: 20 base-pair RNA sequence that complements the 20 DNA nucleotide target sequence on the genome
    Seed sequence: 10-20 base-pair RNA sequence at the beginning of the spacer sequence, extending partway or all the
    way through the sequence, and representing the most important RNA nucleotides for binding to the target sequence .
    Target sequence: 35 base-pair DNA sequence with a 20 base-pair subsequence that matches the 20 RNA nucleotide
    spacer sequence on the sgRNA.
    PAM: 3 base-pair DNA sequence on the guide sequence in the form 5'-NGG-3'
    """
    HsuMatrix = [
        [1.6533, 0.9030, 1.5977, 0.9235, 0.8070, 0.9632, 1.0163, 0.2658, 0.7119, 1.5211, 0.6942, 1.0434, 0.5255, 0.8981,
         0.7164, 0.8399, 0.5376, 0.2821, 0.6898],
        [1.5142, 1.1597, 1.6582, 0.9924, 0.0247, 0.5522, 1.8687, 0.7737, 0.9270, 0.7292, 0.4842, 0.4824, 0.7060, 1.0221,
         0.0181, 0.3496, 0.1811, 0.1362, 0.2700],
        [1.3234, 1.4157, 1.2967, 1.2060, 0.9524, 0.2304, 1.0163, 0.8100, 1.1559, 0.7075, 1.5791, 0.3490, 0.0899, 0.0497,
         0.0045, 0.2267, 0.2153, 0.5250, 0.4965],
        [1.5366, 1.2771, 1.2689, 1.2197, 1.0645, 0.7791, 1.2445, 0.9885, 1.5319, 0.7184, 1.7381, 0.4166, 0.1285, 0.0720,
         0.0549, 0.2261, 0.3119, 0.1343, 0.0601],
        [1.7347, 0.8215, 1.1579, 1.1816, 0.7380, 0.9004, 0.8368, 0.2997, 0.6210, 1.1400, 0.3561, 0.6192, 0.1799, 0.2665,
         0.2793, 0.2613, 0.1152, 0.1680, 0.3372],
        [1.0186, 0.9649, 2.2504, 1.6222, 0.2405, 0.7561, 1.0651, 0.1102, 1.4293, 0.3533, 0.6178, 0.7269, 0.1165, 0.0367,
         0.5013, 0.4147, 0.1786, 0.5315, 0.1664],
        [1.6719, 0.9688, 1.0732, 1.0869, 0.6475, 1.0142, 0.8635, 0.3059, 0.4487, 0.9046, 0.4327, 0.5576, 0.1379, 0.0722,
         0.3279, 0.2420, 0.0433, 0.1351, 0.4403],
        [1.1662, 0.4544, 2.7867, 1.0461, 0.6036, 0.8132, 0.7875, 0.6882, 1.3655, 0.1240, 0.1953, 0.2497, 0.0132, 0.0227,
         0.0478, 0.3682, 0.3175, 0.5621, 0.4588],
        [1.1916, 1.0954, 2.8709, 1.1310, 0.5160, 0.6439, 1.0322, 0.5356, 1.2868, 0.0780, 0.2592, 0.2675, 0.0469, 0.0252,
         0.0052, 0.0218, 0.1718, 0.6970, 0.2720],
        [1.4786, 1.0820, 1.2952, 0.7450, 0.9763, 0.4912, 0.9272, 0.6022, 1.0375, 0.3047, 0.8210, 0.0303, 0.0365, 0.0592,
         0.0253, 0.1553, 0.1006, 0.2175, 0.0275],
        [0.0400, 0.9954, 1.6466, 1.3410, 0.0102, 0.5428, 2.3401, 0.4367, 0.2143, 0.3405, 0.2640, 0.0935, 0.0462, 0.0688,
         0.0165, 0.3659, 0.0546, 0.0857, 0.2543],
        [0.0345, 1.0478, 1.0507, 1.4075, 0.0540, 0.6396, 2.0810, 0.4585, 0.1555, 0.1369, 0.1026, 0.0417, 0.0105, 0.0458,
         0.0099, 0.2114, 0.0552, 0.0253, 0.0596]]
    MismatchToHsuIndexDict = {'G->T': 0, 'A->C': 1, 'G->G': 2, 'T->G': 3, 'T->T': 4, 'C->A': 5, 'C->T': 6, 'G->A': 7,
                              'A->A': 8, 'A->G': 9, 'T->C': 10, 'C->C': 11}
    NucleotideFeaturesDict = {'AA19': -0.0974, 'TT18': -0.0944, 'TT13': -0.0862, 'CT26': -0.0843, 'GC25': -0.0735,
                              'T21': -0.0687, 'TG23': -0.0664, 'AG23': -0.0543,
                              'G30': -0.0463, 'A4': -0.0422, 'AG34': -0.0419, 'GA34': -0.0378, 'A18': -0.0338,
                              'C25': -0.0316, 'C31': -0.0307, 'G1': -0.0297,
                              'C16': -0.0216, 'A14': -0.0185, 'A11': -0.0183, 'T34': -0.0176, 'AA10': -0.0169,
                              'A19': -0.0156, 'G34': -0.0142, 'C30': -0.0132,
                              'GA31': -0.0123, 'T24': -0.0120, 'A15': -0.0106, 'G4': -0.0054, 'GG9': -0.0016,
                              'T23': -0.0014, 'C15': -0.0005, 'C26': -0.0004,
                              'T27': -0.0003, 'A31': 0.0016, 'GT18': 0.0024, 'C9': 0.0024, 'GA20': 0.0097,
                              'A25': 0.0105, 'A12': 0.0116, 'A32': 0.0124, 'T22': 0.0132,
                              'C20': 0.0151, 'G17': 0.0155, 'G18': 0.0165, 'T30': 0.0173, 'A13': 0.0176, 'G19': 0.0179,
                              'A27': 0.0191, 'G11': 0.0209, 'TG3': 0.0229,
                              'GC3': 0.0247, 'G14': 0.0251, 'GG10': 0.0268, 'G12': 0.0276, 'G32': 0.0307, 'A22': 0.0319,
                              'G20': 0.0340, 'C21': 0.0343, 'TT17': 0.0349,
                              'T13': 0.0354, 'G26': 0.0361, 'A24': 0.0375, 'C22': 0.0376, 'G16': 0.0380, 'GG12': 0.0419,
                              'TG18': 0.0459, 'TG31': 0.0481, 'A35': 0.0486,
                              'G15': 0.0511, 'C24': 0.0530, 'TG15': 0.0534, 'GT11': 0.0537, 'GC9': 0.0542,
                              'CA30': 0.0578, 'GT24': 0.0610, 'G13': 0.0614, 'CA24': 0.0622,
                              'AG10': 0.0637, 'G10': 0.0677, 'C13': 0.0695, 'GT31': 0.0734, 'GG13': 0.0744,
                              'C27': 0.0799, 'G27': 0.0852, 'CC21': 0.0889, 'CC23': 0.0951,
                              'G22': 0.1011, 'G24': 0.1055, 'GT23': 0.1067, 'GG25': 0.1116, 'G9': 0.1146}
    PAMDensityScoreMatrix = [[1, 1, 1, 2, 2, 8, 10],
                             [1, 1, 1, 2, 2, 8, 10],
                             [1, 1, 1, 2.5, 2.5, 10, -1],
                             [1, 1, 1.66, 2.5, 2.5, -1, -1],
                             [1.66, 1.66, 3.33, 3.33, -1, -1, -1],
                             [2.5, 2.5, 3.33, -1, -1, -1, -1],
                             [5, 5, -1, -1, -1, -1, -1]]

    def __init__(self, spacerOrTargetSequence, genome):
        """Initialization method that takes in a 23 bp DNA spacer/target or a 20 bp DNA spacer/target sequence,
        and (2) an associated metaGenome class file, and instantiates, checking to ensure the guideSequence is in the
        proper format (6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream), as well as that given genome
        is a MetaGenome, trying to complete the guideSequence using the genome otherwise. """
        self.__spacerSequence = ""  # RNA string of length 23
        self.__onTargetSequences = []  # List of length A with DNA strings of length 35
        self.__onTargetScores = []  # List of length A with doubles
        self.__targets = [] # List of TargetSequence class objects, some on-targets, some off-targets
        if not isinstance(genome, MetaGenome):
            return
        if not isValidTargetSpacerInput(spacerOrTargetSequence):
            print("Input sequence is not a valid spacer sequence: " + spacerOrTargetSequence)
            return
        self.__spacerSequence = convertToDNA(spacerOrTargetSequence[0:20])
        print("DEBUG: Finding targets from spacer for sequence " + self.__spacerSequence)
        targetSequencesInfo = genome.bowtieFindTargetsFromSpacer(self.__spacerSequence)
        print("DEBUG: Found " + str(len(targetSequencesInfo)) + " targets from spacer for sequence " + self.__spacerSequence)
        for targetSequenceInfo in targetSequencesInfo: # Instantiate all target sequences and calculate on-target scores
            targetInstance = TargetSequence(self, targetSequenceInfo[0], targetSequenceInfo[1])
            if targetInstance.sequence_type == "On-Target":
                targetInstance.calculateScores()
                self.__onTargetScores.append(targetInstance.on_target_score)
            self.__targets.append(targetInstance)

        for targetSequence in self.__targets: # Now calculate off-target scores (which require an on-target)
            if targetInstance.sequence_type == "Off-Target":
                targetInstance.calculateScores()

    def calcOnTargetScore(self, targetSequence):
        """Method that takes in a 35 bp guide DNA sequence and returns the calculated on-target score of it, assuming a
        perfectly complementary RNA spacer sequence on the sgRNA."""
        if len(targetSequence) != 35:
            return 0
        CRISPRscanSubscore = 0
        reversePAMs = 0
        forwardPAMs = 0
        penaltyScore = 0
        GCADensityScore = 0
        for c in range(len(targetSequence)):
            if 5 < c < 25:
                if targetSequence[c:c + 1] == "GG":
                    forwardPAMs = forwardPAMs + 1
                if targetSequence[c:c + 1] == "CC":
                    reversePAMs = reversePAMs + 1
                if targetSequence[c] == "G":
                    GCADensityScore = GCADensityScore + 1.0
                elif targetSequence[c] == "C":
                    GCADensityScore = GCADensityScore + 0.5
                elif targetSequence[c] == "A":
                    GCADensityScore = GCADensityScore - 0.1
            if targetSequence[c] + str(c + 1) in self.NucleotideFeaturesDict.keys():
                CRISPRscanSubscore = CRISPRscanSubscore + self.NucleotideFeaturesDict[targetSequence[c] + str(c + 1)]
            if c < len(targetSequence) - 1 and (targetSequence[c:(c + 2)] + str(c + 1)) \
                    in self.NucleotideFeaturesDict.keys():
                CRISPRscanSubscore = CRISPRscanSubscore + self.NucleotideFeaturesDict[targetSequence[c:(c + 2)] +
                                                                                      str(c + 1)]
        CRISPRscanSubscore = (CRISPRscanSubscore + 0.6555) / 1.9495 * 100
        GCADensityScore = GCADensityScore / 20
        if self.PAMDensityScoreMatrix[reversePAMs][forwardPAMs] > 1:
            penaltyScore = self.PAMDensityScoreMatrix[reversePAMs][forwardPAMs] * GCADensityScore
        elif self.PAMDensityScoreMatrix[reversePAMs][forwardPAMs] == 1:
            penaltyScore = self.PAMDensityScoreMatrix[reversePAMs][forwardPAMs] - GCADensityScore / 5
        return CRISPRscanSubscore / penaltyScore

    def calcOffTargetScore(self, offTargetSequence, onTargetScore):
        """Method that takes in a 35 bp string DNA on-target sequence and a 35 bp string DNA potential off-target
        sequence and returns the calculated off-target score of the given off-target sequence with respect to the
        given DNA on-target sequence. """
        if len(offTargetSequence) != 35:
            return 1
        HsuMismatchSubscore = 1
        proximityMismatchSubscore = 0
        steppedProximityMismatchSubscore = 1
        offTargetSpacerSequence = offTargetSequence[6:26]
        for c in range(len(offTargetSpacerSequence)):
            if offTargetSpacerSequence[c] != self.__spacerSequence[c] and offTargetSpacerSequence[c] != "N":
                proximityMismatchSubscore = proximityMismatchSubscore + 1 / (20 - c)
                if c > 0:
                    indexDict = self.MismatchToHsuIndexDict[
                        self.__spacerSequence[c] + '->' + complementaryDNA(offTargetSpacerSequence[c])]
                    HsuMismatchSubscore = HsuMismatchSubscore * self.HsuMatrix[indexDict][c - 1]
                if 20 - c <= 6:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.1
                elif 20 - c <= 12:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.05
                else:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.0125
        proximityMismatchSubscore = (3.5477 - proximityMismatchSubscore) / 3.5477
        activityRatio = self.calcOnTargetScore(offTargetSequence) / onTargetScore
        return 200 * (math.sqrt(HsuMismatchSubscore) + steppedProximityMismatchSubscore) * (activityRatio ** 2) \
               * (proximityMismatchSubscore ** 6) / 4

    def calcOffTargetEstimate(self, offTargetSequence):
        """Method that takes in a 35 bp string DNA potential off-target sequence and returns an estimate of the
        calculated off-target score of the given off-target sequence (without using a reference guide sequence) if
        a potential RNA guide sequence containing the class spacer sequence were to attempt binding."""
        if len(offTargetSequence) != 35:
            return 1
        HsuMismatchSubscore = 1
        proximityMismatchSubscore = 0
        steppedProximityMismatchSubscore = 1
        offTargetSpacerSequence = offTargetSequence[6:26]
        for c in range(len(offTargetSpacerSequence)):
            if convertToDNA(self.__spacerSequence[c]) != offTargetSpacerSequence[c]:
                proximityMismatchSubscore = proximityMismatchSubscore + 1 / (20 - c)
                if c > 0:
                    mismatchIdentity = convertToDNA(self.__spacerSequence[c]) + '->' + \
                                       complementaryDNA(offTargetSpacerSequence[c])
                    if mismatchIdentity in self.HsuMatrix:
                        indexDict = self.MismatchToHsuIndexDict[mismatchIdentity]
                        HsuMismatchSubscore = HsuMismatchSubscore * self.HsuMatrix[indexDict][c - 1]
                if 20 - c <= 6:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.1
                elif 20 - c <= 12:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.05
                else:
                    steppedProximityMismatchSubscore = steppedProximityMismatchSubscore - 0.0125
        proximityMismatchSubscore = (3.5477 - proximityMismatchSubscore) / 3.5477
        return 200 * (math.sqrt(HsuMismatchSubscore) + steppedProximityMismatchSubscore) \
               * (proximityMismatchSubscore ** 6) / 4

    def __calcHeuristics(self):
        """Method that takes in an instance on-target score and an instance list of off-target scores, and calculates a
        general heuristic to measure the effectiveness of the guide sequence. """
        heuristic = 0
        return heuristic

    def getSpacerSequence(self):
        """Method that returns the 20 bp string RNA spacer sequence of the instance."""
        return self.__spacerSequence

    def getOnTargetSequences(self):
        """Method that returns the list of 35 bp string RNA guide sequences of the instance."""
        on_target_sequences = []
        for targetSequence in self.__targets:
            if targetSequence.sequence_type == "On-Target":
                on_target_sequences.append(targetSequence.target_sequence)
        return on_target_sequences

    def getOffTargetSequences(self):
        """Method that returns the 35 bp string DNA off-target sequences of the instance."""
        off_target_sequences = []
        for targetSequence in self.__targets:
            if targetSequence.sequence_type == "Off-Target":
                off_target_sequences.append(targetSequence.target_sequence)
        return off_target_sequences

    def getOnTargetScores(self):
        """Getter method that returns the calculated on target score."""
        on_target_scores = []
        for targetSequence in self.__targets:
            if targetSequence.sequence_type == "On-Target":
                on_target_scores.append(targetSequence.on_target_score)
        return on_target_scores

    def getOffTargetScores(self):
        """Getter method that returns the calculated off target scores."""
        off_target_scores = []
        for targetSequence in self.__targets:
            if targetSequence.sequence_type == "Off-Target":
                off_target_scores.append(targetSequence.off_target_score)
        return off_target_scores

    def getTargets(self):
        return self.__targets

    def getHeuristics(self):
        """Getter method that returns the calculated heuristic of the guideSequence."""
        return self.__heuristics
