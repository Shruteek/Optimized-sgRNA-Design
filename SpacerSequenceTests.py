import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
from GenomeTools import *


class SpacerSequenceTests(unittest.TestCase):
    def testInstantiation(self):
        testGuideSequence = "TTTGGAGAGCCAAGGATTCGGATTCTCGGCTTCCA"
        metaGen = MetaGenome("test.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen)
        print(newSpacer.getOnTargetScore())
        print(newSpacer.getOffTargetScores())
        self.assertTrue(isinstance(newSpacer.getOnTargetScore(), float))
        self.assertTrue(isinstance(newSpacer.getOffTargetScores(), list))
        self.assertEqual(newSpacer.getOffTargetScores(), [])

        metaGen2 = MetaGenome("test2.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen2)
        self.assertTrue(isinstance(newSpacer.getOnTargetScore(), float))
        self.assertTrue(isinstance(newSpacer.getOffTargetScores(), list))
        self.assertEqual(newSpacer.getOffTargetScores(), [])

    def testCalcScores(self):
        testGuideSequence = "TTTGGAGAGCCAAGGATTCGGATTCTCGGCTTCCA"
        knownGuideSequence = "CCGGCGTTTTCCCGGAAGACAAAGTCGGGTTCATA"
        metaGen3 = MetaGenome("pf3_illumina_subset.fasta")
        newSpacer = SpacerSequence(complementaryRNA(knownGuideSequence), metaGen3)
        print(newSpacer.getOffTargetScores())
        self.assertEqual(1, 1)



if __name__ == '__main__':
    unittest.main()
