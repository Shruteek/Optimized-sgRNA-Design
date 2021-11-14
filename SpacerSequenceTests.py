import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
from GenomeTools import *


class SpacerSequenceTests(unittest.TestCase):
    def testInstantiation(self):
        testGuideSequence = "TTTGGAGAGCCAAGGATTCGGATTCTCGGCTTCCA"
        metaGen = MetaGenome("test.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen)
        self.assertTrue(isinstance(newSpacer.getOnTargetScore(), float))
        self.assertTrue(isinstance(newSpacer.getOffTargetScores(), list))
        self.assertEqual(newSpacer.getOffTargetScores(), [])

        metaGen2 = MetaGenome("test2.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen2)
        self.assertTrue(isinstance(newSpacer.getOnTargetScore(), float))
        self.assertTrue(isinstance(newSpacer.getOffTargetScores(), list))
        self.assertEqual(newSpacer.getOffTargetScores(), [])

        metaGen3 = MetaGenome("test.fasta")
        newSpacer = SpacerSequence("THIS SHOULD BE INVALID", metaGen3)
        self.assertEqual(newSpacer.getOnTargetScore(), 0)
        self.assertEqual(newSpacer.getOffTargetScores(), [])
        self.assertEqual(newSpacer.getHeuristic(), 0)

    def testCalcScores(self):
        testGuideSequence = "TTTGGAGAGCCAAGGATTCGGATTCTCGGCTTCCA"
        knownGuideSequence = "CCGGCGTTTTCCCGGAAGACAAAGTCGGGTTCATA"
        metaGen = MetaGenome("test3.fasta")
        newSpacer = SpacerSequence(complementaryRNA(knownGuideSequence), metaGen)
        print(newSpacer.getOnTargetScore())
        print(newSpacer.getOffTargetScores())
        self.assertEqual(newSpacer.getOffTargetScores()[0], 0.5)

        metaGen2 = MetaGenome("pf3_illumina_subset.fasta")
        newSpacer = SpacerSequence(complementaryRNA(knownGuideSequence), metaGen2)
        print(newSpacer.getOnTargetScore())
        print(newSpacer.getOffTargetScores())


if __name__ == '__main__':
    unittest.main()
