import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
from GenomeTools import *


class SpacerSequenceTests(unittest.TestCase):
    def testInstantiation(self):
        testGuideSequence = "TTTGGAGAGCCAAGGATTCGGATTCTCGGCTTCCA"
        metaGen = MetaGenome("datafiles/test.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen)
        onTargetScores = newSpacer.getOnTargetScores()
        offTargetScores = newSpacer.getOffTargetScores()
        self.assertTrue(isinstance(onTargetScores, list))
        for onTargetScore in onTargetScores:
            self.assertTrue(isinstance(onTargetScore, float))
        self.assertTrue(isinstance(offTargetScores, list))
        for offTargetScoreList in offTargetScores:
            self.assertTrue(isinstance(offTargetScoreList, list))
            for offTargetScore in offTargetScoreList:
                self.assertTrue(isinstance(offTargetScore, float))
        for guideSequence in newSpacer.getGuideSequences():
            self.assertEqual(newSpacer.getSpacerSequence(), guideSequence[6:26])
        self.assertEqual(len(onTargetScores), len(offTargetScores))
        self.assertEqual(newSpacer.getOffTargetScores(), [[]])

        metaGen2 = MetaGenome("datafiles/test2.fasta")
        newSpacer = SpacerSequence(complementaryRNA(testGuideSequence), metaGen2)
        onTargetScores = newSpacer.getOnTargetScores()
        offTargetScores = newSpacer.getOffTargetScores()
        self.assertTrue(isinstance(onTargetScores, list))
        for onTargetScore in onTargetScores:
            self.assertTrue(isinstance(onTargetScore, float))
        self.assertTrue(isinstance(offTargetScores, list))
        for offTargetScoreList in offTargetScores:
            self.assertTrue(isinstance(offTargetScoreList, list))
            for offTargetScore in offTargetScoreList:
                self.assertTrue(isinstance(offTargetScore, float))
        for guideSequence in newSpacer.getGuideSequences():
            self.assertEqual(newSpacer.getSpacerSequence(), guideSequence[6:26])
        self.assertEqual(len(onTargetScores), len(offTargetScores))
        self.assertEqual(newSpacer.getOffTargetScores(), [[]])

        metaGen3 = MetaGenome("datafiles/test.fasta")
        newSpacer = SpacerSequence("THIS SHOULD BE INVALID", metaGen3)
        self.assertEqual(newSpacer.getOnTargetScores(), [])
        self.assertEqual(newSpacer.getOffTargetScores(), [])
        self.assertEqual(newSpacer.getHeuristics(), [])
        self.assertEqual(len(newSpacer.getOnTargetScores()), len(newSpacer.getOffTargetScores()))

    def testCalcScores(self):
        knownTargetGuideSequence = "CCGGCGTTTTCGCGGAAGACAAAGTCGGGTTCATA"
        metaGen = MetaGenome("datafiles/test3.fasta")
        newSpacer = SpacerSequence(complementaryRNA(knownTargetGuideSequence), metaGen)
        onTargetScores = newSpacer.getOnTargetScores()
        offTargetScores = newSpacer.getOffTargetScores()
        self.assertTrue(isinstance(onTargetScores, list))
        for onTargetScore in onTargetScores:
            self.assertTrue(isinstance(onTargetScore, float))
        self.assertTrue(isinstance(offTargetScores, list))
        for offTargetScoreList in offTargetScores:
            self.assertTrue(isinstance(offTargetScoreList, list))
            for offTargetScore in offTargetScoreList:
                self.assertTrue(isinstance(offTargetScore, float))
        for guideSequence in newSpacer.getGuideSequences():
            self.assertEqual(newSpacer.getSpacerSequence(), guideSequence[6:26])
        self.assertEqual(len(onTargetScores), len(offTargetScores))

        metaGen2 = MetaGenome("datafiles/pf3_illumina_subset.fasta")
        newSpacer = SpacerSequence(complementaryRNA(knownTargetGuideSequence), metaGen2)
        self.assertEqual(len(onTargetScores), len(offTargetScores))
        print(newSpacer.getGuideSequences())
        print(newSpacer.getOffTargetSequences())
        print(newSpacer.getOnTargetScores())
        print(newSpacer.getOffTargetScores())


if __name__ == '__main__':
    unittest.main()
