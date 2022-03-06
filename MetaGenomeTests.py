import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
from GenomeTools import *
class MetaGenomeTests(unittest.TestCase):
    def testInstantiation(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        self.assertEqual(metaGen.size(), 1)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen.getName(), "Generic MetaGenome")

        metaGen2 = MetaGenome("datafiles/test2.fasta")
        self.assertEqual(metaGen2.size(), 4)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen2.getName(), "Generic MetaGenome")

    def testIsValidFasta(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        self.assertFalse(isValidFasta("test"))
        self.assertFalse(isValidFasta(".fasta"))
        self.assertFalse(isValidFasta(".FASTA"))
        self.assertFalse(isValidFasta("test.fasta"))
        self.assertFalse(isValidFasta("test.FASTA"))
        self.assertFalse(isValidFasta("t.FASTA"))
        self.assertFalse(isValidFasta("t.fasta"))

    def testFindOffTargets(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        keySequence = complementaryRNA("CCATCGT")
        offTargets = metaGen.findTargetsFromSpacer(keySequence)
        self.assertEqual(len(offTargets), 1)
        self.assertEqual(len(offTargets[0]), 0)

        metaGen2 = MetaGenome("datafiles/test.fasta")
        offTargets = metaGen2.findTargetsFromSpacer(keySequence)
        self.assertEqual(len(offTargets), 1)
        self.assertEqual(len(offTargets[0]), 0)

    def testFindSequenceWithSubsequence(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        self.assertEqual(len(metaGen.getSequence("C")), 1)
        self.assertEqual(metaGen.getSequence("C")[0], "CCATCGT")
        self.assertEqual(len(metaGen.getSequence("N")), 0)
        self.assertEqual(metaGen.getSequence("N"), [])
        self.assertEqual(len(metaGen.getSequence(0)), 1)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")

        metaGen2 = MetaGenome("datafiles/test.fasta")
        self.assertEqual(len(metaGen2.getSequence("C")), 1)
        self.assertEqual(metaGen2.getSequence("C")[0], "CCATCGT")
        self.assertEqual(len(metaGen2.getSequence("N")), 0)
        self.assertEqual(metaGen2.getSequence("N"), [])
        self.assertEqual(len(metaGen2.getSequence(0)), 1)
        self.assertEqual(metaGen2.getSequence(0)[0], "CCATCGT")

    def testSetName(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        self.assertEqual(metaGen.getName(), "Generic MetaGenome")
        metaGen.setName("MetaGen1")
        self.assertEqual(metaGen.getName(), "MetaGen1")

        metaGen2 = MetaGenome("datafiles/test2.fasta")
        self.assertEqual(metaGen2.getName(), "Generic MetaGenome")

    def testAddSequences(self):
        metaGen = MetaGenome("datafiles/test.fasta")
        self.assertEqual(metaGen.size(), 1)
        metaGen.addSequences("datafiles/test2.fasta")
        self.assertEqual(metaGen.size(), 5)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen.getSequence(2)[0], "TCCGGTTACGTATTG")
        metaGen.addSequences("datafiles/test.fasta")
        self.assertEqual(metaGen.size(), 6)
        self.assertEqual(metaGen.getSequence(5)[0], "CCATCGT")



if __name__ == '__main__':
    unittest.main()
