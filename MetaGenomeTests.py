import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence

class MetaGenomeTests(unittest.TestCase):
    def testInstantiation(self):
        metaGen = MetaGenome("test.fasta")
        self.assertEqual(metaGen.size(), 1)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen.getName(), "Generic MetaGenome")

        metaGen2 = MetaGenome("test2.fasta")
        self.assertEqual(metaGen2.size(), 4)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen2.getName(), "Generic MetaGenome")

    def testIsValidFasta(self):
        metaGen = MetaGenome("test.fasta")
        self.assertFalse(metaGen.isValidFasta("test"))
        self.assertFalse(metaGen.isValidFasta(".fasta"))
        self.assertFalse(metaGen.isValidFasta(".FASTA"))
        self.assertTrue(metaGen.isValidFasta("test.fasta"))
        self.assertTrue(metaGen.isValidFasta("test.FASTA"))
        self.assertTrue(metaGen.isValidFasta("t.FASTA"))
        self.assertTrue(metaGen.isValidFasta("t.fasta"))

    def testFindOffTargets(self):
        metaGen = MetaGenome("test.fasta")
        keySequence = SpacerSequence.complementaryRNA(None, "CCATCGT")
        offTargets = metaGen.findOffTargets(keySequence)
        self.assertEqual(len(offTargets), 1)
        self.assertEqual(len(offTargets[0]), 1)
        self.assertEqual(str(offTargets[0][0]), "CCATCGT")

        metaGen2 = MetaGenome("test.fasta")
        offTargets = metaGen2.findOffTargets(keySequence)
        self.assertEqual(len(offTargets), 1)
        self.assertEqual(len(offTargets[0]), 1)
        self.assertEqual(str(offTargets[0][0]), "CCATCGT")

    def testFindSequenceWithSubsequence(self):
        metaGen = MetaGenome("test.fasta")
        self.assertEqual(len(metaGen.getSequence("C")), 1)
        self.assertEqual(metaGen.getSequence("C")[0], "CCATCGT")
        self.assertEqual(len(metaGen.getSequence("N")), 1)
        self.assertEqual(metaGen.getSequence("N")[0], None)
        self.assertEqual(len(metaGen.getSequence(0)), 1)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")

        metaGen2 = MetaGenome("test.fasta")
        self.assertEqual(len(metaGen2.getSequence("C")), 1)
        self.assertEqual(metaGen2.getSequence("C")[0], "CCATCGT")
        self.assertEqual(len(metaGen2.getSequence("N")), 1)
        self.assertEqual(metaGen2.getSequence("N")[0], None)
        self.assertEqual(len(metaGen2.getSequence(0)), 1)
        self.assertEqual(metaGen2.getSequence(0)[0], "CCATCGT")

    def testSetName(self):
        metaGen = MetaGenome("test.fasta")
        self.assertEqual(metaGen.getName(), "Generic MetaGenome")
        metaGen.setName("MetaGen1")
        self.assertEqual(metaGen.getName(), "MetaGen1")

        metaGen2 = MetaGenome("test2.fasta")
        self.assertEqual(metaGen2.getName(), "Generic MetaGenome")

    def testAddSequences(self):
        metaGen = MetaGenome("test.fasta")
        self.assertEqual(metaGen.size(), 1)
        metaGen.addSequences("test2.fasta")
        self.assertEqual(metaGen.size(), 5)
        self.assertEqual(metaGen.getSequence(0)[0], "CCATCGT")
        self.assertEqual(metaGen.getSequence(2)[0], "TCCGGTTACGTATTG")
        metaGen.addSequences("test.fasta")
        self.assertEqual(metaGen.size(), 6)
        self.assertEqual(metaGen.getSequence(5)[0], "CCATCGT")



if __name__ == '__main__':
    unittest.main()
