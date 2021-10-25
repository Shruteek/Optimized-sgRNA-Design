import unittest
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
from Sequence import Sequence

class SpacerSequenceTests(unittest.TestCase):
    def testInitialization(self):
        metaGen = MetaGenome("test.fasta")
        newSpacer = SpacerSequence("CTTATATCACGTCCATAACGGGG", metaGen)
        self.assertEqual(1, 1)

    def test_something_else(self):
        self.assertEqual(1, 1)


if __name__ == '__main__':
    unittest.main()
