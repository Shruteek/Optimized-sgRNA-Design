import unittest
from Bio import SeqIO
import Sequence

class SequenceTests(unittest.TestCase):
    def test_setMismatchStrictness(self):
        self.assertEqual(True, False)

    def test_init(self):
        self.assertEqual(True, False)

    def test_correlateStrings(self):
        key = "ACTG"
        genome = "TACTG"
        for record in SeqIO.parse("test.fasta", "fasta"):
            if isinstance(record, SeqIO.SeqRecord):
                print("Successful SeqRecord opened!")
            seqRecord = record
        seq = Sequence.Sequence(seqRecord)
        print(seq.correlateSequences(key, genome))


if __name__ == '__main__':
    unittest.main()
