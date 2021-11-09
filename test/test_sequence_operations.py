import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from CRISPR_comparison_toolkit.cctk import sequence_operations as seqops

class TestRevComp(unittest.TestCase):
	
	def test_rev_comp_simple(self):
		
		self.assertEqual(seqops.rev_comp("ATCGatcg"), "cgatCGAT")

	def test_rev_comp_weird(self):
		self.assertEqual(
			seqops.rev_comp("abracadabraNnN...#!/"),
			"/!#...NnNtrbtdtgtrbt")

class TestHamming(unittest.TestCase):

	def test_hamming_with_sequence(self):

		self.assertEqual(seqops.hamming("ATCG", "ATCG"),
			(0,4,4))

		self.assertEqual(seqops.hamming("ATCG", "GATCGG"),
			(6,6,1))

class TestNeedle(unittest.TestCase):

	def test_needle_with_similar_sequences(self):
		self.assertEqual(seqops.needle(
				"ATTG",
				"ATCG",
				match=1),
			("ATTG", "ATCG"))

	def test_needle_with_different_sequences(self):

		self.assertEqual(seqops.needle(
				"TTTTTTATCGTTTTTT",
				"ATCG",
				match=1),
			("TTTTTTATCGTTTTTT", "------ATCG------"))

		self.assertEqual(seqops.needle(
				"ACTCGCTGCCTGTGA",
				"CGATCACGATGTAGC",
				match=1),
			("-ACTCGCTGCCTGT-GA", "CGATCAC-G-ATGTAGC"))

	def test_needle_with_CRISPR_arrays(self):

		self.assertEqual(seqops.needle(
				[1,2,3,4,5,6,7],
				[8,9,10,11,12,5,7]),
			(['-',1,2,3,4,5,6,7], [8,9,10,11,12,5,'-',7]))

		self.assertEqual(seqops.needle(
				[1,2,3,4,5,6,7],
				[8,9,3,11,12,5,7]),
			([1,2,3,'-',4,5,6,7], [8,9,3,11,12,5,'-',7]))



if __name__ == '__main__':
	unittest.main()
