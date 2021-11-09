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





if __name__ == '__main__':
	unittest.main()
