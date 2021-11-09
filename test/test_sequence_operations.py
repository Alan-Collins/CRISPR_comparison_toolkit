import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from CRISPR_comparison_toolkit.cctk import sequence_operations as seqops

class TestRevComp(unittest.TestCase):
	
	def test_crev_comp(self):
		
		self.assertEqual(seqops.rev_comp("ATCGatcg"), "cgatCGAT")




if __name__ == '__main__':
	unittest.main()
