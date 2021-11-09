import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import dendropy
from CRISPR_comparison_toolkit.cctk import tree_operations as treeops

class TestCreateNodeIds(unittest.TestCase):

	def test_create_internal_node_ids(self):

		self.assertEqual(treeops.create_internal_node_ids(
				5,
				"Pre_",
				"letters",
				), 
			["Pre_a", "Pre_b", "Pre_c", "Pre_d", "Pre_e"])

		self.assertEqual(treeops.create_internal_node_ids(
				30,
				chars="numbers",
				),
			[str(i) for i in range(1,31)])

	def test_create_internal_node_ids_Raise(self):

		with self.assertRaises(ValueError):
			treeops.create_internal_node_ids(1)
			treeops.create_internal_node_ids(5, chars="wrong")
		with self.assertRaises(TypeError):
			treeops.create_internal_node_ids("1")
			treeops.create_internal_node_ids(5.0)
			treeops.create_internal_node_ids(5, prefix=10)


class TestScaleBranches(unittest.TestCase):

	def setUp(self):
		self.tree = dendropy.Tree.get(
			data="((A:1,B:2):3,(C:1,D:5):4);", 
			schema="newick")

	def test_scale_branch_lengths(self):

		self.assertEqual(treeops.scale_branches(self.tree, 10)
			.as_string("newick").strip(), 
			"((A:2,B:4):6,(C:2,D:10):8):0;")

		

if __name__ == '__main__':
	unittest.main()
