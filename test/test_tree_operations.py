import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import dendropy
from CRISPR_comparison_toolkit.cctk import tree_operations

class TestCreateNodeIds(unittest.TestCase):

	def test_create_internal_node_ids(self):

		self.assertEqual(tree_operations.create_internal_node_ids(
				5,
				"Pre_",
				"letters",
				), 
			["Pre_a", "Pre_b", "Pre_c", "Pre_d", "Pre_e"])

		self.assertEqual(tree_operations.create_internal_node_ids(
				30,
				chars="numbers",
				),
			[str(i) for i in range(1,31)])

		with self.assertRaises(ValueError):
			tree_operations.create_internal_node_ids(1)
			tree_operations.create_internal_node_ids(5, chars="wrong")
		with self.assertRaises(TypeError):
			tree_operations.create_internal_node_ids("1")
			tree_operations.create_internal_node_ids(5.0)
			tree_operations.create_internal_node_ids(5, prefix=10)


class TestScaleBranches(unittest.TestCase):

	def setUp(self):
		self.tree = dendropy.Tree.get(
			data="((A:1,B:2):3,(C:1,D:5):4);", 
			schema="newick")

	def test_scale_branch_lengths(self):

		self.assertEqual(tree_operations.scale_branches(self.tree, 10)
			.as_string("newick").strip(), 
			"((A:2,B:4):6,(C:2,D:10):8):0;")

		

if __name__ == '__main__':
	unittest.main()
