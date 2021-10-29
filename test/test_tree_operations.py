import sys
sys.path.append('../cctk')

import unittest
import tree_operations

class TestTreeOps(unittest.TestCase):

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


if __name__ == '__main__':
	unittest.main()