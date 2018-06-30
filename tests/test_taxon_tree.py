import os
import sys
import unittest
from BioSQL import BioSeqDatabase
from biosqlx.taxon_tree import TaxonTree, IntegrityError
from . import connection_parameters

class TaxonTreeTest(unittest.TestCase):
    """Check TaxonTree interface."""

    def setUp(self):
        # drop any old database and create a new one:
        testdb, dbdriver, dbuser, dbpassword, dbhost = connection_parameters(create=True)
        # connect to new database:
        self.server = BioSeqDatabase.open_database(driver=dbdriver,
                                                   user=dbuser, passwd=dbpassword,
                                                   host=dbhost, db=testdb)
        self._create_taxonomy()
        self.taxon_tree = TaxonTree(self.server.adaptor)
        self.testdb = testdb

    def tearDown(self):
        os.remove(self.testdb)

    def _create_taxonomy(self):
        r"""add in dummy values to taxon and taxon_name tables

        level
        1                    1(1)22
                _______________|___________________
               |               |                   |
        2    2(2)5           6(4)11             12(7)21
               |            /     \            /      \
        3    3(3)4       7(5)8   9(6)10    13(8)16   17(10)20
                                              |          |
        4                                  14(9)15   18(11)19



        """
        nodes = [
                (1, -1, 1, 22),
                (2, 1, 2, 5),
                (3, 2, 3, 4),
                (4, 1, 6, 11),
                (5, 4, 7, 8),
                (6, 4, 9, 10),
                (7, 1, 12, 21),
                (8, 7, 13, 16),
                (9, 8, 14, 15),
                (10, 7, 17, 20),
                (11, 10, 18, 19)
                ]
        names = [
                (1, 'one', 'scientific name'),
                (2, 'two', 'scientific name'),
                (3, 'three', 'scientific name'),
                (4, 'four', 'scientific name'),
                (5, 'five', 'scientific name'),
                (6, 'six', 'scientific name'),
                (7, 'seven', 'scientific name'),
                (8, 'eight', 'scientific name'),
                (9, 'nine', 'scientific name'),
                (10, 'ten', 'scientific name'),
                (11, 'eleven', 'scientific name')
                ]
        insert_nodes = 'INSERT INTO taxon(taxon_id, parent_taxon_id, left_value, right_value) VALUES(%s, %s, %s, %s)'
        insert_names = 'INSERT INTO taxon_name(taxon_id, name, name_class) VALUES(%s, %s, %s)'
        self.server.adaptor.executemany(insert_nodes, nodes)
        self.server.adaptor.executemany(insert_names, names)

    def test_add_to_terminal_node(self):
        """test addition of new taxonomy to a terminal (leaf) node."""
        three = self.taxon_tree.find_elements(name='three')[0]
        twelve = self.taxon_tree.add("twelve", "scientific_name", parent=three)
        self.assertEqual(twelve._parent_id, three._id)
        self.assertEqual(twelve._right_val, 5)
        self.assertEqual(twelve._left_val, 4)

    def test_add_to_internal_node(self):
        """test addition of new taxonomy to an internal node."""
        four = self.taxon_tree.find_elements(name='four')[0]
        five = self.taxon_tree.find_elements(name='five')[0]
        six = self.taxon_tree.find_elements(name='six')[0]
        twelve = self.taxon_tree.add('twelve', 'scientific_name', parent=four)

        five_after = self.taxon_tree.find_elements(name='five')[0]
        six_after = self.taxon_tree.find_elements(name='six')[0]
        four_after = self.taxon_tree.find_elements(name='four')[0]

        self.assertEqual(twelve._parent_id, four._id)
        self.assertEqual(twelve._left_val, 7)
        self.assertEqual(twelve._right_val, 8)
        self.assertEqual(four_after._left_val, 6)
        self.assertEqual(four_after._right_val, 13)
        self.assertEqual(five_after._left_val, 9)
        self.assertEqual(five_after._right_val, 10)
        self.assertEqual(six_after._left_val, 11)
        self.assertEqual(six_after._right_val, 12)

    def test_move_terminal_to_other_terminal(self):
        '''test movement of a terminal node to another terminal node'''
        three = self.taxon_tree.find_elements(name='three')[0]
        five = self.taxon_tree.find_elements(name='five')[0]
        self.taxon_tree.move(three, five)

        three_after = self.taxon_tree.find_elements(name='three')[0]
        five_after = self.taxon_tree.find_elements(name='five')[0]
        two_after = self.taxon_tree.find_elements(name='two')[0]

        self.assertEqual(three_after._parent_id, five_after._id)
        self.assertEqual(three_after._left_val, 6)
        self.assertEqual(three_after._right_val, 7)
        self.assertEqual(five_after._left_val, 5)
        self.assertEqual(five_after._right_val, 8)
        self.assertEqual(two_after._left_val + 1, two_after._right_val)

    def test_move_terminal_to_internal(self):
        '''test movement of a terminal node to an internal node'''
        three = self.taxon_tree.find_elements(name='three')[0]
        four = self.taxon_tree.find_elements(name='four')[0]
        self.taxon_tree.move(three, four)

        three_after = self.taxon_tree.find_elements(name='three')[0]
        four_after = self.taxon_tree.find_elements(name='four')[0]
        two_after = self.taxon_tree.find_elements(name='two')[0]

        self.assertEqual(three_after._parent_id, four_after._id)
        self.assertEqual(three_after._left_val, 5)
        self.assertEqual(three_after._right_val, 6)
        self.assertEqual(four_after._left_val, 4)
        self.assertEqual(four_after._right_val, 11)
        self.assertEqual(two_after._left_val + 1, two_after._right_val)

    def test_move_internal_higher_internal(self):
        '''test movement of an internal node to a higher internal node'''
        eight = self.taxon_tree.find_elements(name='eight')[0]
        one = self.taxon_tree.find_elements(name='one')[0]
        self.taxon_tree.move(eight, one)

        eight_after = self.taxon_tree.find_elements(name='eight')[0]
        one_after = self.taxon_tree.find_elements(name='one')[0]

        self.assertEqual(one._left_val, one_after._left_val)
        self.assertEqual(one._right_val, one_after._right_val)
        self.assertEqual(eight_after._parent_id, one_after._id)
        self.assertEqual(eight_after._left_val, 2)
        self.assertEqual(eight_after._right_val, 5)

    def test_move_internal_lower_internal(self):
        '''test movement of an internal node to a lower internal node'''
        four = self.taxon_tree.find_elements(name='four')[0]
        eight = self.taxon_tree.find_elements(name='eight')[0]
        self.taxon_tree.move(four, eight)

        eight_after = self.taxon_tree.find_elements(name='eight')[0]
        four_after = self.taxon_tree.find_elements(name='four')[0]

        self.assertEqual(four_after._parent_id, eight_after._id)
        self.assertEqual(eight_after._left_val, 7)
        self.assertEqual(eight_after._right_val, 16)
        self.assertEqual(four_after._left_val, 8)
        self.assertEqual(four_after._right_val, 13)

    def test_move_internal_to_terminal(self):
        '''test movement of an internal node to a terminal node'''
        four = self.taxon_tree.find_elements(name='four')[0]
        nine = self.taxon_tree.find_elements(name='nine')[0]
        self.taxon_tree.move(four, nine)

        nine_after = self.taxon_tree.find_elements(name='nine')[0]
        four_after = self.taxon_tree.find_elements(name='four')[0]

        self.assertEqual(four_after._parent_id, nine_after._id)
        self.assertEqual(nine_after._left_val, 8)
        self.assertEqual(nine_after._right_val, 15)
        self.assertEqual(four_after._left_val, 9)
        self.assertEqual(four_after._right_val, 14)

    def test_move_parent_under_child_raises_exception(self):
        '''test movement of a node underneath one of it's children raises an exception'''
        four = self.taxon_tree.find_elements(name='four')[0]
        five = self.taxon_tree.find_elements(name='five')[0]
        with self.assertRaises(IntegrityError):
            self.taxon_tree.move(four, five)

    def test_move_node_under_itself_raises_exception(self):
        '''test movement of a node underneath itself raises an exception'''
        four = self.taxon_tree.find_elements(name='four')[0]
        other_four = self.taxon_tree.find_elements(name='four')[0]
        with self.assertRaises(IntegrityError):
            self.taxon_tree.move(four, other_four)

