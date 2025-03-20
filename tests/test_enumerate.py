import unittest
from unittest.mock import patch
from stereofixer.StereoEnumerator import StereoEnumerator 

class TestStereoEnumeration(unittest.TestCase):
    def setUp(self):
        self.stereo_enumerator = StereoEnumerator()

    def test_enumerate_stereoisomers_per_compound(self):
        isomers = self.stereo_enumerator.enumerate_stereoisomers_per_compound("C(C)(CC)(CCC)[O-]")
        isomers.sort()
        self.assertEqual(isomers, ['CCC[C@@](C)([O-])CC', 'CCC[C@](C)([O-])CC'])

    def test_get_all_stereo_combinations(self):
        combinations = self.stereo_enumerator.get_all_stereo_combinations(['C(CCC)(C)([O-])CC', 'C'])
        self.assertEqual(len(combinations), 2)

    def test_remove_stereochemistry_from_atoms(self):
        reaction_string = "[CH3:1][C@H:2]([OH:3])[CH2:4][C@:5]([CH3:6])([Cl:8])[OH:9].[NH3:7]>>[CH3:1][C@H:2]([OH:3])[CH2:4][C@:5]([CH3:6])([NH2:7])[Cl:8].[OH2:9]"
        atom_number = 5
        result = self.stereo_enumerator.remove_stereochemistry_from_atoms(reaction_string, atom_number)
        expected_result = '[CH3:1][C@H:2]([OH:3])[CH2:4][C:5]([CH3:6])([Cl:8])[OH:9].[NH3:7]>>[CH3:1][C@H:2]([OH:3])[CH2:4][C:5]([CH3:6])([NH2:7])[Cl:8].[OH2:9]'
        self.assertEqual(result, expected_result)
    
    def test_enumerate_stereoisomers_per_reaction(self):
        """
        Does not enumerate properly since the hydrogens are not explicitly defined in the atom mapped string
        """
        result, _= self.stereo_enumerator.enumerate_stereoisomers_per_reaction('[C:1][C:2][C:3]([C:4])[O:5].[N:6]>>[C:1][C:2][C:3]([C:4])[N:6].[O:5]',
                                                           mismatched_atoms=[3,])
        expected_result1 = '[C:1][C:2][C:3]([C:4])[O:5].[NH3:6]>>[C:1][C:2][C:3]([C:4])[N:6].[OH2:5]'
        result.sort()
        self.assertEqual(result[0], expected_result1)

    def test_full_one_stereo_center(self):
        """
        It only works with reactions that are atom mapped
        """
        result, _ = self.stereo_enumerator.enumerate_stereoisomers_per_reaction('[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]',
                                                           mismatched_atoms=[3,])
        expected_result1 = '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result2 = '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result3 = '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result4 = '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]'
        result.sort()
        self.assertEqual(result[0], expected_result1)
        self.assertEqual(result[1], expected_result2)
        self.assertEqual(result[2], expected_result3)
        self.assertEqual(result[3], expected_result4)

    # def test_full_two_stereo_centers(self):
    # this test does not work - incorrect input SMILES, have to fix SMILES
    #     result, _ = self.stereo_enumerator.enumerate_reaction('[OH:1][C@@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[NH3:9]'\
    #                                                           '>>[NH2:9][C@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[OH2:1]',
    #                                                        mismatched_atoms=[2,])
    #     print(result['SMILES'])
    #     expected_result1 = '[OH:1][C@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[NH3:9]'\
    #                                                           '>>[NH2:9][C@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[OH2:1]'
    #     expected_result2 = '[OH:1][C@@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[NH3:9]'\
    #                                                           '>>[NH2:9][C@@:2]([CH3:3])([Cl:4])[CH2:5][C@H:6]([CH3:7])[OH:8].[OH2:1]'
    #     result.sort_values(by='SMILES', inplace=True)
    #     result.reset_index(inplace=True)
    #     self.assertEqual(result['SMILES'][0], expected_result1)
    #     self.assertEqual(result['SMILES'][1], expected_result2)
