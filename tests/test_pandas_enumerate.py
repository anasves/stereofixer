import unittest
from unittest.mock import patch
from stereofixer.PandasStereoEnumerator import PandasStereoEnumerator 

class TestPandasStereoEnumeration(unittest.TestCase):
    def setUp(self):
        self.pandas_stereo_enumerator = PandasStereoEnumerator()

    def test_filter_stereo_options_result(self):
        expected_result1 = '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result2 = '[CH3:1][CH2:2][C@@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result3 = '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@@H:3]([CH3:4])[NH2:6].[OH2:5]'
        expected_result4 = '[CH3:1][CH2:2][C@H:3]([CH3:4])[OH:5].[NH3:6]>>[CH3:1][CH2:2][C@H:3]([CH3:4])[NH2:6].[OH2:5]'
        reaction_stero_options = [expected_result1, expected_result2, expected_result3, expected_result4]
        df_filtered = self.pandas_stereo_enumerator.filter_stereo_options_results(reaction_stero_options)
        self.assertEqual(len(df_filtered), 2)