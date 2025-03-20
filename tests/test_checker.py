import unittest
from unittest.mock import patch
from stereofixer.StereoChecker import StereoChecker

class TestStereoChecker(unittest.TestCase):
    def setUp(self):
        self.stereo_checker = StereoChecker()

    def test_1(self):
        """
        
        """
        input_stereo_string1 = 'C(C)(CC)(CCC)[O-]>>[C@](C)(CC)(CCC)[O-]'
        result = self.stereo_checker.smiles_stereo_analysis(input_stereo_string1)

        self.assertEqual(result['case'], 'missing stereo, enumerate')

    def test_2(self):
        """
        
        """
        input_stereo_string2 = '[C@](C)(CC)(CCC)[O-]>>C(C)(CC)(CCC)[O-]'
        result = self.stereo_checker.smiles_stereo_analysis(input_stereo_string2)

        self.assertEqual(result['case'], 'missing stereo, enumerate')

    def test_3(self):
        """
        d
        """
        input_stereo_string3 = '[C@](C)(CC)(CCC)[O-]>>[C@@](C)(CC)(CCC)[O-]'
        result = self.stereo_checker.smiles_stereo_analysis(input_stereo_string3)

        self.assertEqual(result['case'], 'isomerisation')

    def test_4(self):
        """
        
        """
        input_stereo_string4 = 'C(C)(C)(CCC)[O-].C>>[C@@](C)(CC)(CCC)[O-]'
        result = self.stereo_checker.smiles_stereo_analysis(input_stereo_string4)

        self.assertEqual(result['case'], 'missing stereo, enumerate')