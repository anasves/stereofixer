import pandas as pd

from .StereoChecker import AtomMapper, StereoChecker
from .StereoEnumerator import StereoEnumerator

class PandasStereoEnumerator():

    def __init__(self):
        self.rxn_mapper = AtomMapper()
        self.se = StereoEnumerator()
        self.sc = StereoChecker()

    def get_atommapped_reaction(self, rxn_smiles):
        """
        Use RXNMapper to map the given reaction
        """
        rxn_smiles = rxn_smiles.replace('.[H+]', '').replace('[H+].', '').replace('[H][H].', '').replace('.[H][H]', '')
        res = self.se.rxn_mapper.map_one_reaction(rxn_smiles)['mapped_rxn']
        return res

    def exclude_combinations_with_bad_stereo(self, df):
        """
        function to filter out the atoms with the stereo difference
        only reactions with the consistent stereo tag are left
        """
        df['case'] = df['mapped_rxn'].apply(self.sc.get_stereo_case)
        return df[df['case']=='all good']

    def filter_stereo_options_results(self, reaction_stereo_options: list[str]) -> list[str]:
        """
        : param reaction_stereo_options : list of reaction SMILES with varying stereotags to be checked
        : return : list of reactions with only correct stereochemistry
        """
        df = pd.DataFrame(reaction_stereo_options, columns=['SMILES'])
        df['mapped_rxn']=df['SMILES'].apply(self.get_atommapped_reaction)
        df = self.exclude_combinations_with_bad_stereo(df)
        print(df)

        return df['SMILES'].to_list()
