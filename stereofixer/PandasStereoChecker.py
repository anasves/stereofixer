import os
import pandas as pd

from tqdm import tqdm
tqdm.pandas()

from StereoChecker import AtomMapper, StereoChecker

class PandasStereoChecker():
	
    def __init__(self, projectname, no_h_smiles=True, canonicalize_rxns=True):
        """
        This class is implementation of atom mapping of the reaction network of RheaDB to transform it into atom transition network.
        """
        self.projectname = projectname
        self.rxn_mapper = AtomMapper(canonicalize_rxns=canonicalize_rxns)
        self.sc = StereoChecker()
        
        self.no_h_smiles = no_h_smiles
        self.canonicalize_rxns = canonicalize_rxns
        os.makedirs(projectname, exist_ok=True)
    
    def check_all_reactions(self, reaction_smiles_df):
        
        df = self.calculate_atom_mapped_reaction_data(reaction_smiles_df)

        df['stereo_difference']=df['mapped_rxn'].apply(self.sc.get_differences_in_stereo_per_atom)
        df['reaction_stereoisomer_count_total']=df['mapped_rxn'].apply(self.sc.get_reaction_stereoisomer_count)
        df['reaction_stereoisomer_count_unenumerates']=df['mapped_rxn'].apply(self.sc.get_reaction_stereoisomer_count)
        df['current_stereo_atoms']=df['mapped_rxn'].apply(self.sc.get_all_stereo_atoms, args=[self.get_current_stereo_atoms,])
        df['all_theoretically_stereo_atoms']=df['mapped_rxn'].apply(self.sc.get_all_stereo_atoms, args=[self.get_all_theoretically_stereo_atoms,])
        df['disappearing_stereo_atoms']=df['all_theoretically_stereo_atoms'].apply(self.sc.get_disappearing_stereo_atoms)
        df[['unassigned_stereo_atoms', 'num_unassigned_stereo_atoms']]=df.apply(self.sc.get_unassigned_stereo_atoms, axis=1, result_type='expand')
        df[['case','mismatched_atoms']]=df.apply(self.sc.get_mismatched_atoms, axis=1, result_type='expand')

        df.to_csv(os.path.join(self.projectname, f'stereochekcer_checks_done.tsv'), sep='\t', index=False)
        return df

    def get_mismatched_atoms(self, row):
        stereo_difference = row['stereo_difference']
        num_unassigned_stereo_atoms = row['num_unassigned_stereo_atoms']
        return  self.sc.get_mismatched_atoms_one_reaction(self, stereo_difference, num_unassigned_stereo_atoms)

    def get_unassigned_stereo_atoms(self, row):
        all_theoretically_stereo_atoms = row['all_theoretically_stereo_atoms']
        current_stereo_atoms = row['current_stereo_atoms']
        return self.sc.get_unassigned_stereo_atoms_one_reaction(all_theoretically_stereo_atoms, current_stereo_atoms)

    ###########################
    # RXN Mapper: atom mapping
    ###########################
    def calculate_atom_mapped_reaction_data(self, reaction_smiles_df):

        df = reaction_smiles_df.copy()
        
        if self.no_h_smiles:
            df['rxnsmiles_no_h'] = df['reaction_smiles'].apply(
                    lambda x: x.replace('.[H+]', '').replace('[H+].', '').replace('[H][H].', '').replace('.[H][H]', ''))
            
            # atom map reactions
            print("atom map reactions")
            df[['mapped_rxn', 'confidence']] = df['rxnsmiles_no_h'].progress_apply(self.rxn_mapper.map_one_reaction)

        else:
            df[['mapped_rxn', 'confidence']] = df['reaction_smiles'].progress_apply(self.rxn_mapper.map_one_reaction)
        
        # remove the reactions with no atom mapping
        df = df[df['mapped_rxn'] != 'tokenlength_error']

        return df