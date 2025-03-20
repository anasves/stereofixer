# include helper functions for loading and saving reaction networks:

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

def load_reaction_network(project_folder):
    """
    Load a reaction network from TSV files located in the specified project folder.

    Parameters:
        project_folder (str): Path to the project folder containing the TSV files.
        
    Returns:
        dict: A dictionary representing the reaction network, with two keys:
            - "reactions": DataFrame containing reaction data.
            - "compounds": DataFrame containing compound data.
    """
    reactions_file = os.path.join(project_folder, 'reactions.tsv')
    compounds_file = os.path.join(project_folder, 'compounds.tsv')

    reaction_smiles_file = os.path.join(project_folder, 'reaction_smiles.tsv')

    if os.path.exists(reactions_file) and os.path.exists(compounds_file):
        # Load reactions and compounds from TSV files
        reactions_df = pd.read_csv(reactions_file, sep='\t')
        reactions_df = reactions_df.astype({'compound_id':str})
        compounds_df = pd.read_csv(compounds_file, sep='\t')
        compounds_df = compounds_df.astype({'compound_id':str})
        # Create a dictionary to represent the reaction network
        reaction_network = {
            "reactions": reactions_df,
            "compounds": compounds_df
        }
        return reaction_network

    elif os.path.exists(reaction_smiles_file):
        reaction_smiles_df = pd.read_csv(reaction_smiles_file, sep='\t')
        return reaction_smiles_df

    else:
        print('No input. Check your folder structure')
        exit()


def save_reaction_network(fixed_reaction_network, project_folder, reaction_id):
    """
    Save a reaction network to TSV files in the specified project folder.

    Parameters:
        fixed_reaction_network (dict): A dictionary representing the fixed reaction network.
        project_folder (str): Path to the project folder where the TSV files will be saved.
        
    Returns:
        None
    """
    df_chebiid_to_inchi = pd.read_csv('chebiId_inchi.tsv', sep='\t')
    inchi_2_chebiid = dict(zip(df_chebiid_to_inchi['InChI'].to_list(), df_chebiid_to_inchi['CHEBI_ID'].to_list()))
    list_results = list()
    for index, row in fixed_reaction_network.iterrows():
        reaction_id_n = f"{row['reaction_id']}_{index}"
        reaction_equation = reaction_smiles_to_equation(reaction_id_n, row['mapped_rxn'], inchi_2_chebiid, project_folder)
        list_results.append([reaction_id_n, reaction_equation, row['mapped_rxn']])
    df = pd.DataFrame(list_results, columns=['reaction_id', 'equation', 'mapped_reaction_smiles'])
    df.to_csv(os.path.join(project_folder, f'RHEA:{reaction_id}_all_reaction_options.tsv'), sep='\t', index=False)
    

def save_reaction_smiles_as_rdf(rxn_smiles, reaction_id, index):
    """
    Convert a SMILES string representation of a reaction to an RDF file.

    Parameters:
        rxn_smiles (str): The SMILES string of the reaction.
        reaction_id (str): An identifier for the reaction.
        index (int): An index to differentiate multiple reactions.
        
    Returns:
        None
    """
    # Convert the SMILES string to a reaction object
    reaction = AllChem.ReactionFromSmarts(rxn_smiles)

    # Save the reaction to an RDF file
    with open(f"reaction_{reaction_id}_{index}.rdf", "w") as rdf_file:
        rdf_file.write(Chem.rdChemReactions.ReactionToRxnBlock(reaction))

def reaction_smiles_to_equation(reaction_id_n, rxn_smiles, inchi_2_chebiid, project_folder):
    """
    Convert a SMILES string representation of a reaction to equation + .mol files.

    Parameters:
        rxn_smiles (str): The SMILES string of the reaction.
        reaction_id (str): An identifier for the reaction.
        index (int): An index to differentiate multiple reactions.
        
    Returns:
        None
    """
    # Convert the SMILES string to a reaction object
    mol_reactants = [Chem.MolFromSmiles(mol) for mol in rxn_smiles.split('>>')[0].split('.')]
    mol_products = [Chem.MolFromSmiles(mol) for mol in rxn_smiles.split('>>')[1].split('.')]
    # Function to get InChIKeys and export MolFiles
    # Get InChIKeys and export MolFiles for reactants, products, and agents
    Path().joinpath(*[project_folder, 'fixed_reactions', reaction_id_n]).mkdir(parents=True, exist_ok=True)
    reactant_keys = get_inchi_and_export_molfiles(mol_reactants, inchi_2_chebiid, project_folder, reaction_id_n)
    product_keys = get_inchi_and_export_molfiles(mol_products, inchi_2_chebiid, project_folder, reaction_id_n)

    # Create reaction equation using InChIKeys
    reaction_equation = f"{' + '.join(reactant_keys)} >> {' + '.join(product_keys)}"
    return reaction_equation
    
def get_inchi_and_export_molfiles(mol_list, inchi_2_chebiid, project_folder, reaction_id_n):
    stoich_matrix_ids = []
    for mol in mol_list:
        if mol:
            inchi = Chem.MolToInchi(mol)
            chebiid = inchi_2_chebiid.get(inchi)
            if chebiid:
                chebiid = 'CHEBI:'+str(chebiid)
                file_name = os.path.join(project_folder, 'fixed_reactions', reaction_id_n, f"{chebiid}.mol")
                stoich_matrix_ids.append(chebiid)
            else:
                inchi_key = Chem.MolToInchiKey(mol)
                stoich_matrix_ids.append(inchi_key)
                file_name = os.path.join(project_folder, 'fixed_reactions', reaction_id_n, f"{inchi_key}.mol")
            # Export MolFile
            Chem.MolToMolFile(mol, file_name)
    return stoich_matrix_ids

    
def merge_tsvs_into_one_file(tsvs_folder, output_file_name):
    """
    Merge all TSV files in a specified folder into a single TSV file.

    Parameters:
        tsvs_folder (str): The folder containing TSV files to merge.
        output_file_name (str): The name of the output file to save the merged data.
        
    Returns:
        None
    """
    files = [i for i in os.listdir(tsvs_folder) if i.endswith('.tsv')]
    all_dfs = []
    for f in files:
        df = pd.read_csv(os.path.join(tsvs_folder, f), sep='\t')
        df['reaction_id']=f.replace('.tsv', '')
        all_dfs.append(df)
    df_all = pd.concat(all_dfs)
    df_all.to_csv(output_file_name, sep='\t', index=False)
