from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .AtomMapper import AtomMapper
import itertools
        
class StereoEnumerator():

    def __init__(self):
        self.rxn_mapper = AtomMapper()

    def enumerate_stereoisomers_per_compound(self, smiles):
        """
        Accept one SMILES - return list of SMILES with all stereoforms
        """
        m = Chem.MolFromSmiles(smiles)
        opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
        isomers = tuple(EnumerateStereoisomers(m, options=opts))
        return list(set(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers))

    def get_all_stereo_combinations(self, side_smiles_list):
        """
        side_smiles_list - smiles of the compounds on one side of the equation (reactants or products)
        accept smiles with the given stereoform - return all combinations of stereoforms
        """
        return list(itertools.product(*[self.enumerate_stereoisomers_per_compound(smiles) for smiles in side_smiles_list]))

    def remove_stereochemistry_from_atoms(self, mapped_rxn, atom_number):
        """
        mapped_rxn - reaction SMILES with atom map number per atom
        atom_number - int = number of the atoms with the wrong stereochemistry identified
        Since the atom corresponding to atom_number has incorrect stereochemistry, the existing stereochemistry tag will be removed
        """
        # Convert the reaction string into an RDKit reaction object
        rxn = AllChem.ReactionFromSmarts(mapped_rxn)

        # Iterate over reactants and products to remove stereochemistry
        for mol in rxn.GetReactants():
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == atom_number:
                    # Remove stereochemistry
                    atom.SetChiralTag(Chem.CHI_UNSPECIFIED)

        for mol in rxn.GetProducts():
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == atom_number:
                    # Remove stereochemistry
                    atom.SetChiralTag(Chem.CHI_UNSPECIFIED)

        # Convert molecules back to SMILES without the extra '&' characters
        reactants = ".".join(Chem.MolToSmiles(reactant, isomericSmiles=True) for reactant in rxn.GetReactants())
        products = ".".join(Chem.MolToSmiles(product, isomericSmiles=True) for product in rxn.GetProducts())

        # Reconstruct the reaction string
        return f"{reactants}>>{products}"

    def enumerate_stereoisomers_per_reaction(self, mapped_rxn, mismatched_atoms=list()):
        """
        : param mapped_rxn : atom mapped reaction
        : param mismatched_atoms : list of int, identified in StereoChecker
        """
        for atom_number in mismatched_atoms:
            mapped_rxn = self.remove_stereochemistry_from_atoms(mapped_rxn, atom_number)
        reactants, products = (i.split('.') for i in mapped_rxn.split('>>'))
        reactants = ['.'.join(i) for i in self.get_all_stereo_combinations(reactants)]
        products = ['.'.join(i) for i in self.get_all_stereo_combinations(products)]
        print(reactants)
        print(products)
        # logic to exclude cases with too may stereo options
        # This code is to track and correct a few stereo options, if the number is too big, it is better to check for the correct stereo version in the literature directly
        if len(reactants)*len(products)>100:
            return [], len(reactants)*len(products)
        react_prod_combs = list(itertools.product(reactants, products))
        print(react_prod_combs)
        if len(react_prod_combs)>100:
            return [], len(react_prod_combs)
        reaction_stereo_options = ['>>'.join(i) for i in react_prod_combs]

        return reaction_stereo_options, len(react_prod_combs)
    

