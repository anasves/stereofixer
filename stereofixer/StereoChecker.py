from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions, GetStereoisomerCount
from .AtomMapper import AtomMapper

from tqdm import tqdm

tqdm.pandas()

opts_unassigned = StereoEnumerationOptions(onlyUnassigned=False, unique=True) # onlyUnassigned=False : the atoms with stereo will be reenumerated
opts_keep_assigned = StereoEnumerationOptions(onlyUnassigned=True, unique=True) # onlyUnassigned=True : the atoms with stereo will NOT be reenumerated



class StereoChecker():

    def __init__(self, canonicalize_rxns: bool = True) -> None:
        """
        This class is implementation of atom mapping of the reaction network of RheaDB to transform it into atom transition network.
        """
        self.rxn_mapper = AtomMapper(canonicalize_rxns=canonicalize_rxns)

    def get_differences_in_stereo_per_atom(self, mapped_reaction) -> str:
        """
        checks if the stereotag holds on the left and right sides of the equation
        """
        reactants = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[0].split('.')]
        products = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[1].split('.')]
        reactants = [self.assign_atom_stereotags(mol) for mol in reactants]
        products = [self.assign_atom_stereotags(mol) for mol in products]

        reactant_stereo_maps = []
        for mol in reactants:
            reactant_stereo_maps.extend(self.get_map_stereo(mol))

        product_stereo_maps = []
        for mol in products:
            product_stereo_maps.extend(self.get_map_stereo(mol))

        difference = self.get_set_difference(set(reactant_stereo_maps) , set(product_stereo_maps))
        return ';'.join(difference)

    def assign_atom_stereotags(self, mol):
        """
        Assign stereotags to the molecules of the rxn
        """
        chirality = dict(Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True))

        for atom in mol.GetAtoms():
            if atom.GetIdx() in chirality.keys():
                atom.SetProp('stereotag', chirality[atom.GetIdx()])
            else:
                atom.SetProp('stereotag', 'None')
        return mol

    def get_stereotag(self, atom):
        if atom.HasProp('_CIPCode'):
            return atom.GetProp('_CIPCode')
        return 'None'


    def get_set_difference( self, set1 , set2):
        # Find elements in set1 but not in set2
        unmatched_in_set1 = set1 - set2
        # Find elements in set2 but not in set1
        unmatched_in_set2 = set2 - set1
        # If you want a single set of all unmatched elements
        unmatched_elements = unmatched_in_set1.union(unmatched_in_set2)
        return unmatched_elements

    def get_map_stereo(self, mol):
        map_stereo = list()
        for atom in mol.GetAtoms():
            stereotag = atom.GetProp('stereotag')
            map_num = atom.GetAtomMapNum()
            map_stereo.append(str(map_num)+'-'+stereotag)
        return map_stereo

    def get_reaction_stereoisomer_count(self, mapped_reaction, options_enumeration=opts_keep_assigned):
        reactants = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[0].split('.')]
        products = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[1].split('.')]
        reactants = [str(GetStereoisomerCount(mol,  options=options_enumeration)) for mol in reactants]
        products = [str(GetStereoisomerCount(mol,  options=options_enumeration)) for mol in products]
        return '.'.join(reactants)+'>>'+'.'.join(products)

    def get_mismatched_atoms_one_reaction(self, stereo_difference, num_unassigned_stereo_atoms):
        atom_num_to_tags = self.parse_string_to_dict(stereo_difference)
        print(atom_num_to_tags)
        stereoremoval = []
        missing_stereo_am = []
        missing_stereo = []
        stereoremoval_am = []
        isomerisation = []
        for am, tag_pair in atom_num_to_tags.items():
            if any(at==None for at in tag_pair):
                stereoremoval.append((am, set(tag_pair)-{None}))
                stereoremoval_am.append(am)
            elif set(tag_pair) == {'R', 'S'}:
                isomerisation.append((am, tag_pair))
            elif any(at=='?' for at in tag_pair):
                missing_stereo.append((am, set(tag_pair)-{None}))
                missing_stereo_am.append(am)

        print('stereoremoval', stereoremoval)
        print('isomerisation', isomerisation)
        if len(missing_stereo)>0:
            return 'missing stereo, enumerate', missing_stereo
        if num_unassigned_stereo_atoms ==0 and len(isomerisation)==0: # clean case
            return 'all good', None
        if num_unassigned_stereo_atoms==0 and len(isomerisation)>0: #isomerisation only case
            return 'isomerisation', isomerisation
        if len(isomerisation)==0 and num_unassigned_stereo_atoms>0:
            return 'stereo_mismatched_flat_(should_be_3D)_plus_3D', stereoremoval
        if len(isomerisation)>0 and num_unassigned_stereo_atoms>0:
            isomerisation.extend(stereoremoval)
            return 'both isomerisation and stereoremoval happening: check reaction', isomerisation
        return 'other case - should not be returned if logic is correct', isomerisation

    def get_disappearing_stereo_atoms(self, input_string):
        """
        Compares all theoretically possible stereo atoms on the right and on the left of the equation and checks if
        the stereochemistry disappeares naturally in this reaction and not is result of missing annotation
        """
        # Split the input string by '>>'
        left_part, right_part = input_string.split('>>')

        # Split the left part by '-'
        left_numbers = set(left_part.replace('.', '-').split('-'))

        # Split the right part by either '.' or '-' to handle both delimiters
        right_numbers = set(right_part.replace('.', '-').split('-'))

        # Find numbers that are only in the left set
        only_in_left = left_numbers - right_numbers

        # Find numbers that are only in the right set
        only_in_right = right_numbers - left_numbers

        # return atoms with stereo on both sides
        # return common_numbers
        # return atoms for which stereo disappears
        result = only_in_left.union(only_in_right) - {''}
        return [int(i) for i in result]

    def get_unassigned_stereo_atoms_one_reaction(self, all_theoretically_stereo_atoms, current_stereo_atoms):
        # Split the input string by '>>'
        left_part1, right_part1 = all_theoretically_stereo_atoms.split('>>')
        left_part2, right_part2 = current_stereo_atoms.split('>>')

        # Split the left part by '-'
        left_numbers1 = set(left_part1.replace('.', '-').split('-'))
        left_numbers2 = set(left_part2.replace('.', '-').split('-'))

        # Split the right part by either '.' or '-' to handle both delimiters
        right_numbers1 = set(right_part1.replace('.', '-').split('-'))
        right_numbers2 = set(right_part2.replace('.', '-').split('-'))

        left_unassigned = left_numbers1 - left_numbers2
        right_unassigned = right_numbers1 - right_numbers2

        # return atoms with stereo on both sides
        # return common_numbers
        # return atoms for which stereo disappears
        result = left_unassigned.union(right_unassigned) - {''}
        return [int(i) for i in result], len(result)

    def parse_string_to_dict(self, input_string):
        """
        Output:
        {'26-?', '3-None', '14-?', '16-None', '17-None', '32-None', '12-None', '28-?', '1-None', '24-?', '20-?', '13-?', '33-None',
        '8-?', '5-None', '22-None', '35-None', '18-?', '6-None', '10-?', '34-None', '9-None', '23-None', '27-None', '31-None', '2-?',
        '21-None', '4-None', '7-None', '11-None', '25-None', '29-None', '19-None', '30-?', '15-?'}
        """

        if len(input_string)==0:
            return dict()
        # Split the string by semicolons to get the individual parts
        parts = input_string.split(';')
        result_dict = {}
        for part in parts:
            # Split each part by the hyphen to separate the number and the value
            number, value = part.split('-')
            number = int(number)  # Convert the number to an integer
            # Convert 'None' to the actual None value
            if value == 'None':
                value = None
            # If the number is already in the dictionary, append the new value
            if number in result_dict:
                result_dict[number].append(value)
            else:
                # Otherwise, create a new list with the value
                result_dict[number] = [value]
        # Convert lists to tuples
        for key in result_dict:
            result_dict[key] = tuple(result_dict[key])
        return result_dict

    def get_all_stereo_atoms(self, mapped_reaction, func):
        """
        func - get_all_theoretically_stereo_atoms or get_current_stereo_atoms
        """
        reactants = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[0].split('.')]
        products = [Chem.MolFromSmiles(smiles) for smiles in mapped_reaction.split('>>')[1].split('.')]
        reactants = [func(mol) for mol in reactants]
        products = [func(mol) for mol in products]
        return '.'.join(reactants)+'>>'+'.'.join(products)

    def get_all_theoretically_stereo_atoms(self, mol):
        mol=self.assign_atom_stereotags(mol)
        stereo_atoms=[]
        for atom in mol.GetAtoms():
            stereotag = atom.GetProp('stereotag')
            if stereotag!='None':
                stereo_atoms.append(atom.GetAtomMapNum())
        return '-'.join([str(i) for i in stereo_atoms])


    def get_current_stereo_atoms(self, mol):
        mol=self.assign_atom_stereotags(mol)
        stereo_atoms=[]
        for atom in mol.GetAtoms():
            stereotag = atom.GetProp('stereotag')
            if stereotag!='None' and stereotag!='?':
                stereo_atoms.append(atom.GetAtomMapNum())
        return '-'.join([str(i) for i in stereo_atoms])

    def get_stereo_case(self, mapped_rxn):
        """
        Possible cases:
        'all good'
        'isomerisation only'
        'stereo_mismatched_flat_(should_be_3D)_plus_3D'
        'both isomerisation and stereoremoval happening: check reaction'
        """
        stereo_difference = self.get_differences_in_stereo_per_atom(mapped_rxn)
        current_stereo_atoms=self.get_all_stereo_atoms(mapped_rxn, self.get_current_stereo_atoms)
        all_theoretically_stereo_atoms=self.get_all_stereo_atoms(mapped_rxn, self.get_all_theoretically_stereo_atoms)
        _, num_unassigned_stereo_atoms = self.get_unassigned_stereo_atoms_one_reaction(all_theoretically_stereo_atoms, current_stereo_atoms)
        case, _ = self.get_mismatched_atoms_one_reaction(stereo_difference, num_unassigned_stereo_atoms)
        return case


    def smiles_stereo_analysis(self, input_stereo_string: str):

        analysis_result = dict()
        mapped_rxn = self.rxn_mapper.map_one_reaction(input_stereo_string)
        mapped_rxn = mapped_rxn['mapped_rxn']
        analysis_result["mapped_rxn"]=mapped_rxn

        stereo_difference = self.get_differences_in_stereo_per_atom(mapped_rxn)
        reaction_stereoisomer_count_total = self.get_reaction_stereoisomer_count(mapped_rxn, options_enumeration=opts_unassigned)
        reaction_stereoisomer_count_unenumerated = self.get_reaction_stereoisomer_count(mapped_rxn, options_enumeration=opts_keep_assigned)
        current_stereo_atoms = self.get_all_stereo_atoms(mapped_rxn, self.get_current_stereo_atoms)
        all_theoretically_stereo_atoms = self.get_all_stereo_atoms(mapped_rxn, self.get_all_theoretically_stereo_atoms)
        disappearing_stereo_atoms = self.get_disappearing_stereo_atoms(all_theoretically_stereo_atoms)
        unassigned_stereo_atoms, num_unassigned_stereo_atoms = self.get_unassigned_stereo_atoms_one_reaction(all_theoretically_stereo_atoms, current_stereo_atoms)
        case, mismatched_atoms = self.get_mismatched_atoms_one_reaction( stereo_difference, num_unassigned_stereo_atoms)

        analysis_result['stereo_difference'] = stereo_difference
        analysis_result['reaction_stereoisomer_count_total'] = reaction_stereoisomer_count_total
        analysis_result['reaction_stereoisomer_count_unenumerated'] = reaction_stereoisomer_count_unenumerated
        analysis_result['current_stereo_atoms'] = current_stereo_atoms
        analysis_result['all_theoretically_stereo_atoms'] = all_theoretically_stereo_atoms
        analysis_result['disappearing_stereo_atoms'] = disappearing_stereo_atoms
        analysis_result['unassigned_stereo_atoms'] = unassigned_stereo_atoms
        analysis_result['num_unassigned_stereo_atoms'] = num_unassigned_stereo_atoms
        analysis_result['case'] = case
        analysis_result['mismatched_atoms'] = mismatched_atoms

        return analysis_result
