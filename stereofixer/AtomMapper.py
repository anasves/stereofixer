import warnings

# Suppress warnings from RXN mapper - mostly deprecation warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from rxnmapper import RXNMapper
    
class AtomMapper:

    def __init__(self, canonicalize_rxns=True):
        """
        It is a wrapper for RXNMapper, created by Schwaller et al.
        """
        self.rxn_mapper = RXNMapper()
        self.canonicalize_rxns = canonicalize_rxns

    def map_one_reaction(self, rxnsmiles):
        """
        Use function from RXNMapper to get atom map for one reaction.
        Atom mapper has 512 token length limitation. This function handles the exception resulting from this limitation.
        :param rxnsmiles: reaction smiles
        :return: atom map for one reaction
        """
        try:
            res = self.rxn_mapper.get_attention_guided_atom_maps([rxnsmiles], canonicalize_rxns=self.canonicalize_rxns)[0]
            return {'mapped_rxn':res['mapped_rxn'], 'confidence':res['confidence']}
        except Exception as e:
            print(e)
            return {'mapped_rxn':'tokenlength_error', 'confidence':'err'}