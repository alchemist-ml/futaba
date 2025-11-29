import pubchempy as pcp

class Chemical():
    def __init__(self, name, formula, iupac, synonyms, weight, cano_smiles, iso_smiles):
        self.name = name
        self.formula = formula
        self.iupac_name = iupac
        self.synonyms = synonyms
        self.weight = weight
        self.cano_smiles = cano_smiles
        self.iso_smiles = iso_smiles

def chem_info(query, query_type):
    compounds = pcp.get_compounds(query, query_type)

    if not compounds:
        return f"Error: No compound found for '{query}' with type: <{query_type}>"
    
    compound  = compounds[0]

    if query_type == 'name':
        name = query
        formula = compound.molecular_formula
    elif query_type == 'formula':
        formula = query
        name = compound.synonyms[0]

    iupac = compound.iupac_name
    synonyms = compound.synonyms[1:4] if compound.synonyms else []
    weight = compound.molecular_weight
    cano_smiles = compound.connectivity_smiles
    iso_smiles = compound.smiles

    return Chemical(name, formula, iupac, synonyms, weight, cano_smiles, iso_smiles)
