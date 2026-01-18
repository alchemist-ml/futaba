import pubchempy as pcp
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from rdkit.Chem import AllChem, Draw

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


def rxn(smiles, model_name):
    tokenizer = AutoTokenizer.from_pretrained(model_name, return_tensors='pt')
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)

    inp = tokenizer(smiles, return_tensors='pt')
    output = model.generate(**inp, num_beams=4, num_return_sequences=1, return_dict_in_generate=True, output_scores=True)
    output = tokenizer.decode(output['sequences'][0], skip_special_tokens=True).replace(' ', '').rstrip('.')
    return output


def fwd(reactant: list, reagent: list) -> str:
    rc_smiles = reactant[0]
    re_smiles = reagent[0]

    for rc in reactant[1:]:
        rc_smiles += '.'+rc

    for re in reagent[1:]:
        re_smiles += '.'+re
    
    model_fwd = "sagawa/ReactionT5v2-forward"
    return rxn(f"REACTANT:{rc_smiles}REAGENT:{re_smiles}", model_fwd)


def retro(product: list) -> str:
    prod_smiles = product[0]

    for prod in product[1:]:
        prod_smiles += '.'+prod
    
    model_retro = "sagawa/ReactionT5v2-retrosynthesis"
    return rxn(prod_smiles, model_retro)


def visualization(reactant: str, reagent: str, product: str):
    rxn_smiles = f"{reactant}>{reagent}>{product}"
    rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=(800, 800))
    img.save(rxn_smiles+'.png')