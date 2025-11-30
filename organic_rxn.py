from chem_info import chem_info
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from rdkit.Chem import AllChem, Draw

def rxn(smiles, model_name):
    tokenizer = AutoTokenizer.from_pretrained(model_name, return_tensors='pt')
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)

    inp = tokenizer(smiles, return_tensors='pt')
    output = model.generate(**inp, num_beams=4, num_return_sequences=1, return_dict_in_generate=True, output_scores=True)
    output = tokenizer.decode(output['sequences'][0], skip_special_tokens=True).replace(' ', '').rstrip('.')
    return output


def compound_data(compound: str, query_type: str) -> str:
    c = chem_info(compound, query_type)
    return f"""
Name: {c.name}
Formula: {c.formula}
IUPAC name: {c.iupac_name}
Synonyms: {c.synonyms}
Molecular wt.: {c.weight}
Canonical SMILES: {c.cano_smiles}
Isomeric SMILES: {c.iso_smiles}
          """

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