import pubchempy as pcp
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from rdkit.Chem import AllChem, Draw
from flask import Flask, request, jsonify
from flask_cors import CORS
import io
import base64

app = Flask(__name__)
CORS(app)

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
        name = compound.synonyms[0] if compound.synonyms else query
    elif query_type == 'smiles':
        formula = compound.molecular_formula
        name = compound.synonyms[0] if compound.synonyms else query

    iupac = compound.iupac_name
    synonyms = compound.synonyms[1:4] if compound.synonyms else []
    weight = compound.molecular_weight
    cano_smiles = compound.connectivity_smiles
    iso_smiles = compound.smiles

    return Chemical(name, formula, iupac, synonyms, weight, cano_smiles, iso_smiles)

def rxn(smiles, model_name):
    if "retrosynthesis" in model_name.lower():
        from chemformer_helper import chemformer
        return chemformer.predict(smiles).replace(' ', '').rstrip('.')
    
    tokenizer = AutoTokenizer.from_pretrained(model_name, return_tensors='pt')
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)

    inp = tokenizer(smiles, return_tensors='pt')
    output = model.generate(**inp, num_beams=4, num_return_sequences=1, return_dict_in_generate=True, output_scores=True)
    output = tokenizer.decode(output['sequences'][0], skip_special_tokens=True).replace(' ', '').rstrip('.')
    return output

def fwd(reactant: list, reagent: list) -> str:
    if not reactant:
        return "No reactants provided"
    
    rc_smiles = reactant[0]
    re_smiles = reagent[0] if reagent else ""

    for rc in reactant[1:]:
        rc_smiles += '.'+rc

    for re in reagent[1:]:
        re_smiles += '.'+re
    
    model_fwd = "sagawa/ReactionT5v2-forward"
    return rxn(f"REACTANT:{rc_smiles}REAGENT:{re_smiles}", model_fwd)

def create_reaction_image(reactants, products):
    try:
        if not reactants or not products:
            return None
        
        reactant_smiles = '.'.join(reactants)
        product_smiles = '.'.join(products)
        
        rxn_smiles = f"{reactant_smiles}>>{product_smiles}"
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        
        img = Draw.ReactionToImage(rxn, subImgSize=(600, 300))
        
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_base64 = base64.b64encode(img_bytes.getvalue()).decode('utf-8')
        
        return img_base64
    except Exception as e:
        print(f"Image creation error: {e}")
        return None

@app.route('/chem_info', methods=['POST'])
def get_chem_info():
    data = request.json
    query = data['query']
    query_type = data['query_type']
    chem = chem_info(query, query_type)
    if isinstance(chem, str):
        return jsonify({"error": chem})
    return jsonify({
        "name": chem.name,
        "formula": chem.formula,
        "iupac_name": chem.iupac_name,
        "synonyms": chem.synonyms,
        "weight": str(chem.weight),
        "cano_smiles": chem.cano_smiles,
        "iso_smiles": chem.iso_smiles
    })

@app.route('/predict', methods=['POST'])
def predict():
    data = request.json
    reactants = data['reactants']
    reagents = data.get('reagents', [])
    
    reaction_products = fwd(reactants, reagents)
    
    image_base64 = create_reaction_image(reactants, reaction_products.split('.'))
    
    return jsonify({
        "reaction": reaction_products,
        "image_base64": image_base64
    })

if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=5000)