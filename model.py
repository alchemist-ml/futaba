from transformers import AutoTokenizer, AutoModelForSeq2SeqLM

def rxn(smiles, model_name):
    tokenizer = AutoTokenizer.from_pretrained(model_name, return_tensors='pt')
    model = AutoModelForSeq2SeqLM.from_pretrained(model_name)

    inp = tokenizer(smiles, return_tensors='pt')
    output = model.generate(**inp, num_beams=4, num_return_sequences=1, return_dict_in_generate=True, output_scores=True)
    output = tokenizer.decode(output['sequences'][0], skip_special_tokens=True).replace(' ', '').rstrip('.')
    return output

def forward(smiles):
    model_fwd = "sagawa/ReactionT5v2-forward"
    return rxn(smiles, model_fwd)

def retro(smiles):
    model_retro = "sagawa/ReactionT5v2-retrosynthesis"
    return rxn(smiles, model_retro)
