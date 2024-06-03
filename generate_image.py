from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def generate_chemical_image(compound_name, output_path):
    try:
        # Muudab aine nime üheks aine molekuliks
        mol = Chem.MolFromSmiles(compound_name)
        if mol is None:
            return False, None
        
        
        img = Draw.MolToImage(mol) # Loob aine molekulist pildi
        img.save(output_path) # Salvestab pildi static/images kausta
        
        
        molecular_formula = CalcMolFormula(mol) # Arvutab brutovalemi
        
        return True, molecular_formula 
    except Exception as e:
        print(f"Error generating image for {compound_name}: {e}")
        return False, None
