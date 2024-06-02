from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def generate_chemical_image(compound_name, output_path):
    try:
        # Convert the compound name to a molecule
        mol = Chem.MolFromSmiles(compound_name)
        if mol is None:
            return False, None
        
        # Generate the image
        img = Draw.MolToImage(mol)
        img.save(output_path)
        
        # Calculate the molecular formula
        molecular_formula = CalcMolFormula(mol)
        
        return True, molecular_formula
    except Exception as e:
        print(f"Error generating image for {compound_name}: {e}")
        return False, None
