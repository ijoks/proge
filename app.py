from flask import Flask, render_template, request, send_from_directory
import os
import csv
from generate_image import generate_chemical_image

app = Flask(__name__) # Loob Flaski veebilehe

# Teeb kindlaks, kas static ja images kaust on olemas
if not os.path.exists('static/images'):
    os.makedirs('static/images')

# Uus sõnastik, kus on ainete nimed eesti keeles
translations = {}
with open('compound_translations.csv', mode='r', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile) # Loeb faili
    next(reader)  # Jätab esimese rea vahele
    for row in reader:
        if len(row) != 2:
            continue  # Jätab vahele read, millel pole kahte tulpa
        estonian_name, smiles = row 
        translations[estonian_name.strip().lower()] = smiles.strip() # Paneb sõnastikku

 # Funktsioon vanade piltide kustutamiseks   
def delete_old_images(directory='static/images'):
    for filename in os.listdir(directory): # Teeb listi kõikidest failidest kaustas
        file_path = os.path.join(directory, filename) # Loob faili aadressi
        try:
            if os.path.isfile(file_path): # Kontrollib kas faili address viib failini, mitte näiteks kaustani
                os.unlink(file_path) # Kustutab faili
        except Exception as e:
            print(f"Error deleting file {file_path}: {e}") # Prindib errori kui faili ei saa kustutada

@app.route('/', methods=['GET']) # Loob kodulehe
def index():
    return render_template('index.html') # Võtab disaini index.html failist (ja ka styles.css failist)

@app.route('/generate', methods=['POST']) # Loob lehe, kus näidata valemit
def generate():
    compound = request.form['compound'].strip().lower() # Võtab kasutaja sisestatud aine nimetuse
    
    # Kontrollib, kas sisend on sõnastikus olemas, kui on siis tõlgib SMILES formaati
    if compound in translations:
        compound = translations[compound] 

    image_path = os.path.join('static', 'images', f'{compound}.png') # Loob tee aine graafilise valemini

    delete_old_images() # Kustutab vanad pildid kaustast

    # Proovib luua pilti
    success, molecular_formula = generate_chemical_image(compound, image_path)
    
    if success:
        return render_template('index.html', image_url=image_path, molecular_formula=molecular_formula) # Kui pilt õnnestus luua, näitab seda veebilehel koos brutovalemiga
    else:
        error_message = "Vigane sisend. Palun proovi uuesti ;)" 
        return render_template('index.html', error_message=error_message) # Kui pilti ei õnnestunud luua, näitab veebilehel errorit

if __name__ == '__main__': 
    app.run(debug=True) # Käivitab veebilehe ja lülitab sisse debuggingu'u