import xml.etree.ElementTree as ET
import sys
import numpy as np
import pandas as pd
import os

def check_normal_termination(file_path):
    # Ouvrir le fichier et lire tout son contenu
    with open(file_path, 'r') as file:
        content = file.read()

    # Vérifier si la chaîne 'Normal termination' est présente dans le contenu
    if 'Normal termination of Gaussian' in content:
        return True
    else:
        return False

def extract_isotropic_data(gaussian_file_path):
    isotropic_values = []
    
    # Lire le fichier Gaussian ligne par ligne
    with open(gaussian_file_path, 'r') as file:
        for line in file:
            # Vérifier si la ligne contient "Isotropic"
            if "Isotropic" in line:
                # Extraire le 5e champ de la ligne (après avoir splitté sur des espaces multiples)
                fields = line.split()
                if len(fields) >= 5:
                    isotropic_values.append(float(fields[4]))  # Convertir le 5e champ en flottant
    
    return isotropic_values

def add_isotropic_data_to_xml(xml_file_path, isotropic_values, formalism="U"):
    # Charger le fichier XML
    tree = ET.parse(xml_file_path)
    root = tree.getroot()
    
    # Ajouter un nouvel élément pour les données isotropiques
    isotropic_elem = ET.SubElement(root, '{}_IsotropicValues'.format(formalism))
    
    # Ajouter chaque valeur isotropique sous forme d'un sous-élément <Value>
    for value in isotropic_values:
        ET.SubElement(isotropic_elem, 'Value').text = f"{value:.6f}"
    
    # Sauvegarder le fichier XML mis à jour
    tree.write(xml_file_path, encoding='utf-8', xml_declaration=True)


if __name__=="__main__":
    fichier_xyz = sys.argv[1]
    fichier_xml = fichier_xyz.replace(".xyz", ".xml")
    fichier_log = fichier_xyz.replace(".xyz", ".log")
    for formalism in ['U', 'R']:
        prefixed_fichier_log = "{}_{}".format(formalism, fichier_log)
        if check_normal_termination(prefixed_fichier_log):
            data = extract_isotropic_data(prefixed_fichier_log)
            add_isotropic_data_to_xml(fichier_xml, data, formalism)

