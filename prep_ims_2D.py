import xml.etree.ElementTree as ET
import sys
import numpy as np
import pandas as pd
import os

# Étape 1 : Lecture du fichier XYZ
def read_xyz(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[2:]  # On saute les deux premières lignes (commentaires)
    atoms = []
    coordinates = []
    for line in lines:
        parts = line.split()
        atom = parts[0]
        coord = list(map(float, parts[1:]))
        atoms.append(atom)
        coordinates.append(coord)
    return np.array(atoms), np.array(coordinates)

# Étape 2 : Calcul du plan moyen de la molécule
def compute_mean_plane(coordinates):
    centroid = np.mean(coordinates, axis=0)  # Calcul du centre de masse
    centered_coords = coordinates - centroid
    
    # Calcul par SVD pour obtenir le vecteur normal du plan moyen
    _, _, vh = np.linalg.svd(centered_coords)
    normal_vector = vh[-1]  # Le vecteur normal est la dernière colonne de Vh
    
    return centroid, normal_vector

# Étape 3 : Réorientation de la molécule
def reorient_molecule(coordinates, normal_vector):
    z_axis = np.array([0, 0, 1])
    if np.allclose(normal_vector, z_axis):
        return coordinates, np.eye(3)  # Pas de réorientation nécessaire

    rotation_axis = np.cross(normal_vector, z_axis)
    rotation_angle = np.arccos(np.dot(normal_vector, z_axis))
    rotation_matrix = rotation_matrix_from_axis_angle(rotation_axis, rotation_angle)
    
    reoriented_coords = np.dot(coordinates, rotation_matrix.T)
    return reoriented_coords, rotation_matrix

def rotation_matrix_from_axis_angle(axis, angle):
    axis = axis / np.linalg.norm(axis)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    ux, uy, uz = axis
    return np.array([
        [cos_a + ux**2 * (1 - cos_a), ux * uy * (1 - cos_a) - uz * sin_a, ux * uz * (1 - cos_a) + uy * sin_a],
        [uy * ux * (1 - cos_a) + uz * sin_a, cos_a + uy**2 * (1 - cos_a), uy * uz * (1 - cos_a) - ux * sin_a],
        [uz * ux * (1 - cos_a) - uy * sin_a, uz * uy * (1 - cos_a) + ux * sin_a, cos_a + uz**2 * (1 - cos_a)]
    ])

# Étape 4 : Vérification de la planéité de la molécule
def is_molecule_flat(coordinates, normal_vector, threshold=0.5):
    distances = np.dot(coordinates, normal_vector)
    return np.max(np.abs(distances)) < threshold

# Étape 5 : Génération de la grille
def generate_grid(coordinates, margin=0.5, step=0.2, z_value=1.0):
    x_min, y_min = np.min(coordinates[:, :2], axis=0) - margin
    x_max, y_max = np.max(coordinates[:, :2], axis=0) + margin
    
    x_points = np.arange(x_min, x_max, step)
    y_points = np.arange(y_min, y_max, step)
    
    grid_points = np.array([[x, y, z_value] for x in x_points for y in y_points])
    
    return grid_points, x_points, y_points

# Étape 6 : Sauvegarder les résultats dans un fichier xml

def save_grid_to_xml(file_path, grid_points, x_points, y_points, step, reoriented_coords, atoms, rotation_matrix):
    # Créer l'élément racine <Grid>
    root = ET.Element('Grid')

    # Ajouter les paramètres de la grille
    params = ET.SubElement(root, 'Parameters')
    
    # Origine de la grille
    origin = grid_points[0]
    origin_elem = ET.SubElement(params, 'Origin')
    ET.SubElement(origin_elem, 'X').text = f"{origin[0]:.6f}"
    ET.SubElement(origin_elem, 'Y').text = f"{origin[1]:.6f}"
    ET.SubElement(origin_elem, 'Z').text = f"{origin[2]:.6f}"

    # Vecteur directeur selon X
    grid_x_vector = np.array([step, 0, 0])
    grid_x_elem = ET.SubElement(params, 'GridXVector')
    ET.SubElement(grid_x_elem, 'X').text = f"{grid_x_vector[0]:.6f}"
    ET.SubElement(grid_x_elem, 'Y').text = f"{grid_x_vector[1]:.6f}"
    ET.SubElement(grid_x_elem, 'Z').text = f"{grid_x_vector[2]:.6f}"

    # Vecteur directeur selon Y
    grid_y_vector = np.array([0, step, 0])
    grid_y_elem = ET.SubElement(params, 'GridYVector')
    ET.SubElement(grid_y_elem, 'X').text = f"{grid_y_vector[0]:.6f}"
    ET.SubElement(grid_y_elem, 'Y').text = f"{grid_y_vector[1]:.6f}"
    ET.SubElement(grid_y_elem, 'Z').text = f"{grid_y_vector[2]:.6f}"

    # Nombre de points selon X et Y
    ET.SubElement(params, 'NumberOfXPoints').text = str(len(x_points))
    ET.SubElement(params, 'NumberOfYPoints').text = str(len(y_points))

    # Ajouter la matrice de rotation utilisée pour réorienter la molécule
    rotation_elem = ET.SubElement(root, 'RotationMatrix')
    for i in range(3):
        row_elem = ET.SubElement(rotation_elem, f'Row{i+1}')
        for j in range(3):
            ET.SubElement(row_elem, f'Col{j+1}').text = f"{rotation_matrix[i, j]:.6f}"

    # Ajouter la géométrie réorientée de la molécule
    geometry = ET.SubElement(root, 'Geometry')
    for atom, coord in zip(atoms, reoriented_coords):
        atom_elem = ET.SubElement(geometry, 'Atom')
        ET.SubElement(atom_elem, 'Element').text = atom
        ET.SubElement(atom_elem, 'X').text = f"{coord[0]:.6f}"
        ET.SubElement(atom_elem, 'Y').text = f"{coord[1]:.6f}"
        ET.SubElement(atom_elem, 'Z').text = f"{coord[2]:.6f}"

    # Convertir l'arbre XML en une chaîne et l'enregistrer dans un fichier
    tree = ET.ElementTree(root)
    xml_file = os.path.splitext(file_path)[0] + '.xml'
    tree.write(xml_file, encoding='utf-8', xml_declaration=True)

# Étape 7 : Générer un fichier .com pour Gaussian avec grille
def save_gaussian_com(file_path, atoms, reoriented_coords, grid_points):

    charge=0
    # Compter le nombre d'atomes de carbone ('C')
    num_carbon = sum(1 for atom in atoms if atom == 'C')

    # Définir la multiplicité : 1 si le nombre de carbones est pair, 2 sinon
    multiplicity = 1 if num_carbon % 2 == 0 else 2
    multiplicity = 1

    com_file = os.path.splitext(file_path)[0] + '.com'  # Remplacer l'extension par .com
    
    header = f"""%nproc=16
%mem=32GB
#P Ub3lyp/6-311++G(d,p) SCF(Tight) CPHF(Separate) Int(Grid=SuperFine) NMR geom=connectivity

Titre

{charge} {multiplicity}"""
    
    geometry = "\n".join([f"{atom} {x:.6f} {y:.6f} {z:.6f}" for atom, (x, y, z) in zip(atoms, reoriented_coords)])
    
    grid_geometry = "\n".join([f"Bq {x:.6f} {y:.6f} {z:.6f}" for (x, y, z) in grid_points])
    
    total_centers = len(atoms) + len(grid_points)
    
    integer_list = "\n".join([str(i) for i in range(1, total_centers + 1)])
    
    with open(com_file, 'w') as file:
        file.write(header + "\n" + geometry + "\n" + grid_geometry + "\n\n" + integer_list + "\n\n")

# Fonction principale pour traiter un fichier XYZ et générer la grille et le fichier .com
def process_molecule_xyz(file_path, output_xml):
    atoms, coordinates = read_xyz(file_path)
    centroid, normal_vector = compute_mean_plane(coordinates)
    
    if not is_molecule_flat(coordinates, normal_vector):
        raise ValueError("La molécule n'est pas suffisamment plate.")
    
    reoriented_coords, rotation_matrix = reorient_molecule(coordinates - centroid, normal_vector)
    
    grid_points, x_points, y_points = generate_grid(reoriented_coords, z_value=1.0)
    
    save_grid_to_xml(output_xml, grid_points, x_points, y_points, 0.2, reoriented_coords, atoms, rotation_matrix)
    
    save_gaussian_com(file_path, atoms, reoriented_coords, grid_points)
    
    return rotation_matrix

# Exemple d'appel de la fonction principale
# process_molecule_xyz('molecule.xyz', 'grid_output.xml')

if __name__=="__main__":
    fichier_xyz = sys.argv[1]
    fichier_xml = fichier_xyz.replace(".xyz", ".xml")
    process_molecule_xyz(fichier_xyz, fichier_xml)

