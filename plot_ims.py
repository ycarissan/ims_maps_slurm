import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import KDTree

def getBonds(points):
    """
    This function takes a liste of 3D coordiantes as input and returns the three closest atoms of each
    atom at a maximum distance of 1.9 angstrom.
    This allows to draw bonds if required.
    """

# Create a KDTree
    tree = KDTree(points)

# Define max distance
    max_distance = 1.9

# Find the 4 closests neigbourg (including atom itself)
# within maximal distance
    distances, indices = tree.query(points, k=4, distance_upper_bound=max_distance)

# Treat results
    bonds = []
    for i, (dist, idx) in enumerate(zip(distances, indices)):
        for d, j in zip(dist[1:], idx[1:]):  # Start at 1 to exclude atom itself
            if d <= max_distance:
                bonds.append([i,j])
    return bonds

def getCoords(geomfile):
    """
    Returns labels and coordinates from an xyz file
    """
    thrs=1e-6
    with open(geomfile) as fd:
        coords = []
        labels = []
        for line in fd.readlines()[2:]:
            lbl, x, y, z = line.split()
            labels.append(lbl)
            coords.append([float(x), float(y), float(z)])
        coords=np.array(coords)

    ptp = np.ptp(coords, axis=0)
    if (ptp[0]<thrs):
        # 0, 1, 2
        coords[:, [0,2]] = coords[:, [2,0]]
        # 2, 1, 0
    #    coords[:, [0,1]] = coords[:, [1,0]]
        # 1, 2, 0
    elif (ptp[1]<thrs):
        # 0, 1, 2
        coords[:, [0,1]] = coords[:, [1,0]]
        # 1, 0, 2
        coords[:, [0,2]] = coords[:, [2,0]]
        # 2, 0, 1
    return labels, coords

def read_grid_from_xml(file_path):
    # Charger le fichier XML
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    # Lire les paramètres de la grille
    params = root.find('Parameters')
    
    # Origine de la grille
    origin_elem = params.find('Origin')
    origin = np.array([
        float(origin_elem.find('X').text),
        float(origin_elem.find('Y').text),
        float(origin_elem.find('Z').text)
    ])
    
    # Vecteur directeur selon X
    grid_x_elem = params.find('GridXVector')
    grid_x_vector = np.array([
        float(grid_x_elem.find('X').text),
        float(grid_x_elem.find('Y').text),
        float(grid_x_elem.find('Z').text)
    ])
    
    # Vecteur directeur selon Y
    grid_y_elem = params.find('GridYVector')
    grid_y_vector = np.array([
        float(grid_y_elem.find('X').text),
        float(grid_y_elem.find('Y').text),
        float(grid_y_elem.find('Z').text)
    ])
    
    # Nombre de points selon X et Y
    num_x_points = int(params.find('NumberOfXPoints').text)
    num_y_points = int(params.find('NumberOfYPoints').text)
    
    # Lire la matrice de rotation
    rotation_elem = root.find('RotationMatrix')
    rotation_matrix = np.zeros((3, 3))
    for i in range(3):
        row_elem = rotation_elem.find(f'Row{i+1}')
        for j in range(3):
            rotation_matrix[i, j] = float(row_elem.find(f'Col{j+1}').text)
    
    # Lire la géométrie réorientée de la molécule
    geometry_elem = root.find('Geometry')
    atoms = []
    coords = []
    for atom_elem in geometry_elem.findall('Atom'):
        atom = atom_elem.find('Element').text
        x = float(atom_elem.find('X').text)
        y = float(atom_elem.find('Y').text)
        z = float(atom_elem.find('Z').text)
        atoms.append(atom)
        coords.append([x, y, z])
    
    # Convertir la géométrie en un tableau NumPy
    reoriented_coords = np.array(coords)
    
    # Lire les valeurs isotropiques (si elles existent)
    isotropic_values = []
    isotropic_elem = root.find('IsotropicValues')
    if isotropic_elem is not None:
        for value_elem in isotropic_elem.findall('Value'):
            isotropic_values.append(float(value_elem.text))
    
    # Retourner les données lues
    return {
        'origin': origin,
        'grid_x_vector': grid_x_vector,
        'grid_y_vector': grid_y_vector,
        'num_x_points': num_x_points,
        'num_y_points': num_y_points,
        'rotation_matrix': rotation_matrix,
        'atoms': atoms,
        'reoriented_coords': reoriented_coords,
        'isotropic_values': isotropic_values  # Ajout des valeurs isotropiques
    }

# In[3]:


def generer_image(dirname, couples):
    
    for radical, fichiers in couples.items():
        molecule_dir=dirname
        molecule_name=radical
        for a in ["R", "U"]:
            if fichiers[a]:
                print("{} ".format(fichiers[a]))
                try:
                    map_type=a
            
                    ims3d_color_scale=True #True or False
                    
                    show_atoms=True #True or False
                    show_bonds=True #True or False
                    
                    save_as_image=True # True or False
                    image_format="png" # "png" by default
                    
                    geomfile="{}/{}_reopt.xyz".format(molecule_dir,molecule_name)
                    xmlfile="{}/{}_{}.xml".format(molecule_dir, map_type, molecule_name)

                    image_dir="img"
                    
                    imgfile="{}/{}_{}.{}".format(image_dir, map_type, molecule_name, image_format)
                    if not os.path.isfile(imgfile):
                        print("being done")

                        labels, coords = getCoords(geomfile)
                        bonds=getBonds(coords)
                        d = read_grid_from_xml(xmlfile)
                        Origin, VectX, VectY, nX, nY, data = d["origin"], d["grid_x_vector"], d["grid_y_vector"], d["num_x_points"], d["num_y_points"], d["isotropic_values"]    
    
                        X, Y = np.meshgrid(np.linspace(Origin[1], Origin[1]+nY*VectY[1], nY)
                                           , np.linspace(Origin[0], Origin[0]+nX*VectX[0], nX))
                        
                        Z = np.array(data).reshape(nX,nY)
                        
                        #
                        # Definition of the color scale
                        #
                        if ims3d_color_scale:
                            colors = [(0.5, 0.0, 0.0),
                                      (0.8, 0.2, 0.2),
                                      (1.0, 0.6, 0.6),
                                      (1.0, 0.9, 0.8),
                                      (0.6, 0.7, 0.95),
                                      (0.0, 0.4, 0.85),
                                      (0.0, 0.0, 0.4)]
                        
                        levels = [-22, -16.5, -11, -5.5, 5.5, 11, 16.5, 22]
                        
                        # plot
                        fig, ax = plt.subplots()
                        
                        #
                        # Adding contours
                        #
                        if ims3d_color_scale:
                            contours = ax.contourf(X, Y, Z, levels=levels, colors=colors)
                        else:
                            contours = ax.contourf(X, Y, Z)
                            
                        isocontours = ax.contour(X, Y, Z, levels=levels, colors=['black'], linewidths=1)
                        plt.colorbar(contours)
                        ax.autoscale(False)
                        ax.set_xlim(np.min([np.min(X), np.min(Y)]), np.max([np.max(X), np.max(Y)]))
                        
                        #
                        # Adding the bonds if requested
                        #
                        if show_bonds:
                            for bond in bonds:
                                iat = bond[0]
                                jat = bond[1]
                                if jat > iat: # Bonds are present twice in the list, we treat only one of both cases
                                    ax.plot([coords[iat][0], coords[jat][0]],[coords[iat][1], coords[jat][1]],
                                            linewidth=3, color='gray', alpha=1)
                        
                        #
                        # Adding the atoms if requested
                        #
                        if show_atoms:
                            for lbl, at in zip(labels, coords):
                                if 'C' in lbl.upper():
                                    col = 'black'
                                if 'H' in lbl.upper():
                                    col = 'white'
                                circle = Circle((at[0], at[1]), 0.2, edgecolor="white", facecolor=col, fill=True, alpha=.8, zorder=2)
                                ax.add_patch(circle)
                        
                        #
                        # Save as image if requested
                        #
                        if save_as_image:
                            plt.savefig(imgfile)
                    else:
                        print("img exists")
                except:
                    with open("errors","a")as fn:
                        fn.write("{}\n".format(fichiers[a]))
                with open("done","a")as fn:
                    fn.write("{}\n".format(fichiers[a]))
