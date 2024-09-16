#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import sys
import os

def getBonds(points):
    """
    This function takes a list of 3D coordinates as input and returns the three closest atoms of each
    atom at a maximum distance of 1.9 angstrom.
    This allows to draw bonds if required.
    """
    from scipy.spatial import KDTree
    # Create a KDTree
    tree = KDTree(points)

    # Define max distance
    max_distance = 1.9

    # Find the 4 closest neighbours (including atom itself) within maximal distance
    distances, indices = tree.query(points, k=4, distance_upper_bound=max_distance)

    # Treat results
    bonds = []
    for i, (dist, idx) in enumerate(zip(distances, indices)):
        for d, j in zip(dist[1:], idx[1:]):  # Start at 1 to exclude atom itself
            if d <= max_distance:
                bonds.append([i, j])
    return bonds

def getCoordsFromXML(filename):
    """
    Returns labels and coordinates from an XML file where the coordinates are stored.
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    # Read atom labels and coordinates from XML
    atoms = []
    coords = []
    geometry_elem = root.find('Geometry')
    
    for atom_elem in geometry_elem.findall('Atom'):
        atom = atom_elem.find('Element').text
        x = float(atom_elem.find('X').text)
        y = float(atom_elem.find('Y').text)
        z = float(atom_elem.find('Z').text)
        atoms.append(atom)
        coords.append([y, x, z]) # X et Y sont inversés dans le plot :(

    return atoms, np.array(coords)

def getDataFromXML(filename):
    """
    Returns data extracted from the XML file.
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    # Extract origin
    origin_elem = root.find('Parameters/Origin')
    Origin = [float(origin_elem.find('X').text), float(origin_elem.find('Y').text), float(origin_elem.find('Z').text)]

    # Extract grid vectors X and Y
    vectX_elem = root.find('Parameters/GridXVector')
    VectX = [float(vectX_elem.find('X').text), float(vectX_elem.find('Y').text), float(vectX_elem.find('Z').text)]

    vectY_elem = root.find('Parameters/GridYVector')
    VectY = [float(vectY_elem.find('X').text), float(vectY_elem.find('Y').text), float(vectY_elem.find('Z').text)]

    # Extract number of points in X and Y
    nX = int(root.find('Parameters/NumberOfXPoints').text)
    nY = int(root.find('Parameters/NumberOfYPoints').text)

    # Extract isotropic values
    data_U = [float(value_elem.text) for value_elem in root.findall('U_IsotropicValues/Value')]
    data_R = None
    try:
        data_R = [float(value_elem.text) for value_elem in root.findall('R_IsotropicValues/Value')]
    except:
        pass

    return Origin, VectX, VectY, nX, nY, data_U, data_R

def main(molecule_dir, molecule_name):
    """
    Generate an image from the data stored in the XML file.
    """

    ims3d_color_scale = True  # Color scale for 3D maps
    show_atoms = True         # Show atoms in the plot
    show_bonds = True         # Show bonds between atoms in the plot
    save_as_image = True      # Save the generated image
    image_format = "png"      # Save the image as PNG

    # File paths
    xmlfile = "{}/{}.xml".format(molecule_dir, molecule_name)

    # Load coordinates from the XML file (reoriented molecule)
    labels, coords = getCoordsFromXML(xmlfile)
    bonds = getBonds(coords)

    # Load data from the XML file
    Origin, VectX, VectY, nX, nY, data_U, data_R = getDataFromXML(xmlfile)

    formalisms = ["U"]
    if data_R: formalisms.append("R")

    for formalism in formalisms:
        # Generate the grid
        X, Y = np.meshgrid(np.linspace(Origin[1], Origin[1]+nY*VectY[1], nY),
                           np.linspace(Origin[0], Origin[0]+nX*VectX[0], nX))

        if formalism=="U":
            Z = np.array(data_U).reshape(nX, nY)
        elif formalism=="R":
            Z = np.array(data_R).reshape(nX, nY)

        # Define the color scale
        if ims3d_color_scale:
            colors = [(0.5, 0.0, 0.0),
                      (0.8, 0.2, 0.2),
                      (1.0, 0.6, 0.6),
                      (1.0, 0.9, 0.8),
                      (0.6, 0.7, 0.95),
                      (0.0, 0.4, 0.85),
                      (0.0, 0.0, 0.4)]
            levels = [-22, -16.5, -11, -5.5, 5.5, 11, 16.5, 22]

        # Plot the data
        fig, ax = plt.subplots()

        # Adding contours
        if ims3d_color_scale:
            contours = ax.contourf(X, Y, Z, levels=levels, colors=colors)
        else:
            contours = ax.contourf(X, Y, Z)

        isocontours = ax.contour(X, Y, Z, levels=levels, colors=['black'], linewidths=1)
        plt.colorbar(contours)
        ax.autoscale(False)
        ax.set_xlim(np.min([np.min(X), np.min(Y)]), np.max([np.max(X), np.max(Y)]))

        # Adding the bonds if requested
        if show_bonds:
            for bond in bonds:
                iat = bond[0]
                jat = bond[1]
                if jat > iat:  # Bonds are present twice in the list, we treat only one of both cases
                    ax.plot([coords[iat][0], coords[jat][0]], [coords[iat][1], coords[jat][1]],
                            linewidth=3, color='gray', alpha=1)

        # Adding the atoms if requested
        if show_atoms:
            for lbl, at in zip(labels, coords):
                if 'C' in lbl.upper():
                    col = 'black'
                if 'H' in lbl.upper():
                    col = 'white'
                circle = Circle((at[0], at[1]), 0.2, edgecolor="white", facecolor=col, fill=True, alpha=.8, zorder=2)
                ax.add_patch(circle)

        # Save as image if requested
        image_dir = "img"
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)

        if save_as_image:
            imgfile = "{}/{}_{}.{}".format(image_dir, formalism, molecule_name, image_format)
            plt.savefig(imgfile)

        plt.show()

# Vérification des arguments de la ligne de commande
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python main.py <dir> <molecule_name>")
    else:
        main(sys.argv[1], sys.argv[2])

