import os
import re
import argparse
from plot_ims import generer_image

def lister_couples(repertoire):
    """Liste les couples de fichiers .xyz et les fichiers R/U correspondants."""
    fichiers = os.listdir(repertoire)
    couples_fichiers = {}

    # Traiter d'abord tous les fichiers .xyz
    for fichier in fichiers:
        # Vérifie si c'est un fichier .xyz (radical peut contenir des chiffres, des lettres et des tirets)
        match_xyz = re.match(r'([-\w]+)_reopt\.xyz', fichier)
        if match_xyz:
            radical = match_xyz.group(1)
            # Initialise un dictionnaire pour ce radical avec .xyz rempli
            couples_fichiers[radical] = {'xyz': fichier, 'R': None, 'U': None}

    # Traiter ensuite les fichiers R_radical.csv et U_radical.csv
    for fichier in fichiers:
        match_r_u = re.match(r'(R|U)_([-\w]+)\.csv', fichier)
        if match_r_u:
            prefix, radical = match_r_u.groups()
            # Vérifier que le radical a été enregistré dans le premier passage (pour les fichiers .xyz)
            if radical in couples_fichiers:
                couples_fichiers[radical][prefix] = fichier
            else:
                # Si on trouve un fichier R ou U mais pas de fichier .xyz, on initialise tout de même une entrée
                couples_fichiers[radical] = {'xyz': None, 'R': None, 'U': None}
                couples_fichiers[radical][prefix] = fichier

    return couples_fichiers

def afficher_couples(couples_fichiers):
    """Affiche la liste des couples de fichiers .xyz avec leurs correspondants R/U."""
    print("Liste des couples de fichiers :\n")
    for radical, fichiers in couples_fichiers.items():
        print(f"Radical: {radical}")
        print(f"  .xyz: {fichiers['xyz'] if fichiers['xyz'] else 'Aucun fichier .xyz trouvé'}")
        print(f"  R_: {fichiers['R'] if fichiers['R'] else 'Aucun fichier R trouvé'}")
        print(f"  U_: {fichiers['U'] if fichiers['U'] else 'Aucun fichier U trouvé'}")
        print("\n")

def afficher_radicaux_sans_U(couples_fichiers):
    """Affiche les radicaux pour lesquels il manque un fichier U."""
    radicaux_sans_U = [radical for radical, fichiers in couples_fichiers.items() if fichiers['U'] is None]

    if radicaux_sans_U:
        print("Radicaux sans fichier U:")
        for radical in radicaux_sans_U:
            print(f"- {radical}")
    else:
        print("Tous les radicaux ont un fichier U.")

def main():
    parser = argparse.ArgumentParser(description="Lister les couples de fichiers ou vérifier les radicaux sans fichier U.")
    parser.add_argument("repertoire", help="Le répertoire contenant les fichiers.")
    parser.add_argument("-c", action="store_true", help="Vérifier s'il y a des radicaux sans fichier U.")
    parser.add_argument("-l", action="store_true", help="Lister les couples de fichiers.")
    parser.add_argument("-g", action="store_true", help="Génerer les images pour tous les couples")

    args = parser.parse_args()

    # Lister les couples de fichiers dans le répertoire
    couples_fichiers = lister_couples(args.repertoire)

    # Exécuter selon l'option choisie
    if args.c:
        afficher_radicaux_sans_U(couples_fichiers)
    elif args.l:
        afficher_couples(couples_fichiers)
    elif args.g:
        generer_image(args.repertoire, couples_fichiers)
    else:
        print("Veuillez choisir une option : '-c' pour vérifier les radicaux sans fichier U, '-l' pour lister les couples.")

if __name__ == "__main__":
    main()

