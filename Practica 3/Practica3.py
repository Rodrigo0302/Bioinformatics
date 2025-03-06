#  Write code in Python to obtain the file 7FCD.pdb and 
#   calculate the centroid and the center of mass (COM) of each of the chains.
# Hint: Visit: https://tinyurl.com/y4y94zfj

from periodictable import C, H, N, O, S, P
import numpy as np


def calculate_centroid(coordinates):
    """Calculates the centroid of a set of coordinates."""
    return np.mean(coordinates, axis=0)


def calculate_center_of_mass(coordinates, masses):
    """Calculates the center of mass of a set of atoms."""
    weighted_coordinates = coordinates * masses[:, np.newaxis]
    return np.sum(weighted_coordinates, axis=0) / np.sum(masses)


def process_pdb(filename):
    """Processes a PDB file and calculates the centroid and COM for each chain."""

    atom_data = {}  # Dictionary to store atom data for each chain
    f = open(filename, "r")

    text = f.read().split("\n")

    f.close()

    LAtoms = []
    for line in text:
        if line.startswith("ATOM"):
            LAtoms.append(line)
    for line in LAtoms:
        fields = line.split()
        #print(fields)
        
        chain_id = fields[4]
        chain_id = chain_id[0]
        x = float(fields[5])
        y = float(fields[6])
        z = float(fields[7])
        element = fields[-1]
        
        if element == "C":
            mass = C.mass
        elif element == "H":
            mass = H.mass
        elif element == "N":
            mass = N.mass
        elif element == "O":
            mass = O.mass
        elif element == "S":
            mass = S.mass
        elif element == "P":
            mass = P.mass
        else:
            mass = 0.0
            
        if chain_id not in atom_data:
            atom_data[chain_id] = {"coordinates": [], "masses": []}
        atom_data[chain_id]["coordinates"].append(np.array([x, y, z]))
        atom_data[chain_id]["masses"].append(mass)
    
    for chain_id, data in atom_data.items():
        coordinates = np.array(data["coordinates"])
        #print(coordinates)
        masses = np.array(data["masses"])
        #print(masses)
        centroid = calculate_centroid(coordinates)
        com = calculate_center_of_mass(coordinates,masses)
        print(f"Cadena {chain_id}:")
        print(f"Centroide: {centroid}")
        print(f"Centro de Masa: {com}")
        
                        
                        
process_pdb("7FCD.pdb")