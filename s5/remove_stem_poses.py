import schrodinger
from schrodinger import structure
from schrodinger import structutils
from schrodinger.structutils import analyze
import numpy as np
import sys
import os
import glob

base_path = sys.argv[1]
top_ligand = sys.argv[2]

for file in os.listdir(base_path):
    if file[-3:] != ".sd":
        continue
    else:
        file = file[:-3]
    if top_ligand in file:
        continue    
    if '_' not in file: 
        continue
    underscore_divided = file.split('_')
    if len(underscore_divided) != 2:
        continue
    useful_part = underscore_divided[1]
    useful_part_split = useful_part.split('-')
    try:
        useful_part_split[0] = int(useful_part_split[0])
        useful_part_split[1] = int(useful_part_split[1])
        useful_part_split[2] = int(useful_part_split[2])
    except ValueError:
        continue

    # if we made it this far, it's a good file!

    # Acceptable range of residues
    start_acc = 16
    stop_acc = 45

    # Initialize acceptability of pose to false
    is_acceptable = False

    # Initialize dictionary {distance: resnum}
    dis_dict = {}

    # Initialize list of acceptable residue numbers
    acc_list_floats = np.linspace(start_acc, stop_acc, stop_acc-start_acc+1)
    acc_list = [int(i) for i in acc_list_floats]

    # Input files
    DNA_file = f"../s3/{underscore_divided[0]}_{useful_part_split[0]}-{str(useful_part_split[1]).rjust(4,'0')}.mae"
    Ligand_file = f"{base_path}/{file}.sd"

    # Import structures
    DNAst = structure.StructureReader.read(DNA_file)
    Ligst = structure.StructureReader.read(Ligand_file)

    # Get COM of ligand
    COM_L = schrodinger.structutils.analyze.center_of_mass(Ligst)

    for counter, residue in enumerate(DNAst.residue):
        # Get COM of residue
        COM_RES = schrodinger.structutils.analyze.center_of_mass(residue.extractStructure())
        # Get distance between residue COM and ligand COM
        vec_length = np.linalg.norm(COM_L - COM_RES)
        # Add values to dictionary
        #dis_dict[str(residue.resnum())] = vec_length
        dis_dict[vec_length] = residue.resnum

    # Sort the dictionary
    dict_keys = list(dis_dict.keys())
    dict_keys.sort()
    dis_dict_sorted = {i: dis_dict[i] for i in dict_keys}

    closest_residues = list(dis_dict_sorted.values()) # List of residues closest -> furthest
    closest_three = closest_residues[:3]

    # If at least one of the closest 3 is acceptable, pose is acceptable
    for residue_number in closest_three:
        if residue_number in acc_list:
            is_acceptable = True

    #print(f'base_path: {base_path}',flush=True)
    #print(f'top_ligand: {top_ligand}',flush=True)
    #print(is_acceptable,flush=True)
    #print(file,flush=True)
    # Remove unacceptable files
    if not is_acceptable:
        os.remove(f"{base_path}/{file}.sd")
