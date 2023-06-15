import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import sys

other_ligand = sys.argv[1]
top_ligand = sys.argv[2]

ref = mda.Universe(top_ligand)
mobile = mda.Universe(other_ligand)
align.alignto(mobile,ref,select="(element N)")
to_write = mobile.select_atoms('all')
other_ligand_beginning = other_ligand.split(".")[:-1] #discard the extension
if len(other_ligand_beginning) > 2: # if we had multiple . in file name and need to merge them together
    other_ligand_beginning_start = other_ligand_beginning[0]
    for i in range(1,len(other_ligand_beginning)-1):
        other_ligand_beginning_start += other_ligand_beginning[i]
    other_ligand_beginning = other_ligand_beginning_start
else:
    other_ligand_beginning = other_ligand_beginning[0]
to_write.write(f'{other_ligand_beginning}_current.mol2')
