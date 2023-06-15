from schrodinger import structure
from schrodinger import structutils
from schrodinger.structutils import transform,minimize
from schrodinger.forcefield import minimizer
import numpy as np
import sys
import copy
import os

def rotate_receptor(receptor, P_indices,counter,kept_prev,ligand,receptor_name):
    #print('rotate_receptor called')
    ################################################
    ### USER NEEDS TO DEFINE START AND END OF LOOP
    #start = 20
    #end = 41
    ligand_st = structure.StructureReader.read(ligand)
    close_atoms = structutils.measure.get_atoms_close_to_structure(receptor,ligand_st,3)
    for atom_index in close_atoms:
        #print(atom_index)
        for atom_index2 in receptor.atom[atom_index].getResidue().getAtomIndices():
            if atom_index2 not in close_atoms:
                close_atoms.append(atom_index2)
    #print(close_atoms)
    #exit()
    #################################################
 
    to_rotate = np.random.uniform(-.5,.5,counter)
    np.save('to_rotate.npy',to_rotate)
    if kept_prev:
        to_rotate = np.load('to_rotate.npy')
    # roughly corresponding to angle in radians
    # the real angle is probably a little bigger
    #print(len(to_rotate))
    #print(len(P_indices))
    for i in range(2,len(P_indices)):
    #for i in range(1,len(P_indices)):
        if P_indices[i] not in close_atoms:
            continue
        # we aren't going to rotate the first or last residues
        # should change the 1 to a 2 if 5' end is phosphorylated
        # i think this will also avoid rotating the 2nd to last 
        # residue (on the 3' end), assuming the 3' end isn't phosphorylated
        # ACTUALLY, IT IS ROTATING THE 2ND TO LAST RESIDUE
        diff = np.array(receptor.atom[P_indices[i]].xyz) - np.array(receptor.atom[P_indices[i-1]].xyz)
        # now we need an arbitrary vector that is perpendicular to this one
        ortho = [1,1,-(diff[0]+diff[1])/diff[2]]
        ortho = np.array(ortho)/np.linalg.norm(ortho)
        new_vector = np.cross(diff,ortho)
        new_vector = new_vector/np.linalg.norm(new_vector)
        norm_diff = new_vector-ortho
        norm_diff = norm_diff/np.linalg.norm(norm_diff)
        #theta = to_rotate[i-1] 
        theta = to_rotate[i-2]
        factor = theta*2*(np.sqrt(2)/np.pi)
        # formula not actually rotating by the specified angle, but it's close
        rotated_ortho = ortho+factor*norm_diff
        rotated_ortho = rotated_ortho/np.linalg.norm(rotated_ortho)
        matrix = transform.get_alignment_matrix(ortho,rotated_ortho)
        old_coords = copy.deepcopy(receptor.atom[P_indices[i-1]].xyz)
        transform.translate_to_origin(receptor,receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        transform.transform_structure(receptor,matrix,receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        new_coords = receptor.atom[P_indices[i-1]].xyz
        P_vector = [old_coords[index] - new_coords[index] for index in range(len(new_coords))]
        transform.translate_structure(receptor,P_vector[0],P_vector[1],P_vector[2],receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        #print(receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
    
    
    for counter,residue in enumerate(receptor.residue):
        if residue.getAtomIndices()[0] not in close_atoms:
            # choice of index was arbitrary, because all atoms from residue should be included
            continue

        #if counter <start or counter > end:
        #    continue
        for atom1 in residue.getAtomIndices():
            if receptor.atom[atom1].pdbname.strip() == "C1'":
                first = receptor.atom[atom1]
                for atom2 in first.bonded_atoms:
                    if atom2.element == "N":
                        second = atom2
                    else:
                        third = atom2
                for atom3 in second.bonded_atoms:
                    if atom3 != first:
                        fourth = atom3
        adjust_amount = np.random.uniform(-30,30) # degrees
        dihedral = receptor.measure(third,first,second,fourth)
        receptor.adjust(dihedral+adjust_amount,third,first,second,fourth)
    
    for i in range(1,len(P_indices)):
        if P_indices[i] not in close_atoms:
            continue

        #if i<start or i>end:
        #    continue
        diff = np.array(receptor.atom[P_indices[i]].xyz) - np.array(receptor.atom[P_indices[i-1]].xyz)
        # we started like the base rotation loop
        # now, we need a random vector with magnitude 1
        if np.linalg.norm(diff) >= 8:
            minimizer.minimize_substructure(receptor,close_atoms)
            diff = np.array(receptor.atom[P_indices[i]].xyz) - np.array(receptor.atom[P_indices[i-1]].xyz)
        while_loop_counter = 0
        while True:
            while_loop_counter += 1
            unit = np.random.rand(3) - 0.5 
            # random float in [-0.5,0.5)
            unit = unit/np.linalg.norm(unit) * 2 
            # scaling to larger amplitude phosphate movement
            new_diff = diff + unit
            #new_diff = diff
            if np.linalg.norm(new_diff) < 8:
                break # we don't want to stretch out our bases too much
            elif while_loop_counter == 2000: # don't want to get stuck in infinite loop
                break
        old_coords = copy.deepcopy(receptor.atom[P_indices[i-1]].xyz)
        alignment_matrix = transform.get_alignment_matrix(diff,new_diff)
        # we need to scale the result by some constants, so we'll do this
        # by modifying alignment_matrix element-wise
        #scale_factor = np.linalg.norm(new_diff)/np.linalg.norm(alignment_matrix,axis=1)
        #print(scale_factor.shape)
        for row in range(alignment_matrix.shape[1]):
            alignment_matrix[row,:] = alignment_matrix[row,:]/np.linalg.norm(diff)
            alignment_matrix[row,:] = alignment_matrix[row,:]*np.linalg.norm(new_diff)
        transform.translate_to_origin(receptor,receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        origin_coords = copy.deepcopy(receptor.atom[P_indices[i-1]].xyz)
        transform.transform_structure(receptor,alignment_matrix,receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        #new_coords = receptor.atom[P_indices[i-1]].xyz
        P_vector = [old_coords[index] - origin_coords[index] for index in range(len(origin_coords))]
        #print(P_vector)
        
        transform.translate_structure(receptor,P_vector[0]+unit[0],P_vector[1]+unit[1],P_vector[2]+unit[2],receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
        #transform.translate_structure(receptor,unit[0]+old_coords[0],unit[1]+old_coords[1],unit[2]+old_coords[2],receptor.atom[P_indices[i-1]].getResidue().getAtomIndices())
    

    # sometimes, the stereochemistry is messed up here, so we're going to fix it
    for i in range(1,len(receptor.atom)+1):
        if receptor.atom[i].pdbname.strip() == "C1'":
            if receptor.atom[i].chirality != "R":
                for other_atom in receptor.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-receptor.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]
        elif receptor.atom[i].pdbname.strip() == "C3'":
            if receptor.atom[i].chirality != "S":
                for other_atom in receptor.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-receptor.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]
        elif receptor.atom[i].pdbname.strip() == "C4'":
            if receptor.atom[i].chirality != "R":
                for other_atom in receptor.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-receptor.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]

        
    
    #minimize.minimize_structure(receptor,max_steps=200) 
    minimizer.minimize_substructure(receptor,close_atoms)
    minimize.minimize_structure(receptor,max_steps=100)
    #receptor.write('foo2.mae')
    #exit()
    receptor.write(f'{receptor_name[:-4]}_intermediate.mol2') # receptor is a maestro extension file
    receptor.write(f'{receptor_name[:-4]}_intermediate.mae')
    # will evaluate with rbdock, keep if it's better 
    return


run = sys.argv[1] # this is actually the ligand we should be using 
receptor_name = sys.argv[2] 
kept_prev = sys.argv[3]
# top_ligand = sys.argv[4] doesn't matter anymore
#AFFINITY_INFO = sys.argv[5]
# AFFINITY_INFO has the information about the rest of the ligands

'''
# we don't care about other ligand affinities in this script,
# just what their names are so we can dock them
# actually we don't care about that
ligands = AFFINITY_INFO.split(',') 
ligands = [ligand.split(':')[0] for ligand in ligands]
for i in range(len(ligands)):
    if ligands[i][0] == "'":
        ligands[i] = ligands[i][1:]
    if ligands[i][-1] == "'":
        ligands[i] = ligands[i][:-1]
# trimming the literal ' in case they snuck in from the bash script
'''

receptor = structure.StructureReader.read(receptor_name)

counter = 0
for residue in receptor.residue:
    counter += 1
counter -= 2 # we don't rotate the first or last residue
counter2 = 0


P_indices = []
for atom in receptor.atom:
    counter2 += 1
    if atom.element.strip() == 'P':
        P_indices.append(counter2)
rotate_receptor(receptor,P_indices,counter,kept_prev,run,receptor_name)
