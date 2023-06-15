''' going to leav this out of the algorithm for now.
    can come back to it later if we want'''

import os
import subprocess
import MDAnalysis as mda
import numpy as np
import copy
import sys
import json
from glob import glob
from MDAnalysis.analysis.distances import distance_array

def get_inter_score(file,pose_index):   
    # file: relative path, assuming we're in {BASE_PATH}/tertiary_tructures/{GENERATION}
    # pose_index: 0-indexed location of the pose that we want to score in the sd file
    with open(file,'r') as f:
        prev_indices = []
        scores = []
        for counter,line in enumerate(f):
            if '>  <SCORE.INTER>' in line:
                prev_indices.append(counter)
    with open(file,'r') as f:
        for counter,line in enumerate(f):
            if counter-1 in prev_indices:
                scores.append(float(line.strip()))
    score = scores[pose_index]
    return (score,pose_index)

def check_affinity_trend(affinity_list,top_ligand_score,alternate_ligand_scores,strict=False):
    # affinity_list: imported from .json file
    # top_ligand_score: float
    # alternate_ligand_scores: dict where the keys are the alternate ligands
    # values are lists of tuples
    # each tuple in the list represents an alternate ligand pose
    # the first element in the tuple is the score associated with the alternate pose
    # the second element is the 0-indexed location of the pose within the sd file for that ligand
    # strict: if True, clusters that don't contain every ligand in affinity_list
    # will be regarded as breaking the trend (False returned)
    # if False, clusters that don't contain every ligand will only be evaluated based on 
    # ligands that they do contain (may return True or False)
    
    if len(alternate_ligand_scores.keys()) == 0:
        if strict:
            return False  
        else:
            return True


    # first, we'll check that top_ligand_score is lower than all others
    for key in alternate_ligand_scores.keys():
        for item in alternate_ligand_scores[key]:
            if item[0] < top_ligand_score: # top_ligand_score not the best
                return False

    # if we've made it this far, top_ligand_score is the best
    for i in range(1,len(affinity_list)-1): # no need to compare the last element of affinity_list
        for ligand1 in affinity_list[i]:
            for j in range(i+1,len(affinity_list)):
                for ligand2 in affinity_list[j]:
                    # this ligand might not be included within the cluster
                    if ligand2 in alternate_ligand_scores.keys():
                        for item in alternate_ligand_scores[ligand2]:
                            if item[0] < min([score_tuple[0] for score_tuple in alternate_ligand_scores[ligand1]]):
                                return False
                    else: # might want to do some smart counting at some point to check our sampling 
                          # ideally, we would have all ligands in every cluster, but realistically,
                          # this probably won't happen often
                        if strict:
                            return False 
                        else:
                            pass
    # if we've made it this far, our ordering is perfect!
    return True

BASE_PATH = sys.argv[1]
GENERATION = sys.argv[2]
affinity_info = sys.argv[3][1:-1].split(',') # cut out the curly braces, make into list
affinity_info = [element.split(':') for element in affinity_info]
#affinity_info = {key:value for [key,value] in affinity_info}
#os.chdir(f'{BASE_PATH}/tertiary_structures/{GENERATION}')
os.chdir(f'{BASE_PATH}/s4/') # should already be there, but checking

#with open(f'{BASE_PATH}/affinity_list.json','r') as f:
#    affinity_list = json.load(f)
#assert len(affinity_list[0]) == 1 # need a single best ligand

#############################################################3
###############################################################
# now we have to decide how to process affinity info.
# the current method is crude but good enough for now
final_affinity_list = []
affinities = [name_affinity_tuple[1] for name_affinity_tuple in affinity_info]
names =      [name_affinity_tuple[0] for name_affinity_tuple in affinity_info]
bottom_tier = []
keep_indices = []
for i in range(len(names)):
    if affinities[i] >= 1000:
        bottom_tier.append(names[i]) # we'll tack this on at the end
    elif affinity == min(affinities):
        final_affinity_list.append([names[i]]) # require a single optimal ligand
    else: 
        keep_indices.append(i)
affinities = affinities[keep_indices]
names = names[keep_indices]
# group ligands within 50 um of each other (well, not exactly, but kind of)
while len(affinities) > 0:
    affinity_list_updater = []
    keep_indices = []
    affinities_min = min(affinities) # remember this is the highest affinity
    affinities = [element-affinities_min for element in affinities] # shift
    for i in range(len(names)):
        if affinities[i] <= 50:
            affinity_list_updater.append(names[i])
        else: 
            keep_indices.append(i)
    affinities = affinities[keep_indices]
    names = names[keep_indices]
    final_affinity_list.append(affinity_list_updater)
final_affinity_list.append(bottom_tier)
##########################################################################
########################################################################
affinity_list = final_affinity_list
for file in os.listdir('../s3'):   
    if file[-5:] != '.mol2':
        continue
    # otherwise, we're working with the receptor
    receptor = mda.Universe(file)
    receptor_number = file[:-5].split('_')[1] # 2dstructurepredictionnumber-3dstructurepredicitonumber 
    all_receptor_atoms = receptor.select_atoms('element N').positions

    ligands = {}
    for ligand_mol2 in os.listdir('.'):
        if ligand_mol2[-5:] != '.mol2': # not actually a .mol2 file
            continue
        if '_' not in ligand_mol2: # not ligands
            continue
        if ligand_mol2[:-5].split('_')[1] != receptor_number: # not the right ligands
            continue
        ligands.update({ligand_mol2[:-5]:mda.Universe(ligand_mol2)})
    
    # note: should probably identify bad (too far from receptor) poses of the
    # optimal ligand and get rid of those before comparing to other ligands
    # this will be more efficient, and how significant that boost is, I 
    # don't know right now
    # could also move this static definition fo select_string outside the loop
    # but probably insignificant compared to the other operations going on here
    top_ligand_universe = ligands[affinity_list[0][0]]
    select_file = f'{BASE_PATH}/s4/common_select_string.txt'
    with open(select_file,'r') as f:
        for line in f:
            select_string = line
    top_ligand_common_atoms = top_ligand_universe.select_atoms(select_string)
    for top_ligand_pose_counter,top_ligand_pose in enumerate(top_ligand_universe.trajectory):
        distances = distance_array(top_ligand_pose,all_receptor_atoms) 
        # comparing an atomgroup being dynamically modified over the trajectory
        # to a static numpy array
        flag = 0
        for row in distances:
            for item in row: # note that we're double counting over the full matrix
                             # rather than a triangle
                if item < 4: #################PARAMETER################
                    flag += 1
        if int(flag/2) < 4: # this pose not worth exploring
                            # (we're requiring 2 atoms to be within
                            # a reasonable distance)
            continue 
        
        # get docking score
        #receptor_number = file[:-5] we already read this in, not sure why we'd
        # need to do it again. Also, the receptor number should be a different string
        # that file[:-5] (see above)
        ligand_name = affinity_list[0][0]
        sd_filename = glob(f'*{receptor_number}_{ligand_name}.sd')[0]
                 # might have a prefix to the numerical part of the
                 # receptor id, like 'mcsym_'
        top_ligand_score = get_inter_score(sd_filename,top_ligand_pose_counter)[0]

        # see how it compares to other ligands
        alternate_ligand_scores = {}
        for other_ligand_counter in range(1,len(affinity_list)):
            next_ligands = affinity_list[other_ligand_counter]
            for next_ligand in next_ligands:
                next_ligand_universe = ligands[next_ligand]
                # now we need to select the common atoms for the purpose
                # of assessing rmsd. And we should probably do per-atom RMSD
                # because having different numbers of atoms could affect that
                # going to read them in from file
                # relying on select_string getting equal indices for corresponding atoms
                next_ligand_common_atoms = next_ligand_universe.select_atoms(select_string)
                length = len(top_ligand_common_atoms)
                #assert False not in[top_ligand_common_atoms[i]==next_ligand_common_atoms[i] for i in range(length)]
                
                for next_ligand_pose_counter,next_ligand_pose in enumerate(next_ligand_universe.trajectory):
                    distances = distance_array(top_ligand_common_atoms,next_ligand_common_atoms)
                    # so we have all pairwise distances between all corresponding atoms,
                    # now we just want distances between corresponding atoms
                    distances = [distances[i][i] for i in range(len(distances))]
                    normalized_rmsd = np.sqrt(np.average(np.square(distances)))/length
                    if normalized_rmsd < 0.1: # can change this parameter, probably want to experiment a bit
                        # maybe now we read in docking scores from files and compare, then only save the ones
                        # that match our ranking system? but this should be done outside the ligand loop
                        # so maybe we just append ligand score to a list for now, then check at the end
                        # if length of this list is small compared to the number of ligands, that says 
                        # we migth have a sampling problem. If they're equal, then it's more likely that
                        # all ligands are being well-sampled. Then we check for the qualitative docking
                        # score agreement. But first things first, go pu to the ligand loop and get
                        # the docking score after it is verified to be a decent pose
                        #receptor_number = file[:-5]
                        # again, I don't think we need to re-declare receptor number
                        ligand_name = next_ligand
                        sd_filename = f'{receptor_number}_{ligand_name}.sd'
                        if ligand_name not in alternate_ligand_scores.keys():
                            alternate_ligand_scores.update({ligand_name:[get_inter_score(sd_filename,next_ligand_pose_counter)]})
                        else:
                            alternate_ligand_scores[ligand_name].append(get_inter_score(sd_filename,next_ligand_pose_counter))
                        # now we have a dictionary of lists of tuples, where the first element of each tuple
                        # is the docking score, and the second element is its index in the sd file (pose_counter)

                        # so we have an alternate ligand pose within 0.1A per-atom rmsd of the top ligand pose
                        # let's save it in a directory called 
                        # f'{receptor_number}_{next_ligand}_pose{next_ligand_pose_counter}'
                        # we won't bother to save the top ligand pose because it's easy to access that
                        # by looping through the original top ligand sd file
                        pose_write_file_name = f'{receptor_number}_{next_ligand}_pose{next_ligand_pose_counter}'
                        dir_name = f'{receptor_number}_top_ligand_pose{top_ligand_pose_counter}'
                        if dir_name not in os.listdir('.'):
                            os.mkdir(dir_name)
                        to_write = next_ligand_universe.select_atoms('all')
                        to_write.write(f'{dir_name}/{pose_write_file_name}.mol2')
                        # have to do mol2 because mda can't write sd
                        # so then will need to convert to sd to regenerate for ligmapping

                        # that'll be the end of this script
                        # in the future, we'll dock more precisely using ligmap                            
                        # note that we might have more than one pose saved for each
                        # alternative ligand. See the ligmapping script for comments
                        # on how I want to handle this 

        '''
        # I don't expect any of these to agree with the qualitative affinity trends, 
        # but it doesn't hurt to write that code now
        # then I'll probably want to go above this block and save the clusters
        # for input to the 1G ligand mapping round
        
        if check_affinity_trend(affinity_list,top_ligand_score,alternate_ligand_scores):
            # we're going to want to save the poses    
            print(f'pose counter: {pose_counter}')
            print(f'top ligand score: {top_ligand_score}')
            print(f'alternate ligands: {alternate_ligand_scores}')
        else:
            # don't want to do anything because the trend doesn't agree
            pass
        '''
