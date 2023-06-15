import sys
import os
import numpy as np
import subprocess
import json
from schrodinger import structure

def get_score(file,pose_index):
    # file: relative path, assuming we're in {BASE_PATH}/tertiary_tructures/{GENERATION}
    # pose_index: 0-indexed location of the pose that we want to score in the sd file
    with open(file,'r') as f:
        prev_indices = []
        scores = []
        for counter,line in enumerate(f):
            if '>  <SCORE>' in line:
                prev_indices.append(counter)
    with open(file,'r') as f:
        for counter,line in enumerate(f):
            if counter-1 in prev_indices:
                scores.append(float(line.strip()))
    score = scores[pose_index]
    return (score,pose_index)


AFFINITY_INFO = sys.argv[1]
#print(AFFINITY_INFO)
run = sys.argv[2]
top_ligand = sys.argv[3]
receptor = sys.argv[4]
if top_ligand[0] == "'":
    top_ligand = top_ligand[1:]
if top_ligand[-1] == "'":
    top_ligand = top_ligand[:-1]

#build affinity list
affinity_info = AFFINITY_INFO.split(',') 
affinity_info = [ligand.split(':') for ligand in affinity_info]
for i in range(len(affinity_info)):
    affinity_info[i] = [affinity_info[i][0].strip(),affinity_info[i][1].strip()]
    if affinity_info[i][0][0] == "'":
        affinity_info[i][0] = affinity_info[i][0][1:]
    if affinity_info[i][0][-1] == "'":
        affinity_info[i][0] = affinity_info[i][0][:-1]
# trimming the literal ' in case they snuck in from the bash script


final_affinity_list = [[top_ligand]]
affinities = [int(name_affinity_tuple[1]) for name_affinity_tuple in affinity_info]
# it's a list, not a tuple, but doesn't matter
names =      [name_affinity_tuple[0] for name_affinity_tuple in affinity_info]
bottom_tier = []
keep_indices = []
for i in range(len(names)):
    if affinities[i] >= 1000:
        bottom_tier.append(names[i]) # we'll tack this on at the end
    elif affinities[i] == min(affinities):
        continue # we've already got top_ligand in there on its own
    else:
        keep_indices.append(i)
affinities = [affinities[i] for i in range(len(affinities)) if i in keep_indices]
names = [names[i] for i in range(len(names)) if i in keep_indices]
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
    #affinities = affinities[keep_indices]
    #names = names[keep_indices]
    affinities = [affinities[i] for i in range(len(affinities)) if i in keep_indices]
    names = [names[i] for i in range(len(names)) if i in keep_indices]
    final_affinity_list.append(affinity_list_updater)
final_affinity_list.append(bottom_tier)
affinity_list = final_affinity_list


# flatten list for ease of use below
ligands = []
for sublist in affinity_list:
    for item in sublist:
        ligands.append(item)
all_scores = {}
for ligand in ligands:
    #print("for")
    if ligand == top_ligand:
        scores = []
        for index in range(5):    
            try:
                scores.append(get_score(f'{ligand}_intermediate.sd',index)[0])
            except IndexError:
                if index > 0:
                #if False: # actually, I told it to do 5 poses, so it always should, but not the first time round
                    pass # possible that it returned fewer than 5 poses
                else:
                    print('else')
                    print(f'{ligand}_run{run}.sd')
                    raise
        score = min(scores)
        location = scores.index(min(scores))
        all_scores.update({ligand:score})

        for counter, conformer in enumerate(structure.StructureReader(f'{ligand}_intermediate.sd')):
            #print("for")
            if counter == location:
                #print("IF")
                conformer.write(f'temp_{ligand}.sd')
                conformer.write(f'temp_{ligand}.mol2')
                break
        subprocess.run(['cp',f'temp_{ligand}.sd',f'{ligand}_intermediate.sd'])
        subprocess.run(['cp',f'temp_{ligand}.mol2',f'{ligand}_intermediate.mol2'])

    else: # ligand != top_ligand, but we want to get that on the next loop
        pass

for ligand in ligands:
    if ligand == top_ligand:
        continue
    else:
        #ligand_st = structure.StructureReader.read(f'{ligand}_intermediate.sd')
        #ligand_st.write(f'{ligand}_intermediate.mol2')
        # align to top ligand and score
        # align function writes sd file. input: {ligand}_run0.sd, output: write file {ligand}_current.sd
        #subprocess.run(['bash','check_condition_helper.sh',ligand,f'{top_ligand}_intermediate.mol2'])
        ligand_st = structure.StructureReader.read(f'../{ligand}.sd')
        ligand_st.write(f'{ligand}.mol2')
        subprocess.run(['bash','../check_condition_helper.sh',ligand,f'{top_ligand}_intermediate.mol2',receptor])
        #print(ligand)
        #subprocess.run(['python','make_score_prm.py','intermediate.mol2',f'{ligand}_current.sd'])
        #subprocess.run(['rbcavity','-r',f'{ligand}_score.prm','-W'])
        #subprocess.run(['rbdock','-r',f'{ligand}_score.prm','-i',
        #                f'{ligand}_current.sd','-p','score.prm','-o',f'{ligand}_current_score'])
        all_scores.update({ligand:get_score(f'{ligand}_current_score.sd',0)[0]})

condition1 = all_scores[affinity_list[0][0]] == min(all_scores.values())
with open('scores.txt','r') as f:
    for line in f:
        score = float(line.strip())
condition2 = all_scores[affinity_list[0][0]] < score
if condition1 and condition2:
    print(all_scores)
    with open('goodtogo','w') as f:
        f.write('')
    with open('scores.txt','a') as f:
        f.write(f'{all_scores[affinity_list[0][0]]}\n')
elif condition1:
    #param = 20*np.random.rand()
    #param = 0.1/abs(score)
    temperature = 0.001
    # as the score gets lower, it should get harder to jump out of well
    boltzmann = np.exp(-1/temperature)
    probabilities = [boltzmann, 1-boltzmann]
    switch = np.random.choice([True,False],1,p=probabilities)[0]
    if switch:
        #print(all_scores)
        with open('goodtogo','w') as f:
            f.write('')


'''
with open('new_vectors.txt','r') as f:
    for line in f:
        new_value = float(line.strip())
with open('differences.txt','r') as f:
    for line in f:
        old_value = float(line.strip())
    # after this loop, old_value will be the latest entry in the file
if new_value < old_value:
    with open('goodtogo','w') as f:
        f.write('')
    with open('differences.txt','a') as f:
        f.write(f'{new_value}\n')
    with open('new_score.txt','r') as f:
        for line in f:
            score = line.strip()
    with open('scores.txt','a') as f:
        f.write(f'{score}\n')
    with open('prev_used','w') as f:
        f.write('')
else:
    if os.path.exists('prev_used'):
        os.remove('prev_used') # loop modification was a failure, so don't
                               # want to keep the modifications used to do that and repeat
    temperature = old_value/100
    # as the score gets lower, it should get harder to jump out of well
    boltzmann = np.exp(-abs(new_value-old_value)/temperature)
    probabilities = [boltzmann, 1-boltzmann]
    switch = np.random.choice([True,False],1,p=probabilities)[0]
    #switch = True
    if switch: 
        with open('goodtogo','w') as f:
            f.write('')
        with open('differences.txt','a') as f:
            f.write(f'{new_value}\n')
        with open('new_score.txt','r') as f:
            for line in f:
                score = line.strip()
        with open('scores.txt','a') as f:
            f.write(f'{score}\n')



if new_value < threshold:
    with open('threshold','w') as f:
        f.write('')
'''
