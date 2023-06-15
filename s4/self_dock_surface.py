import os
import subprocess
import sys

BASE_PATH = sys.argv[1]
AFFINITY_INFO = sys.argv[2]
my_env = os.environ.copy()
my_env['RBT_HOME'] = f'{BASE_PATH}/s4/'
os.chdir(f'{BASE_PATH}/s4') # should already be here, but we'll do it
                            # again, just in case

top_ligand_name = AFFINITY_INFO[1:-1].split(':')[0][1:-1]
# cut out the curly braces from the string-serialized dictionary,
# then get the name (dictionary key), then cut off the literal quote characters

for file in os.listdir('.'): 
    if file[-4:] != '.prm':
        continue

    subprocess.run(['rbcavity', '-W', '-r', file],env=my_env) 

    subprocess.run(['rbdock', '-i', f'{top_ligand_name}.sd', '-o', f'{file[:-4]}_{ligand[:-3]}', '-r', file, '-p', 'dock.prm','-n','1'],env=my_env)
    
    '''
    for ligand in os.listdir('.'): 
        if ligand[-3:] != '.sd' or top_ligand_name not in ligand:
            continue
        if file[:-4] in ligand: # it's an output pose, not input ligand
            continue
        subprocess.run(['rbdock', '-i', ligand, '-o', f'{file[:-4]}_{ligand[:-3]}', '-r', file, '-p', 'dock.prm','-n','1'],env=my_env)
        #subprocess.run(['/home/clarkf/ADFR/bin/obabel',f'{file[:-4]}_{ligand[:-3]}.sd','-O',f'{file[:-4]}_{ligand[:-3]}.mol2']) I'm not sure this is needed in this new workflow
    '''
