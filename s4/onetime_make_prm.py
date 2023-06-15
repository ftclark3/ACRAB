import os
import subprocess
import sys

BASE_PATH = sys.argv[1]
grid_gen_type = sys.argv[2]

if grid_gen_type == 'surface':
    site_mapper = 'RbtSphereSiteMapper'
    filestring_addition = '' # for compatibility for using this script for onetime ligmapping
    REF_MOL_string = ''
    receptor = sys.argv[5] # expecting .mol2 file
elif grid_gen_type == 'ligmap':
    site_mapper = 'RbtLigandSiteMapper'
    filestring_addition = '_ligmap'
    REF_MOL = f'{sys.argv[3]}'
    REF_MOL_string = f'\n    REF_MOL {REF_MOL}'
    #REF_MOL_string = f'REF_MOL {REF_MOL}'
    GENERATION = sys.argv[4]
    receptor = sys.argv[5]
else:
    print('unknown mapping method')
    raise SystemExit()

if grid_gen_type == 'ligmap': # code too different from surface to have it within the same block
    prm_string = f"""RBT_PARAMETER_FILE_V1.00
TITLE {receptor[:-5]}
    
RECEPTOR_FILE {receptor}
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: LIGAND OR SPHERE MAPPING METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER {site_mapper}{REF_MOL_string}
    RADIUS 6.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    MIN_VOLUME 100
    GRIDSTEP 0.5
END_SECTION

############################
## CAVITY RESTRAINT PENALTY
############################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION"""
    with open(f'{BASE_PATH}/tertiary_structures/{GENERATION}/{REF_MOL.strip()[:-3]}.prm','w') as f:
        f.write(prm_string)

else:
    #receptors = [file for file in os.listdir('tertiary_structures/1G') if '_' not in file and file[-5:] == '.mol2']
    #for receptor in receptors:
    prm_string = f"""RBT_PARAMETER_FILE_V1.00
TITLE {receptor[:-5]}

RECEPTOR_FILE {receptor}
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: LIGAND OR SPHERE MAPPING METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER {site_mapper}{REF_MOL_string}
    RADIUS 200.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 10
    VOL_INCR 0.0
    MIN_VOLUME 100
    GRIDSTEP 0.5
END_SECTION

############################
## CAVITY RESTRAINT PENALTY
############################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION"""
    with open(f'{BASE_PATH}/s4/{receptor[:-5].split("/")[-1]}.prm','w') as f:
        f.write(prm_string)
