import sys

receptor = sys.argv[1]
run = sys.argv[2]
affinity_info = sys.argv[3] # string of format 'key':value,'key':value,etc
ligands = affinity_info.split(',')
ligands = [ligand.split(':')[0][1:-1] for ligand in ligands]
for ligand in ligands:
    prm_string = f"""RBT_PARAMETER_FILE_V1.00
TITLE basic_dock

RECEPTOR_FILE {receptor} 
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {ligand}_{run}
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
    with open(f'{ligand}_current.prm','w') as f:
        f.write(prm_string)
