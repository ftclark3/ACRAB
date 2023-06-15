#!/bin/bash
#echo "check_condition_helper called"
source /home/clarkf/anaconda3/etc/profile.d/conda.sh
unset PYTHONHOME
conda activate mda
python ../align.py "$1.mol2" $2 # writes mol2 file
conda deactivate
"/opt/schrodinger2022-4/run" ../mol2tosd.py "$1_current.mol2" "$1_current.sd"
#python make_score_prm.py intermediate.mol2 "$1_current.sd"
#receptor=$( ls * )
echo "RBT_PARAMETER_FILE_V1.00
TITLE basic_dock

RECEPTOR_FILE $3 
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL $1_current.sd
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
END_SECTION" > $1_score.prm

rbcavity -r $1_score.prm -W >> foo
rbdock -r $1_score.prm -i $1_current.sd -p score.prm -o $1_current_score >> foo
