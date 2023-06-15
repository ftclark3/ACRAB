#!/bin/bash
source $1 # the conda shell
BASE_NAME=$3
BASE_PATH="$2/${BASE_NAME}"
AFFINITY_INFO=${10} # this will be a string serialization of a python dictionary

# assuming we're in BASE_PATH/s4, we want to go through each receptor in s3 
for file in $( ls ../s3 | grep ".mol2" ); do
  python onetime_make_prm.py ${BASE_PATH} "surface" "_" "_" "../s3/${file}.mol2"
  # need to pass in "_" as dummy variables to onetime_make_prm.py script
  # so that the onetime_make_prm.py script remains compatible with some 
  # earlier workflows I was building
done

# still in BASE_PATH/s4, so we want to go through each prm file and do
# a surface docking. We'll find ligands in the s4 section of acrab_path
python self_dock_surface.py ${BASE_PATH} 
# this time, the loop is inside the python script

conda activate mda
python onetime_cluster_finder.py ${BASE_PATH} "_" $AFFINITY_INFO
# dummy "_" instead of GENERATION parameter from an older workflow
conda deactivate








# have to move the below command to a call within self_dock_ligmap.py
# and include the name of the reference file as sys.argv[3]
#python onetime_make_prm.py /home/clarkf/aptamer/pipeline/leighton2 ligmap
    
# need to move 1G receptors to 2G because self_dock_ligmap.py
# expects them to already be there 
#mkdir /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/2G
#cp /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/1G/*.mol2 /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/2G/
#rm /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/2G/*_*.mol2

# these *_ligmap.prm files are actually useless because ligand mapping requires the ligand pose
#cp /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/1G/*_ligmap.prm /home/clarkf/aptamer/pipeline/leighton2/tertiary_structures/2G/

#conda activate mda
#python self_dock_ligmap.py /home/clarkf/aptamer/pipeline/leighton2 2G
#conda deactivate
