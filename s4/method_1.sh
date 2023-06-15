#!/bin/bash
source $1 # the conda shell
BASE_NAME=$3
BASE_PATH="$2/${BASE_NAME}"
AFFINITY_INFO=${10} # this will be a string serialization of a python dictionary
SCHRODINGER_PATH=$5
#echo "AFFINITY_IFNO"
#echo $AFFINITY_INFO

# function to facilitate multiprocessing below
rbcavity_rbdock() {
  # rbcavity argument -> $1
  # rbdock argument 1 (-i) -> $2
  # rbdock argument 2(-o) -> $3
  # rbdock argument 3 (-r) equal to $1
  rbcavity -W -r $1
  rbdock -i $2 -o $3 -r $1 -p dock.prm -n 5
}



# assuming we're in BASE_PATH/s4, we want to go through each receptor in s3 
for file in $( ls ../s3 | grep ".mol2" ); do
  python onetime_make_prm.py ${BASE_PATH} "surface" "_" "_" "../s3/${file}"
  # don't need the .mol2 extension here, it's included in the ls return value
  # need to pass in "_" as dummy variables to onetime_make_prm.py script
  # so that the onetime_make_prm.py script remains compatible with some 
  # earlier workflows I was building
done

# still in BASE_PATH/s4, so we want to go through each prm file and do
# a surface docking. We'll find ligands in the s4 section of acrab_path
#python self_dock_surface.py ${BASE_PATH} $AFFINITY_INFO
# this time, the loop is inside the python script
# or do we even need a python script?
AFFINITY_INFO_TRIM_1="${AFFINITY_INFO#{\'}" # delete shortest match of literal {' from begining
# vim doing incorrect syntax highlighting here, the ' is escaped correctly
top_ligand="${AFFINITY_INFO_TRIM_1%%\':*}" # delete longest match of expansion ':* from end
for file in $( ls | grep .prm ); do
  #rbcavity -W -r $file
  #rbdock -i ${top_ligand}.sd -o "${file%.prm}"_"${top_ligand}" -r $file -p dock.prm -n 100 
  # remember that the program automatically appends the .sd extension to the -o switch
  ((i=i%3)); ((i++==0)) && wait # apparently, this command initializes i
  rbcavity_rbdock $file ${top_ligand}.sd "${file%.prm}"_"${top_ligand}" & 
done

conda activate mda
for file1 in $( ls ../s3/*-*.mol2 ); do
    # unlike in the above make_prm case, the whole relative path is returned by ls
    file1_name=${file1%.mol2} 
    for file2 in $( ls "${file1_name##*\/}"_*.sd ); do 
      # delete longer match of anything ending in / 
      # from beginning (so we're just working with the file name)
      # first, we need to convert to mol2 format
      "$SCHRODINGER_PATH"/run mol2sd_converter.py $file2 "${file2%.sd}".mol2
      python multi_lig_sd_stage_5_prep.py $file1 "${file2%.sd}".mol2
    done
done
conda deactivate






#conda activate mda
#python onetime_cluster_finder.py ${BASE_PATH} "_" $AFFINITY_INFO
# dummy "_" instead of GENERATION parameter from an older workflow
#conda deactivate








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
