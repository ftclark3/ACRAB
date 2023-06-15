#!/bin/bash
source $1 # conda shell
SCHRODINGER_PATH=$5 # schrodinger python interpreter
AFFINITY_INFO=${10}
 
AFFINITY_INFO_TRIM_1="${AFFINITY_INFO#{\'}" # delete shortest match of literal {' from begining
AFFINITY_INFO_TRIM_2="${AFFINITY_INFO_TRIM_1%\}}" # delete shortest match of literal } from end
# vim doing incorrect syntax highlighting here, the ' is escaped correctly
AFFINITY_INFO=$AFFINITY_INFO_TRIM_2
top_ligand="${AFFINITY_INFO_TRIM_1%%\':*}" # delete longest match of expansion ':* from end
top_ligand_array=("$top_ligand")
# going to pass AFFINITY_INFO to the python programs
# and make them parse it themselves
# at this point, its just a comma-separated list of key:value
# with keys likely in literal single quotes

for file in $( ls ../s4/*-*-*.sd ); do
  file_name="${file##*\/}"
  if echo "$file_name" | grep "$top_ligand" ; then continue; fi
  # shortest match of expression */ from beginning 
  # (so the file path excluding the name)

  # we're going to put each ligand in its own directory
  mkdir "${file_name%.sd}"
  cd "${file_name%.sd}"
  echo "100" > scores.txt # an initial threshold
  # the first accepted docking score must be lower than 100
  # pretty much anything should be lower, so this isn't very restrictive
  cp "../$file" ./"${file_name%.sd}"_run0.sd
  # just getting putting the result from previous run
  # here under a different name, to make it easier 
  # in the future
  # we grep $top_ligand because the ligand name could 
  # have a hyphen in it, but we don't want to include files
  # with the ligand name because they're redundant and in a
  # form that I'm not prepared to deal with (they contain 
  # multiple structures per file, and I want one structure
  # per file)
  # now we need to make a parameter file for each of our first ligands
  # remove the last - and everything after it to get the receptor name



  #receptor_name="${file_name%\-*}".mol2
  #${SCHRODINGER_PATH}/run ../mol2tomae.py ../../s3/$receptor_name
  receptor_name="${file_name%\-*}".mae 
  cp "../../s3/${receptor_name}" .
  cp "../../s3/${receptor_name%.mae}.mol2" . 
  echo "RBT_PARAMETER_FILE_V1.00
TITLE basic_dock

RECEPTOR_FILE ${receptor_name%.*}.mol2
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL ${file_name%.sd}_run0.sd
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
END_SECTION" > "${file_name%.sd}"_run0.prm
  
  flag=0
  run=0
  while [ $flag -eq 0 ]; do
    use_prev=0
    #if [ $flag -eq 1 ]; then exit; fi don't see why we need this again
    if test -f prev_used; then use_prev=1; fi
    # modify receptor
    "$SCHRODINGER_PATH"/run python3 ../rotate_loop.py ${file_name%.sd}_run${run}.sd "${receptor_name%.*}".mae $use_prev $top_ligand
    exit
    # now, we need to make a parameter file for the next round of docking
    # on the very first run, this will overwrite the previous receptor,
    # but those receptors are stored in s4, so it's fine
    echo "RBT_PARAMETER_FILE_V1.00
TITLE basic_dock

RECEPTOR_FILE ${receptor_name%.*}_intermediate.mol2
RECEPTOR_FLEX 3.0

##############################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
##############################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL ${file_name%.sd}_run${run}.sd
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
END_SECTION" > "${file_name%.sd}"_run${run}.prm
    rbcavity -r "${file_name%.sd}"_run${run}.prm -W >> foo
    #rbdock -r "${file_name%.sd}"_run${run}.prm -i ${file_name%.sd}_run${run}.sd -p dock.prm -n 5 -o "${file_name%.sd}"_intermediate >> foo
    rbdock -r "${file_name%.sd}"_run${run}.prm -i ${file_name%.sd}_run${run}.sd -p dock.prm -n 5 -o "${top_ligand}"_intermediate >> foo
    "$SCHRODINGER_PATH"/run ../check_condition.py "$AFFINITY_INFO" ${run} ${top_ligand_array} ${receptor_name%.*}_intermediate.mol2
    if test -f goodtogo
    then
      rm goodtogo
      cp intermediate.mol2 run$(( $run + 1 )).mol2
      cp ${receptor_name%.mol2}_intermediate.mae "${receptor_name%.mol2}".mae # we're recycling names, but whatever
    # this isn't quite right, so take some time to compare to the previous working version and see how the names
    # differ and how we can make this better
      #cp 1-MX_intermediate.sd 1-MX_run$(( $run + 1 )).sd
      #cp 3-MX_intermediate.sd 3-MX_run$(( $run + 1 )).sd
      #cp 7-MX_intermediate.sd 7-MX_run$(( $run + 1 )).sd
      #cp 13-DMX_intermediate.sd 13-DMX_run$(( $run + 1 )).sd
      #cp 17-DMX_intermediate.sd 17-DMX_run$(( $run + 1 )).sd
      #cp 37-DMX_intermediate.sd 37-DMX_run$(( $run + 1 )).sd
      #cp 13-DMU_intermediate.sd 13-DMU_run$(( $run + 1 )).sd
      #cp 17-DMU_intermediate.sd 17-DMU_run$(( $run + 1 )).sd
      #cp 37-DMU_intermediate.sd 37-DMU_run$(( $run + 1 )).sd
      #cp 137-TMX_intermediate.sd 137-TMX_run$(( $run + 1 )).sd
      #cp 137-TMU_intermediate.sd 137-TMU_run$(( $run + 1 )).sd
      cp "${file_name%.sd}"_intermediate.sd "${file_name%run*.sd}"run$(( $run + 1 )).sd 
      cp ${receptor_name%.mol2}_intermediate.mae "${receptor_name%.mol2}".mae # we're recycling names, but whatever
      run=$(( $run + 1 ))    
      continue
  else
    exit #temporary
    #continue
  fi

  done
done





exit



flag=0
run=0
#declare -a ligand_string=("1-MX" "3-MX" "7-MX" "13-DMX" "17-DMX" "37-DMX" "137-TMX" "13-DMU" "17-DMU" "37-DMU" "137-TMU")
#actual_ligand_string="1-MX 3-MX 7-MX 13-DMX 17-DMX 37-DMX 137-TMX 13-DMU 17-DMU 37-DMU 137-TMU"
#top_ligand="17-DMU"
#declare -a top_ligand_array=("17-DMU")

while [ $flag -eq 0 ]
do
  use_prev=0
  if [ $flag -eq 1 ]; then exit; fi
  if test -f prev_used; then use_prev=1; fi
  #echo "before rotate_loop call"
  "$SCHRODINGER_PATH"/run python3 rotate_loop.py run${run}.sd run${run}.mae $use_prev $top_ligand $AFFINITY_INFO
  #ligand is passed to the program but appears to be unnecessary, so
  # I will not worry about changing it now that we have multiple ligands
  # writes out 'intermediate.mol2' and intermediate.mae
  #echo "after rotate_loop call"
  #declare -a ligand_string=("1-MX" "3-MX" "7-MX" "13-DMX" "17-DMX" "37-DMX" "137-TMX" "13-DMU" "17-DMU" "37-DMU" "137-TMU")
  #actual_ligand_string="1-MX 3-MX 7-MX 13-DMX 17-DMX 37-DMX 137-TMX 13-DMU 17-DMU 37-DMU 137-TMU"
  #top_ligand="17-DMU"
  #declare -a top_ligand_array=("17-DMU")
  #set -- $actual_ligand_string
  # not actually a string but you get the idea
  #python make_current_prm.py intermediate.mol2 run${run}.sd ${actual_ligand_string}
  python make_current_prm.py intermediate.mol2 run${run}.sd $AFFINITY_INFO

  for i in "${ligand_string[@]}"
  do  
      if [ $i == ${top_ligand} ]
      then
          rbcavity -r "${i}_current.prm" -W >> foo 
          rbdock -r "${i}_current.prm" -i "${i}_run${run}.sd" -p dock.prm -n 5 -o "${i}_intermediate" >> foo
          continue
      fi
  done
  #"/opt/schrodinger2022-4/run" check_condition.py ${top_ligand_array} ${run}
  
  "$SCHRODINGER_PATH"/run check_condition.py ${actual_ligand_string} ${run} ${top_ligand_array}
  #if test -f test_file; then exit; fi #temporary 
  if test -f goodtogo 
    then
    rm goodtogo
    cp intermediate.mol2 run$(( $run + 1 )).mol2
    #cp 1-MX_intermediate.sd 1-MX_run$(( $run + 1 )).sd
    #cp 3-MX_intermediate.sd 3-MX_run$(( $run + 1 )).sd
    #cp 7-MX_intermediate.sd 7-MX_run$(( $run + 1 )).sd
    #cp 13-DMX_intermediate.sd 13-DMX_run$(( $run + 1 )).sd
    #cp 17-DMX_intermediate.sd 17-DMX_run$(( $run + 1 )).sd
    #cp 37-DMX_intermediate.sd 37-DMX_run$(( $run + 1 )).sd
    #cp 13-DMU_intermediate.sd 13-DMU_run$(( $run + 1 )).sd
    #cp 17-DMU_intermediate.sd 17-DMU_run$(( $run + 1 )).sd
    #cp 37-DMU_intermediate.sd 37-DMU_run$(( $run + 1 )).sd
    #cp 137-TMX_intermediate.sd 137-TMX_run$(( $run + 1 )).sd
    #cp 137-TMU_intermediate.sd 137-TMU_run$(( $run + 1 )).sd
    cp ${top_ligand}_intermediate.sd ${top_ligand}_run$(( $run + 1 )).sd

    cp intermediate.mae run$(( $run + 1 )).mae
    run=$(( $run + 1 ))
    #if test -f threshold; then flag=1; rm threshold; fi    
    continue  
  else
    #exit #temporary
    continue
  fi
  # now calculate deviation of intermediate.mol2-intermediate.sd
  # now compare to previous deviation and decide whether to accept or reject
  # we'll have python create the file goodtogo to signal to bash that
  # can continue with the loop
  # if accept:
    # check if file exists using if test -f goodtogo; then
    # rm goodtogo
    # cp intermediate.mol2 run$(( $run + 1 )).mol2
    # cp intermediate.sd run$(( $run + 1 )).sd
    # run=$(( $run + 1 ))
    # continue  
    #fi
  # else # (reject):
    # don't increment run
    # continue
done
