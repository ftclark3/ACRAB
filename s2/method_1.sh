#!/bin/bash
#source /home/clarkf/anaconda3/etc/profile.d/conda.sh
source $1
#BASE_NAME='leighton3'
BASE_NAME=$3
#BASE_PATH="/home/clarkf/aptamer/pipeline/${BASE_NAME}"
BASE_PATH="$2/${BASE_NAME}"

# $4 for MDA path
# $5 for schrodinger path
# $6 for SPOT-RNA path
# $7 for MCFold path
# mfold binary expected to already be in path

#start with SPOT-RNA
conda activate SPOT-RNA
mkdir ${BASE_PATH}/s2/SPOT-RNA
mv ${BASE_PATH}/s2/bpseq2dotb.py ${BASE_PATH}/s2/SPOT-RNA/
cd $6
python SPOT-RNA.py --inputs ${BASE_PATH}/${BASE_NAME}.fasta --outputs ${BASE_PATH}/s2/SPOT-RNA
conda deactivate
# current working directory at this point int he program: 
# ~/clarkf/SPOT-RNA
cd ${BASE_PATH}/s2/SPOT-RNA
rm *.ct
rm *.prob
python bpseq2dotb.py ${BASE_NAME}.bpseq
# just keeping .dbn files
cd ..
for file in $( ls ./SPOT-RNA/*.dbn ); do 
  sed -n '1p' $file >> dot_brackets.txt
  echo '' >> dot_brackets.txt
done

# now MCFold 
#export MCSYM_DB="/usr/local/share/mcsymdb-4.2.1.bin.gz"
TOP=15
#cd /usr/local/share
export QUERY_STRING="pass=lucy&sequence=$(cat ${BASE_PATH}/${BASE_NAME}.fasta | grep -v '>')&pseudoknotted=true&top=${TOP}"
# we'll let it predict pseudoknots, because it's reasonable about not always predicting them (based on a single sample though)
cd $7
mkdir ${BASE_PATH}/s2/MCFold
#/usr/local/share/mcfold.static.exe > ${BASE_PATH}/secondary_structures/MCFold/${BASE_NAME}_results.txt 
./mcfold.static.exe > ${BASE_PATH}/s2/MCFold/${BASE_NAME}_results.txt 
cd ${BASE_PATH}
sed -n '10,$p' ${BASE_PATH}/s2/MCFold/${BASE_NAME}_results.txt > tmp 
# above command is removing the first few lines of the file which are bad for some reason
cp tmp ${BASE_PATH}/s2/MCFold/${BASE_NAME}_results.txt
rm tmp
echo "import os
import sys
number = int(sys.argv[1]) # number of predictions
os.chdir('./s2/MCFold')

for file in os.listdir('.'):
    
    if file[-4:] == '.txt':
        lines = []
        with open(file,'r') as f:
            for counter3,line in enumerate(f):
                lines.append(line)
        
        
        end = False
        for i in range(1,number+10): # the 10 here is kind of arbitrary
                                     # just want to make sure we have enough iterations
                                     # to hit the break statement
            if end:
                if lines[-i][0] in '.(': #haven't hit the first one yet
                    continue
                else:
                    start = -i+1
                    break
            else:
                if lines[-i][0] in '.(':
                    end = True
                else:
                    continue
        dotbs = lines[start:-4] #5th to last line always the last structure
        for counter1,dotb in enumerate(dotbs):
            for counter2,character in enumerate(dotb):
                if character == ' ':
                    dotb = dotb[:counter2]
                    with open(f'{file[:-4]}_{counter1+1}.dotb','w') as f:
                        f.write(dotb)" > dot_bracket_extract_mcfold.py

python dot_bracket_extract_mcfold.py ${TOP} #only needs to be run once
cd ${BASE_PATH}/s2
for file in $( ls ./MCFold/*.dotb ); do
  if [[ "$file" == *".dotb" ]]; then
    # .dotb is just another suffix im using, equivalent to .dbn
    sed -n '1p' $file >> dot_brackets.txt
    echo '' >> dot_brackets.txt
  fi
done

## now mfold
### UNFORTUNATELY, MFOLD DOESN'T USE THE SAME BPSEQ FORMAT
### AS MCFOLD, SO MY CONVERSION THING WON'T WORK
### GOING TO SKIP THIS FOR NOW, SO I CAN GET THE WHOLE
### PIPELINE IMPLEMENTED 
#mkdir ${BASE_PATH}/s2/mfold
#cd ${BASE_PATH}/s2/mfold
#mfold SEQ=${BASE_PATH}/${BASE_NAME}.fasta NA=DNA MAX=5
#mkdir tmp
##mv *.ct tmp
##rm * 
#mv tmp/* .
#rmdir tmp 
