"${SCHRODINGER}/utilities/prepwizard" $1.mae $1_prepwizard.mae -noepik -nometaltreat -nohtreat -antibody_cdr_scheme Kabat -propka_pH 7.4 -f S-OPLS -fix -rmsd 0.3 -watdist 5.0 -JOBNAME $1 -HOST localhost:4
# trying more alkaline (pHs of 8 and 9) than the default 7.4, so that backbone phosphates are deprotonated, but didn't work, even after removing the noepik flag
