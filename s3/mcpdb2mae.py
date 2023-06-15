import sys
from schrodinger import structure
from schrodinger.structutils import build
from schrodinger.forcefield import minimizer
#from schrodinger.forcefield.common import AtomTypingFailure
#from schrodinger.structutils import assignbondorders
import subprocess
import time
import os
import shutil

def mcpdb2mae(in_path_out_path_tuple):
    in_path = in_path_out_path_tuple[0]
    out_path = in_path_out_path_tuple[1]
    lines = []
    codes = []
    with open(in_path,'r') as f:
        for line in f:
            if "ATOM" not in line:
                continue
            if "O2*" in line or "NA" in line: # NA added for minimization
                continue # ignoring RNA stuff
            lines.append(line)
            # starting with backbone
            if "O1P" in line:
                codes.append("OP1")            
            elif "O2P" in line:
                codes.append("OP2")
            elif "O3*" in line:
                codes.append("O3'")
            elif "O5*" in line:
                codes.append("O5'")
            elif "C1*" in line:
                codes.append("C1'")
            elif "C2*" in line:
                codes.append("C2'")
            elif "C3*" in line:
                codes.append("C3'")
            elif "C4*" in line:
                codes.append("C4'")
            elif "O4*" in line:
                codes.append("O4'")
            elif "C5*" in line:
                codes.append("C5'")
            elif "H1*" in line:
                codes.append("H1'")
            elif "H2*" in line:
                codes.append("H2'")
            elif "H3*" in line:
                codes.append("H3'")
            elif "H4*" in line:
                codes.append("H4'")
            elif "H5*" in line:
                codes.append("H5'")
            
            # side chain adjustments
            # G
            elif "1H2" in line:
                codes.append("H21")
            elif "2H2" in line:
                codes.append("H22")
          
	    #C
            elif "1H4" in line:
                codes.append("H41")
            elif "2H4" in line:
                codes.append("H42")
          
	   #A
            elif "1H6" in line:
                codes.append("H61")
            elif "2H6" in line:
                codes.append("H62")
   
            # 5' terminal phosphate is missing an oxygen
            # maybe we just make it hydroxyl capped, maybe
            # we add it in, but we'll take care of that later

  
  
            # lastly, we need to just read the code from the
            # string in the old pdb file, if we have nothing
            # else to do
            else:
                potential_codes = [part for part in line.split(' ') if part != '']
                code = potential_codes[2]
                codes.append(code)
          
    for counter in range(len(lines)):
        line = lines[counter]
        code = codes[counter]
        #print(code)
        if len(code) < 4:
            name = f" {code.ljust(3,' ')}"
        elif len(code) == 4:
            name = code
        else:
            return f'weird atom name !!!: {code}'
            #print('weird atom name!!!')
            #print(code)
            #exit()
        line_info = [part for part in line.split(' ') if part != '']
        new_line = f"{line[:11]} {name} D{line_info[3].ljust(2,' ')}{line[20:]}"
        #print(new_line)
        lines[counter] = new_line
   
    all_lines = []
    #with open(sys.argv[1],'r') as f:
    with open(in_path, 'r') as f:
        for line in f: 
            if "O2*" not in line and "NA" not in line:
                all_lines.append(line)
   
    counter1 = 0
    for counter2,line in enumerate(all_lines):
        if "ATOM" in line:
            all_lines[counter2] = lines[counter1]
            counter1 += 1
   
    all_lines_string = f'{all_lines[0]}'
    for i in range(1,len(all_lines)):
        all_lines_string += f'{all_lines[i]}'
   
    # all_lines_string now represents a pdb format file
    # with accurate atom names
   
    # get it into schrodinger format
    with open(f'{in_path}_temp.pdb','w') as f:
        f.write(all_lines_string)
    try:
        st = structure.StructureReader.read(f'{in_path}_temp.pdb')
    except StopIteration:
        #print(sys.argv[2])
        #print(all_lines_string)
        #raise
        return f'no structure in file: {in_path}_temp.pdb'
    os.remove(f'{in_path}_temp.pdb')
    
    #for residue in st.residue:
    #    print(residue.pdbres)
    #    print(residue.resname)
    #exit()
    # get lowest resnum
    lowest_resnum = 9999
    for residue in st.residue:
        if residue.resnum < lowest_resnum:
            lowest_resnum = residue.resnum
   
   
   
    #########################################################3
    ### MUTATING DU TO DT
    for residue in st.residue:
        if residue.pdbres == "DU  ": 
            residue.pdbres = "DT  "
            # note that this is inconsistent with pdb syntax expectations,
            # but maestro supposedly handles this okay, so I'm deferring to
            # maestro's formatting expectations
            for index in residue.getAtomIndices():
                if st.atom[index].pdbname == " H5 ":
                    #print(f'......{st.atom[index].pdbname}.....')
                    to_change = index
            st.atom[to_change].element = 'C'
            st.atom[to_change].pdbname = " C7 "
    ############################################################
   
    # adding hydrogens probably helped schrodinger figure out
    # which atoms should be bonded, but now we have to manually
    # correct a few backbone bonds, so it makes sense to delete 
    # hydrogens, then add them back after these corrections
    build.delete_hydrogens(st)
   
    #############################################################
    #### FIXING ALL THE BONDS IN OUR STRUCTURE
   
   
    # make sure every OP1-P bond is double
    for i in range(1,len(st.atom)+1):
        to_delete = []
        P = None
        if st.atom[i].pdbname.strip() == "OP1":
            for atom in st.atom[i].bonded_atoms:
                if atom.pdbname.strip() == "P":
                    P = atom
                else:
                    to_delete.append(atom)
            assert P, f'st.atom[i].bonded_atoms: {st.atom[i].bonded_atoms}'
            #assert len(to_delete) == 0, to_delete
            for thing in to_delete:
                st.deleteBond(st.atom[i],thing)
            st.addBond(st.atom[i],P,2)
   
    # fix backbone bonds
    for i in range(1,len(st.atom)+1):
        atom = st.atom[i]
        current_res = atom.resnum
   
        if atom.pdbname.strip() == "O3'":
            C3p_flag = False
            to_delete_indices = [] # delete nothing by default
            for other_atom in atom.bonded_atoms:
                if other_atom.pdbname.strip() == "C3'":
                    C3p_flag = True
                if other_atom.pdbname.strip() != "P" and other_atom.pdbname.strip() != "C3'":
                    to_delete_indices.append(other_atom.index)
            if len(to_delete_indices) != 0:
                for to_delete_index in to_delete_indices:
                    st.deleteBond(atom,st.atom[to_delete_index])
                #atom.deleteBond(st.atom[to_delete_index])
            if C3p_flag == False:
                to_add_list = [C3p for C3p in st.atom if C3p.pdbname.strip() == "C3'" and C3p.resnum == current_res]
                # should only be 1 atom in that list
                to_add = to_add_list[0]
                st.addBond(atom,st.atom[to_add],1)
 
        elif atom.pdbname.strip() == "O5'":
            C5p_flag = False
            to_delete_indices = []
            for other_atom in atom.bonded_atoms:
                if other_atom.pdbname.strip() == "C5'":
                    C5p_flag = True
                if other_atom.pdbname.strip() != "P" and other_atom.pdbname.strip() != "C5'":
                    to_delete_indices.append(other_atom.index)
            if len(to_delete_indices) != 0:
                for to_delete_index in to_delete_indices:
                    st.deleteBond(atom,st.atom[to_delete_index])
                    #atom.deleteBond(st.atom[to_delete_index])
            if C5p_flag == False:
                to_add_list = [C5p for C5p in st.atom if C5p.pdbname.strip() == "C5'" and C5p.resnum == current_res]
                #print(to_add_list)
                # should only be 1 atom in that list
                to_add = to_add_list[0]
    
                #atom.addBond(to_add,1)
                st.addBond(atom,st.atom[to_add],1)
           
        elif atom.pdbname.strip() == "P": 
            # make sure we aren't bonded to any carbons
            to_delete_indices = []
            for other_atom in atom.bonded_atoms:
                if "C" in other_atom.pdbname:
                    #print("C found")
                    to_delete_indices.append(other_atom.index)
            for to_delete_index in to_delete_indices:
                st.deleteBond(atom,st.atom[to_delete_index])
                #print('bond deleted')
                #atom.deleteBond(st.atom[to_delete_index])
    
            # handle O5', assuming adding a bond to an atom that already exists
            # won't raise an error
            to_add_list = [O5p for O5p in st.atom if O5p.pdbname.strip() == "O5'" and O5p.resnum == current_res]
            #print(to_add_list)        
            to_add = to_add_list[0]
            assert len(to_add_list) > 0, f"O5' not found! index {i}"
            #atom.addBond(to_add,1)
            st.addBond(atom,to_add,1)
   
            # now doing the same for O3'
            to_add_list = [O3p for O3p in st.atom if O3p.pdbname.strip() == "O3'" and O3p.resnum == current_res-1]
            if current_res != lowest_resnum: #if we're dealing with the first residue, won't be a 3' to connect to
                to_add = to_add_list[0]
                #atom.addBond(to_add,1)
                st.addBond(atom,to_add,1)
          
    # somehow, we still have bonds between carbon and phosphorus atoms
    # need to get rid of those
    # also, a second objective: delete bonds between P and O3' occurring
    # within the same residue
    for i in range(1,len(st.atom)+1):
        to_delete = []
        if "P" == st.atom[i].pdbname.strip():
            for atom in st.atom[i].bonded_atoms:
                #if "C" in atom.pdbname:
                if atom.pdbname.strip() not in ["OP1","OP2","OP3","O5'","O3'"]:
                    to_delete.append(atom)
                if atom.pdbname.strip() == "O3'" and atom.resnum == st.atom[i].resnum:
                    to_delete.append(atom)
                    # this is the second objective of the comment above
   
            for atom in to_delete:
                st.deleteBond(st.atom[i],atom)
 
    # make sure OP2 is only bonded to phosphorus
    for i in range(1,len(st.atom)+1):
        found_P = False
        to_delete = []
        if st.atom[i].pdbname.strip() == "OP2":
            for atom in st.atom[i].bonded_atoms:
                if atom.pdbname.strip() != "P":
                    found_P = True
                    to_delete.append(atom)
            st.deleteBond(st.atom[i],atom)
            if found_P == False:
                #for other_atom in st.atom[i].getResidue().atom:
                #    if other_atom.pdbname.strip() == "P":
                #        st.addBond(st.atom[i],other_atom,1)
                for j in range(1,len(st.atom)+1):
                    if st.atom[j].pdbname.strip() == "P":
                         if st.atom[j].resnum == st.atom[i].resnum:
                             st.addBond(st.atom[i],st.atom[j],1)

    # make sure there are no bonds between residues that shouldn't exist
    for i in range(1,len(st.atom)+1):
        if st.atom[i].pdbname.strip() not in ["O3'","P"]:
            to_delete = []
            for other_atom in st.atom[i].bonded_atoms:
                if other_atom.resnum != st.atom[i].resnum:
                    to_delete.append(other_atom)
            for other_atom in to_delete:
                st.deleteBond(st.atom[i],other_atom)
            if st.atom[i].pdbname.strip() == "O4":
                other_atom = [other_atom for other_atom in st.atom[i].bonded_atoms][0]
                assert other_atom.pdbname.strip() == "C4", (other_atom.pdbname.strip(), other_atom.index)
                #if other_atom.pdbname.strip() == "C4":
                #    st.addBond(st.atom[i],other_atom,2) # thymine amide oxygen
                #else:
                #    st.deleteBond(st.atom[i],other_atom)
    
    # add bonds between some atoms in G and C
    for residue in st.residue:
        if residue.pdbres.strip() in ["DG","DA"]:
            for index in residue.getAtomIndices():
                if st.atom[index].pdbname.strip() == "C4":
                    C4 = index
                elif st.atom[index].pdbname.strip() == "C5":
                    C5 = index
                elif st.atom[index].pdbname.strip() == "C2":
                    C2 = index
                elif st.atom[index].pdbname.strip() == "N3":
                    N3 = index
                elif st.atom[index].pdbname.strip() == "N7":
                    N7 = index
                elif st.atom[index].pdbname.strip() == "C8":
                    C8 = index
                elif st.atom[index].pdbname.strip() == "O6":
                    O6 = index
                elif st.atom[index].pdbname.strip() == "C6":
                    C6 = index
                elif st.atom[index].pdbname.strip() == "N1":
                    N1 = index
            st.addBond(st.atom[C4],st.atom[C5],2)
            st.addBond(st.atom[C2],st.atom[N3],2)
            st.addBond(st.atom[N7],st.atom[C8],2)
            if residue.pdbres.strip() == "DG":
                st.addBond(st.atom[C6],st.atom[O6],2)
            if residue.pdbres.strip() == "DA":
                st.addBond(st.atom[C6],st.atom[N1],2)
        elif residue.pdbres.strip() == "DC":
            for index in residue.getAtomIndices():
                if st.atom[index].pdbname.strip() == "O2":
                    O2 = index
                elif st.atom[index].pdbname.strip() == "C2":
                    C2 = index
                elif st.atom[index].pdbname.strip() == "C5":
                    C5 = index
                elif st.atom[index].pdbname.strip() == "C6":
                    C6 = index
                elif st.atom[index].pdbname.strip() == "C4":
                    C4 = index
                elif st.atom[index].pdbname.strip() == "N3":
                    N3 = index
            st.addBond(st.atom[O2],st.atom[C2],2)
            st.addBond(st.atom[C5],st.atom[C6],2)
            st.addBond(st.atom[C4],st.atom[N3],2)
        elif residue.pdbres.strip() == "DT":
            for index in residue.getAtomIndices():
                if st.atom[index].pdbname.strip() == "O2":
                    O2 = index
                elif st.atom[index].pdbname.strip() == "C2":
                    C2 = index
                elif st.atom[index].pdbname.strip() == "O4":
                    O4 = index
                elif st.atom[index].pdbname.strip() == "C4":
                    C4 = index
                elif st.atom[index].pdbname.strip() == "C5":
                    C5 = index
                elif st.atom[index].pdbname.strip() == "C6":
                    C6 = index
            st.addBond(st.atom[O2],st.atom[C2],2)
            st.addBond(st.atom[O4],st.atom[C4],2)
            st.addBond(st.atom[C5],st.atom[C6],2)
    
    #########################################################################
    ### THE PROTEIN PREP, RUN BELOW, REQUIRES A HYDROGEN ON N1 IN DG,
    ### SO WE NEED TO ADD HYDROGENS TO THAT AND THYMINE C6, BUT NOT OTHERS
    build.add_hydrogens(st)

    # add OP3 to 5' terminal phosphorus
    resnums = []
    for residue in st.residue:
        resnums.append(residue.resnum)
    residue_5p = min(resnums)
    for residue in st.residue:
        if residue.resnum == residue_5p:
            for atom in residue.atom:
                if atom.pdbname.strip() == "P":
                    for other_atom in atom.bonded_atoms:
                        if other_atom.element == "H":
                            other_atom.element = "O"
                            other_atom.pdbname = "OP3 "
                            other_atom.formal_charge = -1

    to_delete = []
    for i in range(1,len(st.atom)+1):
        if st.atom[i].element == "H":
            for other_atom in st.atom[i].bonded_atoms:
                if other_atom.getResidue().pdbres.strip() in ["DC","DA","DT"]:
                     to_delete.append(i)
                #elif other_atom.getResidue().pdbres.strip() == "DT":
                #    if other_atom.pdbname.strip() != "C6":
                #        to_delete.append(i)
                elif other_atom.getResidue().pdbres.strip() == "DG":
                    if other_atom.pdbname.strip() != "N1":
                        to_delete.append(i)
                else:
                    #print('weird pdbres!')
                    #exit()
                    return 'weird pdbres!'
            
                #if (other_atom.element != 'N' or \
                #       other_atom.getResidue().pdbres.strip() != "DG")\
                #       and (other_atom.pdbname.strip() != "C6" or \
                #            other_atom.getResidue().pdbres.strip() != "DT"):
                #    to_delete.append(i)
    st.deleteAtoms(to_delete)
    ###########################################################################
   
    ############################################################################
    '''
    # so now we have our bonds figured out
    # but we should prepare (which includes adding hydrogens) before minimizing
    #st.write(sys.argv[2])
    st.write(out_path)
    #subprocess.run(['bash','proteinprep.sh',sys.argv[2].split('.')[0]])
    subprocess.run(['bash','proteinprep.sh',out_path.split('.')[0]])
    # need to wait for the job to finish
    #while not os.path.exists(f'{sys.argv[2].split(".")[0]}_prepwizard.mae'):
    #    if os.path.exists(f'{sys.argv[2].split(".")[0]}-001/{sys.argv[2].split(".")[0]}-001.log.failed.1'):
    #        print(f'{sys.argv[2].split(".")[0]} proteinprep job failed')
    #	exit()
    while not os.path.exists(f'{out_path.split(".")[0]}_prepwizard.mae'):
        if os.path.exists(f'{out_path.split(".")[0]}-001/{out_path.split(".")[0]}-001.log.failed.1'):
            return f'{out_path.split(".")[0]} proteinprep job failed'
        time.sleep(0.01)
    '''
    ###########################################################################
    # protein prep is silly, so we need to remove hydrogens from phosphates
    #st = structure.StructureReader.read(f'{sys.argv[2].split(".")[0]}_prepwizard.mae')
    '''st = structure.StructureReader.read(f'{out_path.split(".")[0]}_prepwizard.mae')'''
    # don't need to load structure in again if we aren't doing proteinprep
    to_delete = []
    for i in range(1,len(st.atom)+1):
        if st.atom[i].element == "H":
            for other_atom in st.atom[i].bonded_atoms:
                if other_atom.element == "P":
                    st.atom[i].element = "O"
                    st.atom[i].formal_charge = -1
                    st.atom[i].pdbname = " OP3" 
                    # " O3P" should be used according to
                    # https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
                    # version 3.3 released in 2011, which appears to be the most recent update
      
                elif other_atom.pdbname.strip() == "OP2":
                    to_delete.append(st.atom[i])
                    other_atom.formal_charge = -1
                    st.deleteBond(st.atom[i],other_atom)
    st.deleteAtoms(to_delete)
    ########################################################################
    ### NOW WE NEED TO MAKE SURE ALL OP2 AND OP3 ARE NEGATIVE (FOR SOME REASON,
    ### IT SEEMS THAT PROTEINPREP DOESN'T ALWAYS ADD HYDROGENS TO PHOSPHATES,
    ### RESULTING IN SOME CHARGE ADJUSTMENTS THAT MIGHT BE MISSED BY THE ABOVE
    ### LOOP)
    for i in range(1,len(st.atom)+1):
        if st.atom[i].pdbname.strip() in ["OP2","OP3"]:
            st.atom[i].formal_charge = -1
    ########################################################################
    #########################################################################
    ### CHECK CHIRALITY
    for i in range(1,len(st.atom)+1):
        if st.atom[i].pdbname.strip() == "C1'":
            if st.atom[i].chirality != "R":
                for other_atom in st.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-st.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]
        elif st.atom[i].pdbname.strip() == "C3'":
            if st.atom[i].chirality != "S":
                for other_atom in st.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-st.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]
        elif st.atom[i].pdbname.strip() == "C4'":
            if st.atom[i].chirality != "R":
                for other_atom in st.atom[i].bonded_atoms:
                    if other_atom.element == "H":
                        old_vec = [other_atom.xyz[j]-st.atom[i].xyz[j] for j in range(3)]
                        other_atom.xyz = [other_atom.xyz[j]-2*old_vec[j] for j in range(3)]
    ########################################################################
    ### MINIMIZE AND CLEAN UP
    build.add_hydrogens(st) # '''seeing if we can make it work without protein prep'''
    minimization_options = minimizer.MinimizationOptions(max_step=100)
    try:
        minimizer.minimize_structure(st,minimization_options)
    except Exception as e: #AtomTypingFailure
        lines_list = str(e).split('\n')
        for line_counter,line in enumerate(lines_list):
            if line.find('mmlewis error:') != -1:
                #print(f'post-proteinprep minimization failure: {line}')
                return f'minimization failure: {line}\n    {lines_list[line_counter+1]}\nout_path: {out_path}'
    #finally:
        #st.write(sys.argv[2])
        #os.remove(f'{sys.argv[2].split(".")[0]}_prepwizard.mae')
        #os.remove(f'{sys.argv[2].split(".")[0]}.log')
        #shutil.rmtree(f'{sys.argv[2].split(".")[0]}-001')
        #st.write(out_path) # don't want to write bad structures now
                       # (that was just for debugging)
        '''
        # seeing if we can make it work without proteinprep
        #os.remove(f'{out_path.split(".")[0]}_prepwizard.mae')
        #os.remove(f'{out_path.split(".")[0]}.log')
        #shutil.rmtree(f'{out_path.split(".")[0]}-001')
        '''   
    st.write(out_path)
    st.write(f'{out_path[:-5]}.mae') # also need maestro format
    return "OK" # if we've made it to the end, everything should be alright
	    
if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    in_out_tuple = (in_path,out_path)
    result = mcpdb2mae(in_out_tuple)
    print(result,flush=True)
