import sys
import os 
import subprocess
import multiprocessing
#from mcpdb2mae import mcpdb2mae
from glob import glob

# to modify the number of predictions that we direct
# mcsym to make in the .mcc script
# (default is 1000) 
def modify_num_predictions_only(script,number):
    lines = []
    with open(script,'r') as f:
        for line in f:
            if "model_limit = 1000" in line:
                location = line.index("model_limit = 1000")
                #foo = f"{line[:location]}model_limit = {number}{line[location+(18+2-len(str(number))):]}"
                foo = f"{line[:location]}model_limit = {number},{line[location+(19-2+len(str(number))):]}"
                lines.append(foo)
                continue
            else:
                lines.append(line)
    to_write = ""
    for line in lines:
        to_write += line
    with open(script,'w') as f:
        f.write(to_write)


# need to run schrodinger-dependent parts as subprocess
# because multiprocessing is disable in schrodinger, preventing
# the schrodinger-installed python from being able to run this
# script efficiently
def mcpdb2mae_runner(in_out_tuple):
    in_file = in_out_tuple[0]
    #print(f'in_file: {in_file}',flush=True)
    out_file = in_out_tuple[1]
    #print(f'out file: {out_file}',flush=True)
    in_file = in_file.strip()
    out_file = out_file.strip()
    #print(in_file,flush=True)
    #print(out_file,flush=True)
    completed_process = subprocess.run([f'{schrodinger_path}/run',\
        'python3','mcpdb2mae.py',in_file,out_file],capture_output=True)
    return (completed_process.stdout,completed_process.stderr)
    # a string representing the result of the script

def make_tertiaries(counter,seq_dbn_tuple):
    # have kind of a hacky way to get a functional 
    # mcsym script, based on mcsymize.py followed by
    # some edits to the output script
    
    #print(counter)
    #print(seq_dbn_tuple) 
    sequence = seq_dbn_tuple[0]
    dbn = seq_dbn_tuple[1]
    with open(f'dbn_{counter}.mcc','w') as f:
        subprocess.run(['python',f'{mcsymizer_path}/mcsymizer.py',\
            '--sequence1',sequence,'--structure1',dbn,'--db_path',\
            f'{mcsym_path}/MCSYM-DB','--no_dangling','--name',f'mcsym_{counter}'],\
            stdout=f,text=True)
            #'>',f'dbn_{counter}.mcc'])
    fix_mcc_script(sequence,dbn,f'dbn_{counter}.mcc')
    with open(f'dbn_{counter}.log','w') as f:    
        subprocess.run([f'{mcsym_path}/mcsym',f'dbn_{counter}.mcc',\
            '-D',f'{mcsym_path}/MCSYM-DB'],stdout=f,text=True)
     
''' 
def cleanup_tertiaries():
    # then need to get into maestro format and clean up
    # will copy back to pdb at the beginning of the next stage
    # for surface docking
    in_files = glob('*.pdb.gz') 
    out_files = []
    for file in in_files:
        out_files.append(f'{file[:-7]}.pdb')
    multiprocessing.set_start_method('fork')
    with multiprocessing.pool.Pool as pool:
        for result in pool.map(mcpdb2mae,zip(in_files,out_files)):
            if result != "OK":
                print(result) 
                # shouldn't need to flush because of fork
'''

def fix_mcc_script(sequence,dbn,script):
    sequence = sequence.replace("T","U")
    first_backward_flag = False
    dot_counts = 0
    last_dot = -100
    for i in range(1,len(dbn)+1):
        if dbn[-i] in '[(':
            break # crude way to make sure we aren't getting hairpins
        if dbn[-i] in '])':
            first_backward_flag = True
            dot_counts = 0 # need to reset the count
            continue
        if dbn[-i] == '.' and first_backward_flag == True: 
            dot_counts += 1
            #if dot_counts == 4:
            if dot_counts == 4 and dbn[-i-1] not in '([': # hacky way to deal with 4-nt hairpins 
                last_dot = i-3
                break
    print(f'last_dot:{last_dot}',flush=True)
    if last_dot != -100: 
        locations = []
        context = []
        beginning = []
        end = []
        index = len(dbn) - last_dot + 1 
        # residues are 1-indexed in MCSYM chain identification
        with open(script,'r') as f:
            for counter1,line in enumerate(f):
                #if f'A{last_dot}' in line and '<-' in line:
                if f'A{index}' in line and '<-' in line:
                    #print(line,flush=True)
                    locations.append((counter1,line))
            #assert len(locations)==1,f'len(locations){len(locations)},last_Dot:{last_dot}'
            if len(locations) > 1: #mcsym broke up the 3 nt thing into 2 2nt thing,
                                   # so maybe we're okay
                                   # still, we need to modify the number of predictions
                modify_num_predictions_only(script,100)
                return
                
            location = locations[0]
        with open(script,'r') as f:
            for counter2,line in enumerate(f):
                if counter2 < location[0]-2:
                    beginning.append(line)
                elif location[0]-2 <= counter2 <= location[0]+2:
                    context.append(line)
                elif location[0]+2 < counter2:
                    end.append(line)
                else:
                    print("weird!",flush=True)
                    raise SystemExit()
        print(f'first context: {context}',flush=True)
        assert 'ss2' in context[1], context # hopefully it's the error I think it is
        temp_1 = context[1].split('/')
        for counter,element in enumerate(temp_1):
            if element == "ss2":
                temp_1[counter] = "ss3"
                #temp_1[counter+1] = sequence[last_dot-1:last_dot+2]
                temp_1[counter+1] = sequence[index-1:index+2]
        new_1 = ''
        for item in temp_1:
            new_1 += item
            if item != temp_1[-1]: # don't need trailing slash
                new_1 += '/'
        context[1] = new_1
        
        hashtag_counter = 0
        A_counter = 0
        context[2] = list(context[2])
        # need to make it into a list of characters to modify specific indices
        for counter,element in enumerate(context[2]):
            if element == "#":
                hashtag_counter += 1
                if hashtag_counter == 2:
                    hashtag_counter = 0 # don't care about this anymore
                    assert int(context[2][counter+1]) == 2, context[2][counter+1]
                    context[2][counter+1] = '3' # need string for ''.join(context[2]) later
            if element == "A":
                A_counter += 1
                if A_counter == 2:
                    A_counter = 0
                    #assert int(''.join(context[2][counter+1:counter+3]))==last_dot,\
                    assert int(''.join(context[2][counter+1:counter+3]))==index,\
                        f'counter:{"".join(context[2][counter+1:counter+3])},last_dot:{last_dot}'
                        #f'counter:{"".join(context[2][counter+1:counter+3])},last_dot:{last_dot}'
                    #context[2] = f'{"".join(context[2][:counter+1])}{last_dot+1}{"".join(context[2][counter+3:])}'
                    context[2] = f'{"".join(context[2][:counter+1])}{index+1}{"".join(context[2][counter+3:])}'
        print(f'last context:{context}',flush=True)    
        to_write = ''
        for line in beginning:
            to_write += line
        for line in context:
            to_write += line
        for line in end:
            # similar to content of function modify_num_predictions_only, but not exactly the same 
            # because we already have the file loaded in and don't need to open it
            if "model_limit = 1000" in line:
                number = 100 # can change this
                location = line.index("model_limit = 1000")
                #foo = f"{line[:location]}model_limit = 10{line[location+18:]}"
                #foo = f"{line[:location]}model_limit = 100{line[location+18:]}"
                #foo = f"{line[:location]}model_limit = {number}{line[location+(18+2-len(str(number))):]}"
                foo = f"{line[:location]}model_limit = {number},{line[location+(19-2+len(str(number))):]}"
                to_write += foo
                continue
            to_write += line
        with open(script,'w') as f:
            f.write(to_write)
    else:
        # corresponding to the if statement currently on line 81
        # this is the scenario where no other modification of mcsym
        # script is needed
        # still, we might only want to sample 100 structures instead of 1000
        modify_num_predictions_only(script,100)
        # content of function below
        '''
        lines = []
        with open(script,'r') as f:
            for line in f:
                if "model_limit = 1000" in line:
                    location = line.index("model_limit = 1000")
                    foo = f"{line[:location]}model_limit = 10{line[location+18:]}"
                    lines.append(foo)
                    continue
                else:
                    lines.append(line)
        to_write = ""
        for line in lines:
            to_write += line
        with open(script,'w') as f:
            f.write(to_write)
        '''

if __name__ == "__main__":
    #print(sys.argv)
    acrab_path = sys.argv[2]
    jobname = sys.argv[3]
    schrodinger_path = sys.argv[5] 
    mcsymizer_path = sys.argv[8]
    mcsym_path = sys.argv[9]
    os.chdir(acrab_path)
    dbns = []
    with open(f'{acrab_path}/{jobname}/{jobname}.fasta','r') as f:
        for line in f:
            if line[0] != '>':
                sequence = line.strip()
    with open(f'{acrab_path}/{jobname}/s2/dot_brackets.txt','r') as f:
        for line in f:
            dbns.append(line.strip())
    #print(sequence)
    #print(dbns)
    os.chdir(f'{jobname}/s3')
    #with multiprocessing.Pool() as pool:
    #    pool.map(make_tertiaries,enumerate(dbns)) 
        # function followed by list of different inputs to use
    #map(make_tertiaries,enumerate(dbns))
    # this below worked, so apparently I don't understand the map command
    #for counter, dbn in enumerate(dbns):
    #    make_tertiaries((counter,dbn))
    sequences = []
    while len(sequences) < len(dbns):
        sequences.append(sequence)
    #for counter, seqence_dbn_tuple in enumerate(zip(sequences,dbns)):
    #    make_tertiaries(counter,sequence_dbn_tuple)
    multiprocessing.set_start_method('fork')
    with multiprocessing.Pool() as pool:
        pool.starmap(make_tertiaries,enumerate(zip(sequences,dbns)))
        #for result in pool.starmap(make_tertiaries,enumerate(dbns)):
            #print(result)
            
    # need unzipped pdb files
    pdbgzs = glob(f'{acrab_path}/{jobname}/s3/*.pdb.gz')
    for pdbgz in pdbgzs:
        subprocess.run(['gunzip',pdbgz])

    # now need to get into maestro format and clean up
    # will copy back to pdb at the beginning of the next stage
    # for surface docking
    in_files = glob('*.pdb')
    out_files = []
    for file in in_files:
        out_files.append(f'{file[:-4]}.mol2') # need mol2 format for docking
    with multiprocessing.pool.Pool() as pool:
        for result in pool.map(mcpdb2mae_runner,zip(in_files,out_files)):
            if str(result[0],'UTF-8').strip() != "OK":
                # result is a tuple of stdout,stderr
                print(result,flush=True)
                # flushing because output is 1 level up (in acrab.py)
                # not sure if this really makes sense, but just to be safe
