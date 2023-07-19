# example usage: python acrab.py config/sequence_name.toml
def bash_or_python_run(script,f,program):
    if program != 'schrodinger':
        subprocess.run([program,script,\
            config['conda_shell'],config["acrab_path"],config["jobname"],\
            config['mda_shell'],config['schrodinger_path'],config['SPOT-RNA_path'],\
            config['MCFold_path'],config['MC-Symizer_path'],config['MC-Sym_path'],str(config['affinity_info'])],\
            stdout=f,stderr=f,text=True)
    elif program == 'schrodinger':
        subprocess.run([f'{config["schrodinger_path"]}/run','python3',script,\
            config['conda_shell'],config["acrab_path"],config["jobname"],\
            config['mda_shell'],config['schrodinger_path'],config['SPOT-RNA_path'],\
            config['MCFold_path'],config['MC-Symizer_path'],config['MC-Sym_path'],str(config['affinity_info'])],\
            stdout=f,stderr=f,text=True)

###########################################################
### STAGE 0: GET CONFIG INFORMATION
print("STAGE 0: GET CONFIG INFORMATION")
import toml
import sys
import os
import subprocess
config = toml.load(sys.argv[1])
print("STAGE 0: COMPLETE!")

#############################################################
### STAGE 1: SET UP JOB
# stage 1: set up expected directory structure
# and check that jobname doesn't already exist
print("STAGE 1: SET UP JOB")
os.chdir(config['acrab_path'])
if config['jobname'] in os.listdir('.'):
    print(f'{config["jobname"]} directory already exists!\
               \n Aborting run')
    raise SystemExit()
else:
    os.mkdir(config['jobname'])
    os.chdir(config['jobname'])
    os.mkdir('s2')
    os.mkdir('s3')
    os.mkdir('s4')
    os.mkdir('s5')
    # write fasta file corresponding to sequence
    with open(f'{config["jobname"]}.fasta','w') as f:
        f.write(f">{config['jobname']}\n{config['sequence']}")
    print("STAGE 1: COPYING CONFIG FILE TO JOB RUN DIRECTORY \
              FOR FUTURE REFERENCE")
    subprocess.run(['cp',f'../{sys.argv[1]}','./backup_config.toml'])
    print("STAGE 1: COMPLETE!")
##############################################################
### STAGE 2: PREDICT SECONDARY STRUCTURES
print("STAGE 2: PREDICT SECONDARY STRUCTURES")
print(f"STAGE 2: COPYING STAGE RUN SCRIPT(s)\n\
       {config['acrab_path']}/{config['s2_protocol']}\
       AND ASSOCIATED HELPER SCRIPTS\
       \nto current job run directory \n{os.getcwd()}/s2")
os.chdir('s2')
subprocess.run(['cp',f'../../{config["s2_protocol"]}','.'])
s2_helpers = config["s2_helpers"].split(',')
for helper in s2_helpers:
    subprocess.run(['cp',f'../../{helper}','.'])
print("STAGE 2: RUNNING PREDICTION")
extension = config["s2_protocol"].split('.')[-1]
name_part = config["s2_protocol"].split('/')[-1]
beginning_list = name_part.split(".")[:-1]
beginning = ""
for item in beginning_list:
    beginning += item
if extension in ["sh","bash"]:
    with open(f'{beginning}.log','w') as f:
        # for broad compatibility, we'll pass all relevant
        # config information to bash script as arguments
        # whether it uses that information is up to the
        # particulars of the custom script
        bash_or_python_run(name_part,f,'bash') 
elif extension == "py":
    with open(f'{beginning}.log','w') as f:
        #subprocess.run([config['base_python_path'],config['s2_protocol']])
        bash_or_python_run(name_part,f,'python')
else:
    print("STAGE 2: PREDICTION SCRIPT NOT RECOGNIZED.\
            IT IS EXPECTED TO BE A BASH OR PYTHON SCRIPT")
    raise SystemExit()
os.chdir('..')
print("STAGE 2: COMPLETE!")
#print(f'current directory: {os.getcwd()}')
# at this point, we're in the root job directory
#############################################################
### STAGE 3: PREDICT TERTIARY STRUCTURES
print("STAGE 3: PREDICT TERTIARY STRUCTURES")
print(f"STAGE 3: COPYING STAGE RUN SCRIPT\n\
       {config['acrab_path']}/{config['s3_protocol']}\
       AND ASSOCIATED HELPER SCRIPTS\
       \nto current job run directory \n{os.getcwd()}/s3")
os.chdir('s3')
subprocess.run(['cp',f'../../{config["s3_protocol"]}','.'])
s3_helpers = config["s3_helpers"].split(',')
if s3_helpers != [""]:
    for helper in s3_helpers:
        subprocess.run(['cp',f'../../{helper}','.'])
print("STAGE 3: RUNNING PREDICTION")
extension = config["s3_protocol"].split('.')[-1]
name_part = config["s3_protocol"].split('/')[-1]
beginning_list = name_part.split(".")[:-1]
beginning = ""
for item in beginning_list:
    beginning += item
if config["s3_schrodinger"]==True:
    with open(f'{beginning}.log','w') as f:
        bash_or_python_run(name_part,f,'schrodinger')
elif extension in ["sh","bash"]:
    with open(f'{beginning}.log','w') as f:
        # for broad compatibility, we'll pass all relevant
        # config information to bash script as arguments
        # whether it uses that information is up to the
        # particulars of the custom script
        bash_or_python_run(name_part,f,'bash')
elif extension == "py":
    with open(f'{beginning}.log','w') as f:
        #subprocess.run([config['base_python_path'],config['s2_protocol']])
        bash_or_python_run(name_part,f,'python')
else:
    print("STAGE 3: PREDICTION SCRIPT NOT RECOGNIZED.\
            IT IS EXPECTED TO BE A BASH OR PYTHON SCRIPT")
    raise SystemExit()
os.chdir('..')
print("STAGE 3: COMPLETE!")
#############################################################
### STAGE 4: PREDICT TERTIARY STRUCTURES
print("STAGE 4: DOCK TO STRUCTURES")
print(f"STAGE 4: COPYING STAGE RUN SCRIPT\n\
       {config['acrab_path']}/{config['s4_protocol']}\
       AND ASSOCIATED HELPER SCRIPTS\
       \nto current job run directory \n{os.getcwd()}/s4")
os.chdir('s4')
subprocess.run(['cp',f'../../{config["s4_protocol"]}','.'])
s4_helpers = config["s4_helpers"].split(',')
if s4_helpers != [""]:
    for helper in s4_helpers:
        subprocess.run(['cp',f'../../{helper}','.'])
print("STAGE 4: DOCKING")
extension = config["s4_protocol"].split('.')[-1]
name_part = config["s4_protocol"].split('/')[-1]
beginning_list = name_part.split(".")[:-1]
beginning = ""
for item in beginning_list:
    beginning += item
if config["s4_schrodinger"]==True:
    with open(f'{beginning}.log','w') as f:
        bash_or_python_run(name_part,f,'schrodinger')
elif extension in ["sh","bash"]:
    with open(f'{beginning}.log','w') as f:
        # for broad compatibility, we'll pass all relevant
        # config information to bash script as arguments
        # whether it uses that information is up to the
        # particulars of the custom script
        bash_or_python_run(name_part,f,'bash')
elif extension == "py":
    with open(f'{beginning}.log','w') as f:
        #subprocess.run([config['base_python_path'],config['s2_protocol']])
        bash_or_python_run(name_part,f,'python')
else:
    print("STAGE 4: PREDICTION SCRIPT NOT RECOGNIZED.\
            IT IS EXPECTED TO BE A BASH OR PYTHON SCRIPT")
    raise SystemExit()
os.chdir('..')
print("STAGE 4: COMPLETE!")
#########################################################3
### STAGE 5: LOOP MODELING
print("STAGE 5: LOOP MODELING")
print(f"STAGE 5: COPYING STAGE RUN SCRIPT\n\
       {config['acrab_path']}/{config['s5_protocol']}\
       AND ASSOCIATED HELPER SCRIPTS\
       \nto current job run directory \n{os.getcwd()}/s5")
os.chdir('s5')
subprocess.run(['cp',f'../../{config["s5_protocol"]}','.'])
s5_helpers = config["s5_helpers"].split(',')
if s5_helpers != [""]:
    for helper in s5_helpers:
        subprocess.run(['cp',f'../../{helper}','.'])
print("STAGE 5: MODELING LOOPS")
extension = config["s5_protocol"].split('.')[-1]
name_part = config["s5_protocol"].split('/')[-1]
beginning_list = name_part.split(".")[:-1]
beginning = ""
for item in beginning_list:
    beginning += item
if config["s5_schrodinger"]==True:
    with open(f'{beginning}.log','w') as f:
        bash_or_python_run(name_part,f,'schrodinger')
elif extension in ["sh","bash"]:
    with open(f'{beginning}.log','w') as f:
        # for broad compatibility, we'll pass all relevant
        # config information to bash script as arguments
        # whether it uses that information is up to the
        # particulars of the custom script
        bash_or_python_run(name_part,f,'bash')
elif extension == "py":
    with open(f'{beginning}.log','w') as f:
        #subprocess.run([config['base_python_path'],config['s2_protocol']])
        bash_or_python_run(name_part,f,'python')
else:
    print("STAGE 5: PREDICTION SCRIPT NOT RECOGNIZED.\
            IT IS EXPECTED TO BE A BASH OR PYTHON SCRIPT")
    raise SystemExit()
os.chdir('..')
print("STAGE 5: COMPLETE!")

'''
###########################################################3
### STAGE 6: EVALUATE RESULTS
print("STAGE 6: EVALUATE RESULTS")
print(f"STAGE 6: COPYING STAGE RUN SCRIPT\n\
       {config['acrab_path']}/{config['s6_protocol']}\
       AND ASSOCIATED HELPER SCRIPTS\
       \nto current job run directory \n{os.getcwd()}/s6")
os.chdir('s6')
subprocess.run(['cp',f'../../{config["s6_protocol"]}','.'])
s4_helpers = config["s6_helpers"].split(',')
if s5_helpers != [""]:
    for helper in s5_helpers:
        subprocess.run(['cp',f'../../{helper}','.'])
print("STAGE 6: EVALUATING RESULTS")
extension = config["s6_protocol"].split('.')[-1]
name_part = config["s6_protocol"].split('/')[-1]
beginning_list = name_part.split(".")[:-1]
beginning = ""
for item in beginning_list:
    beginning += item
if config["s6_schrodinger"]:
    with open(f'{beginning}.log','w') as f:
        bash_or_python_run(name_part,f,'schrodinger')
elif extension in ["sh","bash"]:
    with open(f'{beginning}.log','w') as f:
        # for broad compatibility, we'll pass all relevant
        # config information to bash script as arguments
        # whether it uses that information is up to the
        # particulars of the custom script
        bash_or_python_run(name_part,f,'bash')
elif extension == "py":
    with open(f'{beginning}.log','w') as f:
        #subprocess.run([config['base_python_path'],config['s2_protocol']])
        bash_or_python_run(name_part,f,'python')
else:
    print("STAGE 6: PREDICTION SCRIPT NOT RECOGNIZED.\
            IT IS EXPECTED TO BE A BASH OR PYTHON SCRIPT")
    raise SystemExit()
os.chdir('..')
print("STAGE 6: COMPLETE!")
'''
