import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

def identify_good_poses(receptor_file,mol2_file):
    # mda can't read in an sd file, so we have to 
    # first convert it to mol2
    good_poses = []
    receptor = mda.Universe(receptor_file)
    receptor_atoms = receptor.select_atoms('all')
    ligand_poses = mda.Universe(mol2_file)
    for counter,ligand_pose in enumerate(ligand_poses.trajectory):
        distances = distance_array(ligand_pose,receptor_atoms)
        # comparing an atomgroup being dynamically modified over the trajectory
        # to a static numpy array
        flag = 0
        for row in distances:
            for item in row: # note that we're double counting over the full matrix
                             # rather than a triangle
                if item < 4: #################PARAMETER################
                    flag += 1
        if int(flag/2) < 4: # this pose not worth exploring
                            # (we're requiring 2 atoms to be within
                            # a reasonable distance)
            continue
        else:
            good_poses.append(counter)
    return good_poses

def write_sd_strings(sd_file,good_poses):
    # sd_file: file with a bunch of ligands,
    #          will be used to derive base string to which "-#.sd" 
    #          will be appended when writing new pose
    #          ex. mcsym_0-0102-3.sd is
    #          secondary structure 0, tertiary structure 102,
    #          bound pose 3
    # good_poses: list of integers corresponding to  
    #             0-indexed ligand poses in sd_file
    #             which are deemed "good" poses
   
    base_name_list = sd_file.split("_")[:-1]
    base_name = ""
    for item in base_name_list:
        base_name += item
        base_name += "_" 
        # we'll put it back
    # then remove the last one
    base_name = base_name[:-1] 
    pose_strings = []
    pose_string = ''
    with open(sd_file,'r') as f:
        # need there to be an identical sd file to the mol2 file
        for line in f:
            pose_string += line
            if line.strip() == "$$$$":
                pose_strings.append(pose_string)
                pose_string = ''
            else:
                pass
    pose_strings = [pose_string for counter,pose_string in enumerate(pose_strings)\
                    if counter in good_poses] 
    for counter,pose_string in enumerate(pose_strings):
        with open(f'{base_name}-{counter}.sd','w') as f:
            # note that we're not preserving the previous order
            # for example, if pose 5 is found to be bad, 
            # the pose formerly known as pose 6 will
            # be written as pose 5, and so on
            f.write(pose_string)

def prep_method_1(receptor_file,extensionless_filename):
    # to facilitate stage 4 workflow in method_1.sh
    good_poses = identify_good_poses(receptor_file,f'{extensionless_filename}.mol2')
    write_sd_strings(f'{extensionless_filename}.sd',good_poses)
    
if __name__ == "__main__":
    import sys
    receptor_file = sys.argv[1]
    mol2_file = sys.argv[2]
    extensionless_filename = mol2_file[:-5] # holds the ligands
    prep_method_1(receptor_file,extensionless_filename)
