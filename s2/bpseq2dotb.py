# for convert bpseq file output by SPOT-RNA to dot bracket notation
import sys
import os
import warnings
#os.chdir(f'/home/clarkf/aptamer/pipeline/{sys.argv[1]}/secondary_structures/SPOT-RNA/')
for file in os.listdir('.'):
    if 'container' in file:
        continue
    if file[-6:] != '.bpseq':
        continue
    # otherwise, we're working with the kind of file that we want
    with open(file,'r') as f:
        pairs = {}
        for line_number,line in enumerate(f):
            if line_number == 0:
                continue # just a comment line
            # otherwise, we're working with lines that show pairing info
            info = line.split(' ')
            pairs.update({info[0].strip():info[2].strip()})
    dotb = []
    for dummy in range(len(pairs.keys())):
        dotb.append('.')
    for start_base in range(1,len(pairs.keys())+1): #start_base is an int
        if int(pairs[str(start_base)]) > int(start_base):
            dotb[start_base-1] = '(' # because of 1-indexing
    # we might need to go back and change some of the # '(' to higher-order ones because of pseudoknots, but this is a good start
    # now, we'll close all with a regular ')'
    for start_base in range(1,len(pairs.keys())+1):
        if dotb[start_base-1] == '.' and pairs[str(start_base)] != '0':
            dotb[start_base-1] = ')'
    

    # now we need to identify pseudoknots
    certain_backwards_parentheses = []
    for i in range(len(dotb)):
        if dotb[i] == '(':
            pair_location = int(pairs[str(i+1)])-1 # pairs is one-indexed,
             # so we add one to the dotb index and subtract one from the pairs index
            for index in certain_backwards_parentheses:
                if i<index<pair_location:
                    dotb[i] = '['
                    dotb[pair_location] = ']' 
            if dotb[i] == '(': # if it wasn't changed
                certain_backwards_parentheses.append(pair_location) 
            # first one is always (, so the one that pairs with it must be )

    '''
    num_pairs = int(len([base for base in dotb if base in '()'])/2)
    addressed_pairs = 0
    while addressed_pairs < num_pairs:

        right_pairing_counter = 0
        left_pairing_counter = 0
        for start_base in range(1,len(pairs.keys())+1):
            if dotb[-start_base] in '([':
                right_pairing_counter += 1
                if right_pairing_counter > addressed_pairs:
                    start_index = -start_base
                    break
        for end_base in range(start_index+1,0): #remember, start_index is negative
            if dotb[end_base] in ')]':
                left_pairing_counter += 1
                if left_pairing_counter > addressed_pairs: 
                    end_index = end_base
                    break
        #print(pairs[str(len(pairs.keys())+start_index+1)]) #changing from 0-indexed dotb to 1-indexed dict
        #print(str(len(pairs.keys())+start_index+1))

        #print(str(len(pairs.keys())+end_index))
        #print('\n')
        if pairs[str(len(pairs.keys())+start_index+1)] == str(len(pairs.keys())+end_index+1): # going from 0-indexed to 1-indexed
            #dotb[end_index] = ')'
            # the above assignment should have already happened
            assert dotb[end_index] == ')', dotb[end_index]
        else: 
            # mark as part of pseudoknot
            open_right_counter = 0
            open_left_counter = 0
            for base in dotb[:end_index]:
                #print(start_index,end_index)
                if base == '(':
                    open_right_counter += 1
                elif base in ')]': # asymmetry here because we might have already changed
                                   # open left ) to ]
                    open_left_counter += 1
            if open_right_counter > open_left_counter:
                dotb[end_index] = ']'
                dotb[start_index] = '['
            #pass
        addressed_pairs += 1
    '''
    dotb_string = ''
    for item in dotb:
        dotb_string += item
    if '[' in dotb_string or ']' in dotb_string:
        warnings.warn(f'WARNING: Pseudoknot in predicted secondary structure {file}. My attempt to convert the secondary structure to dot-bracket notation may be incorrect. Please manually verify with raw output secondary structure format from the core software')
    outfile_split = file.split('.')[:-1]
    outfile = ""
    for part in outfile_split:
        outfile += part
    with open(f'{outfile}.dbn','w') as f:
        f.write(dotb_string)
    
