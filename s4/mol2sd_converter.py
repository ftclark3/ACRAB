from schrodinger import structure
def mol2tosd(in_file,out_file):
    structures = []
    for mol2 in structure.StructureReader(in_file):
        structures.append(mol2)
    writer = structure.StructureWriter(out_file)
    writer.extend(structures)
    writer.close()
def sdtomol2(in_file,out_file):
    structures = [] 
    for sd in structure.StructureReader(in_file):
        structures.append(sd)
    writer = structure.StructureWriter(out_file)
    writer.extend(structures)
    writer.close()

if __name__ == "__main__":
    import sys 
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    if in_file.split('.')[-1] == 'sd':
        sdtomol2(in_file,out_file)
    elif in_file.split('.')[-1] == 'mol2':
        mol2tosd(in_file,out_file)
    else:
        raise SystemExit("file format not recognized!")
