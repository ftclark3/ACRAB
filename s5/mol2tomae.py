from schrodinger import structure
def mol2tomae(mol2_path):
    st = structure.StructureReader.read(mol2_path)

    # fix weird residue numbering thing
    #indices = []
    #for counter,residue in enumerate(st.residue):
    #    indices.append(residue.getAtomIndices())
    #    one_indexed = counter + 1
    #for counter,residue in enumerate(st.residue):
    #    one_indexed = counter +1
    #    residue.resnum = one_indexed

    name_part = mol2_path.split('/')[-1][:-5]
    st.write(f'{name_part}.mae')

if __name__ == "__main__":
    import sys
    mol2_path = sys.argv[1]
    mol2tomae(mol2_path)
