from schrodinger import structure
import sys
mol2 = sys.argv[1]
sd = sys.argv[2]
st = structure.StructureReader.read(mol2)
st.write(sd)
