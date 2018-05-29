from Bio.PDB.Chain import Chain


p = PDBParser()
structure=p.get_structure('cytb5',mypdb.pdb)
model=structure[0]

my_chain = Chain("C")
model.add(my_chain)

print(model)
