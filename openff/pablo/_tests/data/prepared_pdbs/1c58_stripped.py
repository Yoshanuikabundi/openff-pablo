import MDAnalysis as mda
import requests

url = "https://files.rcsb.org/download/1c58.pdb.gz"
response = requests.get(url)

with open("1c58.pdb.gz", "wb") as file:
    file.write(response.content)
u = mda.Universe("1c58.pdb.gz", topology_format="PDB")
sel = u.select_atoms('(not resname HOH MG) and not(altloc B)')
for atom in sel.atoms:
    atom.altLoc = ''
sel.write("1c58_stripped.pdb")




