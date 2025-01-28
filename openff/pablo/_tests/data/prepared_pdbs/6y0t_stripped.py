import MDAnalysis as mda
import requests

url = "https://files.rcsb.org/download/6y0t.pdb.gz"
response = requests.get(url)

with open("6y0t.pdb.gz", "wb") as file:
    file.write(response.content)
u = mda.Universe("6y0t.pdb.gz", topology_format="PDB")
sel = u.select_atoms('(not resname HOH MG K) and not(altloc B)')
for atom in sel.atoms:
    atom.altLoc = ''
sel.write("6y0t_stripped.pdb")




