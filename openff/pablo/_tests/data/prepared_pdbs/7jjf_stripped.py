import MDAnalysis as mda
import requests

url = "https://files.rcsb.org/download/7jjf.pdb.gz"
response = requests.get(url)

with open("7jjf.pdb.gz", "wb") as file:
    file.write(response.content)

u = mda.Universe("7jjf.pdb.gz", topology_format="PDB")
sel = u.select_atoms('not resname HOH MG')
sel.write("7jjf_stripped.pdb")




