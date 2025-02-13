import inspect
from io import StringIO
from pathlib import Path

from openff.toolkit import Topology


def main():
    in_fn = Path("2MUM_neutralized.json")
    out_fn_template = "2MUM_{}.pdb"

    # top = topology_from_pdb(in_fn)
    top = Topology.from_json(in_fn.read_text())

    reuse_serial(top, out_fn_template)
    reuse_resseq(top, out_fn_template)
    # letters_in_serial(top, out_fn_template)
    # letters_in_resseq(top, out_fn_template)
    # extended_width_serial(top, out_fn_template)
    # extended_width_resseq(top, out_fn_template)
    # icode(top, out_fn_template)
    # discontiguous_serial(top, out_fn_template)
    # discontiguous_resseq(top, out_fn_template)


def write_pdb(top: Topology, out_fn_template: str):
    with StringIO() as f:
        top.to_file(f, file_format="PDB", keep_ids=True)
        pdbtext = f.getvalue()

    lines = []
    i = 0
    for line in pdbtext.splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            line = (
                line[:6]
                + f"{top.atom(i).metadata['atom_serial']: >5}"
                + line[11:22]
                + f"{top.atom(i).metadata['residue_number']: >4}"
                + line[26:]
            )
            i += 1
        lines.append(line)

    out_fn = Path(out_fn_template.format(inspect.stack()[1][3]))
    out_fn.write_text("\n".join(lines))


def reuse_serial(top: Topology, out_fn_template: str):
    top = Topology(top)

    res_iter = iter(top.residues)
    for i in range(5):
        reused_serial = next(res_iter).atom(-1).metadata["atom_serial"]
    for res in res_iter:
        for atom in res.atoms:
            atom.metadata["atom_serial"] = reused_serial

    write_pdb(top, out_fn_template)


def reuse_resseq(top: Topology, out_fn_template: str):
    top = Topology(top)

    res_iter = iter(top.residues)
    for i in range(5):
        reused_resseq = next(res_iter).identifier[1]
    for res in res_iter:
        for atom in res.atoms:
            atom.metadata["residue_number"] = reused_resseq

    write_pdb(top, out_fn_template)


def letters_in_serial(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def letters_in_resseq(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def extended_width_serial(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def extended_width_resseq(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def icode(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def discontiguous_serial(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


def discontiguous_resseq(top: Topology, out_fn_template: str):
    top = Topology(top)

    write_pdb(top, out_fn_template)


if __name__ == "__main__":
    main()
