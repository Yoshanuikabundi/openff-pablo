from typing import Union, Optional, List, Dict, Iterable, Tuple
import uuid
from io import StringIO
from numpy.typing import NDArray
import numpy as np

import rdkit.Chem
import rdkit.Chem.Draw
from IPython.display import SVG
import nglview
from nglview.base_adaptor import Structure, Trajectory

from openff.toolkit import Molecule, Topology
from openff.units import unit

import openmm

DEFAULT_WIDTH = 400
DEFAULT_HEIGHT = 300

Color = Iterable[float]
BondIndices = Tuple[int, int]

def unicode_to_ascii(s):
    import unicodedata
    
    return unicodedata.normalize('NFKD', s).encode('ascii','ignore').decode()

def draw_molecule(
    molecule: Union[Molecule, rdkit.Chem.rdchem.Mol],
    *,
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
    highlight_atoms: Optional[Union[List[int], Dict[int, Color]]] = None,
    highlight_bonds: Optional[
        Union[List[BondIndices], Dict[BondIndices, Color]]
    ] = None,
    atom_notes: Optional[Dict[int, str]] = None,
    bond_notes: Optional[Dict[BondIndices, str]] = None,
    emphasize_atoms: Optional[List[int]] = None,
    explicit_hydrogens: Optional[bool] = None,
    color_by_element: Optional[bool] = None,
) -> SVG:
    """Draw a molecule

    Parameters
    ==========

    molecule
        The molecule to draw.
    image_width
        The width of the resulting image in pixels.
    image_height
        The height of the resulting image in pixels.
    highlight_atoms
        A list of atom indices to highlight, or a map from indices to colors.
        Colors should be given as triplets of floats between 0.0 and 1.0.
    highlight_bonds
        A list of pairs of atom indices indicating bonds to highlight, or a map
        from index pairs to colors. Colors should be given as triplets of floats
        between 0.0 and 1.0.
    atom_notes
        A map from atom indices to a string that should be printed near the
        atom.
    bond_notes
        A map from atom index pairs to a string that should be printed near the
        bond.
    emphasize_atoms
        A list of atom indices to emphasize by drawing other atoms (and their
        bonds) in light grey.
    explicit_hydrogens
        If ``False``, allow uncharged monovalent hydrogens to be hidden. If
        ``True``, make all hydrogens explicit. If ``None``, defer to the
        provided molecule.
    color_by_element
        If True, color heteroatoms according to their element; if False, color
        atoms and bonds monochromatically. By default, uses black and white when
        highlight_atoms or highlight_bonds is provided, and color otherwise.

    Raises
    ======

    KeyError
        When an atom or bond in highlight_atoms or highlight_bonds is missing
        from the image, including when it is present in the molecule but hidden.
    """
    # We're working in RDKit
    try:
        rdmol = molecule.to_rdkit()
    except AttributeError:
        rdmol = rdkit.Chem.rdchem.Mol(molecule)

    # Process color_by_element argument
    if color_by_element is None:
        color_by_element = highlight_atoms is None and highlight_bonds is None

    if color_by_element:
        set_atom_palette = lambda draw_options: draw_options.useDefaultAtomPalette()
    else:
        set_atom_palette = lambda draw_options: draw_options.useBWAtomPalette()

    # Process explicit_hydrogens argument
    # If we need to remove atoms, create a map from the original indices to the
    # new ones.
    if explicit_hydrogens is None:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
    elif explicit_hydrogens:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
        rdmol = rdkit.Chem.AddHs(rdmol, explicitOnly=True)
    else:
        idx_map = {
            old: new
            for new, old in enumerate(
                a.GetIdx()
                for a in rdmol.GetAtoms()
                if a.GetAtomicNum() != 1 and a.GetMass() != 1
            )
        }
        rdmol = rdkit.Chem.RemoveHs(rdmol, updateExplicitCount=True)

    # Process highlight_atoms argument for highlightAtoms and highlightAtomColors
    # highlightAtoms takes a list of atom indices
    # highlightAtomColors takes a mapping from atom indices to colors
    if highlight_atoms is None:
        highlight_atoms = []
        highlight_atom_colors = None
    elif isinstance(highlight_atoms, dict):
        highlight_atom_colors = {
            idx_map[i]: tuple(c) for i, c in highlight_atoms.items() if i in idx_map
        }
        highlight_atoms = list(highlight_atoms.keys())
    else:
        highlight_atoms = [idx_map[i] for i in highlight_atoms if i in idx_map]
        highlight_atom_colors = None

    # Process highlight_bonds argument for highlightBonds and highlightBondColors
    # highlightBonds takes a list of bond indices
    # highlightBondColors takes a mapping from bond indices to colors
    if highlight_bonds is None:
        highlight_bonds = []
        highlight_bond_colors = None
    elif isinstance(highlight_bonds, dict):
        highlight_bond_colors = {
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b]).GetIdx(): tuple(v)
            for (i_a, i_b), v in highlight_bonds.items()
            if i_a in idx_map and i_b in idx_map
        }

        highlight_bonds = list(highlight_bond_colors.keys())
    else:
        highlight_bonds = [
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b]).GetIdx()
            for i_a, i_b in highlight_bonds
            if i_a in idx_map and i_b in idx_map
        ]
        highlight_bond_colors = None

    # Process bond_notes argument and place notes in the molecule
    if bond_notes is not None:
        for (i_a, i_b), note in bond_notes.items():
            if i_a not in idx_map or i_b not in idx_map:
                continue
            rdbond = rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            rdbond.SetProp("bondNote", unicode_to_ascii(str(note)))

    # Process atom_notes argument and place notes in the molecule
    if atom_notes is not None:
        for i, note in atom_notes.items():
            if i not in idx_map:
                continue
            rdatom = rdmol.GetAtomWithIdx(idx_map[i])
            rdatom.SetProp("atomNote", unicode_to_ascii(str(note)))

    # Resolve kekulization so it is the same for all drawn molecules
    rdkit.Chem.rdmolops.Kekulize(rdmol, clearAromaticFlags=True)

    # Compute 2D coordinates
    rdkit.Chem.rdDepictor.Compute2DCoords(rdmol)

    # Construct the drawing object and get a handle to its options
    drawer = rdkit.Chem.Draw.MolDraw2DSVG(width, height)
    draw_options = drawer.drawOptions()

    # Specify the scale to fit all atoms
    # This is important for emphasize_atoms
    coords_2d = next(rdmol.GetConformers()).GetPositions()[..., (0, 1)]
    drawer.SetScale(
        width,
        height,
        rdkit.Geometry.rdGeometry.Point2D(*(coords_2d.min(axis=0) - 1.0)),
        rdkit.Geometry.rdGeometry.Point2D(*(coords_2d.max(axis=0) + 1.0)),
    )

    # Set the colors used for each element according to the emphasize_atoms and
    # color_by_element arguments
    if emphasize_atoms:
        draw_options.setAtomPalette(
            {i: (0.8, 0.8, 0.8) for i in range(rdmol.GetNumAtoms())}
        )
    else:
        set_atom_palette(draw_options)

    # Draw the molecule
    # Note that if emphasize_atoms is used, this will be the un-emphasized parts
    # of the molecule
    drawer.DrawMolecule(
        rdmol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
    )

    # Draw an overlapping molecule for the emphasized atoms
    if emphasize_atoms:
        # Set the atom palette according to the color_by_element argument
        set_atom_palette(draw_options)

        # Create a copy of the molecule that removes atoms that aren't emphasized
        emphasized = rdkit.Chem.rdchem.RWMol(rdmol)
        emphasized.BeginBatchEdit()
        for i in set(idx_map) - set(emphasize_atoms):
            emphasized.RemoveAtom(idx_map[i])
        emphasized.CommitBatchEdit()

        # Draw the molecule. The scale has been fixed and we're re-using the
        # same coordinates, so this will overlap the background molecule
        drawer.DrawMolecule(emphasized)

    # Finalize the SVG
    drawer.FinishDrawing()

    # Return an SVG object that we can view in notebook
    svg_contents = drawer.GetDrawingText()
    return SVG(svg_contents)

from openff.units import ensure_quantity
import mdtraj
import numpy
import nglview
from pathlib import Path

def nglview_show_openmm(topology, positions, image_molecules=False):
    top = mdtraj.Topology.from_openmm(topology)
    
    if isinstance(positions, str) or isinstance(positions, Path):
        traj = mdtraj.load(positions, top=top)
        if image_molecules:
            traj.image_molecules(inplace=True)
    else:
        positions = ensure_quantity(positions, "openmm").value_in_unit(openmm.unit.nanometer)
        xyz = np.asarray([positions])
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            l1, l2, l3, alpha, beta, gamma = mdtraj.utils.box_vectors_to_lengths_and_angles(
                *np.asarray(box_vectors.value_in_unit(openmm.unit.nanometer))
            )
            unitcell_angles, unitcell_lengths = [alpha, beta, gamma], [l1, l2, l3]
        else:
            unitcell_angles, unitcell_lengths = None, None
        traj = mdtraj.Trajectory(xyz, top, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)
    return nglview.show_mdtraj(traj)
