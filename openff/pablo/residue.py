"""
Classes for defining custom residues.
"""

from copy import deepcopy
from dataclasses import InitVar, dataclass
from functools import cached_property
from typing import (
    Collection,
    Iterator,
    Literal,
    Mapping,
    Self,
)

from openff.toolkit import Molecule
from openff.units import elements

__all__ = [
    "AtomDefinition",
    "BondDefinition",
    "ResidueDefinition",
]


@dataclass(frozen=True)
class AtomDefinition:
    """
    Description of an atom in a residue from the Chemical Component Dictionary (CCD).
    """

    name: str
    synonyms: tuple[str, ...]
    symbol: str
    leaving: bool
    charge: int
    aromatic: bool
    stereo: Literal["S", "R"] | None


@dataclass(frozen=True)
class BondDefinition:
    """
    Description of a bond in a residue from the Chemical Component Dictionary (CCD).
    """

    atom1: str
    atom2: str
    order: int
    aromatic: bool
    stereo: Literal["E", "Z"] | None


@dataclass(frozen=True)
class ResidueDefinition:
    """
    Description of a residue from the Chemical Component Dictionary (CCD).
    """

    residue_name: str
    parent_residue_name: str | None
    description: str
    linking_bond: BondDefinition | None
    atoms: tuple[AtomDefinition, ...]
    bonds: tuple[BondDefinition, ...]
    _skip_post_init_validation: InitVar[bool] = False

    def __post_init__(self, _skip_post_init_validation: bool):
        if _skip_post_init_validation:
            return

        self._validate()

    def _validate(self):
        if self.linking_bond is None and True in {atom.leaving for atom in self.atoms}:
            raise ValueError(
                f"{self.residue_name}: Leaving atoms were specified, but there is no linking bond",
                self,
            )
        if len(set(atom.name for atom in self.atoms)) != len(self.atoms):
            raise ValueError(
                f"{self.residue_name}: All atoms must have unique canonical names"
            )

        all_leaving_atoms = {atom.name for atom in self.atoms if atom.leaving}
        assigned_leaving_atoms = self.prior_bond_leaving_atoms.union(
            self.posterior_bond_leaving_atoms
        )
        unassigned_leaving_atoms = all_leaving_atoms.difference(assigned_leaving_atoms)
        if len(unassigned_leaving_atoms) != 0:
            raise ValueError(
                f"{self.residue_name}: Leaving atoms could not be assigned to a bond: {unassigned_leaving_atoms}"
            )

    @classmethod
    def from_molecule(
        cls,
        name: str,
        molecule: Molecule,
        linking_bond: BondDefinition | None = None,
        description: str = "",
        parent_residue_name: str | None = None,
    ) -> Self:
        atoms: list[AtomDefinition] = []
        for i, atom in enumerate(molecule.atoms):
            atoms.append(
                AtomDefinition(
                    name=atom.name,
                    synonyms=(),
                    symbol=atom.symbol,
                    leaving=bool(atom.metadata.get("leaving_atom")),
                    charge=atom.formal_charge,
                    stereo=atom.stereochemistry,
                    aromatic=atom.is_aromatic,
                )
            )
        bonds: list[BondDefinition] = []
        for bond in molecule.bonds:
            bonds.append(
                BondDefinition(
                    atom1=bond.atom1.name,
                    atom2=bond.atom2.name,
                    order=bond.bond_order,
                    aromatic=bond.is_aromatic,
                    stereo=bond.stereochemistry,
                )
            )

        return cls(
            residue_name=name,
            parent_residue_name=parent_residue_name,
            description=description,
            linking_bond=linking_bond,
            atoms=tuple(atoms),
            bonds=tuple(bonds),
        )

    @classmethod
    def from_capped_molecule(
        cls,
        name: str,
        molecule: Molecule,
        leaving_atom_indices: Collection[int],
        linking_bond: BondDefinition,
        description: str = "",
    ) -> Self:
        molecule = deepcopy(molecule)
        for i in leaving_atom_indices:
            molecule.atom(i).metadata["leaving_atom"] = True
        return cls.from_molecule(
            name=name,
            molecule=molecule,
            linking_bond=linking_bond,
            description=description,
        )

    @classmethod
    def from_smiles(
        cls,
        name: str,
        mapped_smiles: str,
        atom_names: Mapping[int, str],
        leaving_atoms: Collection[int] = (),
        linking_bond: BondDefinition | None = None,
        description: str = "",
    ) -> Self:
        molecule = Molecule.from_mapped_smiles(mapped_smiles)
        leaving_atom_indices = set(leaving_atoms)
        for i, atom in enumerate(molecule.atoms, start=1):
            if i in leaving_atom_indices:
                atom.metadata["leaving_atom"] = True
            atom.name = atom_names[i]

        return cls.from_molecule(
            name=name,
            molecule=molecule,
            linking_bond=linking_bond,
            description=description,
        )

    def to_openff_molecule(self) -> Molecule:
        molecule = Molecule()
        atoms = {}
        for atom in self.atoms:
            atoms[atom.name] = molecule.add_atom(
                atomic_number=elements.NUMBERS[atom.symbol],
                formal_charge=atom.charge,
                is_aromatic=atom.aromatic,
                stereochemistry=atom.stereo,
                name=atom.name,
                metadata={
                    "residue_name": self.residue_name,
                    "leaving": atom.leaving,
                },
            )

        for bond in self.bonds:
            molecule.add_bond(
                atom1=atoms[bond.atom1],
                atom2=atoms[bond.atom2],
                bond_order=bond.order,
                is_aromatic=bond.aromatic,
                stereochemistry=bond.stereo,
            )

        molecule.properties.update(
            {
                "linking_bond": self.linking_bond,
            }
        )

        return molecule

    @cached_property
    def name_to_atom(self) -> dict[str, AtomDefinition]:
        """Map from each atoms' name and synonyms to the value of a field."""
        mapping = {atom.name: atom for atom in self.atoms}
        canonical_names = set(mapping)
        for atom in self.atoms:
            for synonym in atom.synonyms:
                if synonym in mapping and mapping[synonym] != atom:
                    raise ValueError(
                        f"synonym {synonym} degenerately defined for canonical"
                        + f" names {mapping[synonym]} and {atom.name} in"
                        + f" residue {self.residue_name}"
                    )
                if synonym in canonical_names:
                    raise ValueError(
                        f"synonym {synonym} of atom {atom.name} clashes with"
                        + f" another canonical name in residue {self.residue_name}"
                    )
                mapping[synonym] = atom
        return mapping

    def atoms_bonded_to(self, atom_name: str) -> Iterator[str]:
        for bond in self.bonds:
            if bond.atom1 == atom_name:
                yield bond.atom2
            if bond.atom2 == atom_name:
                yield bond.atom1

    def _leaving_fragment_of(self, linking_atom: str) -> Iterator[str]:
        atoms_to_check = list(self.atoms_bonded_to(linking_atom))
        checked_atoms: set[str] = set()
        while atoms_to_check:
            atom_name = atoms_to_check.pop()
            if self.name_to_atom[atom_name].leaving:
                yield atom_name
                atoms_to_check.extend(
                    filter(
                        lambda x: x not in checked_atoms,
                        self.atoms_bonded_to(atom_name),
                    )
                )
            checked_atoms.add(atom_name)

    @cached_property
    def posterior_bond_leaving_atoms(self) -> set[str]:
        return (
            set()
            if self.linking_bond is None
            else set(self._leaving_fragment_of(self.posterior_bond_linking_atom))
        )

    @cached_property
    def prior_bond_leaving_atoms(self) -> set[str]:
        return (
            set()
            if self.linking_bond is None
            else set(self._leaving_fragment_of(self.prior_bond_linking_atom))
        )

    @property
    def prior_bond_linking_atom(self) -> str:
        if self.linking_bond is None:
            raise ValueError("not a linking residue")
        return self.linking_bond.atom2

    @property
    def posterior_bond_linking_atom(self) -> str:
        if self.linking_bond is None:
            raise ValueError("not a linking residue")
        return self.linking_bond.atom1
