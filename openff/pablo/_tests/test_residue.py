import pytest
from openff.toolkit import Molecule

from openff.pablo.chem import DISULFIDE_BOND, PEPTIDE_BOND
from openff.pablo.residue import AtomDefinition, BondDefinition, ResidueDefinition


class TestBondDefinition:
    def test_flipped(self):
        assert (
            BondDefinition.with_defaults("H1", "H2")
            == BondDefinition.with_defaults("H2", "H1").flipped()
        )


class TestResidueDefinition:
    @pytest.fixture
    def cys_mapped_smiles(self) -> str:
        return r"[N:1]([C@@:2]([C:5]([S:6][H:13])([H:11])[H:12])([C:3]([O:7][H:14])=[O:4])[H:10])([H:8])[H:9]"

    @pytest.fixture
    def cys_mapped_atomnames(self) -> dict[int, str]:
        return {
            1: "N",
            2: "CA",
            3: "C",
            4: "O",
            5: "CB",
            6: "SG",
            7: "OXT",
            8: "H",
            9: "H2",
            10: "HA",
            11: "HB2",
            12: "HB3",
            13: "HG",
            14: "HXT",
        }

    def test_from_molecule(
        self,
        cys_def: ResidueDefinition,
        cys_mapped_smiles: str,
        cys_mapped_atomnames: dict[int, str],
    ):
        cysteine = Molecule.from_mapped_smiles(cys_mapped_smiles)
        for i, atom in enumerate(cysteine.atoms, start=1):
            atom.name = cys_mapped_atomnames[i]
            if i in {7, 14, 13, 9}:
                atom.metadata["leaving_atom"] = True

        from_molecule = ResidueDefinition.from_molecule(
            cysteine,
            "CYS",
            linking_bond=PEPTIDE_BOND,
            crosslink=DISULFIDE_BOND,
            description="CYSTEINE",
        )
        assert from_molecule == cys_def

    def test_from_capped_molecule(
        self,
        cys_def: ResidueDefinition,
        cys_mapped_smiles: str,
        cys_mapped_atomnames: dict[int, str],
    ):
        cysteine = Molecule.from_mapped_smiles(cys_mapped_smiles)
        for i, atom in enumerate(cysteine.atoms, start=1):
            atom.name = cys_mapped_atomnames[i]

        from_molecule = ResidueDefinition.from_capped_molecule(
            cysteine,
            residue_name="CYS",
            leaving_atom_indices={6, 13, 12, 8},
            linking_bond=PEPTIDE_BOND,
            crosslink=DISULFIDE_BOND,
            description="CYSTEINE",
        )
        assert from_molecule == cys_def

    def test_from_smiles(
        self,
        cys_def: ResidueDefinition,
        cys_mapped_smiles: str,
        cys_mapped_atomnames: dict[int, str],
    ):
        cys_from_smiles = ResidueDefinition.from_smiles(
            residue_name="CYS",
            mapped_smiles=cys_mapped_smiles,
            atom_names=cys_mapped_atomnames,
            leaving_atoms={7, 14, 13, 9},
            crosslink=DISULFIDE_BOND,
            linking_bond=PEPTIDE_BOND,
            description="CYSTEINE",
        )
        assert cys_from_smiles == cys_def

    def test_from_smiles_raises_on_bad_atom_names(
        self,
        cys_def: ResidueDefinition,
        cys_mapped_smiles: str,
        cys_mapped_atomnames: dict[int, str],
    ):
        # Missing atom name
        atom_names = {k: v for k, v in cys_mapped_atomnames.items() if k != 5}
        with pytest.raises(ValueError):
            ResidueDefinition.from_smiles(
                residue_name="CYS",
                mapped_smiles=cys_mapped_smiles,
                atom_names=atom_names,
                leaving_atoms={7, 14, 13, 9},
                crosslink=DISULFIDE_BOND,
                linking_bond=PEPTIDE_BOND,
                description="CYSTEINE",
            )

    def test_leaving_atoms_must_be_associated_with_bond(
        self,
        cys_mapped_smiles: str,
        cys_mapped_atomnames: dict[int, str],
    ):
        kwargs = dict(
            residue_name="CCL",
            mapped_smiles=cys_mapped_smiles,
            atom_names=cys_mapped_atomnames,
            crosslink=DISULFIDE_BOND,
            linking_bond=PEPTIDE_BOND,
            description="CYSTEINE",
        )
        ResidueDefinition.from_smiles(**kwargs)  # type:ignore[arg]
        ResidueDefinition.from_smiles(**kwargs, leaving_atoms={7, 14, 13, 9})  # type:ignore[arg]
        with pytest.raises(ValueError):
            ResidueDefinition.from_smiles(**kwargs, leaving_atoms=[14])  # type:ignore[arg]

    def test_clashing_names_forbidden(self):
        with pytest.raises(ValueError):
            ResidueDefinition(
                residue_name="HHG",
                atoms=(
                    AtomDefinition.with_defaults("H", "H"),
                    AtomDefinition.with_defaults("H", "H"),
                ),
                bonds=(),
                linking_bond=None,
                crosslink=None,
                description="",
            )

    def test_to_openff_molecule(
        self,
        cys_def: ResidueDefinition,
        cys_mapped_smiles: str,
    ):
        cys_molecule = Molecule.from_mapped_smiles(cys_mapped_smiles)
        assert cys_def.to_openff_molecule() == cys_molecule

    def test_to_openff_molecule_roundtrip(
        self,
        hoh_def: ResidueDefinition,
    ):
        hoh_molecule = hoh_def.to_openff_molecule()
        hoh_def_roundtripped = ResidueDefinition.from_molecule(
            hoh_molecule,
            residue_name="HOH",
        )
        assert hoh_def_roundtripped == hoh_def

    def test_name_to_atom(self, hoh_def_with_synonyms: ResidueDefinition):
        assert (
            hoh_def_with_synonyms.name_to_atom["O"]
            == hoh_def_with_synonyms.name_to_atom["O1"]
        )
        assert (
            hoh_def_with_synonyms.name_to_atom["H1"]
            == hoh_def_with_synonyms.name_to_atom["HA"]
        )
        assert (
            hoh_def_with_synonyms.name_to_atom["H2"]
            == hoh_def_with_synonyms.name_to_atom["HB"]
        )
        assert hoh_def_with_synonyms.name_to_atom["O1"].symbol == "O"
        assert hoh_def_with_synonyms.name_to_atom["HA"].symbol == "H"
        assert hoh_def_with_synonyms.name_to_atom["HB"].symbol == "H"
        assert hoh_def_with_synonyms.name_to_atom["O1"].name == "O"
        assert hoh_def_with_synonyms.name_to_atom["HA"].name == "H1"
        assert hoh_def_with_synonyms.name_to_atom["HB"].name == "H2"

    def test_atom_bonded_to(self, hoh_def: ResidueDefinition):
        assert set(hoh_def.atoms_bonded_to("O")) == {"H1", "H2"}
        assert set(hoh_def.atoms_bonded_to("H1")) == {"O"}
        assert set(hoh_def.atoms_bonded_to("H2")) == {"O"}

    @pytest.mark.parametrize(
        ("bond_name", "leaving_atoms"),
        [
            ("prior_bond", {"H2"}),
            ("posterior_bond", {"OXT", "HXT"}),
            ("crosslink", {"HG"}),
        ],
    )
    def test_leaving_atom_props(
        self,
        cys_def: ResidueDefinition,
        bond_name: str,
        leaving_atoms: set[str],
    ):
        assert getattr(cys_def, bond_name + "_leaving_atoms") == leaving_atoms

    def test_prior_bond_linking_atom(
        self,
        cys_def: ResidueDefinition,
    ):
        assert cys_def.prior_bond_linking_atom == "N"

    def test_posterior_bond_linking_atom(
        self,
        cys_def: ResidueDefinition,
    ):
        assert cys_def.posterior_bond_linking_atom == "C"
