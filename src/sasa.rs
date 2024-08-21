use freesasa_rs::{set_verbosity, structure::Structure, FreesasaVerbosity};
use pdbtbx::{Atom, Residue};
use std::process;

/// Represents an atom with additional solvent accessible surface area (SASA) information.
#[derive(Debug)]
pub struct ExtendedAtom {
    /// The underlying atom.
    pub atom: Atom,
    /// The solvent accessible surface area of the atom.
    pub sasa: f64,
}

impl ExtendedAtom {
    /// Creates a new `ExtendedAtom` instance.
    ///
    /// # Arguments
    ///
    /// * `atom` - The underlying `Atom` object.
    /// * `sasa` - The solvent accessible surface area of the atom.
    ///
    /// # Returns
    ///
    /// A new `ExtendedAtom` instance.
    pub fn new(atom: Atom, sasa: f64) -> Self {
        ExtendedAtom { atom, sasa }
    }
}

/// Represents a residue with extended SASA information.
pub struct ExtendedRes {
    /// The underlying residue.
    pub residue: Residue,
    /// The chain identifier of the residue.
    pub chain: String,
    /// The SASA of the backbone atoms.
    pub sasa_bb: f64,
    /// The SASA of the side chain atoms.
    pub sasa_sc: f64,
    /// The relative SASA of the backbone atoms.
    pub rel_sasa_bb: f64,
    /// The relative SASA of the side chain atoms.
    pub rel_sasa_sc: f64,
    /// The total relative SASA of the residue.
    pub rel_sasa_total: f64,
}

impl ExtendedRes {
    /// Creates a new `ExtendedRes` instance.
    ///
    /// # Arguments
    ///
    /// * `residue` - The underlying `Residue` object.
    /// * `chain` - The chain identifier of the residue.
    /// * `sasa_bb` - The SASA of the backbone atoms.
    /// * `sasa_sc` - The SASA of the side chain atoms.
    /// * `rel_sasa_bb` - The relative SASA of the backbone atoms.
    /// * `rel_sasa_sc` - The relative SASA of the side chain atoms.
    /// * `rel_sasa_total` - The total relative SASA of the residue.
    ///
    /// # Returns
    ///
    /// A new `ExtendedRes` instance.
    pub fn new(
        residue: Residue,
        chain: String,
        sasa_bb: f64,
        sasa_sc: f64,
        rel_sasa_bb: f64,
        rel_sasa_sc: f64,
        rel_sasa_total: f64,
    ) -> Self {
        ExtendedRes {
            residue,
            chain,
            sasa_bb,
            sasa_sc,
            rel_sasa_bb,
            rel_sasa_sc,
            rel_sasa_total,
        }
    }
}

/// Calculates the Solvent Accessible Surface Area (SASA) for a given PDB structure.
///
/// This function performs SASA calculations for each atom in the structure and aggregates
/// the results at the residue level. It also calculates relative SASA values based on
/// standard accessible surface areas for each residue type.
///
/// # Arguments
///
/// * `pdbtbx_struct` - A mutable reference to a `pdbtbx::PDB` structure.
///
/// # Returns
///
/// A `Vec<ExtendedRes>` containing SASA information for each residue in the structure.
///
/// # Panics
///
/// This function will panic if:
/// - The number of atoms in the pdbtbx structure doesn't match the number of atoms in the SASA calculation.
/// - A residue name is encountered that doesn't have a corresponding relative accessibility value.
///
/// # Steps
///
/// 1. Calculates SASA for each atom using the freesasa library.
/// 2. Creates `ExtendedAtom` instances for each atom, storing SASA values.
/// 3. Aggregates atom SASA values into `ExtendedRes` instances for each residue.
/// 4. Calculates relative SASA values for backbone, side chain, and total residue.
///
/// # Notes
///
/// - The function modifies the input `pdbtbx_struct` by setting B-factors to SASA values.
/// - SASA calculations are performed using the freesasa library with error-only verbosity.
/// - Relative SASA values are calculated based on standard accessible surface areas defined in `REL_ASA`.
///
pub fn calculate_sasa(mut pdbtbx_struct: pdbtbx::PDB) -> Vec<ExtendedRes> {
    // ================================================================================
    // Calculate the SASA of each atom
    set_verbosity(FreesasaVerbosity::Error);
    let structure = Structure::from_pdbtbx(&pdbtbx_struct).unwrap();

    let result = structure.calculate_sasa().unwrap();
    let atom_sasa = result.atom_sasa();

    let pdbtbx_atom_count = pdbtbx_struct.atoms().count();
    let structure_atom_count = atom_sasa.len();

    if pdbtbx_atom_count != structure_atom_count {
        println!("Warning: The number of atoms in the pdbtbx structure ({}) does not match the number of atoms in the structure ({})", pdbtbx_atom_count, structure_atom_count);
        panic!("Aborting");
    }

    // Create a vector of ExtendedAtoms containing the SASA of each atom
    let mut extended_atoms: Vec<ExtendedAtom> = Vec::new();

    // Use the Atom from pdbtbx as base for the ExtendedAtom
    // let mut rng = rand::thread_rng();
    // for atom in pdbtbx_struct.atoms() {

    // Create a vector as big as pdbtbx_struct.atoms()
    // let atom_sasa: Vec<f64> = (0..pdbtbx_atom_count).map(|_| rng.gen()).collect();

    for (atom, sasa) in pdbtbx_struct.atoms_mut().zip(atom_sasa.iter()) {
        // Create a random sasa value for development
        // let n: f64 = rng.gen();

        // Find the atom in the atom_sasa vector
        // let n = atom_sasa.iter().find(|(a, _)| a == atom).unwrap();
        // let mut atom = atom;
        let _ = atom.set_b_factor(*sasa);

        let extended_atom = ExtendedAtom::new(atom.clone(), *sasa);
        extended_atoms.push(extended_atom);
        // Add a new property to the `atom` and a value
    }

    // Loop over the extended atoms and figure out to which residue they belong to
    let mut extended_residues: Vec<ExtendedRes> = Vec::new();
    for chain in pdbtbx_struct.chains() {
        for residue in chain.residues() {
            let mut sasa_bb = 0.0;
            let mut sasa_sc = 0.0;
            for atom in residue.atoms() {
                // Find the atom in the extended atoms vector
                let extended_atom = extended_atoms.iter().find(|&x| x.atom == *atom).unwrap();
                // Add the SASA of the atom to the corresponding residue
                if atom.is_backbone() {
                    sasa_bb += extended_atom.sasa;
                } else {
                    sasa_sc += extended_atom.sasa;
                }
            }
            let extended_residue = ExtendedRes::new(
                residue.clone(),
                chain.id().to_string(),
                sasa_bb,
                sasa_sc,
                0.0,
                0.0,
                0.0,
            );
            extended_residues.push(extended_residue);
        }
    }

    // Calculate the relative SASA of each residue, based on `asa::REL_ASA`
    for (e_res, _res) in extended_residues
        .iter_mut()
        .zip(pdbtbx_struct.residues_mut())
    {
        let resname = e_res.residue.name().unwrap();

        // if let Some(rel_asa)  =
        match REL_ASA.iter().find(|(name, _)| *name == resname) {
            Some(rel_asa) => {
                let rel_sasa_bb = (e_res.sasa_bb / rel_asa.1.bb) * 100_f64;
                let rel_sasa_sc = (e_res.sasa_sc / rel_asa.1.sc) * 100_f64;
                let rel_total = (e_res.sasa_bb + e_res.sasa_sc / rel_asa.1.total) * 100_f64;

                e_res.rel_sasa_bb = rel_sasa_bb;
                e_res.rel_sasa_sc = rel_sasa_sc;
                e_res.rel_sasa_total = rel_total;
            }
            None => {
                eprintln!("\n### ERROR CALCULATING SASA ###");
                eprintln!(
                    "# No relative accessibility found for residue `{}`",
                    resname
                );
                process::exit(1);
            }
        }
    }
    extended_residues
}

/// Represents the Accessible Surface Area (ASA) values for a residue.
///
/// This struct holds the total ASA and its breakdown into backbone (bb) and side chain (sc) components.
pub struct AsaValues {
    /// The total Accessible Surface Area.
    pub total: f64,

    /// The Accessible Surface Area of the backbone atoms.
    pub bb: f64,

    /// The Accessible Surface Area of the side chain atoms.
    pub sc: f64,
}

pub const REL_ASA: &[(&str, AsaValues)] = &[
    (
        "ALA",
        AsaValues {
            total: 107.95,
            bb: 38.54,
            sc: 69.41,
        },
    ),
    (
        "CYS",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "ASP",
        AsaValues {
            total: 140.39,
            bb: 37.70,
            sc: 102.69,
        },
    ),
    (
        "GLU",
        AsaValues {
            total: 172.25,
            bb: 37.51,
            sc: 134.74,
        },
    ),
    (
        "PHE",
        AsaValues {
            total: 199.48,
            bb: 35.37,
            sc: 164.11,
        },
    ),
    (
        "GLY",
        AsaValues {
            total: 80.10,
            bb: 47.77,
            sc: 32.33,
        },
    ),
    (
        "HIS",
        AsaValues {
            total: 182.88,
            bb: 35.80,
            sc: 147.08,
        },
    ),
    (
        "ILE",
        AsaValues {
            total: 175.12,
            bb: 37.16,
            sc: 137.96,
        },
    ),
    (
        "LYS",
        AsaValues {
            total: 200.81,
            bb: 37.51,
            sc: 163.30,
        },
    ),
    (
        "LEU",
        AsaValues {
            total: 178.63,
            bb: 37.51,
            sc: 141.12,
        },
    ),
    (
        "MET",
        AsaValues {
            total: 194.15,
            bb: 37.51,
            sc: 156.64,
        },
    ),
    (
        "ASN",
        AsaValues {
            total: 143.94,
            bb: 37.70,
            sc: 106.24,
        },
    ),
    (
        "PRO",
        AsaValues {
            total: 136.13,
            bb: 16.23,
            sc: 119.90,
        },
    ),
    (
        "GLN",
        AsaValues {
            total: 178.50,
            bb: 37.51,
            sc: 140.99,
        },
    ),
    (
        "ARG",
        AsaValues {
            total: 238.76,
            bb: 37.51,
            sc: 201.25,
        },
    ),
    (
        "SER",
        AsaValues {
            total: 116.50,
            bb: 38.40,
            sc: 78.11,
        },
    ),
    (
        "THR",
        AsaValues {
            total: 139.27,
            bb: 37.57,
            sc: 101.70,
        },
    ),
    (
        "VAL",
        AsaValues {
            total: 151.44,
            bb: 37.16,
            sc: 114.28,
        },
    ),
    (
        "TRP",
        AsaValues {
            total: 249.36,
            bb: 38.10,
            sc: 211.26,
        },
    ),
    (
        "TYR",
        AsaValues {
            total: 212.76,
            bb: 35.38,
            sc: 177.38,
        },
    ),
    (
        "ASH",
        AsaValues {
            total: 140.39,
            bb: 37.70,
            sc: 102.69,
        },
    ),
    (
        "DDZ",
        AsaValues {
            total: 107.95,
            bb: 38.54,
            sc: 69.41,
        },
    ),
    (
        "GLH",
        AsaValues {
            total: 172.25,
            bb: 37.51,
            sc: 134.74,
        },
    ),
    (
        "CYM",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "CSP",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "CYF",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "CYC",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "CFE",
        AsaValues {
            total: 134.28,
            bb: 37.53,
            sc: 96.75,
        },
    ),
    (
        "NEP",
        AsaValues {
            total: 182.88,
            bb: 35.80,
            sc: 147.08,
        },
    ),
    (
        "ALY",
        AsaValues {
            total: 200.81,
            bb: 37.51,
            sc: 163.30,
        },
    ),
    (
        "MLZ",
        AsaValues {
            total: 200.81,
            bb: 37.51,
            sc: 163.30,
        },
    ),
    (
        "MLY",
        AsaValues {
            total: 200.81,
            bb: 37.51,
            sc: 163.30,
        },
    ),
    (
        "M3L",
        AsaValues {
            total: 200.81,
            bb: 37.51,
            sc: 163.30,
        },
    ),
    (
        "HYP",
        AsaValues {
            total: 136.13,
            bb: 16.23,
            sc: 119.90,
        },
    ),
    (
        "SEP",
        AsaValues {
            total: 116.50,
            bb: 38.40,
            sc: 78.11,
        },
    ),
    (
        "TOP",
        AsaValues {
            total: 139.27,
            bb: 37.57,
            sc: 101.70,
        },
    ),
    (
        "TYP",
        AsaValues {
            total: 212.76,
            bb: 35.38,
            sc: 177.38,
        },
    ),
    (
        "PTR",
        AsaValues {
            total: 212.76,
            bb: 35.38,
            sc: 177.38,
        },
    ),
    (
        "TYS",
        AsaValues {
            total: 212.76,
            bb: 35.38,
            sc: 177.38,
        },
    ),
    (
        "PNS",
        AsaValues {
            total: 116.50,
            bb: 38.40,
            sc: 78.11,
        },
    ),
    (
        "QSR",
        AsaValues {
            total: 390.53,
            bb: 26.64,
            sc: 363.88,
        },
    ),
];

#[cfg(test)]
mod test {

    use super::*;
    use std::env;

    // TODO: Add more tests

    #[test]
    fn test_calculate_sasa() {
        let pdb_path = env::current_dir().unwrap().join("tests/data/complex.pdb");
        let pdb = pdbtbx::open_pdb(pdb_path.to_str().unwrap(), pdbtbx::StrictnessLevel::Loose)
            .unwrap()
            .0;
        let extended_residues = calculate_sasa(pdb);

        assert_eq!(extended_residues.len(), 116);
    }
}
