use std::collections::HashSet;
use std::fmt::Write as _;
use std::path::Path;

use crate::core::pml;
use crate::core::tbl_parser::parse_tbl;

/// Generates a PyMOL (`.pml`) visualization directly from a `.tbl` restraints
/// file and the PDB(s) it refers to — without needing a `config.json`.
///
/// Mirrors the format of `Air::gen_pml` (`src/core/air.rs`), reusing its
/// pml-string helpers (`crate::core::pml`), but is built straight from the
/// parsed `.tbl` blocks instead of an `Interactor`/`Air` graph, since pml
/// rendering only ever needs `(resid, chain)` pairs.
pub fn tbl2pml(tbl_path: &str, pdb_paths: &[String], output: &str) -> Result<(), String> {
    let content = std::fs::read_to_string(tbl_path)
        .map_err(|e| format!("Could not read {}: {}", tbl_path, e))?;

    let restraints = parse_tbl(&content)?;

    for pdb in pdb_paths {
        if !Path::new(pdb).is_file() {
            return Err(format!("PDB file not found: {}", pdb));
        }
    }

    let mut pml_script = String::new();

    for pdb in pdb_paths {
        let _ = writeln!(pml_script, "load \"{}\"", pdb);
    }

    pml_script.push_str(&pml::header());

    let mut active: HashSet<(i16, String)> = HashSet::new();
    let mut passive: HashSet<(i16, String)> = HashSet::new();

    for restraint in &restraints {
        // Mirrors `Air::gen_pml` skipping interactors with no partners:
        // an active residue with nothing restrained to it isn't a
        // restraint worth drawing or coloring.
        if restraint.partners.is_empty() {
            continue;
        }

        active.insert((restraint.active.resid, restraint.active.chain.clone()));

        let active_sel = pml::atom_selector(restraint.active.resid, &restraint.active.chain);
        let identifier = format!("{}-{}", restraint.active.resid, restraint.active.chain);

        for partner in &restraint.partners {
            passive.insert((partner.resid, partner.chain.clone()));

            let partner_sel = pml::atom_selector(partner.resid, &partner.chain);

            let _ = writeln!(
                pml_script,
                "distance {}, ({}), ({})",
                identifier, active_sel, partner_sel
            );
        }
    }

    pml_script.push_str(&pml::color_footer(
        passive.iter().map(|(r, c)| (*r, c.as_str())),
        active.iter().map(|(r, c)| (*r, c.as_str())),
    ));

    std::fs::write(output, &pml_script).map_err(|e| format!("Could not write {}: {}", output, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_test_pdb(path: &str) {
        std::fs::write(
            path,
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        )
        .unwrap();
    }

    #[test]
    fn test_tbl2pml_writes_expected_pml() {
        let tbl_content = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let tbl_path = "test_tbl2pml_input.tbl";
        let pdb_path = "test_tbl2pml_input.pdb";
        let pml_path = "test_tbl2pml_output.pml";

        std::fs::write(tbl_path, tbl_content).unwrap();
        write_test_pdb(pdb_path);

        let result = tbl2pml(tbl_path, &[pdb_path.to_string()], pml_path);
        assert!(result.is_ok());

        let pml = std::fs::read_to_string(pml_path).unwrap();

        assert!(pml.contains(&format!("load \"{}\"\n", pdb_path)));
        assert!(pml.contains(
            "distance 1-A, (resi 1 and (name CA or name C1') and chain A), (resi 2 and (name CA or name C1') and chain B)\n"
        ));
        assert!(pml.contains("color green, (resi 2 and chain B)\n"));
        assert!(pml.contains("color red, (resi 1 and chain A)\n"));

        std::fs::remove_file(tbl_path).unwrap();
        std::fs::remove_file(pdb_path).unwrap();
        std::fs::remove_file(pml_path).unwrap();
    }

    #[test]
    fn test_tbl2pml_skips_partnerless_restraint() {
        let tbl_content = "assign ( resid 1 and segid A ) 2.0 2.0 0.0\n\n";
        let tbl_path = "test_tbl2pml_partnerless.tbl";
        let pdb_path = "test_tbl2pml_partnerless.pdb";
        let pml_path = "test_tbl2pml_partnerless.pml";

        std::fs::write(tbl_path, tbl_content).unwrap();
        write_test_pdb(pdb_path);

        let result = tbl2pml(tbl_path, &[pdb_path.to_string()], pml_path);
        assert!(result.is_ok());

        let pml = std::fs::read_to_string(pml_path).unwrap();

        assert!(!pml.contains("distance"));
        assert!(!pml.contains("color red"));

        std::fs::remove_file(tbl_path).unwrap();
        std::fs::remove_file(pdb_path).unwrap();
        std::fs::remove_file(pml_path).unwrap();
    }

    #[test]
    fn test_tbl2pml_missing_tbl_file_errors() {
        let result = tbl2pml(
            "does_not_exist.tbl",
            &["complex.pdb".to_string()],
            "out.pml",
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_tbl2pml_missing_pdb_file_errors() {
        let tbl_content = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let tbl_path = "test_tbl2pml_missing_pdb.tbl";
        std::fs::write(tbl_path, tbl_content).unwrap();

        let result = tbl2pml(
            tbl_path,
            &["does_not_exist.pdb".to_string()],
            "test_tbl2pml_missing_pdb.pml",
        );
        assert!(result.is_err());
        assert!(!Path::new("test_tbl2pml_missing_pdb.pml").exists());

        std::fs::remove_file(tbl_path).unwrap();
    }

    #[test]
    fn test_tbl2pml_unwritable_output_errors_without_panicking() {
        let tbl_content = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let tbl_path = "test_tbl2pml_unwritable_input.tbl";
        let pdb_path = "test_tbl2pml_unwritable_input.pdb";
        std::fs::write(tbl_path, tbl_content).unwrap();
        write_test_pdb(pdb_path);

        // Directory doesn't exist, so the write should fail as `Err`
        // rather than panicking.
        let result = tbl2pml(tbl_path, &[pdb_path.to_string()], "nonexistent_dir/out.pml");
        assert!(result.is_err());

        std::fs::remove_file(tbl_path).unwrap();
        std::fs::remove_file(pdb_path).unwrap();
    }
}
