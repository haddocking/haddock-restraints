use std::collections::HashSet;

use crate::core::tbl_parser::parse_tbl;
use crate::core::utils;

/// Generates a PyMOL (`.pml`) visualization directly from a `.tbl` restraints
/// file and the PDB(s) it refers to — without needing a `config.json`.
///
/// Mirrors the format of `Air::gen_pml` (`src/core/air.rs`), but is built
/// straight from the parsed `.tbl` blocks instead of an `Interactor`/`Air`
/// graph, since pml rendering only ever needs `(resid, chain)` pairs.
pub fn tbl2pml(tbl_path: &str, pdb_paths: &[String], output: &str) -> Result<(), String> {
    let content =
        std::fs::read_to_string(tbl_path).map_err(|e| format!("Could not read {}: {}", tbl_path, e))?;

    let restraints = parse_tbl(&content)?;

    let mut pml = String::new();

    for pdb in pdb_paths {
        pml.push_str(format!("load {}\n", pdb).as_str());
    }

    pml.push_str("set label_size, 0\n");
    pml.push_str("set dash_gap, 0\n");
    pml.push_str("set dash_color, yellow\n");

    let mut active: HashSet<(i16, String)> = HashSet::new();
    let mut passive: HashSet<(i16, String)> = HashSet::new();

    for restraint in &restraints {
        active.insert((restraint.active.resid, restraint.active.chain.clone()));

        let active_sel = format!(
            "resi {} and (name CA or name C1') and chain {}",
            restraint.active.resid, restraint.active.chain
        );
        let identifier = format!("{}-{}", restraint.active.resid, restraint.active.chain);

        for partner in &restraint.partners {
            passive.insert((partner.resid, partner.chain.clone()));

            let partner_sel = format!(
                "resi {} and (name CA or name C1') and chain {}",
                partner.resid, partner.chain
            );

            pml.push_str(
                format!(
                    "distance {}, ({}), ({})\n",
                    identifier, active_sel, partner_sel
                )
                .as_str(),
            );
        }
    }

    pml.push_str("color white\n");
    passive.iter().for_each(|(resnum, chain)| {
        pml.push_str(format!("color green, (resi {} and chain {})\n", resnum, chain).as_str())
    });
    active.iter().for_each(|(resnum, chain)| {
        pml.push_str(format!("color red, (resi {} and chain {})\n", resnum, chain).as_str())
    });

    utils::write_string_to_file(&pml, output).map_err(|e| format!("Could not write {}: {}", output, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tbl2pml_writes_expected_pml() {
        let tbl_content = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let tbl_path = "test_tbl2pml_input.tbl";
        let pml_path = "test_tbl2pml_output.pml";

        std::fs::write(tbl_path, tbl_content).unwrap();

        let result = tbl2pml(tbl_path, &["complex.pdb".to_string()], pml_path);
        assert!(result.is_ok());

        let pml = std::fs::read_to_string(pml_path).unwrap();

        assert!(pml.contains("load complex.pdb\n"));
        assert!(pml.contains(
            "distance 1-A, (resi 1 and (name CA or name C1') and chain A), (resi 2 and (name CA or name C1') and chain B)\n"
        ));
        assert!(pml.contains("color green, (resi 2 and chain B)\n"));
        assert!(pml.contains("color red, (resi 1 and chain A)\n"));

        std::fs::remove_file(tbl_path).unwrap();
        std::fs::remove_file(pml_path).unwrap();
    }

    #[test]
    fn test_tbl2pml_missing_tbl_file_errors() {
        let result = tbl2pml("does_not_exist.tbl", &["complex.pdb".to_string()], "out.pml");
        assert!(result.is_err());
    }
}
