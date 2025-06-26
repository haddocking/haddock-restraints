pub fn handle_ti(input: &str, cutoff: &f64, pml: &Option<String>) {
    let _ = haddock_restraints::true_interface(input, cutoff, pml);
}

pub fn handle_gen_tbl(input: &str, pml: &Option<String>) {
    haddock_restraints::gen_tbl(input, pml)
}
pub fn handle_unambig_ti(input: &str, cutoff: &f64, pml: &Option<String>) {
    let _ = haddock_restraints::unambig_ti(input, cutoff, pml);
}
pub fn handle_restraint_bodies(input: &str, pml: &Option<String>) {
    let _ = haddock_restraints::restraint_bodies(input, pml);
}

pub fn handle_list_interface(input: &str, cutoff: &f64) {
    let _ = haddock_restraints::list_interface(input, cutoff);
}

pub fn handle_generate_z_restraints(
    input: &str,
    output: &str,
    residues: &[Vec<isize>],
    grid_size: &usize,
    grid_spacing: &f64,
) {
    let _ =
        haddock_restraints::generate_z_restraints(input, output, residues, grid_size, grid_spacing);
}
