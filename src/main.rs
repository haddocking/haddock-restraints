mod air;
mod interactor;
mod io;
mod pdb;
// mod sasa;

use air::Air;
use interactor::Interactor;
use pdbtbx::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let restraints_file = io::read_args()?;

    let interactors = io::read_json_file(&restraints_file).unwrap();

    let air = Air::new(interactors);

    let _tbl = air.gen_tbl()?;

    let (pdbtbx_struct, _errors) = pdbtbx::open(
        "/home/rodrigo/repos/haddock-restraints/src/data/1crn.pdb",
        StrictnessLevel::Loose,
    )
    .unwrap();

    // Find what are the neighbours of the active residues

    let target_res = pdb::get_residues(&pdbtbx_struct, vec![36]);
    let neighbors = pdb::neighbor_search(pdbtbx_struct.clone(), target_res, 5.0);

    println!("{:?}", neighbors);

    // let res_w_sasa = sasa::calculate_sasa(pdbtbx_struct);

    Ok(())
}
