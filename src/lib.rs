//! Main library crate for haddock-restraints

// Internal module organization
mod core;

// Public API exports
pub use core::air::*;
pub use core::commands::{
    restraint::restraint_bodies,
    tbl::gen_tbl,
    ti::{list_interface, true_interface, unambig_ti},
    z::generate_z_restraints,
};
pub use core::input::*;
pub use core::interactor::*;
pub use core::sasa::*;
pub use core::structure::*;
pub use core::utils::*;

// Optional: Add any library-level functionality here
