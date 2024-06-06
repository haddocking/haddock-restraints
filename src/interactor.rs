// use crate::sasa;
use crate::structure;
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize, Debug, Clone)]
pub struct Interactor {
    id: u16,
    chain: String,
    active: HashSet<i16>,
    pub passive: HashSet<i16>,
    target: HashSet<u16>,
    structure: Option<String>,
    passive_from_active: Option<bool>,
    surface_as_passive: Option<bool>,
}

#[allow(clippy::too_many_arguments)]
impl Interactor {
    pub fn new(id: u16) -> Self {
        Interactor {
            id,
            chain: String::new(),
            active: HashSet::new(),
            passive: HashSet::new(),
            target: HashSet::new(),
            structure: None,
            passive_from_active: None,
            surface_as_passive: None,
        }
    }

    pub fn is_valid(&self) -> Result<bool, &str> {
        if self.target.is_empty() {
            return Err("Target residues are empty");
        }
        if self.active.intersection(&self.passive).next().is_some() {
            return Err("Active/Passive selections overlap");
        }
        Ok(true)
    }

    pub fn set_passive_from_active(&mut self) {
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let residues = structure::get_residues(
                    &pdb,
                    self.active.iter().map(|x| *x as isize).collect(),
                );

                let neighbors = structure::neighbor_search(pdb.clone(), residues, 5.0);

                // Add these neighbors to the passive set
                neighbors.iter().for_each(|x| {
                    self.passive.insert(*x as i16);
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
    }

    pub fn set_surface_as_passive(&mut self) {
        // match pdbtbx::open_pdb(
        //     self.structure.clone().unwrap(),
        //     pdbtbx::StrictnessLevel::Loose,
        // ) {
        //     Ok((pdb, _warnings)) => {
        //         let sasa = sasa::calculate_sasa(pdb.clone());

        //         // Add these neighbors to the passive set
        //         sasa.iter().for_each(|r| {
        //             // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
        //             if r.rel_sasa_total > 0.7 && r.chain == self.chain {
        //                 self.passive.insert(r.residue.serial_number() as i16);
        //             }
        //         });
        //     }
        //     Err(e) => {
        //         panic!("Error opening PDB file: {:?}", e);
        //     }
        // }
    }

    pub fn id(&self) -> u16 {
        self.id
    }

    pub fn chain(&self) -> &str {
        &self.chain
    }

    pub fn active(&self) -> &HashSet<i16> {
        &self.active
    }

    pub fn passive(&self) -> &HashSet<i16> {
        &self.passive
    }

    pub fn target(&self) -> &HashSet<u16> {
        &self.target
    }

    pub fn structure(&self) -> &str {
        match &self.structure {
            Some(structure) => structure,
            None => "",
        }
    }

    pub fn set_structure(&mut self, structure: &str) {
        self.structure = Some(structure.to_string());
    }

    pub fn set_chain(&mut self, chain: &str) {
        self.chain = chain.to_string();
    }

    pub fn set_active(&mut self, active: Vec<i16>) {
        self.active = active.into_iter().collect();
    }

    pub fn set_passive(&mut self, passive: Vec<i16>) {
        self.passive = passive.into_iter().collect();
    }

    pub fn passive_from_active(&self) -> bool {
        self.passive_from_active.unwrap_or(false)
    }

    pub fn surface_as_passive(&self) -> bool {
        self.surface_as_passive.unwrap_or(false)
    }

    pub fn add_target(&mut self, target: u16) {
        self.target.insert(target);
    }

    pub fn create_block(&self, target_res: Vec<(&str, &i16)>) -> String {
        let mut block = String::new();
        let mut _active: Vec<i16> = self.active().iter().cloned().collect();
        _active.sort();

        // Sort the target residues by residue number
        let mut target_res: Vec<(&str, &i16)> = target_res.clone();
        target_res.sort_by(|a, b| a.1.cmp(b.1));

        for resnum in _active {
            let assign_line = format!(
                "assign ( resid {} and segid {} )\n       (\n",
                resnum,
                self.chain()
            );

            let res_lines: Vec<String> = target_res
                .iter()
                .enumerate()
                .map(|(index, res)| {
                    let mut res_line = format!("        ( resid {} and segid {} )\n", res.1, res.0);
                    if index != target_res.len() - 1 {
                        res_line.push_str("     or\n");
                    }
                    res_line
                })
                .collect();

            block.push_str(&assign_line);
            block.push_str(&res_lines.join(""));
            block.push_str("       ) 2.0 2.0 0.0\n\n");
        }
        block
    }
}

pub fn collect_resnums(interactors: Vec<&Interactor>) -> Vec<(&str, &i16)> {
    let mut resnums = Vec::<(&str, &i16)>::new();
    for interactor in interactors {
        resnums.extend(interactor.active().iter().map(|x| (interactor.chain(), x)));
        resnums.extend(interactor.passive().iter().map(|x| (interactor.chain(), x)));
    }
    resnums
}
