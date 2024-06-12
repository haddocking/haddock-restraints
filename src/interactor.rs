use crate::sasa;
use crate::structure;
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize, Debug, Clone)]
pub struct Interactor {
    id: u16,
    chain: String,
    active: HashSet<i16>,
    active_atoms: Option<Vec<String>>,
    pub passive: HashSet<i16>,
    passive_atoms: Option<Vec<String>>,
    target: HashSet<u16>,
    target_distance: Option<f64>,
    lower_margin: Option<f64>,
    upper_margin: Option<f64>,
    structure: Option<String>,
    passive_from_active: Option<bool>,
    surface_as_passive: Option<bool>,
    filter_buried: Option<bool>,
    filter_buried_cutoff: Option<f64>,
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
            filter_buried: None,
            filter_buried_cutoff: None,
            active_atoms: None,
            passive_atoms: None,
            target_distance: None,
            lower_margin: None,
            upper_margin: None,
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
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let sasa = sasa::calculate_sasa(pdb.clone());

                // Add these neighbors to the passive set
                sasa.iter().for_each(|r| {
                    // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
                    if r.rel_sasa_total > 0.7 && r.chain == self.chain {
                        self.passive.insert(r.residue.serial_number() as i16);
                    }
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
    }

    pub fn remove_buried_residues(&mut self) {
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let sasa = sasa::calculate_sasa(pdb.clone());

                let sasa_cutoff = self.filter_buried_cutoff.unwrap_or(0.7);

                sasa.iter().for_each(|r| {
                    // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
                    if r.rel_sasa_total < sasa_cutoff && r.chain == self.chain {
                        // This residue is not accessible, remove it from the passive and active sets
                        self.passive.remove(&(r.residue.serial_number() as i16));
                        self.active.remove(&(r.residue.serial_number() as i16));
                    }
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
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

    pub fn set_target_distance(&mut self, distance: f64) {
        self.target_distance = Some(distance);
    }

    pub fn set_lower_margin(&mut self, margin: f64) {
        self.lower_margin = Some(margin);
    }

    pub fn set_upper_margin(&mut self, margin: f64) {
        self.upper_margin = Some(margin);
    }

    pub fn passive_from_active(&self) -> bool {
        self.passive_from_active.unwrap_or(false)
    }

    pub fn surface_as_passive(&self) -> bool {
        self.surface_as_passive.unwrap_or(false)
    }

    pub fn filter_buried(&self) -> bool {
        self.filter_buried.unwrap_or(false)
    }

    pub fn add_target(&mut self, target: u16) {
        self.target.insert(target);
    }

    pub fn set_active_atoms(&mut self, atoms: Vec<String>) {
        self.active_atoms = Some(atoms);
    }

    pub fn set_passive_atoms(&mut self, atoms: Vec<String>) {
        self.passive_atoms = Some(atoms);
    }

    pub fn create_block(&self, target_res: Vec<(&str, &i16)>) -> String {
        let mut block = String::new();
        let mut _active: Vec<i16> = self.active().iter().cloned().collect();
        _active.sort();

        // Sort the target residues by residue number
        let mut target_res: Vec<(&str, &i16)> = target_res.clone();
        target_res.sort_by(|a, b| a.1.cmp(b.1));

        // Check if need to use multiline separation
        let multiline = target_res.len() > 1;

        for resnum in _active {
            let atom_str = format_atom_string(&self.active_atoms);
            let mut assign_str = format!(
                "assign ( resid {} and segid {}{} )",
                resnum,
                self.chain(),
                atom_str
            );

            if multiline {
                assign_str += "\n       (\n";
            }

            block.push_str(assign_str.as_str());

            let res_lines: Vec<String> = target_res
                .iter()
                .enumerate()
                .map(|(index, res)| {
                    let passive_atom_str = format_atom_string(&self.passive_atoms);

                    let mut res_line = String::new();
                    if multiline {
                        res_line.push_str(
                            format!(
                                "        ( resid {} and segid {}{} )\n",
                                res.1, res.0, passive_atom_str
                            )
                            .as_str(),
                        );
                    } else {
                        res_line.push_str(
                            format!(
                                " ( resid {} and segid {}{} )",
                                res.1, res.0, passive_atom_str
                            )
                            .as_str(),
                        );
                    }

                    if index != target_res.len() - 1 {
                        res_line.push_str("     or\n");
                    }
                    res_line
                })
                .collect();

            block.push_str(&res_lines.join(""));

            let distance_string = format_distance_string(
                &self.target_distance,
                &self.lower_margin,
                &self.upper_margin,
            );
            if multiline {
                block.push_str(format!("       ) {}\n\n", distance_string).as_str());
            } else {
                block.push_str(format!(" {}\n\n", distance_string).as_str())
            }
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

pub fn format_distance_string(
    target: &Option<f64>,
    lower: &Option<f64>,
    upper: &Option<f64>,
) -> String {
    let target = match target {
        Some(target) => target,
        None => &2.0,
    };

    let lower = match lower {
        Some(lower) => lower,
        None => &2.0,
    };

    let upper = match upper {
        Some(upper) => upper,
        None => &0.0,
    };

    format!("{:.1} {:.1} {:.1}", target, lower, upper)
}

pub fn format_atom_string(atoms: &Option<Vec<String>>) -> String {
    match atoms {
        Some(atoms) => {
            let atoms: Vec<String> = atoms.iter().map(|x| format!(" and name {}", x)).collect();
            atoms.join("")
        }
        None => "".to_string(),
    }
}

// pub fn format_oneline_assign_string()

#[cfg(test)]
mod tests {

    use crate::interactor::Interactor;

    #[test]
    fn test_create_block_multiline() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");

        let observed = interactor.create_block(vec![("B", &2), ("B", &3)]);

        let block = "assign ( resid 1 and segid A )\n       (\n        ( resid 2 and segid B )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_oneline() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");

        let observed = interactor.create_block(vec![("B", &2)]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_active_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);

        let observed = interactor.create_block(vec![("B", &2)]);

        let block =
            "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_passive_atoms(vec!["CA".to_string()]);

        let observed = interactor.create_block(vec![("B", &2)]);

        let block =
            "assign ( resid 1 and segid A ) ( resid 2 and segid B and name CA ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_active_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);
        interactor.set_passive_atoms(vec!["CB".to_string()]);

        let observed = interactor.create_block(vec![("B", &2)]);

        let block =
            "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B and name CB ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_multiline_block_active_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);
        interactor.set_passive_atoms(vec!["CB".to_string()]);

        let observed = interactor.create_block(vec![("B", &2), ("B", &3)]);

        let block = "assign ( resid 1 and segid A and name CA )\n       (\n        ( resid 2 and segid B and name CB )\n     or\n        ( resid 3 and segid B and name CB )\n       ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_with_distance() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_target_distance(5.0);
        interactor.set_lower_margin(0.0);

        let observed = interactor.create_block(vec![("B", &2)]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 5.0 0.0 0.0\n\n";

        assert_eq!(observed, block);
    }
}
