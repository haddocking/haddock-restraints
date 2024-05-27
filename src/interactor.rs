use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize, Debug)]
pub struct Interactor {
    id: u16,
    chain: String,
    active: HashSet<i16>,
    passive: HashSet<i16>,
    target: HashSet<u16>,
}

impl Interactor {
    pub fn is_valid(&self) -> Result<bool, &str> {
        if self.target.is_empty() {
            return Err("Target residues are empty");
        }
        if self.active.intersection(&self.passive).next().is_some() {
            return Err("Active/Passive selections overlap");
        }
        Ok(true)
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

    pub fn create_block(&self, target_res: Vec<(&str, &i16)>) -> String {
        let mut block = String::new();
        for resnum in self.active() {
            let assign_line = format!(
                "assign ( resid {} and segid {} )\n       (\n",
                resnum,
                self.chain()
            );

            let res_lines: Vec<String> = target_res
                .iter()
                .enumerate()
                .map(|(index, res)| {
                    let mut res_line = format!("        ( resid {}  and segid {})\n", res.1, res.0);
                    if index != target_res.len() - 1 {
                        res_line.push_str("     or\n");
                    }
                    res_line
                })
                .collect();

            block.push_str(&assign_line);
            block.push_str(&res_lines.join(""));
            block.push_str("       )\n");
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
