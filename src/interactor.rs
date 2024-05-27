use std::collections::HashSet;

#[derive(Debug)]
pub struct Interactor {
    id: u16,
    chain: String,
    active: HashSet<i16>,
    passive: HashSet<i16>,
    target: Vec<u16>,
}

impl Interactor {
    pub fn new(
        id: u16,
        chain: &str,
        active: Vec<i16>,
        passive: Vec<i16>,
        target: Vec<u16>,
    ) -> Self {
        Interactor {
            id,
            chain: chain.to_string(),
            active: active.into_iter().collect(),
            passive: passive.into_iter().collect(),
            target,
        }
    }

    pub fn is_valid(&self) -> Result<bool, &str> {
        if self.active.is_empty() {
            return Err("Active residues are empty");
        }
        if self.passive.is_empty() {
            return Err("Passive residues are empty");
        }
        if self.target.is_empty() {
            return Err("Target residues are empty");
        }
        if self.active.intersection(&self.passive).next().is_some() {
            return Err("Active and passive selections contain one or more of the same residues");
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

    pub fn target(&self) -> &Vec<u16> {
        &self.target
    }
}
