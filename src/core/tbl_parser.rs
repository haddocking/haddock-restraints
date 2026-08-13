//! Parses HADDOCK `.tbl` AIR restraint files (the *inverse* of
//! `Interactor::create_block`, see `src/core/interactor.rs`).
//!
//! Only `resid <N> ... segid <chain>` selections are extracted from each
//! `assign` block — that's all pml rendering (`tbl2pml`) needs. Extra
//! clauses like `and name CA` or `and attr z gt 42.00` are ignored.

/// A single residue selection extracted from a `.tbl` restraint block.
#[derive(Debug, Clone, PartialEq)]
pub struct ResSelection {
    pub resid: i16,
    pub chain: String,
}

/// One `assign` block: one active residue restrained against a list of
/// (possibly OR'd) partner residues.
///
/// ```text
/// assign ( resid 1 and segid A )
///        (
///         ( resid 2 and segid B )
///      or
///         ( resid 3 and segid B )
///        ) 2.0 2.0 0.0
/// ```
/// parses into `active = {1, A}`, `partners = [{2, B}, {3, B}]`.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedRestraint {
    pub active: ResSelection,
    pub partners: Vec<ResSelection>,
}

/// Parses a HADDOCK `.tbl` restraints file into `ParsedRestraint`s.
///
/// `assign` blocks that don't resolve to any `resid`/`segid` selection
/// (dihedral, hbond, and other non-AIR CNS restraint types share the same
/// `assign` keyword but a different body) are skipped rather than treated
/// as errors, so one such block doesn't discard every restraint already
/// parsed from the rest of the file.
///
/// # Errors
/// Returns `Err` if a `resid`/`segid` clause is present but malformed (a
/// residue number that doesn't parse, or a `segid` with no chain after it),
/// or if the file contains no usable restraint at all.
pub fn parse_tbl(content: &str) -> Result<Vec<ParsedRestraint>, String> {
    let cleaned = strip_comments(content);
    let mut restraints = Vec::new();
    let mut saw_assign_block = false;

    for chunk in cleaned.split("assign").skip(1) {
        saw_assign_block = true;

        let mut selections = extract_selections(chunk)?.into_iter();

        let Some(active) = selections.next() else {
            continue;
        };

        restraints.push(ParsedRestraint {
            active,
            partners: selections.collect(),
        });
    }

    if saw_assign_block && restraints.is_empty() {
        return Err(
            "no resid/segid restraints found in .tbl (only unsupported restraint types?)"
                .to_string(),
        );
    }

    Ok(restraints)
}

/// Strips `!`-prefixed CNS comments (rest-of-line) from tbl content.
fn strip_comments(content: &str) -> String {
    content
        .lines()
        .map(|line| line.split('!').next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Splits `chunk` into whitespace-separated tokens, treating `(` and `)` as
/// their own single-character tokens regardless of surrounding whitespace
/// (so `(resid` and `A)` still tokenize as `(`, `resid` / `A`, `)`).
fn tokenize(chunk: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut chars = chunk.char_indices().peekable();

    while let Some(&(start, c)) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
            continue;
        }
        if c == '(' || c == ')' {
            tokens.push(&chunk[start..start + c.len_utf8()]);
            chars.next();
            continue;
        }
        let mut end = start + c.len_utf8();
        chars.next();
        while let Some(&(idx, c2)) = chars.peek() {
            if c2.is_whitespace() || c2 == '(' || c2 == ')' {
                break;
            }
            end = idx + c2.len_utf8();
            chars.next();
        }
        tokens.push(&chunk[start..end]);
    }

    tokens
}

/// Extracts every residue selection found in `chunk`'s innermost `( ... )`
/// groups — the paren-balanced units with no further parens nested inside
/// them, which is exactly what a single `resid`/`segid` selection is (the
/// OR-wrapping around a partner list is always one paren level up). Within
/// each such group, `resid` and `segid` are located independently, so
/// either order (`resid ... segid ...` or `segid ... resid ...`) and any
/// interleaved clause (`and name CA`, `and attr ...`) works the same way.
///
/// Groups that don't contain both a `resid` and a `segid` (dihedral/hbond
/// clauses, the outer OR-wrapper itself, ...) are skipped, not errored —
/// they're just not residue selections. A `resid`/`segid` keyword that
/// *is* present but missing its value is a genuine parse error, though.
fn extract_selections(chunk: &str) -> Result<Vec<ResSelection>, String> {
    let tokens = tokenize(chunk);

    // Stack of (index just after this group's `(`, has a nested `(` inside).
    let mut stack: Vec<(usize, bool)> = Vec::new();
    let mut selections = Vec::new();

    for (idx, &tok) in tokens.iter().enumerate() {
        match tok {
            "(" => {
                if let Some(parent) = stack.last_mut() {
                    parent.1 = true;
                }
                stack.push((idx + 1, false));
            }
            ")" => {
                if let Some((start, has_nested)) = stack.pop()
                    && !has_nested
                    && let Some(selection) = parse_selection(&tokens[start..idx])?
                {
                    selections.push(selection);
                }
            }
            _ => {}
        }
    }

    Ok(selections)
}

/// Parses a single innermost selection's tokens (with the surrounding
/// parens already stripped) into a `ResSelection`, if it has both a
/// `resid` and a `segid` clause.
fn parse_selection(tokens: &[&str]) -> Result<Option<ResSelection>, String> {
    let mut resid = None;
    let mut chain = None;

    let mut i = 0;
    while i < tokens.len() {
        match tokens[i] {
            "resid" => {
                let value = tokens
                    .get(i + 1)
                    .ok_or_else(|| "\"resid\" keyword with no residue number after it".to_string())?;
                resid = Some(
                    value
                        .parse::<i16>()
                        .map_err(|_| format!("invalid residue number \"{}\"", value))?,
                );
                i += 1;
            }
            "segid" => {
                let value = tokens
                    .get(i + 1)
                    .ok_or_else(|| "\"segid\" keyword with no chain identifier after it".to_string())?;
                chain = Some((*value).to_string());
                i += 1;
            }
            _ => {}
        }
        i += 1;
    }

    Ok(resid.zip(chain).map(|(resid, chain)| ResSelection { resid, chain }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_oneline() {
        let tbl = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(
            restraints,
            vec![ParsedRestraint {
                active: ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                partners: vec![ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                }],
            }]
        );
    }

    #[test]
    fn test_parse_multiline_or() {
        let tbl = "assign ( resid 1 and segid A )\n       (\n        ( resid 2 and segid B )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(
            restraints,
            vec![ParsedRestraint {
                active: ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                partners: vec![
                    ResSelection {
                        resid: 2,
                        chain: "B".to_string()
                    },
                    ResSelection {
                        resid: 3,
                        chain: "B".to_string()
                    },
                ],
            }]
        );
    }

    #[test]
    fn test_parse_ignores_atom_clauses() {
        let tbl = "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B and name CB ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
        assert_eq!(restraints[0].active.chain, "A");
        assert_eq!(restraints[0].partners[0].resid, 2);
        assert_eq!(restraints[0].partners[0].chain, "B");
    }

    #[test]
    fn test_parse_multiple_blocks() {
        let tbl = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n\
                   assign ( resid 5 and segid A ) ( resid 9 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints.len(), 2);
        assert_eq!(restraints[1].active.resid, 5);
    }

    #[test]
    fn test_parse_ignores_comments() {
        let tbl = "! this is a comment\nassign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0 ! inline comment\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
    }

    #[test]
    fn test_parse_atom_clause_before_segid() {
        // `and name CA` sits between `resid` and `segid` here.
        let tbl = "assign ( resid 1 and name CA and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
        assert_eq!(restraints[0].active.chain, "A");
    }

    #[test]
    fn test_parse_segid_before_resid() {
        let tbl = "assign ( segid A and resid 1 ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
        assert_eq!(restraints[0].active.chain, "A");
    }

    #[test]
    fn test_parse_skips_non_air_block_without_erroring() {
        // Not every `assign` block is a resid/segid AIR distance
        // restraint — some CNS restraint types select purely by `attr`.
        // One such block shouldn't discard the valid restraint that
        // follows it.
        let tbl = "assign ( attr store1 ) ( attr store2 ) 5.0 0.0 0.0\n\n\
                   assign ( resid 5 and segid A ) ( resid 9 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].active.resid, 5);
    }

    #[test]
    fn test_parse_invalid_resid_errors() {
        let tbl = "assign ( resid 99999999 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        assert!(parse_tbl(tbl).is_err());
    }

    #[test]
    fn test_parse_segid_without_chain_errors() {
        let tbl = "assign ( resid 1 and segid ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        assert!(parse_tbl(tbl).is_err());
    }

    #[test]
    fn test_parse_empty_file_yields_no_restraints() {
        assert!(parse_tbl("").unwrap().is_empty());
        assert!(parse_tbl("! just a comment\n").unwrap().is_empty());
    }

    #[test]
    fn test_parse_only_non_air_blocks_errors() {
        // Every assign block present resolves to zero selections, so
        // there's nothing to visualize and this should be reported rather
        // than silently returning an empty restraint list.
        let tbl = "assign ( attr store1 ) ( attr store2 ) 5.0 0.0 0.0\n\n";
        assert!(parse_tbl(tbl).is_err());
    }

    #[test]
    fn test_extract_selections_no_surrounding_whitespace() {
        let selections = extract_selections("(resid 1 and segid A)(resid 2 and segid B)").unwrap();
        assert_eq!(
            selections,
            vec![
                ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                },
            ]
        );
    }

    #[test]
    fn test_extract_selections() {
        let selections =
            extract_selections("( resid 1 and segid A ) ( resid 2 and segid B )").unwrap();
        assert_eq!(
            selections,
            vec![
                ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                },
            ]
        );
    }
}
