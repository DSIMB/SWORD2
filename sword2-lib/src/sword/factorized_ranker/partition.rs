#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Segment {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Domain {
    pub segments: Vec<Segment>,
    pub residues: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ParsedPartition {
    pub domains: Vec<Domain>,
    pub residue_to_domain: Vec<usize>,
    pub canonical: String,
}

#[allow(dead_code)]
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub(crate) enum FeatureError {
    #[error("empty or malformed delineation")]
    Malformed,
    #[error("segment {start}-{end} exceeds chain length {chain_len}")]
    OutOfRange {
        start: usize,
        end: usize,
        chain_len: usize,
    },
    #[error("residue {0} occurs in more than one domain")]
    Overlap(usize),
    #[error("partition does not cover every residue")]
    IncompleteCoverage,
    #[error("declared domain count does not match parsed partition")]
    DomainCountMismatch,
    #[error("required structural evidence is unavailable: {0}")]
    MissingContext(&'static str),
    #[error("non-finite feature {0}")]
    NonFinite(&'static str),
    #[error("feature schema or vector length mismatch")]
    SchemaMismatch,
}

pub(crate) fn parse_partition(
    text: &str,
    chain_len: usize,
) -> Result<ParsedPartition, FeatureError> {
    if chain_len == 0 || text.trim().is_empty() {
        return Err(FeatureError::Malformed);
    }

    let mut domains = Vec::new();
    for raw_domain in text.split_whitespace() {
        let mut segments = Vec::new();
        for raw_segment in raw_domain.split(';') {
            let (start, end) = parse_segment(raw_segment)?;
            if end >= chain_len {
                return Err(FeatureError::OutOfRange {
                    start,
                    end,
                    chain_len,
                });
            }
            segments.push(Segment { start, end });
        }
        if segments.is_empty() {
            return Err(FeatureError::Malformed);
        }
        segments.sort();
        domains.push(Domain {
            segments,
            residues: Vec::new(),
        });
    }
    domains.sort_by_key(|domain| domain.segments[0].start);

    let mut residue_to_domain = vec![usize::MAX; chain_len];
    for (domain_index, domain) in domains.iter_mut().enumerate() {
        for segment in &domain.segments {
            for residue in segment.start..=segment.end {
                if residue_to_domain[residue] != usize::MAX {
                    return Err(FeatureError::Overlap(residue));
                }
                residue_to_domain[residue] = domain_index;
                domain.residues.push(residue);
            }
        }
    }
    if residue_to_domain.contains(&usize::MAX) {
        return Err(FeatureError::IncompleteCoverage);
    }

    for domain in &mut domains {
        domain.segments = coalesce(&domain.residues);
    }
    let canonical = domains
        .iter()
        .map(|domain| {
            domain
                .segments
                .iter()
                .map(render_segment)
                .collect::<Vec<_>>()
                .join(";")
        })
        .collect::<Vec<_>>()
        .join(" ");

    Ok(ParsedPartition {
        domains,
        residue_to_domain,
        canonical,
    })
}

fn parse_segment(raw: &str) -> Result<(usize, usize), FeatureError> {
    let pieces: Vec<&str> = raw.split('-').collect();
    if pieces.is_empty() || pieces.len() > 2 || pieces.iter().any(|piece| piece.is_empty()) {
        return Err(FeatureError::Malformed);
    }
    let start = parse_ascii_usize(pieces[0])?;
    let end = if pieces.len() == 2 {
        parse_ascii_usize(pieces[1])?
    } else {
        start
    };
    if end < start {
        return Err(FeatureError::Malformed);
    }
    Ok((start, end))
}

fn parse_ascii_usize(raw: &str) -> Result<usize, FeatureError> {
    if raw.is_empty() || !raw.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(FeatureError::Malformed);
    }
    raw.parse().map_err(|_| FeatureError::Malformed)
}

fn coalesce(residues: &[usize]) -> Vec<Segment> {
    let mut segments = Vec::new();
    let Some((&first, rest)) = residues.split_first() else {
        return segments;
    };
    let mut start = first;
    let mut end = first;
    for &residue in rest {
        if residue == end + 1 {
            end = residue;
        } else {
            segments.push(Segment { start, end });
            start = residue;
            end = residue;
        }
    }
    segments.push(Segment { start, end });
    segments
}

fn render_segment(segment: &Segment) -> String {
    if segment.start == segment.end {
        segment.start.to_string()
    } else {
        format!("{}-{}", segment.start, segment.end)
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_partition, FeatureError};

    #[test]
    fn rejects_overlap_gap_and_out_of_range() {
        assert!(matches!(
            parse_partition("0-3 3-5", 6),
            Err(FeatureError::Overlap(3))
        ));
        assert!(matches!(
            parse_partition("0-2 4-5", 6),
            Err(FeatureError::IncompleteCoverage)
        ));
        assert!(matches!(
            parse_partition("0-2 3-6", 6),
            Err(FeatureError::OutOfRange { .. })
        ));
    }

    #[test]
    fn canonicalizes_domain_and_segment_order() {
        let p = parse_partition("4-5 3;0-1 2", 6).unwrap();
        assert_eq!(p.canonical, "0-1;3 2 4-5");
        assert_eq!(p.residue_to_domain, vec![0, 0, 1, 0, 2, 2]);
    }
}
