//! Seed‑and‑extend local alignments (EMBOSS `supermatcher`‑like).
//!
//! This implements a pragmatic approximation of EMBOSS **supermatcher**: first
//! locate **word matches** of length `wordlen` between the two sequences, then
//! **extend** around each seed by running a local Smith–Waterman (our `matcher`
//! engine with `alternatives=1`) in a window surrounding the seed’s diagonal.
//!
//! Compared with exhaustive SW on full sequences, this is much faster on long
//! sequences, while typically recovering the same dominant local matches.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{supermatcher, SuperMatcherParams, WaterMatrix};
//! let p = SuperMatcherParams{ wordlen: 11, ..Default::default() };
//! let hits = supermatcher("ACGTTTGCA", "TTACGTTTGCAT", &p).unwrap();
//! assert!(hits.len() >= 1);
//! ```
//!
use std::collections::{HashMap, HashSet};
use crate::common::{EmbossersError, WaterMatrix};
use crate::matcher::{matcher, MatcherParams, MatcherHit};

/// Parameters for [`supermatcher`].
#[derive(Clone, Debug)]
pub struct SuperMatcherParams {
    /// Scoring matrix (DNA/BLOSUM62) used by the local extension step.
    pub matrix: WaterMatrix,
    /// Gap open used by the extension step.
    pub gap_open: f32,
    /// Gap extend used by the extension step.
    pub gap_extend: f32,
    /// Integer scale (preserve fractional penalties).
    pub scale: f32,
    /// Word length for seeding (k‑mer size).
    pub wordlen: usize,
    /// Half‑window in residues extracted around a seed for local SW.
    pub flank: usize,
    /// Minimum alignment score to keep.
    pub min_score: i32,
    /// Maximum number of seed windows to extend (to bound runtime).
    pub max_seeds: usize,
}

impl Default for SuperMatcherParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            wordlen: 11,
            flank: 200,
            min_score: 1,
            max_seeds: 5000,
        }
    }
}

/// Run seed‑and‑extend local alignment for a pair of sequences.
pub fn supermatcher(a: &str, b: &str, params: &SuperMatcherParams) -> Result<Vec<MatcherHit>, EmbossersError> {
    if a.len() < params.wordlen || b.len() < params.wordlen {
        return Ok(Vec::new());
    }
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();

    // Build k‑mer index on the shorter sequence to minimize memory.
    let (s_chars, l_chars, s_is_a) = if a_chars.len() <= b_chars.len() {
        (a_chars.clone(), b_chars.clone(), true)
    } else {
        (b_chars.clone(), a_chars.clone(), false)
    };
    let mut index: HashMap<String, Vec<usize>> = HashMap::new();
    for i in 0..=s_chars.len() - params.wordlen {
        let k: String = s_chars[i..i+params.wordlen].iter().collect();
        index.entry(k).or_default().push(i);
    }

    // Find seeds on the longer sequence and map to diagonals.
    let mut diagonals: HashSet<isize> = HashSet::new();
    let mut seeds: Vec<(usize, usize)> = Vec::new(); // (i_long, i_short)
    for j in 0..=l_chars.len() - params.wordlen {
        let k: String = l_chars[j..j+params.wordlen].iter().collect();
        if let Some(pos_list) = index.get(&k) {
            for &i in pos_list {
                // diagonal d = j - i (in long-short coordinates)
                let d = (j as isize) - (i as isize);
                if diagonals.insert(d) {
                    seeds.push((j, i));
                    if seeds.len() >= params.max_seeds { break; }
                }
            }
        }
        if seeds.len() >= params.max_seeds { break; }
    }

    // Prepare matcher parameters for extension
    let mparams = MatcherParams{
        matrix: params.matrix.clone(),
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        scale: params.scale,
        alternatives: 1,
        min_score: params.min_score,
    };

    // Extend around each seed by extracting a window on the original a/b.
    let mut hits: Vec<MatcherHit> = Vec::new();
    for (j_long, i_short) in seeds {
        // Map back to coordinates of a/b
        let (a_start, b_start) = if s_is_a {
            (i_short, j_long)
        } else {
            (j_long, i_short)
        };
        // Window in each sequence
        let a0 = a_start.saturating_sub(params.flank);
        let a1 = (a_start + params.wordlen + params.flank).min(a.len());
        let b0 = b_start.saturating_sub(params.flank);
        let b1 = (b_start + params.wordlen + params.flank).min(b.len());
        let sub_a = &a[a0..a1];
        let sub_b = &b[b0..b1];
        let mut hset = matcher(sub_a, sub_b, &mparams)?;
        // translate local coords to global coords
        for h in hset.drain(..) {
            let ha0 = a0 + h.a_start;
            let hb0 = b0 + h.b_start;
            let hit = MatcherHit{
                score: h.score,
                a_start: ha0, a_end: a0 + h.a_end,
                b_start: hb0, b_end: b0 + h.b_end,
                align_a: h.align_a, align_b: h.align_b,
                cigar: h.cigar,
                pct_identity: h.pct_identity, pct_gaps: h.pct_gaps,
            };
            hits.push(hit);
        }
    }
    // Optional: deduplicate hits by (a_start..a_end,b_start..b_end)
    hits.sort_by_key(|h| (h.a_start, h.b_start, -(h.score as isize)));
    hits.dedup_by(|x,y| x.a_start==y.a_start && x.a_end==y.a_end && x.b_start==y.b_start && x.b_end==y.b_end);
    Ok(hits)
}
