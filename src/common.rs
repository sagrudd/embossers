//! Common helpers shared by the algorithms: minimal FASTA parsing, scoring
//! configuration, and a compact BLOSUM62 implementation.
//!
//! ## FASTA
//! The parser is intentionally permissive and suitable for small/medium files
//! and tests. It supports multi-record inputs and keeps all non‑alphabetic
//! symbols as‑is (conversion to uppercase only).
//!
//! ## Scoring
//! The [`WaterMatrix`] enum describes either a simple DNA match/mismatch
//! scheme or the built‑in BLOSUM62 table for proteins.
//!
//! ## Examples
//! ```rust,no_run
//! use embossers::parse_fasta;
//! let recs = parse_fasta(r#">seq
//! ACGT
//! >p
//! PAWHEAE
//! "#);
//! assert_eq!(recs.len(), 2);
//! assert_eq!(recs[0].seq, "ACGT");
//! ```
//!

/// Errors that can be returned by the algorithms in this crate.
#[derive(thiserror::Error, Debug)]
pub enum EmbossersError {
    /// Returned if `jmin` > `jmax` or either is less than 1.
    #[error("invalid j-range: jmin={jmin}, jmax={jmax}")]
    InvalidJRange { jmin: usize, jmax: usize },
    /// Returned if `lwin` or `step` is zero.
    #[error("window length and step must be > 0 (lwin={lwin}, step={step})")]
    InvalidWindow { lwin: usize, step: usize },
    /// Returned when sequence input is empty or otherwise invalid.
    #[error("invalid sequence input: {0}")]
    InvalidSequence(&'static str),
}

/// A single FASTA sequence (identifier and uppercase sequence letters).
#[derive(Clone, Debug)]
/// A simple in-memory FASTA record parsed by [`parse_fasta`].
pub struct FastaRecord {
    /// Identifier from the FASTA header (text after '>').
    pub id: String,
    /// Raw sequence (uppercase). Non-ACGT or non-amino-acid symbols are kept as-is by the parser.
    pub seq: String,
}

/// Parse a minimal FASTA string into records (tolerant of whitespace).

/// Parse a minimal FASTA string into a vector of [`FastaRecord`].
///
/// *Lines starting with `>` start a new record.* All other lines are appended
/// (without spaces) to the current sequence. Sequences are uppercased.
///
/// ## Panics
/// This function does not panic.
///
/// ## Examples
/// ```rust,no_run
/// use embossers::parse_fasta;
/// let recs = parse_fasta(r#">id\nAC\nGT\n"#);
/// assert_eq!(recs[0].id, "id");
/// assert_eq!(recs[0].seq, "ACGT");
/// ```
pub fn parse_fasta(text: &str) -> Vec<FastaRecord> {
    let mut out: Vec<FastaRecord> = vec![];
    let mut id = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !id.is_empty() { out.push(FastaRecord{ id: id.clone(), seq: seq.to_ascii_uppercase() }); seq.clear(); }
            id = rest.trim().split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !id.is_empty() { out.push(FastaRecord{ id, seq: seq.to_ascii_uppercase() }); }
    out
}

/// Scoring scheme enum shared by `water` and `needle`.
#[derive(Clone, Debug)]

/// Scoring matrix selection used by both local and global aligners.
///
/// - For DNA, supply integer `match_score` and `mismatch` values.
/// - For proteins, use the built‑in BLOSUM62.
pub enum WaterMatrix {
    /// DNA scoring: match/mismatch integer scoring.
    Dna { match_score: i32, mismatch: i32 },
    /// Protein scoring with a static BLOSUM62 matrix.
    Blosum62,
}

/// BLOSUM62 look-up using a fixed 20x20 matrix (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V).

/// Return the BLOSUM62 score for the given amino acids. Unknowns score −4.
///
/// This is a compact, fixed lookup (no external files).
pub fn blosum62_score(x: char, y: char) -> i32 {
    fn idx(c: char) -> Option<usize> {
        match c.to_ascii_uppercase() {
            'A'=>Some(0),'R'=>Some(1),'N'=>Some(2),'D'=>Some(3),'C'=>Some(4),
            'Q'=>Some(5),'E'=>Some(6),'G'=>Some(7),'H'=>Some(8),'I'=>Some(9),
            'L'=>Some(10),'K'=>Some(11),'M'=>Some(12),'F'=>Some(13),'P'=>Some(14),
            'S'=>Some(15),'T'=>Some(16),'W'=>Some(17),'Y'=>Some(18),'V'=>Some(19),
            _=>None
        }
    }
    const M: [[i32;20];20] = [
        [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0], // A
        [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3], // R
        [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3], // N
        [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3], // D
        [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1], // C
        [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2], // Q
        [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2], // E
        [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3], // G
        [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3], // H
        [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3], // I
        [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1], // L
        [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2], // K
        [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1], // M
        [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1], // F
        [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2], // P
        [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2], // S
        [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0], // T
        [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3], // W
        [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1], // Y
        [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4], // V
    ];
    if let (Some(ix), Some(iy)) = (idx(x), idx(y)) {
        M[ix][iy]
    } else {
        -4 // default for unknowns
    }
}

/// Count matches case-insensitively; used for identity %. 

/// Case‑insensitive equality used when computing identity % in alignments.
pub fn equals_case_insensitive(a: char, b: char) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}
