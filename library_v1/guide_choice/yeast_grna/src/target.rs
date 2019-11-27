extern crate bio;
extern crate rust_htslib;

use std::iter;
use std::path;

use bio::io::fasta;
use bio::pattern_matching::shift_and;

use myerr::{MyErr,Annotable};
use seq_loc::*;

/// Describes how a target is specified: the length of the target
/// sequence overall, the PAM sequence and position, and the guide
/// sequence and position.
#[derive(Debug)]
pub struct TargetSpec {
    /// Length of the target seuqence overall
    pub target_len: usize,
    /// Sequence of the PAM
    pub pam: Vec<u8>,
    /// Starting position of the PAM in the target
    pub pam_start: usize,
    /// Starting position of the guide in the target
    pub guide_start: usize,
    /// Length of the guide in the target
    pub guide_len: usize
}

impl TargetSpec {
    /// Extract the sequence of a target.
    ///
    /// The target must be on the forward strand of the
    /// sequence. Guide sequences are copied from the reference
    /// sequence, the PAM is copied from `self`, and other bases
    /// (neither target nor PAM) are given as `N`.
    ///
    /// # Arguments
    /// * `sequ` is the overall reference sequence
    /// * `target_sequ_start` is the 0-based start of the target
    pub fn _get_target_sequ(&self, sequ: &[u8], target_sequ_start: usize) -> Vec<u8> {
        let mut target_sequ: Vec<u8> = iter::repeat(b'N').take(self.target_len).collect();
        
        for pam_pos in 0..self.pam.len() {
            target_sequ[pam_pos + self.pam_start] = self.pam[pam_pos];
        }
        
        for guide_pos in 0..self.guide_len {
            target_sequ[self.guide_start + guide_pos] = sequ[target_sequ_start + self.guide_start + guide_pos];
        }
        
        target_sequ
    }

    pub fn guide_sequence(&self, target: &Target) -> Vec<u8> {
        target.sequ[self.guide_start..(self.guide_start + self.guide_len)].to_owned()
    }
    
    // ZZZ UNTESTED
    // pub fn get_target_id_sequ(&self, sequ: &[u8], loc: &SeqContigLoc) -> Vec<u8> {
    //     let mut target_sequ: Vec<u8> = iter::repeat(b'N').take(self.target_len).collect();
        
    //     for pam_pos in 0..self.pam.len() {
    //         target_sequ[pam_pos + self.pam_start] = self.pam[pam_pos];
    //     }
        
    //     for guide_pos in 0..self.guide_len {
    //         let fwdpos = guide_pos + self.guide_start;
    //         target_sequ[fwdpos] = if loc.strand_fwd().unwrap_or(true) {
    //             sequ[loc.start() + fwdpos]
    //         } else {
    //             dna::complement(sequ[loc.start() + self.target_len - (fwdpos + 1)])
    //         }
    //     }
        
    //     target_sequ        
    // }
}

/// Target identifier plus the sequence of the target
#[derive(Debug,Clone)]
pub struct Target {
    loc: SeqContigLoc,
    sequ: Vec<u8>,
}

impl Target {
    /// Create a target by location and sequence
    pub fn new(loc: SeqContigLoc, sequ: Vec<u8>) -> Self { Target { loc: loc, sequ: sequ } }

    /// Target location
    pub fn loc(&self) -> &SeqContigLoc { &self.loc }

    /// Target position
    pub fn pos(&self) -> SeqPos { 
        SeqPos::new(self.loc.refid(), 
                    self.loc.start() + self.loc.length() / 2,
                    self.loc.strand())
    }

    /// Target sequence
    pub fn sequ(&self) -> &[u8] { &self.sequ }

    /// Target as a Fasta-format record, named according to location
    pub fn record(&self) -> fasta::Record {
        fasta::Record::with_attrs(&self.loc().to_string(), None, &self.sequ)
    }

    /// Construct a target from a BAM record, using the query name and sequence
    ///
    /// The query name is converted into the target location
    ///
    /// # Error
    /// Errors arise in parsing the query name as per `SeqContigLoc::parse()`.
    pub fn from_bam(rec: &rust_htslib::bam::record::Record) -> Result<Self, MyErr> {
        let qname = String::from_utf8(rec.qname().to_vec())?;
        let loc = qname.parse::<SeqContigLoc>()?;
        Ok( Target{ loc: loc, sequ: rec.seq().as_bytes() } )
    }

    /// Construct a TSV-format text record for the target.
    ///
    /// The first field is the target location and the second field is
    /// the target sequence.
    pub fn to_tsv(&self) -> Vec<u8> {
        let mut tsv = (self.loc.to_string() + "\t").into_bytes();
        tsv.extend_from_slice(&self.sequ);
        tsv
    }

    /// Parse a TSV-format text record.
    ///
    /// The first field is parsed as the target location and the
    /// second is used as the sequence.
    /// 
    /// # Error Errors arise in parsing the query name as per
    /// `SeqContigLoc::parse()`, or if there are not exactly two
    /// tab-separated fields.
    pub fn _from_tsv(tsv: &str) -> Result<Target, MyErr> {
        let tabpt = tsv.find('\t')
            .map_or(Err(MyErr::Str("No tab in target entry ".to_string() + tsv)), Ok)?;
        let loc = tsv.get(0..tabpt).
            map_or(Err(MyErr::Str("No target id in target entry ".to_string() + tsv)), Ok)?.
            parse::<SeqContigLoc>()?;
        let sequ_str = tsv.get((tabpt + 1)..tsv.len())
            .map_or(Err(MyErr::Str("No sequence string in target entry ".to_string() + tsv)), Ok)?;
        Ok( Target{ loc: loc, sequ: Vec::from(sequ_str.as_bytes()) } )
    }
}

/// Find all targets in a genome
///
/// # Arguments
/// * `target_spec` specifies what a target looks like
/// * `genome_fasta` is the filename of a genome fasta file
pub fn genomic_targets(target_spec: &TargetSpec, 
                       genome_fasta: &path::Path)
                       -> Result<Vec<Target>,MyErr> {
    let genome_reader = fasta::Reader::from_file(genome_fasta)
        .annot_err(format!("Opening genome fasta file \"{:?}\"", genome_fasta))?;

    let contig_records = genome_reader.records();

    let mut targets: Vec<Target> = Vec::new();
    
    for contig_result in contig_records {
        let contig = contig_result?;
        let mut contig_targets = contig_targets(target_spec, &contig)?;
        targets.append(&mut contig_targets);
    }

    Ok( targets )
}

/// Find all targets in a query sequence
///
/// # Arguments
/// * `spec` specifies what a target looks like
/// * `contig_record` is a Fasta-format record of the query sequence
pub fn contig_targets(spec: &TargetSpec, contig_record: &fasta::Record) -> Result<Vec<Target>, MyErr> {
    let mut fwd = forward_target_records(spec, contig_record);
    let mut rev = reverse_target_records(spec, contig_record);
        
    fwd.append(&mut rev);
    
    Ok(fwd)
}

fn forward_target_records(spec: &TargetSpec, contig_record: &fasta::Record) -> Vec<Target>
{
    let contig_sequ = contig_record.seq();
    let contig_id = contig_record.id();
    strand_targets(spec, &contig_id, Strand::Forward, contig_sequ)
}

fn reverse_target_records(spec: &TargetSpec, contig_record: &fasta::Record) -> Vec<Target>
{
    let rc_contig_sequ = bio::alphabets::dna::revcomp(contig_record.seq());
    let contig_id = contig_record.id();
    strand_targets(spec, &contig_id, Strand::Reverse, &rc_contig_sequ)
}

fn target_get_sequ(spec: &TargetSpec, sequ: &[u8], target_sequ_start: usize) -> Vec<u8> {
    let mut target_sequ: Vec<u8> = iter::repeat(b'N').take(spec.target_len).collect();
    
    for pam_pos in 0..spec.pam.len() {
        target_sequ[pam_pos + spec.pam_start] = spec.pam[pam_pos];
    }
    
    for guide_pos in 0..spec.guide_len {
        target_sequ[spec.guide_start + guide_pos] = sequ[target_sequ_start + spec.guide_start + guide_pos];
    }

    target_sequ
}

fn strand_targets<'a>(spec: &TargetSpec, refid: &str, strand: Strand, sequ: &[u8]) -> Vec<Target> {
    let mut strand_pam_targets = Vec::new();

    let pam_matcher = shift_and::ShiftAnd::new(&spec.pam);
    let pam_starts = pam_matcher.find_all(sequ);

    for pam_start in pam_starts {
        if pam_start >= spec.pam_start {
            let target_sequ_start = pam_start - spec.pam_start;
            let target_sequ_end = target_sequ_start + spec.target_len;

            if target_sequ_end <= sequ.len() {
                let target_sequ = target_get_sequ(spec, sequ, target_sequ_start);

                let target_start = match strand {
                    Strand::Forward => target_sequ_start,
                    Strand::Reverse => (sequ.len() - spec.target_len) - target_sequ_start,
                    _ => target_sequ_start, // ZZZ
                };
                
                let loc = SeqContigLoc::new( refid, target_start, spec.target_len, strand );
                let target = Target::new( loc, target_sequ );
                strand_pam_targets.push(target);
            }
        }
    }
    
    strand_pam_targets
}

