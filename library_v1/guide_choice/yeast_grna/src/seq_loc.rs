extern crate bio;

use std::collections::HashMap;
use std::fmt;
use std::num::ParseIntError;
use std::ops::Deref;
use std::ops::Range;
use std::rc::Rc;
use std::str::FromStr;

use bio::data_structures::interval_tree::{IntervalTree};
use bio::io::bed;
use bio::utils::Interval;
use bio::utils;

#[derive(Debug,Clone,Copy,Hash,PartialEq,Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown
}

#[derive(Debug,Clone)]
pub enum ParseError {
    NoAt(String),
    NoRefid(String),
    NoPosString(String),
    PosParse(String, ParseIntError),
    BadStrand(String),
    Backwards(String),
}

fn no_at(s: &str) -> ParseError { ParseError::NoAt(s.to_owned()) }
fn no_refid(s: &str) -> ParseError { ParseError::NoRefid(s.to_owned()) }
fn no_pos_string(s: &str) -> ParseError { ParseError::NoPosString(s.to_owned()) }
fn pos_parse(s: &str, e: ParseIntError) -> ParseError { ParseError::PosParse(s.to_owned(), e) }
fn bad_strand(s: &str) -> ParseError { ParseError::BadStrand(s.to_owned()) }
fn backwards(s: &str) -> ParseError { ParseError::Backwards(s.to_owned()) }

#[derive(Debug,Clone,PartialEq)]
pub enum LocError {
    NoStrand(SeqContigLoc, String),
    _NegativeStart(SeqContigLoc, String),
    _NegativeLength(SeqContigLoc, String),
}

/// Position on a named sequence, e.g., chromosome and position.
#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct SeqPos {
    refid: String,
    pos: usize,
    strand: Strand,
}

fn break_refid(s: &str) -> Result<(&str, &str),ParseError> {
    let breakpt = s.find(':').ok_or_else(|| no_at(s))?;
    let refid = s.get(0..breakpt).ok_or_else(|| no_refid(s))?;
    let rest = s.get((breakpt+1)..(s.len())).ok_or_else(|| no_pos_string(s))?;
    Ok( (refid, rest) )
}

fn parse_pos_strand(s: &str) -> Result<(&str, Strand),ParseError> {
    if s.get(s.len()-3..s.len()-2) == Some("(") {
        match s.get((s.len()-3)..s.len()) {
            Some("(+)") => Ok( (s.get(0..(s.len()-3)).unwrap_or(""), Strand::Forward) ),
            Some("(-)") => Ok( (s.get(0..(s.len()-3)).unwrap_or(""), Strand::Reverse) ),
            Some("(.)") => Ok( (s.get(0..(s.len()-3)).unwrap_or(""), Strand::Unknown) ),
            res => Err(bad_strand(&format!("Bad strand {:?} in {}", res, s))),
        }
    } else {
        Ok( (s, Strand::Unknown) )
    }
}

impl SeqPos {
    pub fn new(refid: &str, pos: usize, strand: Strand) -> SeqPos {
        SeqPos{ refid: refid.to_owned(), pos: pos, strand: strand }
    }

    pub fn refid(&self) -> &str { &self.refid.as_str() }
    pub fn pos(&self) -> usize { self.pos }
    pub fn strand(&self) -> Strand { self.strand }

    pub fn _to_bed(&self) -> bed::Record {
        let mut bed = bed::Record::new();
        bed.set_chrom(&self.refid);
        bed.set_start(self.pos as u64);
        bed.set_end((self.pos + 1) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.push_aux(match self.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        });
        bed
    }
}

impl fmt::Display for SeqPos {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}{}", self.refid, self.pos,
               match self.strand {
                   Strand::Forward => "(+)",
                   Strand::Reverse => "(-)",
                   Strand::Unknown => "(.)",
               })
    }
}

impl FromStr for SeqPos {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (refid, rest) = break_refid(s)?;
        let (posstr, strand) = parse_pos_strand(rest)?;
        let pos = posstr.parse::<usize>()
            .map_err(|e| pos_parse(&format!("{} => ({}, {} => ({}))", s, refid, rest, posstr), e))?;
        Ok( SeqPos { refid: refid.to_owned(), pos: pos, strand: strand } )
    }
}

impl SeqLocFeature for SeqPos {
    fn refid(&self) -> &str { self.refid() }
    fn start(&self) -> usize { self.pos() }
    fn length(&self) -> usize { 1 }
    fn strand(&self) -> Strand { self.strand }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct SeqContigLoc {
    refid: String,
    start: usize,
    length: usize,
    strand: Strand,
}

impl SeqContigLoc {
    pub fn new(refid: &str, start: usize, length: usize, strand: Strand) -> Self {
        SeqContigLoc{ refid: refid.to_owned(), strand: strand, start: start, length: length }
    }

    pub fn refid(&self) -> &str { &self.refid.as_str() }
    pub fn start(&self) -> usize { self.start }
    pub fn length(&self) -> usize { self.length }
    pub fn strand(&self) -> Strand { self.strand }

    pub fn range(&self) -> Range<usize> { self.start..(self.start + self.length) }

    pub fn set_strand(&mut self, strand: Strand) { self.strand = strand }
    
    pub fn extend_upstream_upto(&mut self, dist: usize) -> Result<usize, LocError> {
        match self.strand {
            Strand::Forward => {
                let actual = if self.start >= dist { dist } else { self.start };
                self.start -= actual;
                self.length += actual;
                Ok(actual)
            },
            Strand::Reverse => {
                self.length += dist;
                Ok(dist)
            },
            Strand::Unknown => Err( LocError::NoStrand(self.clone(), format!("extend_upstream_upto({})", dist)) ),
        }
    }

    pub fn extend_downstream_upto(&mut self, dist: usize) -> Result<usize, LocError> {
        match self.strand {
            Strand::Forward => {
                self.length += dist;
                Ok(dist)
            },
            Strand::Reverse => {
                let actual = if self.start >= dist { dist } else { self.start };
                self.start -= actual;
                self.length += actual;
                Ok(actual)
            },
            Strand::Unknown => Err( LocError::NoStrand(self.clone(), format!("extend_downstream_upto({})", dist)) ),
        }
    }

    pub fn to_bed(&self) -> bed::Record {
        let mut bed = bed::Record::new();
        bed.set_chrom(&self.refid);
        bed.set_start(self.start as u64);
        bed.set_end((self.start + self.length) as u64);
        bed.set_name("");
        bed.set_score("0");
        bed.push_aux(match self.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        });
        bed
    }
}

impl fmt::Display for SeqContigLoc {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}{}", self.refid, self.start, self.start + self.length, 
               match self.strand {
                   Strand::Forward => "(+)",
                   Strand::Reverse => "(-)",
                   Strand::Unknown => "(.)",
               })
    }
}

impl FromStr for SeqContigLoc {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (refid, rest) = break_refid(s)?;
        let (rangestr, strand) = parse_pos_strand(rest)?;
        let dashpt = rangestr.find('-').ok_or(no_pos_string(s))?;
        let start = rangestr.get(0..dashpt).ok_or(no_pos_string(s))?.
            parse::<usize>().map_err(|e| pos_parse(s, e))?;
        let end = rangestr.get((dashpt+1)..rangestr.len()).ok_or(no_pos_string(s))?.
            parse::<usize>().map_err(|e| pos_parse(s, e))?;
        if end < start {
            Err( backwards(s) )
        } else {
            Ok( SeqContigLoc { refid: refid.to_owned(), 
                               start: start, length: end - start,
                               strand: strand } )
        }
    }
}

impl SeqLocFeature for SeqContigLoc {
    fn refid(&self) -> &str { &self.refid }
    fn start(&self) -> usize { self.start }
    fn length(&self) -> usize { self.length }
    fn strand(&self) -> Strand { self.strand }
    fn loc(&self) -> SeqContigLoc { self.clone() }
}

impl <'a> From<&'a bed::Record> for SeqContigLoc {
    fn from(bed: &bed::Record) -> Self {
        Self::new(bed.chrom(), bed.start() as usize, 
                  (bed.end() - bed.start()) as usize,
                  match bed.strand() {
                      Some(utils::Strand::Forward) => Strand::Forward,
                      Some(utils::Strand::Reverse) => Strand::Reverse,
                      Some(utils::Strand::Unknown) => Strand::Unknown,
                      None => Strand::Unknown,
                  })
    }    
}

impl <'a> From<&'a SeqPos> for SeqContigLoc {
    fn from(pos: &SeqPos) -> Self { pos.loc() }
}

pub trait SeqLocFeature {
    fn refid(&self) -> &str;
    fn start(&self) -> usize;
    fn length(&self) -> usize;
    fn end(&self) -> usize { self.start() + self.length() - 1 }
    fn strand(&self) -> Strand;

    fn loc(&self) -> SeqContigLoc {
        SeqContigLoc::new(self.refid(), self.start(), 
                          self.length(), self.strand())
    }

    fn interval(&self) -> Interval<usize> {
        Interval::new(self.start()..(self.start() + self.length())).unwrap()
    }

    fn first(&self) -> Option<usize> {
        match self.strand() {
            Strand::Forward => if self.length() > 0 { Some(self.start()) } else { None },
            Strand::Reverse => if self.length() > 0 { Some(self.start() + self.length() - 1) } else { None },
            Strand::Unknown => None,
        }
    }

    fn last(&self) -> Option<usize> {
        match self.strand() {
            Strand::Forward => if self.length() > 0 { Some(self.start() + self.length() - 1) } else { None },
            Strand::Reverse => if self.length() > 0 { Some(self.start()) } else { None },
            Strand::Unknown => None,
        }
    }

    fn first_seqpos(&self) -> Option<SeqPos> {
        self.first().map(|pos| SeqPos::new(self.refid(), pos, self.strand()))
    }

    fn last_seqpos(&self) -> Option<SeqPos> {
        self.last().map(|pos| SeqPos::new(self.refid(), pos, self.strand()))
    }        
}

impl <F> SeqLocFeature for Rc<F>
    where F: SeqLocFeature
{
    fn refid(&self) -> &str { self.deref().refid() }
    fn start(&self) -> usize { self.deref().start() }
    fn length(&self) -> usize { self.deref().length() }
    fn strand(&self) -> Strand { self.deref().strand() }
}

#[derive(Debug)]
pub struct SeqLocMap<F> 
    where F: SeqLocFeature
{
    map: HashMap<String,IntervalTree<usize,F>>
}

impl<F> SeqLocMap<F>
    where F: SeqLocFeature
{
    pub fn new() -> Self { Self { map: HashMap::new() } }

    pub fn insert(&mut self, f: F) -> ()
    {
        let refid = f.refid().to_owned();
        let loc = f.start()..(f.start() + f.length());
        let ref_itree = self.map.entry(refid).or_insert_with(|| IntervalTree::new());
        ref_itree.insert(loc, f)
    }

    pub fn _build<I>(features: I) -> Self
        where I: Iterator<Item=F>
    {
        let mut map = Self::new();
        for f in features {
            map.insert(f);
        }
        map
    }

    fn strand_match(s1: Strand, s2: Strand) -> bool {
        !( (s1 == Strand::Forward && s2 == Strand::Reverse)
            || (s1 == Strand::Reverse && s2 == Strand::Forward) )
    }

    pub fn features_by_loc<G>(&self, loc: &G) -> Vec<&F>
        where G: SeqLocFeature
    {
        if let Some(itree) = self.map.get(loc.refid()) {
            let entries = itree.find(loc.interval());
            entries.filter_map(|f| if Self::strand_match(f.data().strand(), loc.strand()) { Some(f.data()) } else { None }).collect()
        } else {
            Vec::new()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug,Clone)]
    struct TestFeature {
        name: String,
        loc: SeqContigLoc,
    }

    impl TestFeature {
        fn new(name: &str, loc: SeqContigLoc) -> Self {
            TestFeature { name: name.to_owned(), loc: loc }
        }
        fn name(&self) -> &str { &self.name }
        fn loc(&self) -> &SeqContigLoc { &self.loc }
    }

    impl SeqLocFeature for TestFeature {
        fn refid(&self) -> &str { &self.loc.refid }
        fn start(&self) -> usize { self.loc.start }
        fn length(&self) -> usize { self.loc.length }
        fn strand(&self) -> Strand { self.loc.strand }
        fn loc(&self) -> SeqContigLoc { self.loc.clone() }
    }

    fn feature_names<'a>(xs: Vec<&'a TestFeature>) -> Vec<&'a str> {
        xs.into_iter().map(TestFeature::name).collect()
    }

    #[test]
    fn seq_pos_string() {
        let pos = SeqPos::new("chrV", 166885, Strand::Reverse);
        assert_eq!(pos.to_string(), "chrV:166885(-)");
        assert_eq!(pos.to_string().parse::<SeqPos>().unwrap(), pos);

        let pos = SeqPos::new("chrX", 461829, Strand::Forward);
        assert_eq!(pos.to_string(), "chrX:461829(+)");
        assert_eq!(pos.to_string().parse::<SeqPos>().unwrap(), pos);
    }

    #[test]
    fn seq_contig_loc() {
        let tma20 = SeqContigLoc::new("chrV", 166237, (166886 - 166237), Strand::Reverse);
        let tma22 = SeqContigLoc::new("chrX", 461829, (462426 - 461829), Strand::Forward);

        assert_eq!(tma20.refid(), "chrV");
        assert_eq!(tma20.start(), 166237);
        assert_eq!(tma20.length(), 649);
        assert_eq!(tma20.strand(), Strand::Reverse);

        assert_eq!(tma20.to_string(), "chrV:166237-166886(-)");
        assert_eq!("chrV:166237-166886(-)".parse::<SeqContigLoc>().unwrap(), tma20);

        assert_eq!(tma20.first_seqpos().unwrap().to_string(), "chrV:166885(-)");
        assert_eq!(tma20.last_seqpos().unwrap().to_string(), "chrV:166237(-)");

        assert_eq!(tma22.refid(), "chrX");
        assert_eq!(tma22.start(), 461829);
        assert_eq!(tma22.length(), 597);
        assert_eq!(tma22.strand(), Strand::Forward);

        assert_eq!(tma22.to_string(), "chrX:461829-462426(+)");
        assert_eq!("chrX:461829-462426(+)".parse::<SeqContigLoc>().unwrap(), tma22);

        assert_eq!(tma22.first_seqpos().unwrap().to_string(), "chrX:461829(+)");
        assert_eq!(tma22.last_seqpos().unwrap().to_string(), "chrX:462425(+)");
    }

    #[test]
    fn extend_upto() {
        let mut scl = SeqContigLoc::from_str("chrI:200-300(+)").unwrap();
        assert_eq!(scl.extend_upstream_upto(50), Ok(50));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:150-300(+)").unwrap());
        assert_eq!(scl.extend_downstream_upto(75), Ok(75));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:150-375(+)").unwrap());

        let mut scl = SeqContigLoc::from_str("chrI:500-600(-)").unwrap();
        assert_eq!(scl.extend_upstream_upto(50), Ok(50));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:500-650(-)").unwrap());
        assert_eq!(scl.extend_downstream_upto(75), Ok(75));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:425-650(-)").unwrap());

        let mut scl = SeqContigLoc::from_str("chrI:1000-2000(+)").unwrap();
        assert_eq!(scl.extend_upstream_upto(1500), Ok(1000));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:0-2000(+)").unwrap());
        assert_eq!(scl.extend_downstream_upto(2500), Ok(2500));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:0-4500(+)").unwrap());

        let mut scl = SeqContigLoc::from_str("chrI:40-80(-)").unwrap();
        assert_eq!(scl.extend_upstream_upto(120), Ok(120));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:40-200(-)").unwrap());
        assert_eq!(scl.extend_downstream_upto(80), Ok(40));
        assert_eq!(scl, SeqContigLoc::from_str("chrI:0-200(-)").unwrap());
    }

    #[test]
    fn seq_loc_map() {
        let mut slm = SeqLocMap::new();

        slm.insert(TestFeature::new("TMA20", SeqContigLoc::new("chrV", 166237, (166886 - 166237), Strand::Reverse)));
        slm.insert(TestFeature::new("TMA22", SeqContigLoc::new("chrX", 461829, (462426 - 461829), Strand::Forward)));

        let none: Vec<&str> = Vec::new();

        let qy = SeqContigLoc::new("chrV", 166300, 100, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), vec!["TMA20"]);

        let qy = SeqContigLoc::new("chrV", 166100, 100, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), none);

        let qy = SeqContigLoc::new("chrV", 166900, 100, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), none);

        let qy = SeqContigLoc::new("chrVI", 166300, 100, Strand::Forward);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), none);

        let qy = SeqContigLoc::new("chrV", 166300, 100, Strand::Unknown);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), vec!["TMA20"]);

        let qy = SeqContigLoc::new("chrV", 166230, 8, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), vec!["TMA20"]);

        let qy = SeqContigLoc::new("chrV", 166230, 7, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), none);

        let qy = SeqContigLoc::new("chrV", 166884, 10, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), vec!["TMA20"]);

        let qy = SeqContigLoc::new("chrV", 166885, 10, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), vec!["TMA20"]);

        let qy = SeqContigLoc::new("chrV", 166886, 10, Strand::Reverse);
        assert_eq!(feature_names(slm.features_by_loc(&qy)), none);
    }
}
