extern crate bio;
extern crate rust_htslib;
extern crate itertools;

use std::cmp::*;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::ops::{Range};
use std::path;
use std::str;

use gene::{GeneStart,GenomeAnnot};
use myerr::MyErr;
use seq_loc::*;
use target::{Target};

#[derive(Debug,Clone)]
pub struct TargetChoiceConfig {
    pub tss_upstream: usize,
    pub tss_downstream: usize,
    pub tss_zones: Vec<Range<isize>>,
    pub cds_upstream: usize,
    pub cds_downstream: usize,
    pub cds_zones: Vec<Range<isize>>,
    pub num_per_gene: usize,
    pub work_dir: path::PathBuf,
}

impl TargetChoiceConfig {
    pub fn target_loc_for_start(&self, gene_start: &GeneStart) -> SeqContigLoc {
        let mut loc = gene_start.first_pos().loc();
        loc.extend_upstream_upto(if gene_start.is_cds() { self.cds_upstream } else { self.tss_upstream })
            .unwrap_or_else(|e| panic!("Unexpected failure extending target loc upstream: {:?}", e));
        loc.extend_downstream_upto(if gene_start.is_cds() { self.cds_downstream } else { self.tss_downstream })
            .unwrap_or_else(|e| panic!("Unexpected failure extending target loc downstream: {:?}", e));
        loc.set_strand(Strand::Unknown);
        loc
    }
}

struct GeneStartTargetLoc {
    target_loc: SeqContigLoc,
    gene_start: GeneStart
}

impl GeneStartTargetLoc {
    fn new(target_loc: SeqContigLoc, gene_start: GeneStart) -> Self {
        GeneStartTargetLoc { target_loc: target_loc, gene_start: gene_start }
    }
    fn yorf(&self) -> &str { self.gene_start.yorf() }
}

impl SeqLocFeature for GeneStartTargetLoc {
    fn refid(&self) -> &str { self.target_loc.refid() }
    fn start(&self) -> usize { self.target_loc.start() }
    fn length(&self) -> usize { self.target_loc.length() }
    fn strand(&self) -> Strand { self.target_loc.strand() }
}

fn gene_start_target_locs(config: &TargetChoiceConfig,
                          annot: &GenomeAnnot)
                          -> Result<SeqLocMap<GeneStartTargetLoc>,MyErr> {
    let mut map = SeqLocMap::new();
    for gene in annot.genes().into_iter() {
        match annot.gene_start(gene) {
            Some(gene_start) =>
                Ok(map.insert(GeneStartTargetLoc::new(config.target_loc_for_start(&gene_start), gene_start))),
            None => Err(MyErr::Str(format!("gene_start_target_locs: no GeneStart for {}", gene))),
        }?;
    }
    Ok( map )
}

#[derive(Debug)]
pub struct ScoringInfo {    
    uniqueness: HashMap<SeqContigLoc,bool>,
    atacs: HashMap<SeqContigLoc,f64>,
    annot: GenomeAnnot,
    loc_to_target: HashMap<SeqContigLoc,Target>,
    gene_to_target: HashMap<String,Vec<SeqContigLoc>>,
    target_to_gene: HashMap<SeqContigLoc,Vec<String>>,
}

impl ScoringInfo {
    pub fn new(config: &TargetChoiceConfig,
               annot: GenomeAnnot,
               uniqueness: HashMap<SeqContigLoc,bool>,
               atacs: HashMap<SeqContigLoc,f64>,
               targets: Vec<Target>)
        -> Result<Self,MyErr>
    {
        let gene_start_loc_map = gene_start_target_locs(config, &annot)?;

        let mut loc_to_target = HashMap::new();
        let mut gene_to_target = HashMap::new();
        let mut target_to_gene = HashMap::new();

        let mut table = fs::File::create(config.work_dir.join("target-and-gene-all.txt"))?;

        for target in targets.into_iter() {
            let targeted = gene_start_loc_map.features_by_loc(&target.pos().loc());

            if !targeted.is_empty() {
                for gs_target in targeted {
                    let yorf_entry = gene_to_target.entry(gs_target.yorf().to_owned()).or_insert_with(|| Vec::new());
                    yorf_entry.push(target.loc().clone());

                    let target_entry = target_to_gene.entry(target.loc().clone()).or_insert_with(|| Vec::new());
                    target_entry.push(gs_target.yorf().to_owned());

                    write!(table, "{}\t{}\t{}\t{}\n", target.loc(), gs_target.yorf(),
                           gs_target.gene_start.target_relative(&target),
                           if gs_target.gene_start.is_tss() { "TSS" } else if gs_target.gene_start.is_cds() { "CDS" } else { "Unknown" })?;
                }
                loc_to_target.insert(target.loc().clone(), target);
            }
        }

        let info = ScoringInfo { uniqueness: uniqueness, atacs: atacs, annot: annot,
                                 loc_to_target: loc_to_target,
                                 gene_to_target: gene_to_target,
                                 target_to_gene: target_to_gene };

        info.write_target_table(config)?;
        info.write_gene_table(config)?;
        
        Ok( info )
    }

    pub fn annot(&self) -> &GenomeAnnot { &self.annot }

    pub fn loc_target(&self, loc: &SeqContigLoc) -> Result<&Target,MyErr> {
        self.loc_to_target.get(loc)
            .ok_or_else(|| MyErr::Str(format!("Cannot find target for {}", loc)))
    }

    pub fn gene_targets(&self, yorf: &str) -> Vec<&SeqContigLoc> {
        match self.gene_to_target.get(yorf) {
            Some(targets) => targets.iter().collect(),
            None => Vec::new(),
        }
    }
    
    pub fn is_loc_unique(&self, loc: &SeqContigLoc) -> Result<bool,MyErr> {
        self.uniqueness.get(loc)
            .ok_or_else(|| MyErr::Str(format!("No uniqueness entry for {}", loc)))
            .map(|u| *u )
    }

    pub fn is_loc_specific(&self, loc: &SeqContigLoc, yorf: &str) -> Result<bool,MyErr> {
        let start_yorfs = self.target_to_gene.get(loc).
            ok_or_else(|| MyErr::Str(format!("No targets entry for {}", loc)))?;
        if !start_yorfs.iter().any(|y| y == yorf) {
            Err(MyErr::Str(format!("No target-to-gene entry for {} targeting {}", loc, yorf)))
        } else {
            Ok( start_yorfs.iter().all(|y| y == yorf || self.annot.dubious().contains(y)) )
        }
    }

    pub fn loc_atac(&self, loc: &SeqContigLoc) -> Result<f64,MyErr> {
        self.atacs.get(loc)
            .ok_or_else(|| MyErr::Str(format!("No ATAC-Seq entry for {}", loc)))
            .map(|a| *a )
    }
    
    fn write_target_table(&self, config: &TargetChoiceConfig) -> Result<(),MyErr> {
        let mut table = fs::File::create(config.work_dir.join("target-to-genes.txt"))?;
        for (target, genes) in self.target_to_gene.iter() {
            write!(table, "{}", target)?;
            for gene in genes.iter() {
                write!(table, "\t{}", gene)?;
            }
            write!(table, "\n")?;
        }
        Ok( () )
    }

    fn write_gene_table(&self, config: &TargetChoiceConfig) -> Result<(),MyErr> {
        let mut table = fs::File::create(config.work_dir.join("gene-to-targets.txt"))?;
        for (gene, targets) in self.gene_to_target.iter() {
            write!(table, "{}", gene)?;
            for target in targets.iter() {
                write!(table, "\t{}", target)?;
            }
            write!(table, "\n")?;
        }
        Ok( () )
    }
}

#[derive(Debug,Clone)]
pub struct TargetChoice {
    yorf: String,
    target_score: TargetScore,
    guide_no: usize,
    rank: usize,
    zone: Option<usize>,
}

impl TargetChoice {
    pub fn new(yorf: &str, target_score: TargetScore, guide_no: usize, rank: usize, zone: Option<usize>) -> Self
    {
        TargetChoice{ yorf: yorf.to_owned(), target_score: target_score,
                      guide_no: guide_no, rank: rank, zone: zone }
    }

    pub fn yorf(&self) -> &str { &self.yorf }
    pub fn target_score(&self) -> &TargetScore { &self.target_score }
    pub fn guide_no(&self) -> usize { self.guide_no }
    pub fn rank(&self) -> usize { self.rank }
    pub fn zone(&self) -> Option<usize> { self.zone }

    pub fn _guide_name(&self) -> String {
        format!("{}_g{:02}", self.yorf(), self.guide_no())
    }

    pub fn guide_desc(&self) -> String {
        format!("{}_g{:02}_{:+}_{}_{}_{:0.3}",
                self.yorf, self.guide_no, self.target_score.offset(),
                if self.target_score.unique() { "U" } else { "D" },
                if self.target_score.specific() { "S" } else { "N" },
                self.target_score.atac())
    }
}

#[derive(Debug,Clone)]
pub struct TargetScore {
    target: Target,
    unique: bool,
    specific: bool,
    atac: f64,
    offset: isize,
}

impl TargetScore {
    pub fn target(&self) -> &Target { &self.target }
    pub fn unique(&self) -> bool { self.unique }
    pub fn specific(&self) -> bool { self.specific }
    pub fn atac(&self) -> f64 { self.atac }
    pub fn offset(&self) -> isize { self.offset }
}

pub fn choose_targets(config: &TargetChoiceConfig,
                      info: &ScoringInfo)
                      -> Result<HashMap<String,Vec<TargetChoice>>,MyErr>
{
    let genes = info.annot.genes();

    let mut choices = HashMap::new();

    let mut all_score_file = fs::File::create(config.work_dir.join("all-scores.txt"))?;

    for gene in genes {
        let targets = choose_gene_targets(config, &mut all_score_file, info, gene)?;
        choices.insert(gene.to_owned(), targets);
    }

    Ok(choices)
}

fn choose_gene_targets<W>(config: &TargetChoiceConfig,
                          out: &mut W,
                          info: &ScoringInfo,
                          yorf: &str)
                          -> Result<Vec<TargetChoice>,MyErr>
    where W: Write
{
    let gene_start = info.annot.gene_start(yorf)
        .ok_or_else(|| MyErr::Str(format!("No gene start for {}", yorf)))?;

    let scores = score_all_gene_targets(info, &gene_start)?;
    write_scores(out, &gene_start, &scores)?;
    let mut ranked_scores: Vec<(usize,TargetScore)> = scores.into_iter().enumerate().collect();
    
    let zones = if gene_start.is_tss() {
        Ok(&config.tss_zones)
    } else if gene_start.is_cds() {
        Ok(&config.cds_zones)
    } else {
        Err(MyErr::Str(format!("Unknown gene annotation {:?}", gene_start)))
    }?;

    let mut choices: Vec<(Option<usize>,(usize,TargetScore))> = Vec::new();
    
    for (zoneno, ref zone) in zones.iter().enumerate() {
        match extract_best_in_zone(&mut ranked_scores, &zone) {
            Some(zone_score) => choices.push( (Some(zoneno), zone_score) ),
            None => ()
        };
    }

    for score in ranked_scores {
        choices.push( (None, score) );
    }

    let mut target_choices = Vec::new();
    for (guide_no, (zone, (rank, score))) in choices.into_iter().enumerate() {
        target_choices.push(TargetChoice::new(yorf, score, guide_no, rank, zone));
    }
    
    Ok( target_choices.into_iter().take(config.num_per_gene).collect() )
}

fn write_scores<W>(out: &mut W,
                   gene_start: &GeneStart,
                   scores: &Vec<TargetScore>)
                   -> Result<(),MyErr>
    where W: Write
{
    if scores.is_empty() {
        write!(out, "{}\tNONE\n", gene_start.yorf())?;
    } else {
        for score in scores {
            let sequ_str = str::from_utf8(score.target().sequ())?;
            write!(out, "{}\t{}\t{}\t{:+}\t{}\t{}\n",
                   gene_start.yorf(), sequ_str, score.target().loc(), score.offset(),
                   gene_start, gene_start.first_pos())?;
        }
    }
    Ok( () )
}
   

fn extract_best_in_zone(ranked_target_scores: &mut Vec<(usize, TargetScore)>, 
                        zone: &Range<isize>)
                        -> Option<(usize,TargetScore)>
{
    for idx in 0..ranked_target_scores.len() {
        let offset = ranked_target_scores[idx].1.offset;
        if offset >= zone.start && offset < zone.end {
            let best = ranked_target_scores.remove(idx);
            return Some(best);
        }
    }

    None
}

fn score_all_gene_targets(info: &ScoringInfo,
                          gene_start: &GeneStart)
                          -> Result<Vec<TargetScore>,MyErr>
{
    let mut target_scores = Vec::new();
    let target_locs = info.gene_targets(gene_start.yorf());

    for loc in target_locs {
        let target = info.loc_target(loc)?;
        let score = score_target_for_gene(info, gene_start, target)?;
        target_scores.push(score);
    }

    target_scores.sort_unstable_by(compare_target_scores);
    
    Ok( target_scores )
}

fn score_target_for_gene(info: &ScoringInfo,
                         gene_start: &GeneStart,
                         target: &Target)
                         -> Result<TargetScore,MyErr>
{
    let loc = target.loc();
    let unique = info.is_loc_unique(loc)?;
    let specific = info.is_loc_specific(loc, gene_start.yorf())?;
    let atac = info.loc_atac(loc)?;
    
    let offset = gene_start.target_relative(target);

    Ok( TargetScore{ target: target.clone(), unique: unique,
                     specific: specific, atac: atac, 
                     offset: offset } )
}

fn compare_target_scores(ts1: &TargetScore, ts2: &TargetScore) -> Ordering {
    if ts1.unique && !ts2.unique {
        return Ordering::Less;
    } else if !ts1.unique && ts2.unique {
        return Ordering::Greater;
    }

    if ts1.specific && !ts2.specific {
        return Ordering::Less;
    } else if !ts1.specific && ts2.specific {
        return Ordering::Greater;
    }

    if ts1.atac > ts2.atac {
        return Ordering::Less;
    } else if ts2.atac > ts1.atac {
        return Ordering::Greater;
    }

    return Ordering::Equal;
}
