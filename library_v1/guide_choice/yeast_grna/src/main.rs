extern crate bio;
extern crate rust_htslib;
extern crate itertools;

use std::collections::{HashMap};
use std::env;
use std::fs;
use std::io::{Write};
use std::path;
use std::str;

use bio::io::bed;

use atac_seq::{AtacConfig,target_average_atacs};
use gene::{GeneConfig,GenomeAnnot};
use myerr::{MyErr,Annotable};
use seq_loc::*;
use target::{TargetSpec,Target,genomic_targets};
use target_choice::{TargetChoiceConfig,ScoringInfo,choose_targets};
use unique::{UniqueConfig,target_uniqueness};

mod atac_seq;
mod bed_graph;
mod gene;
mod myerr;
mod seq_loc;
mod target;
mod target_choice;
mod tifs;
mod unique;

#[derive(Debug)]
struct Config {
    genome_fasta: path::PathBuf,
    work_dir: path::PathBuf,
    unique_config: UniqueConfig,
    atac_config: AtacConfig,
    gene_config: GeneConfig,
    target_choice_config: TargetChoiceConfig,
    output_dir: path::PathBuf,
    target_spec: TargetSpec,
    oligo_5ext: String,
    oligo_3ext: String,
    quick_test: Option<usize>,
}

fn main() {
    let argv: Vec<String> = env::args().collect();

    if argv.len() != 2 {
        panic!("Exactly one argument, the work directory: {:?}", argv);
    }

    let base_dir = path::PathBuf::from(&argv[1]);
    let data_dir = base_dir.join("data");
    let work_dir = base_dir.join("work");
    let config = Config {
        genome_fasta: data_dir.join("sacCer3.fa"),
        work_dir: work_dir.clone(),
        unique_config: UniqueConfig {
            genome_bt2: data_dir.join("sacCer3").to_string_lossy().into_owned(),
            work_dir: work_dir.clone(),
        },
        atac_config: AtacConfig {
            atac_wigs: vec![ data_dir.join("GSE66386/A.occ.wig"),
                             data_dir.join("GSE66386/B.occ.wig"),],
            work_dir: work_dir.clone(),
        },
        gene_config: GeneConfig {
            mtif_filename: data_dir.join("S2_tcd_mTIFAnno2.txt"),
            cds_filename: data_dir.join("sacCer3.bed"),
            dubious_filename: data_dir.join("dubious.txt"),
            work_dir: work_dir.clone(),
        },
        target_choice_config: TargetChoiceConfig {
            tss_upstream: 220,
            tss_downstream: 20,
            tss_zones: vec![-220..-141,-140..-61,-60..20],
            cds_upstream: 350,
            cds_downstream: 0,
            cds_zones: vec![-350..-271,-270..-191,-190..-111,-110..-30],
            num_per_gene: 10,
            work_dir: work_dir.clone(),
        },
        output_dir: base_dir.join("output"),
        target_spec: TargetSpec {
            pam: b"GG".to_vec(),
            pam_start: 21,
            target_len: 23,
            guide_start: 0,
            guide_len: 20,
        },
        oligo_5ext: "tctgggagctgcgattggca".to_owned(),
        oligo_3ext: "gttttagagctagaaatagc".to_owned(),
        quick_test: None,
    };

    match run(config) {
        Ok(_) => (),
        Err(e) => panic!("Fatal error: {:?}", e),
    }
}

fn run(config: Config) -> Result<(), MyErr> {
    fs::create_dir_all(&config.output_dir)
        .annot_err_by(|| format!("Creating output directory {:?}", &config.output_dir))?;

    fs::create_dir_all(&config.work_dir)
        .annot_err_by(|| format!("Creating work directory {:?}", &config.work_dir))?;

    let targets = genomic_targets(&config.target_spec, &config.genome_fasta)
        .annot_err_by(|| format!("Generating target list"))?;

    let annot = GenomeAnnot::from_config(&config.gene_config)
        .annot_err_by(|| format!("Reading genome annotation"))?;

    let uniqueness = target_uniqueness(&config.unique_config, targets.iter())
        .annot_err_by(|| "Analyzing uniqueness by genomic alignment".to_string())?;
    
    let atacs = target_average_atacs(&config.atac_config, targets.iter())
        .annot_err_by(|| "Annotating targets with average ATAC values".to_string())?;
    
    write_bed_files(&config, &targets, &uniqueness, &atacs)?;

    let info = ScoringInfo::new(&config.target_choice_config, 
                                annot, uniqueness, atacs,
                                targets)?;

    let gene_choices = choose_targets(&config.target_choice_config, &info)?;

    let mut choice_writer = fs::File::create(config.output_dir.join("target-choices.txt"))?;
    let mut oligo_writer = fs::File::create(config.output_dir.join("target-oligos.txt"))?;
    let mut choice_bed = bed::Writer::to_file(config.output_dir.join("target-choices.bed"))?;

    write!(choice_writer, "Yorf\tGuideNo\tTargetSeq\tTargetLoc\tOffset\tRank\tZone\tUnique\tSpecific\tATAC\tYorfStart\tYorfLoc\n")?;
    write!(oligo_writer, "Yorf_GNo\tTargetLoc\tOligo\n")?;
    for choices in gene_choices.values() {
        for choice in choices {
            let sequ_str = str::from_utf8(choice.target_score().target().sequ())?;
            let start_str = info.annot().gene_start(choice.yorf()).map_or_else(|| "N/A".to_owned(), |gs| gs.to_string());
            let pos_str = info.annot().gene_start(choice.yorf()).map_or_else(|| "N/A".to_owned(), |gs| gs.first_pos().to_string());
            write!(choice_writer, "{}\t{:02}\t{}\t{}\t{:+}\t{}\t{}\t{}\t{}\t{:0.3}\t{}\t{}\n",
                   choice.yorf(), choice.guide_no(), sequ_str, choice.target_score().target().loc(), choice.target_score().offset(),
                   choice.rank(), choice.zone().map_or_else(|| "N/A".to_owned(), |z| z.to_string()),
                   choice.target_score().unique(), choice.target_score().specific(), choice.target_score().atac(),
                   start_str, pos_str)?;

            let mut bed = choice.target_score().target().loc().to_bed();
            bed.set_name(&choice.guide_desc());
            choice_bed.write(&bed)?;

            let mut oligo = config.oligo_5ext.clone();
            oligo += str::from_utf8(&config.target_spec.guide_sequence(choice.target_score().target()))?;
            oligo += config.oligo_3ext.as_str();

            write!(oligo_writer, "{}_{:02}\t{}\t{}\n",
                   choice.yorf(), choice.guide_no(),
                   choice.target_score().target().loc(),
                   oligo)?;
        }
    }
    
    Ok(())
}

fn write_bed_files(config: &Config,
                   targets: &Vec<Target>,
                   uniques: &HashMap<SeqContigLoc,bool>,
                   _atacs: &HashMap<SeqContigLoc,f64>)
                   -> Result<(),MyErr>
{
    let mut target_bed_writer = bed::Writer::to_file(config.work_dir.join("targets.bed"))?;
    let mut unique_bed_writer = bed::Writer::to_file(config.work_dir.join("targets-uniqueness.bed"))?;
    for target in targets.iter() {
        let mut bed = target.loc().to_bed();
        target_bed_writer.write(&bed)?;

        bed.set_score(match uniques.get(target.loc()) {
            Some(&true)  => "900",
            Some(&false) => "500",
            None        => "100"
        });
        unique_bed_writer.write(&bed)?;
    }    
    Ok( () )
}


