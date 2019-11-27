extern crate bio;
extern crate rust_htslib;
extern crate itertools;

use std::collections::HashMap;
use std::fs;                                                                                                                                                                      
use std::io::{Write};
use std::path;
use std::process;

use itertools::Itertools;

use bio::io::fasta;

use rust_htslib::bam::Read;

use myerr::{MyErr,Annotable};
use seq_loc::{SeqContigLoc,Strand};
use target::Target;

#[derive(Debug)]
pub struct UniqueConfig {
    pub genome_bt2: String,
    pub work_dir: path::PathBuf,
}

pub fn target_uniqueness<'a, I>(config: &UniqueConfig, targets: I) -> Result<HashMap<SeqContigLoc,bool>,MyErr>
    where I: Iterator<Item=&'a Target>
{
    fs::create_dir_all(&config.work_dir)?;

    let target_path = create_target_fasta(&config.work_dir, targets)
        .annot_err_by(|| format!("Creating Fasta file of target sequences"))?;

    let bam_path = align_targets(&config, &target_path)
        .annot_err_by(|| format!("Aligning target sequences to genome using bowtie"))?;

    bam_uniqueness(&config.work_dir, &bam_path)
        .annot_err("Picking unique targets from genomic alignments".to_string())
}

fn bam_uniqueness(work_dir: &path::PathBuf, bam_path: &path::Path) -> Result<HashMap<SeqContigLoc,bool>, MyErr> {
    let mut bam_reader = rust_htslib::bam::Reader::from_path(bam_path)
        .map_err(|e| MyErr::Str(format!("{}", e)))?;

    let results: Result<Vec<String>, _> = bam_reader.header().target_names().iter()
        .map(|u8| String::from_utf8(u8.to_vec())).collect();
    let bam_target_names = results?;

    let query_groups = bam_reader.records().group_by(|res| {
        match *res {
            Ok(ref rec) => Ok(rec.qname().to_vec()),
            Err(ref e) => Err(format!("{}", e)),
        }
    });

    let hit_count = fs::File::create(work_dir.join("targets-hit-counts.txt"))?;

    let mut uniq_map = HashMap::new();
    
    for (qres, qgroup) in query_groups.into_iter() {
        let qname_u8 = qres?;
        let qname = String::from_utf8(qname_u8.to_vec())?;

        let qrecords_res: Result<Vec<rust_htslib::bam::Record>,rust_htslib::bam::ReadError>
            = qgroup.collect();
        let qrecords = qrecords_res?;

        writeln!(&hit_count, "{}\t{}", qname, qrecords.len())?;

        if qrecords.len() == 0 {
            return Err( MyErr::Str(format!("Zero hits for {}", qname)) );
        } else {            
            let record = qrecords.get(0)
                .map_or(Err( MyErr::Str("Unexpected failure getting singleton BAM record".to_string()) ), Ok )?;

            let qtarg = Target::from_bam(record)?;
            
            if qrecords.len() == 1 && !record.is_unmapped() {
                let bam_target_name = bam_target_names.get(record.tid() as usize)
                    .map_or(Err(MyErr::Str("No target for TID".to_string())), Ok)?;
                
                if (qtarg.loc().refid() != bam_target_name)
                    || (qtarg.loc().start() != record.pos() as usize)
                    || (record.is_reverse() && qtarg.loc().strand() != Strand::Reverse)
                    || (!record.is_reverse() && qtarg.loc().strand() != Strand::Forward) {
                        return Err( MyErr::Str("Singleton hit not good ".to_string() + &qname) );
                    }
            }
            
            uniq_map.insert(qtarg.loc().clone(), qrecords.len() == 1 && !record.is_unmapped());
        }
        
    }
    
    Ok( uniq_map )
}

fn align_targets(config: &UniqueConfig, target_path: &path::Path) -> Result<path::PathBuf, MyErr> {
    let bowtie_err = fs::File::create(config.work_dir.join("bowtie-output.txt"))?;
    let bowtie_sam_path = target_path.with_extension("sam");

    let target_fa = target_path.to_str().map_or(Err(MyErr::Str("Bad target .fa path".to_string())), Ok)?;
    let bowtie_sam = bowtie_sam_path.to_str().map_or(Err(MyErr::Str("Bad .sam path".to_string())), Ok)?;
    
    let mut bowtie = process::Command::new("bowtie2")
        .args(&["-p36", "-a", "-f", "--np", "0", "-L", "20",
                "-x", &config.genome_bt2,
                "-U", target_fa,
                "-S", bowtie_sam])
        .stderr(bowtie_err)
        .spawn()?;

    let bowtie_exit = bowtie.wait()?;

    if !bowtie_exit.success() {
        return Err(MyErr::Str("bowtie exited with error".to_string()));
    }

    let bowtie_bam_path = target_path.with_extension("bam");
    let bowtie_bam_path_return = bowtie_bam_path.clone();
    let bowtie_bam = bowtie_bam_path.to_str().map_or(Err(MyErr::Str("Bad .bam path".to_string())), Ok)?;
    let mut samtools = process::Command::new("samtools")
        .args(&["view", "-b", "-S", "-o", bowtie_bam, bowtie_sam])
        .spawn()?;
    let samtools_exit = samtools.wait()?;

    if !samtools_exit.success() {
        return Err(MyErr::Str("samtools exited with error".to_string()));
    }

    Ok(bowtie_bam_path_return)
}

fn create_target_fasta<'a, I>(work_dir: &path::PathBuf, targets: I) -> Result<path::PathBuf, MyErr>
    where I: Iterator<Item=&'a Target>
{
    let target_pathbuf = work_dir.join("all-pam-targets.fa");
    let mut target_writer = fasta::Writer::to_file(target_pathbuf.as_path())
        .annot_err(format!("Opening {:?} to write", target_pathbuf))?;
    
    for target in targets {
        target_writer.write_record(&target.record())?;
    }

    Ok(target_pathbuf)
}

