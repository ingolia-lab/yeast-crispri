extern crate bio;

use std::collections::HashMap;
use std::fs;                                                                                                                                                                      
use std::io;
use std::io::{Write};
use std::path;

use bed_graph::BedGraph;
use myerr::MyErr;
use seq_loc::*;
use target::Target;

#[derive(Debug)]
pub struct AtacConfig {
    pub atac_wigs: Vec<path::PathBuf>,
    pub work_dir: path::PathBuf
}

pub fn target_average_atacs<'a, I>(config: &AtacConfig, 
                                   targets: I) 
                                   -> Result<HashMap<SeqContigLoc,f64>, MyErr>
    where I: Iterator<Item=&'a Target>    
{
    fs::create_dir_all(&config.work_dir)?;

    let bedgraphs = read_atacs(&config.atac_wigs)?;

    let mut atac_writer = fs::File::create(config.work_dir.join("target-atacs.txt"))?;
    let mut atac_map = HashMap::new();
    
    for target in targets {
        let mut atacs = Vec::new();

        for bedgraph in bedgraphs.iter() {
            atacs.push(target_average(&bedgraph, &target));
        }

        let avg = atacs_average(&atacs);
        
        atac_writer.write(&target.to_tsv())?;
        write!(atac_writer, "\t{:0.3}", avg)?;
        
        for atac in atacs.iter() {
            match *atac {
                Some(ref x) => write!(atac_writer, "\t{:0.3}", x),
                None => write!(atac_writer, "\tN/A"),
            }?;
        }
        atac_writer.write(b"\n")?;

        atac_map.insert(target.loc().clone(), avg);
    }

    Ok( atac_map )
}

fn read_atacs(atac_wigs: &Vec<path::PathBuf>) -> Result<Vec<BedGraph<f64>>, MyErr> {
    let mut atacs = Vec::new();

    for atac_wig in atac_wigs {
        let atac_file = fs::File::open(atac_wig)?;
        let atac = BedGraph::read(io::BufReader::new(atac_file))?;
        atacs.push(atac);
    }
    
    Ok( atacs )
}

fn target_average(track: &BedGraph<f64>, target: &Target) -> Option<f64> {
    let vals = track.get_value_slice(&target.loc().refid(), target.loc().range());

    vals.map_or(None, average_somes)
}

fn average_somes(opt_values: &[Option<f64>]) -> Option<f64> {
    let nttl = opt_values.len();
    let nsome = opt_values.iter().filter_map(|x| *x).count();
    if (nsome * 2 >= nttl) && (nsome > 0) {
        let ttl: f64 = opt_values.iter().filter_map(|x| *x).sum();
        Some( ttl / (nsome as f64) )
    } else {
        None
    }
}

fn atacs_average(atacs: &Vec<Option<f64>>) -> f64 {
    let n = atacs.len() as f64;
    let ttl: f64 = atacs.iter().map(|x| x.unwrap_or(0.0)).sum();
    ttl / n
}

