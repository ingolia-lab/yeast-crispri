use std::collections::{BTreeMap,HashMap};
use std::fs::File;
use std::io;
use std::io::{BufRead,BufReader,Read,Write};

use bio_types::strand::ReqStrand;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use clap::{App, Arg};
use csv;
use failure;
use failure::bail;

#[derive(Debug)]
pub struct Config {
    pub odm_bedgraph: String,
    pub seqloc_table: String,
    pub output: String,
    pub flank: usize,
    pub max_flank: usize,
}

fn main() {
    let matches = App::new("odm-seq")
        .version("0.1")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("ODM occupancy for target locations")
        .arg(
            Arg::with_name("bedgraph")
                .short("b")
                .long("bedgraph")
                .value_name("ODM.BEDGRAPH")
                .help("BedGraph of ODM occupancy values")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("locations")
                .short("l")
                .long("locations")
                .value_name("SEQLOC.TXT")
                .help("Text file of locations for analysis")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT.TXT")
                .help("Output filename")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let config = Config {
        odm_bedgraph: matches.value_of("bedgraph").unwrap().to_string(),
        seqloc_table: matches.value_of("locations").unwrap().to_string(),
        output: matches.value_of("output").unwrap().to_string(),
        flank: 5,
        max_flank: 40,
    };

    match odm_seq(&config) {
        Ok(_) => (),
        Err(e) => write!(std::io::stderr(), "ERROR {:?}", e).unwrap(),
    }
}

fn odm_seq(config: &Config) -> Result<(), failure::Error>
{
    let odm = ODM::new(&config.odm_bedgraph)?;

    let seqloc_reader: Box<dyn Read> = if config.seqloc_table == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.seqloc_table)?)
    };
    
    let mut writer: Box<dyn Write> = if config.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.output)?)
    };

    for loc_line_res in BufReader::new(seqloc_reader).lines() {
        let loc_str = loc_line_res?;
        let loc: Contig<String, ReqStrand> = loc_str.parse()?;

        let stats = if let Some(nearby_odm) = NearbyODM::nearby_nonempty(&odm, &loc, config.flank, config.max_flank)? {
            let nvals = nearby_odm.vals.len();
            let avg = nearby_odm.vals.iter().cloned().fold(0.0, std::ops::Add::add) / (nvals as f64);
            let low = nearby_odm.vals.iter().cloned().fold(std::f64::MAX, f64::min);
            format!("{}\t{}\t{}\t{:.02}\t{:.02}", nearby_odm.odm_loc, nearby_odm.flank, nvals, avg, low)
        } else {
            format!("NA\tNA\t0\tNA\tNA")
        };

        write!(writer, "{}\t{}\n", loc_str, stats)?;
    }
    
    Ok(())
}

struct NearbyODM {
    odm_loc: Contig<String, ReqStrand>,
    flank: usize,
    vals: Vec<f64>,
}

impl NearbyODM {
    fn new(odm: &ODM, loc: &Contig<String, ReqStrand>, flank: usize) -> Result<NearbyODM, failure::Error> {
        let odm_loc = Contig::new(loc.refid().clone(),
                                  loc.start() - (flank as isize),
                                  loc.length() + flank * 2,
                                  loc.strand());
        let vals = odm.get(&odm_loc)?;
        Ok(NearbyODM { odm_loc: odm_loc, flank: flank, vals: vals })
    }

    fn nearby_nonempty(odm: &ODM, loc0: &Contig<String, ReqStrand>, flank0: usize, max_flank: usize) -> Result<Option<NearbyODM>, failure::Error> {
        let mut flank = flank0;

        while flank <= max_flank {
            let nearby_odm = Self::new(odm, loc0, flank)?;
            if nearby_odm.vals.len() > 0 {
                return Ok(Some(nearby_odm))
            }

            flank *= 2;
        }

        Ok(None)
    }
}

type BedGraphRec = (String, u64, u64, f64);

struct ODM {
    seqs: HashMap<String, SeqODM>,
}

impl ODM {
    fn new(filename: &str) -> Result<Self, failure::Error> {
        let mut odm = ODM { seqs: HashMap::new() };
        
        let mut bgin = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(File::open(filename)?);

        for recres in bgin.deserialize() {
            let rec: BedGraphRec = recres?;
            odm.add_rec(rec)?;
        }

        Ok(odm)
    }

    fn add_rec(&mut self, rec: BedGraphRec) -> Result<(), failure::Error> {
        self.seqs.entry(rec.0).or_insert_with(|| SeqODM::new()).add_at(rec.1, rec.3)
    }

    fn get<S>(&self, loc: &Contig<String, S>) -> Result<Vec<f64>, failure::Error> {
        match self.seqs.get(loc.refid()) {
            Some(seqodm) => seqodm.get(loc.start(), loc.length()),
            None => bail!("No sequence {} in bedgraph", loc.refid()),
        }
    }
}

struct SeqODM {
    odms: BTreeMap<isize, f64>,
}

impl SeqODM {
    fn new() -> Self {
        SeqODM { odms: BTreeMap::new() }
    }

    fn add_at(&mut self, pos: u64, val: f64) -> Result<(), failure::Error> {
        self.odms.insert(pos as isize, val);
        Ok(())
    }

    fn get(&self, start: isize, length: usize) -> Result<Vec<f64>, failure::Error> {
        let vals: Vec<f64> = self.odms.range(start..(start + (length as isize))).map(|(_pos, val)| *val).collect();
        Ok(vals)
    }
}
