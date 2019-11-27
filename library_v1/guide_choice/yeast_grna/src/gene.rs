use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::From;
use std::fmt;
use std::fs;
use std::io;
use std::io::{BufRead,BufReader,Write};
use std::ops::Deref;
use std::path;
use std::rc::Rc;

use bio::io::bed;

use myerr::{MyErr,Annotable};
use seq_loc::*;
use target::Target;
use tifs::{GenomeTssMap,ModalTss};

#[derive(Debug)]
pub struct GeneConfig {
    pub mtif_filename: path::PathBuf,
    pub cds_filename: path::PathBuf,
    pub dubious_filename: path::PathBuf,
    pub work_dir: path::PathBuf,
}

#[derive(Debug)]
pub struct GenomeAnnot {
    cdses: GenomeCdsMap,
    tsses: GenomeTssMap,
    dubious: HashSet<String>,
}

impl GenomeAnnot {
    pub fn from_config(config: &GeneConfig) -> Result<Self,MyErr>
    {
        let tsses = GenomeTssMap::read_from_mtifs(&config.mtif_filename)
            .annot_err_by(|| format!("Reading TIFs from {:?}", &config.mtif_filename))?;
        tsses.write_table(&config.work_dir.join("genome-tss-map.txt"))?;

        let cdses = GenomeCdsMap::read_from_bed(&config.cds_filename)
            .annot_err_by(|| format!("Reading CDSes from {:?}", &config.cds_filename))?;
        cdses.write_table(&config.work_dir.join("genome-cds-map.txt"))?;
        
        let dubious_file = fs::File::open(&config.dubious_filename)?;
        let dubious_result: Result<HashSet<String>,io::Error>
            = BufReader::new(dubious_file).lines().collect();
        let dubious = dubious_result?;

        let annot = GenomeAnnot { cdses: cdses, tsses: tsses, dubious: dubious };
        annot.write_table(&config.work_dir.join("genome-annotation.txt"))?;
        Ok( annot )
    }

    pub fn cdses(&self) -> &GenomeCdsMap { &self.cdses }
    pub fn tsses(&self) -> &GenomeTssMap { &self.tsses }
    pub fn dubious(&self) -> &HashSet<String> { &self.dubious }

    pub fn gene_start(&self, yorf: &str) -> Option<GeneStart> { 
        match self.tsses().tss_by_yorf(yorf) {
            Some(ref tss) => Some( GeneStart::from_tss(tss) ),
            None => match self.cdses().cds_by_yorf(yorf) {
                Some(ref cds) => Some( GeneStart::from_cds(cds) ),
                None => None,
            }
        }
    }
    
    pub fn genes<'a>(&'a self) -> Vec<&'a str> {
        let mut gene_set = HashSet::new();
        for tss in self.tsses.tsses() {
            gene_set.insert(tss.yorf());
        }
        for cds in self.cdses.cdses() {
            gene_set.insert(cds.yorf());
        }
        gene_set.into_iter().collect()
    }

    pub fn write_table(&self, table_filename: &path::Path) -> Result<(),MyErr> {
        let mut table = fs::File::create(table_filename)?;
        for gene in self.genes().iter() {
            match self.gene_start(gene) {
                Some(ref gs) => write!(table, "{}\t{}\n", gene, gs)?,
                None => Err(MyErr::Str(format!("No start for listed gene {}", gene)))?,
            }
        }
        Ok(())
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub enum GeneStart {
    TssStart(ModalTss),
    CdsStart(Cds),
}

impl GeneStart {
    pub fn from_cds(cds: &Cds) -> Self { GeneStart::CdsStart(cds.clone()) }
    pub fn from_tss(tss: &ModalTss) -> Self { GeneStart::TssStart(tss.clone()) }

    pub fn yorf(&self) -> &str {
        match *self {
            GeneStart::TssStart(ref tss) => tss.yorf(),
            GeneStart::CdsStart(ref cds) => cds.yorf(),
        }
    }

    pub fn first_pos(&self) -> SeqPos {
        match *self {
            GeneStart::TssStart(ref tss) => tss.site().clone(),
            GeneStart::CdsStart(ref cds) => cds.first_pos().clone(),
        }
    }

    pub fn target_relative(&self, target: &Target) -> isize {
        let gene_pos = self.first_pos();
        let chr_offset = (target.pos().pos() as isize) - (gene_pos.pos() as isize);

        if gene_pos.strand() == Strand::Reverse {
            -chr_offset
        } else {
            chr_offset
        }
    }

    pub fn is_cds(&self) -> bool {
        match *self {
            GeneStart::TssStart(_) => false,
            GeneStart::CdsStart(_) => true,
        }
    }

    pub fn is_tss(&self) -> bool {
        match *self {
            GeneStart::TssStart(_) => true,
            GeneStart::CdsStart(_) => false,
        }
    }
}

impl fmt::Display for GeneStart {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            GeneStart::CdsStart(ref cds) => write!(f, "CDS={}", cds),
            GeneStart::TssStart(ref tss) => write!(f, "TSS={}", tss),
        }
    }    
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct Cds {
    loc: SeqContigLoc,
    first: SeqPos,
    yorf: String,
}

impl Cds {
    pub fn loc(&self) -> &SeqContigLoc { &self.loc }
    pub fn yorf(&self) -> &str { &self.yorf }
    pub fn first_pos(&self) -> &SeqPos { &self.first }

    pub fn from_bed(bed: &bed::Record) -> Result<Self,MyErr> {
        let loc: SeqContigLoc = From::from(bed);
        if loc.strand() == Strand::Unknown {
            return Err( MyErr::Str(format!("Cds::from_bed no strand")) )
        };
        let yorf = match &bed.name() {
            &Some(ref name) => Ok( name.to_owned() ),
            &None => Err( MyErr::Str(format!("Cds::from_bed no name")) ),
        }?;
        let first = loc.first_seqpos().ok_or_else(|| MyErr::Str(format!("No first_seqpos for {} at {}", yorf, loc)))?;
        Ok( Cds { loc: loc, yorf: yorf.to_owned(), first: first } )
    }
}

impl SeqLocFeature for Cds {
    fn refid(&self) -> &str { self.loc.refid() }
    fn start(&self) -> usize { self.loc.start() }
    fn length(&self) -> usize { self.loc.length() }
    fn strand(&self) -> Strand { self.loc.strand() }
}

impl fmt::Display for Cds {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}@{}", self.yorf, self.loc)
    }    
}

#[derive(Debug)]
pub struct GenomeCdsMap {
    yorf_map: HashMap<String,Rc<Cds>>,
    loc_map: SeqLocMap<Rc<Cds>>
}

impl GenomeCdsMap {
    pub fn new() -> Self {
        GenomeCdsMap{ yorf_map: HashMap::new(),
                      loc_map: SeqLocMap::new(), }
    }

    pub fn insert(&mut self, cds: Cds) -> Result<(),MyErr> {
        let yorf = String::from(cds.yorf.as_str());
        let rc = Rc::new(cds);
        self.yorf_map.insert(yorf, rc.clone())
                .map_or( Ok(()), |_| Err( MyErr::Str(format!("Duplicate CDS entries for {}", rc.deref().yorf)) ) )?;
        self.loc_map.insert(rc.clone());

        Ok( () )
    }

    pub fn cdses<'a>(&'a self) -> Vec<&Cds>
    {
        let values = self.yorf_map.values();
        values.map(Deref::deref).collect()
    }

    pub fn cds_by_yorf(&self, yorf: &str) -> Option<&Cds> {
        self.yorf_map.get(yorf).map(Deref::deref)
    }
    
    pub fn cdses_by_loc(&self, loc: &SeqContigLoc) -> Vec<&Cds>
    {
        self.loc_map.features_by_loc(loc).into_iter().map(|e| e.deref()).collect()
    }

    pub fn read_from_bed(bed_filename: &path::Path) -> Result<Self,MyErr> {
        let mut map = Self::new();
        for result_record in bed::Reader::from_file(bed_filename)?.records() {
            let record = result_record?;
            let cds = Cds::from_bed(&record)?;
            map.insert(cds)?;
        }
        Ok( map )
    }

    pub fn write_table(&self, table_filename: &path::Path) -> Result<(),MyErr> {
        let mut table = fs::File::create(table_filename)?;
        for (yorf, cds_rc) in self.yorf_map.iter() {
            let cds = cds_rc.deref();
            write!(table, "{}\t{}\n", yorf, cds)?;
            let cdses = self.cdses_by_loc(&cds.loc());
            if !cdses.iter().any(|t| t.yorf() == yorf) {
                Err(MyErr::Str(format!("Failed reverse lookup, {} not in {:?}", cds, cdses)))?;
            }
        }
        Ok(())
    }
}
