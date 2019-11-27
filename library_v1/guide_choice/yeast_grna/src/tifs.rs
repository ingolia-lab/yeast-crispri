use std::collections::HashMap;
use std::convert::From;
use std::fmt;
use std::fs;
use std::io::{BufRead,BufReader,Write};
use std::ops::Deref;
use std::path;
use std::rc::Rc;
use std::vec;

use myerr::MyErr;
use seq_loc::*;

 //   4140 Covering >=2 ORFs
 // 184473 Covering one intact ORF
 //   1675 CUT/XUT
 //  16824 Intergenic transcripts
 //   3681 Internal transcripts
 //   4096 Overlap >=2 ORFs
 // 135806 Overlap 3' of one ORF
 //  16761 Overlap 5' of one ORF
 //   3631 SUT
 //      1 type

#[derive(Debug,Clone,PartialEq,Eq)]
pub enum TifType {
    CoveringMulti,
    CoveringOne,
    CutXut,
    Intergenic,
    Internal,
    OverlapMulti,
    Overlap3,
    Overlap5,
    Sut
}

impl TifType {

    pub fn from_str(name: &str) -> Option<Self> {
        match name {
            "Covering >=2 ORFs"       => Some(TifType::CoveringMulti),
            "Covering one intact ORF" => Some(TifType::CoveringOne),
            "CUT/XUT"                 => Some(TifType::CutXut),
            "Intergenic transcripts"  => Some(TifType::Intergenic),
            "Internal transcripts"    => Some(TifType::Internal),
            "Overlap >=2 ORFs"        => Some(TifType::OverlapMulti),
            "Overlap 3' of one ORF"   => Some(TifType::Overlap3),
            "Overlap 5' of one ORF"   => Some(TifType::Overlap5),
            "SUT"                     => Some(TifType::Sut),
            _ => None
        }
    }    
}

#[derive(Debug)]
pub struct GenomeTssMap {
    yorf_map: HashMap<String,Rc<ModalTss>>,
    loc_map: SeqLocMap<Rc<ModalTss>>,
}

impl GenomeTssMap {
    pub fn new() -> Self {
        GenomeTssMap{ yorf_map: HashMap::new(),
                      loc_map: SeqLocMap::new(), }
    }

    pub fn build<I>(tsses: I) -> Result<Self,MyErr>
        where I: Iterator<Item = ModalTss>
    {
        let mut tss_map = Self::new();
        for tss in tsses {
            tss_map.insert(tss)?;
        }
        Ok( tss_map )
    }
    
    pub fn insert(&mut self, tss: ModalTss) -> Result<(),MyErr> {
        let yorf = String::from(tss.yorf.as_str());
        let rc = Rc::new(tss);

        self.yorf_map.insert(yorf, rc.clone())
                .map_or( Ok(()), |_| Err( MyErr::Str(format!("Duplicate TSS entries for {}", rc.deref().yorf)) ) )?;
        self.loc_map.insert(rc.clone());

        Ok( () )
    }

    pub fn tsses<'a>(&'a self) -> Vec<&ModalTss>
    {
        let values = self.yorf_map.values();
        values.map(Deref::deref).collect()
    }

    pub fn tss_by_yorf(&self, yorf: &str) -> Option<&ModalTss> {
        self.yorf_map.get(yorf).map(Deref::deref)
    }
    
    pub fn tsses_by_loc(&self, loc: &SeqContigLoc) -> Vec<&ModalTss>
    {
        self.loc_map.features_by_loc(loc)
            .into_iter().map(|e| e.deref()).collect()
    }

    pub fn read_from_mtifs(mtif_filename: &path::Path) -> Result<Self,MyErr> {
        let mtifs = Tif::read_tif_file(mtif_filename)?;
        let tsses = modal_tsses(mtifs.iter())?;
        Self::build(tsses.into_iter())
    }

    pub fn write_table(&self, table_filename: &path::Path) -> Result<(),MyErr> {
        let mut table = fs::File::create(table_filename)?;
        for (yorf, tss_rc) in self.yorf_map.iter() {
            let tss = tss_rc.deref();
            write!(table, "{}\t{}\n", yorf, tss)?;
            let tsses = self.tsses_by_loc(&tss.loc());
            if !tsses.iter().any(|t| t.yorf() == yorf) {
                Err(MyErr::Str(format!("Failed reverse lookup, {} not in {:?}", tss, tsses)))?;
            }
        }
        Ok(())
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct ModalTss {
    site: SeqPos,
    tss_count: i32,
    total_count: i32,
    yorf: String,
}

impl ModalTss {
    pub fn site(&self) -> &SeqPos { &self.site }
    pub fn yorf(&self) -> &str { &self.yorf }

    // pub fn to_tsv(&self) -> String {
    //     format!("{}\t{}\t{}\t{}", self.yorf, 
    //             self.site, self.tss_count, self.total_count)
    // }

    fn from_tifs<'a, I>(tifs_in: I) -> Result<ModalTss,MyErr>
        where I: Iterator<Item = &'a Tif>
    {
        let mut tifs = tifs_in.peekable();
        let peek_res: Result<(String, &str, bool),MyErr> = {
            let tif0 = tifs.peek()
                .map_or( Err( MyErr::Str(format!("Zero-length TIF list"))), Ok)?;
            let refid = tif0.refid.as_str();
            let reffwd = tif0.reffwd;
            let yorf = match tif0.yorf {
                Some(ref y) => Ok(String::from(y.as_str())),
                None => Err(MyErr::Str(format!("TIF list with no YORF"))), 
            }?;
            Ok( (yorf, refid, reffwd) )
        };

        let (yorf, refid, reffwd) = peek_res?;
        
        let ( tss, tss_count, total_count ) = Self::modal_tss(tifs)?;

        Ok( ModalTss{ site: SeqPos::new(refid, tss, if reffwd { Strand::Forward } else { Strand::Reverse }),
                      tss_count: tss_count, total_count: total_count,
                      yorf: yorf } )
    }
    
    fn modal_tss<'a, I>(tifs: I) -> Result<(usize,i32,i32),MyErr>
        where I: Iterator<Item = &'a Tif>
    {
        let mut tss_counts = HashMap::new();
        let mut total = 0;
        let mut tiffwd_m = None;
        for tif in tifs {
            let tss_count = tss_counts.entry(tif.t5).or_insert(0);
            *tss_count += tif.ypd + tif.gal;
            total += tif.ypd + tif.gal;
            tiffwd_m = match tiffwd_m {
                None => Some(tif.reffwd),
                Some(fwd) => { assert!(fwd == tif.reffwd); tiffwd_m },
            };
        }

        let tiffwd = tiffwd_m.ok_or_else(|| MyErr::Str(format!("No TIFs found")))?;
        let max_count = tss_counts.iter().map(|(_,count)| count).max().ok_or_else(|| MyErr::Str(format!("No modal TSS found")))?;
        
        let max_tsses: Vec<(&usize, &i32)> = tss_counts.iter().filter(|&(_,count)| count == max_count).collect();

        let earlier_tss = |&(tssleft, _): &(&usize, &i32), &(tssright, _): &(&usize, &i32)| {
            if tiffwd {
                tssright.cmp(tssleft)
            } else {
                tssleft.cmp(tssright)
            }
        };

        if max_tsses.len() > 1 {
            print!("Degenerate TSS choices:");
            for &(&tss, &count) in max_tsses.iter() {
                print!(" {};{}", tss, count);
            }
            print!("\n");
        }
        
        // Pick longest 
        match max_tsses.into_iter().max_by(earlier_tss) {
            Some( (&tssmax, &countmax) ) => Ok( (tssmax, countmax, total) ),
            None => Err( MyErr::Str(format!("No modal TSS found")) ),
        }
    }
}

impl SeqLocFeature for ModalTss {
    fn refid(&self) -> &str { self.site.refid() }
    fn start(&self) -> usize { self.site.pos() }
    fn length(&self) -> usize { 1 }
    fn strand(&self) -> Strand { self.site.strand() }    
}

impl fmt::Display for ModalTss {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}@{}", self.yorf, self.site)
    }    
}

pub fn modal_tsses<'a, I>(tifs: I) -> Result<Vec<ModalTss>,MyErr>
    where I: Iterator<Item = &'a Tif>
{
    let mut tifs_by_yorf: HashMap<String,Vec<&Tif>> = HashMap::new();
    
    for tif in tifs.filter(|tif| tif.tif_type == TifType::CoveringOne) {
        let myorf = match tif.yorf {
            Some(ref yorf) => Some(String::from(yorf.as_str())),
            None => None
        };
        
        if let Some(yorf) = myorf {
            let yorf_tifs = tifs_by_yorf.entry(String::from(yorf.as_str())).or_insert_with(|| Vec::new());
            yorf_tifs.push(tif);
        }
    }

    let mut tsses = Vec::new();
    for (_,yorf_tifs) in tifs_by_yorf.drain() {
        let foo: Vec<&Tif> = yorf_tifs;
        let bar: vec::IntoIter<&Tif> = foo.into_iter();
        let tss = ModalTss::from_tifs(bar)?;
        tsses.push(tss);
    }
    
    Ok( tsses )
}

#[derive(Debug)]
pub struct Tif {
    refid: String,
    reffwd: bool,
    t5: usize,
    t3: usize,
    ypd: i32,
    gal: i32,
    tif_type: TifType,
    yorf: Option<String>,
}

impl Tif {
    pub fn read_tif_file(filename: &path::Path) -> Result<Vec<Self>,MyErr> {
        let tif_file = fs::File::open(filename)?;
        let mut tif_lines = BufReader::new(tif_file).lines();

        let header = tif_lines.next()
            .map_or( Err(MyErr::Str("No header in mTIF file".to_string())), |x| x.map_err(MyErr::IO) )?;
        if header != "chr\tstrand\tt5\tt3\typd\tgal\ttype\tname" {
            return Err(MyErr::Str("Bad header in mTIF file: ".to_string() + &header));
        }

        let mut tifs = Vec::new();
        
        for line_res in tif_lines {
            let line = line_res?;
            let tif = Tif::from_tsv(&line)?;
            tifs.push(tif);
        }

        Ok( tifs )
    }

    pub fn from_tsv(tsv: &str) -> Result<Tif, MyErr> {
        let fields : Vec<&str> = tsv.split('\t').collect();

        if fields.len() != 8 {
            return Err( MyErr::Str("Malformed mTIF line: ".to_string() + tsv) );
        }
        
        let refid = Self::refid_remap(fields[0]).map(String::from)?;
        let reffwd = match fields[1] {
            "-" => Ok(false),
            "+" => Ok(true),
            _ => Err( MyErr::Str("Bad mTIF strand: ".to_string() + fields[1]) ),
        }?;
        let t5 = fields[2].parse::<usize>()?;
        let t3 = fields[3].parse::<usize>()?;
        let ypd = fields[4].parse::<i32>()?;
        let gal = fields[5].parse::<i32>()?;
        let tif_type = TifType::from_str(fields[6])
            .map_or( Err( MyErr::Str("Bad TIF type".to_string() + fields[6]) ), Ok )?;
        let yorf = if fields[7] == "NA" {
            None
        } else {
            Some( String::from(fields[7]) )
        };

        Ok( Tif{ refid: refid, reffwd: reffwd, t5: t5, t3: t3,
                 ypd: ypd, gal: gal, tif_type: tif_type, yorf: yorf } )
    }

    fn refid_remap(refid: &str) -> Result<&str, MyErr> {
        match refid {
            "1" => Ok("chrI"),
            "2" => Ok("chrII"),
            "3" => Ok("chrIII"),
            "4" => Ok("chrIV"),
            "5" => Ok("chrV"),
            "6" => Ok("chrVI"),
            "7" => Ok("chrVII"),
            "8" => Ok("chrVIII"),
            "9" => Ok("chrIX"),
            "10" => Ok("chrX"),
            "11" => Ok("chrXI"),
            "12" => Ok("chrXII"),
            "13" => Ok("chrXIII"),
            "14" => Ok("chrXIV"),
            "15" => Ok("chrXV"),
            "16" => Ok("chrXVI"),
            _ => Err( MyErr::Str("Unknown chromosome ".to_string() + refid) )
        }
    }
}
