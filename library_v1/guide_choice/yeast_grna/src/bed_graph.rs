use std::collections::HashMap;
use std::collections::hash_map::Keys;
use std::convert::From;
use std::io::{BufRead};
use std::ops::Range;
use std::str::FromStr;

use myerr::MyErr;

pub struct BedGraph<T> {
    chrom_graphs: HashMap<String, Vec<Option<T>>>,
}

impl <T> BedGraph<T>
    where T: FromStr + Clone,
          MyErr: From<<T as FromStr>::Err>
{
    pub fn read<R: BufRead>(reader: R) -> Result<Self,MyErr> {
        let mut chrom_graphs = HashMap::new();
        
        for line_res in reader.lines() {
            let line = line_res?;

            if !line.starts_with("#") {
                let fields: Vec<&str> = line.split('\t').collect();
                
                let chrom = fields.get(0).map_or( Err( MyErr::myerr("bed_graph: malformed line") ), Ok )?;
                let startstr = fields.get(1).map_or( Err( MyErr::myerr("bed_graph: malformed line") ), Ok )?;
                let start = startstr.parse::<usize>()?;
                let endstr = fields.get(2).map_or( Err( MyErr::myerr("bed_graph: malformed line") ), Ok )?;
                let end = endstr.parse::<usize>()?;

                let valstr = fields.get(3).map_or( Err( MyErr::myerr("bed_graph: malformed line") ), Ok )?;
                let val = valstr.parse::<T>()?;

                let chrom_graph = chrom_graphs.entry(String::from(*chrom)).or_insert_with(|| Vec::new());

                chrom_graph.resize(end, None);
                
                for pos in start..end {
                    if let Some(elem) = chrom_graph.get_mut(pos) {
                        *elem = Some(val.clone());
                    }
                }
            }
        }

        Ok( BedGraph{ chrom_graphs: chrom_graphs } )
    }

    pub fn _chroms(&self) -> Keys<String, Vec<Option<T>>> {
        self.chrom_graphs.keys()
    }

    pub fn get_chrom(&self, chrom: &str) -> Option<&Vec<Option<T>>> {
        self.chrom_graphs.get(chrom)
    }
    
    pub fn _get_value(&self, chrom: &str, index: usize) -> Option<&T>
    {
        self.get_chrom(chrom)
            .map_or(None, |v| v.get(index).map_or(None, |x| x.as_ref()))
    }

    pub fn get_value_slice(&self, chrom: &str, index: Range<usize>) -> Option<&[Option<T>]>
    {
        self.get_chrom(chrom).map_or(None, |v| v.get(index))
    }
}

