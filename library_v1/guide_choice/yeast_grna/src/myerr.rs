extern crate rust_htslib;
extern crate csv;

use std::convert;
use std::io;
use std::num;
use std::str;
use std::string;

use seq_loc;

#[derive(Debug)]
pub enum MyErr {
    IO(io::Error),
    Utf8(str::Utf8Error),
    Str(String),
    ParseInt(num::ParseIntError),
    ParseFloat(num::ParseFloatError),
    BamRead(rust_htslib::bam::ReadError),
    Annot(Box<MyErr>, String),
    CSV(csv::Error),
    ParseSeqLoc(seq_loc::ParseError),
}

impl MyErr {
    pub fn myerr(error: &str) -> MyErr {
        MyErr::Str( String::from(error) )
    }
}

impl convert::From<io::Error> for MyErr {
    fn from(error: io::Error) -> Self { MyErr::IO(error) }
}

impl convert::From<str::Utf8Error> for MyErr {
    fn from(error: str::Utf8Error) -> Self {
        MyErr::Utf8(error)
    }
}

impl convert::From<string::FromUtf8Error> for MyErr {
    fn from(error: string::FromUtf8Error) -> Self {
        MyErr::Utf8(error.utf8_error())
    }
}

impl convert::From<String> for MyErr {
    fn from(error: String) -> Self { MyErr::Str(error) }
}

impl convert::From<num::ParseIntError> for MyErr {
    fn from(error: num::ParseIntError) -> Self {
        MyErr::ParseInt(error)
    }
}

impl convert::From<num::ParseFloatError> for MyErr {
    fn from(error: num::ParseFloatError) -> Self {
        MyErr::ParseFloat(error)
    }
}

impl convert::From<rust_htslib::bam::ReadError> for MyErr {
    fn from(error: rust_htslib::bam::ReadError) -> Self {
        MyErr::BamRead(error)
    }
}

impl convert::From<csv::Error> for MyErr {
    fn from(error: csv::Error) -> Self {
        MyErr::CSV(error)
    }
}

impl convert::From<seq_loc::ParseError> for MyErr {
    fn from(error: seq_loc::ParseError) -> Self {
        MyErr::ParseSeqLoc(error)
    }
}

pub trait Annotable<T,E: Into<MyErr>> {
    fn annot_err(self, annot: String) -> Result<T, MyErr>;
    fn annot_err_by<F>(self, annot: F) -> Result<T, MyErr>
        where F: FnOnce() -> String;
}

impl <T, E: Into<MyErr>> Annotable<T,E> for Result<T,E> {
    fn annot_err(self, annot:String) -> Result<T, MyErr> {
        self.map_err(|e| MyErr::Annot(Box::new(e.into()), annot))
    }

    fn annot_err_by<F>(self, annot: F) -> Result<T, MyErr>
        where F: FnOnce() -> String
    {
        self.map_err(|e| MyErr::Annot(Box::new(e.into()), annot()))
    }
}
