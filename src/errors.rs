use std::io;
use std::str;

error_chain! {
    errors {
        /// Raised in some cases if panics are caught
        Panic(message: String) {
            description("panic")
            display("panic: {}", message)
        }
        /// Raised for internal errors in the libraries.  Should not happen.
        Internal(message: &'static str) {
            description("internal error")
            display("internal error: {}", message)
        }
        MismatchKSizes {
            description("different ksizes cannot be compared")
        }
        MismatchDNAProt {
            description("DNA/prot minhashes cannot be compared")
        }
        MismatchMaxHash {
            description("mismatch in max_hash; comparison fail")
        }
        MismatchSeed {
            description("mismatch in seed; comparison fail")
        }
        InvalidDNA(message: String) {
            description("Invalid DNA character in input")
            display("invalid DNA character in input: {}", message)
        }
        InvalidProt(message: String) {
            description("Invalid protein character in input")
            display("invalid protein character in input: {}", message)
        }
    }

    foreign_links {
        Io(io::Error);
        Utf8Error(str::Utf8Error);
        ParseInt(::std::num::ParseIntError);
    }
}

#[repr(u32)]
pub enum SourmashErrorCode {
    // no error
    NoError = 0,
    // panics and internals
    Panic = 1,
    Internal = 2,
    Msg = 3,
    Unknown = 4,
    // Compatibility errors
    MismatchKSizes = 101,
    MismatchDNAProt = 102,
    MismatchMaxHash = 103,
    MismatchSeed = 104,
    // Input sequence errors
    InvalidDNA = 1101,
    InvalidProt = 1102,
    // external errors
    Io = 100001,
    Utf8Error = 100002,
    ParseInt = 100003,
}

impl SourmashErrorCode {
    pub fn from_kind(kind: &ErrorKind) -> SourmashErrorCode {
        match *kind {
            ErrorKind::Panic(..) => SourmashErrorCode::Panic,
            ErrorKind::Msg(..) => SourmashErrorCode::Msg,
            ErrorKind::Internal(..) => SourmashErrorCode::Internal,
            ErrorKind::MismatchKSizes => SourmashErrorCode::MismatchKSizes,
            ErrorKind::MismatchDNAProt => SourmashErrorCode::MismatchDNAProt,
            ErrorKind::MismatchMaxHash => SourmashErrorCode::MismatchMaxHash,
            ErrorKind::MismatchSeed => SourmashErrorCode::MismatchSeed,
            ErrorKind::InvalidDNA(..) => SourmashErrorCode::InvalidDNA,
            ErrorKind::InvalidProt(..) => SourmashErrorCode::InvalidProt,
            ErrorKind::Io(..) => SourmashErrorCode::Io,
            ErrorKind::Utf8Error(..) => SourmashErrorCode::Utf8Error,
            ErrorKind::ParseInt(..) => SourmashErrorCode::ParseInt,
            ErrorKind::__Nonexhaustive { .. } => unreachable!(),
        }
    }
}
