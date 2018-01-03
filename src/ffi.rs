use std::ffi::CStr;
use std::mem;
use std::os::raw::c_char;
use ::{_hash_murmur, KmerMinHash};

use std::panic;
use std::thread;
use std::cell::RefCell;

use backtrace::Backtrace;

use errors::{ErrorKind, Error, Result, SourmashErrorCode};

thread_local! {
    pub static LAST_ERROR: RefCell<Option<Error>> = RefCell::new(None);
    pub static LAST_BACKTRACE: RefCell<Option<(Option<String>, Backtrace)>> = RefCell::new(None);
}

macro_rules! ffi_fn (
    // a function that catches panics and returns a result (err goes to tls)
    (
        $(#[$attr:meta])*
        unsafe fn $name:ident($($aname:ident: $aty:ty),* $(,)*) -> Result<$rv:ty> $body:block
    ) => (
        #[no_mangle]
        $(#[$attr])*
        pub unsafe extern "C" fn $name($($aname: $aty,)*) -> $rv
        {
            $crate::ffi::landingpad(|| $body)
        }
    );

    // a function that catches panics and returns nothing (err goes to tls)
    (
        $(#[$attr:meta])*
        unsafe fn $name:ident($($aname:ident: $aty:ty),* $(,)*) $body:block
    ) => {
        #[no_mangle]
        $(#[$attr])*
        pub unsafe extern "C" fn $name($($aname: $aty,)*)
        {
            // this silences panics and stuff
            $crate::ffi::landingpad(|| { $body; Ok(0 as ::std::os::raw::c_int) });
        }
    }
);

#[no_mangle]
pub extern "C" fn hash_murmur(kmer: *const c_char, seed: u64) -> u64 {
    let c_str = unsafe {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    _hash_murmur(c_str.to_bytes(), seed)
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_new(n: u32, k: u32, prot: bool,
                                         seed: u64, mx: u64, track_abundance: bool)
                                         -> *mut KmerMinHash {
    mem::transmute(Box::new(KmerMinHash::new(n, k, prot, seed, mx, track_abundance)))
}

#[no_mangle]
pub extern "C" fn kmerminhash_free(ptr: *mut KmerMinHash) {
    if ptr.is_null() { return }
    unsafe { Box::from_raw(ptr); }
}

ffi_fn! {
unsafe fn kmerminhash_add_sequence(ptr: *mut KmerMinHash, sequence: *const c_char, force: bool) ->
    Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    mh.add_sequence(c_str.to_bytes(), force)
}
}

#[no_mangle]
pub extern "C" fn kmerminhash_add_hash(ptr: *mut KmerMinHash, h: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    mh.add_hash(h);
}

#[no_mangle]
pub extern "C" fn kmerminhash_add_word(ptr: *mut KmerMinHash, word: *const c_char) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = unsafe {
        assert!(!word.is_null());

        CStr::from_ptr(word)
    };

    mh.add_word(c_str.to_bytes());
}

#[no_mangle]
pub extern "C" fn kmerminhash_get_min_idx(ptr: *mut KmerMinHash, idx: u64) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins[idx as usize]
}

#[no_mangle]
pub extern "C" fn kmerminhash_get_mins_size(ptr: *mut KmerMinHash) -> usize {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins.len()
}

#[no_mangle]
pub extern "C" fn kmerminhash_is_protein(ptr: *mut KmerMinHash) -> bool {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.is_protein
}

#[no_mangle]
pub extern "C" fn kmerminhash_seed(ptr: *mut KmerMinHash) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.seed
}

#[no_mangle]
pub extern "C" fn kmerminhash_num(ptr: *mut KmerMinHash) -> u32 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.num
}

#[no_mangle]
pub extern "C" fn kmerminhash_ksize(ptr: *mut KmerMinHash) -> u32 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.ksize
}

#[no_mangle]
pub extern "C" fn kmerminhash_max_hash(ptr: *mut KmerMinHash) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.max_hash
}

#[no_mangle]
pub extern "C" fn kmerminhash_merge(ptr: *mut KmerMinHash, other: *const KmerMinHash) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = unsafe {
       assert!(!other.is_null());
       &*other
    };

    let merged = mh.merge(other_mh);
    mh.mins = merged
}

#[no_mangle]
pub extern "C" fn kmerminhash_count_common(ptr: *mut KmerMinHash, other: *const KmerMinHash) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = unsafe {
       assert!(!other.is_null());
       &*other
    };

    mh.count_common(other_mh)
}

fn notify_err(err: Error) {
    if let Some(backtrace) = err.backtrace() {
        LAST_BACKTRACE.with(|e| {
            *e.borrow_mut() = Some((None, backtrace.clone()));
        });
    }
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = Some(err);
    });
}

/// Clears the last error.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_clear() {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = None;
    });
    LAST_BACKTRACE.with(|e| {
        *e.borrow_mut() = None;
    });
}

/// Initializes the library
#[no_mangle]
pub unsafe extern "C" fn sourmash_init() {
    set_panic_hook();
}

/// Returns the last error code.
///
/// If there is no error, 0 is returned.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_get_last_code() -> SourmashErrorCode {
    LAST_ERROR.with(|e| {
        if let Some(ref err) = *e.borrow() {
            SourmashErrorCode::from_kind(err.kind())
        } else {
            SourmashErrorCode::NoError
        }
    })
}

pub unsafe fn set_panic_hook() {
    panic::set_hook(Box::new(|info| {
        let backtrace = Backtrace::new();
        let thread = thread::current();
        let thread = thread.name().unwrap_or("unnamed");

        let msg = match info.payload().downcast_ref::<&str>() {
            Some(s) => *s,
            None => {
                match info.payload().downcast_ref::<String>() {
                    Some(s) => &**s,
                    None => "Box<Any>",
                }
            }
        };

        let panic_info = match info.location() {
            Some(location) => {
                format!("thread '{}' panicked with '{}' at {}:{}",
                                     thread, msg, location.file(),
                                     location.line())
            }
            None => {
                format!("thread '{}' panicked with '{}'", thread, msg)
            }
        };

        LAST_BACKTRACE.with(|e| {
            *e.borrow_mut() = Some((Some(panic_info), backtrace));
        });
    }));
}

pub unsafe fn landingpad<F: FnOnce() -> Result<T> + panic::UnwindSafe, T>(
    f: F) -> T
{
    match panic::catch_unwind(f) {
        Ok(rv) => rv.map_err(|err| notify_err(err)).unwrap_or(mem::zeroed()),
        Err(err) => {
            use std::any::Any;
            let err = &*err as &Any;
            let msg = match err.downcast_ref::<&str>() {
                Some(s) => *s,
                None => {
                    match err.downcast_ref::<String>() {
                        Some(s) => &**s,
                        None => "Box<Any>",
                    }
                }
            };
            notify_err(ErrorKind::Panic(msg.to_string()).into());
            mem::zeroed()
        }
    }
}
