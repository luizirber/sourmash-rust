use std::ffi::CStr;
use std::mem;
use std::os::raw::c_char;
use {_hash_murmur, KmerMinHash};

#[no_mangle]
pub extern "C" fn hash_murmur(kmer: *const c_char, seed: u64) -> u64 {
    let c_str = unsafe {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    _hash_murmur(c_str.to_bytes(), seed)
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_new(
    n: u32,
    k: u32,
    prot: bool,
    seed: u64,
    mx: u64,
    track_abundance: bool,
) -> *mut KmerMinHash {
    mem::transmute(Box::new(KmerMinHash::new(
        n,
        k,
        prot,
        seed,
        mx,
        track_abundance,
    )))
}

#[no_mangle]
pub extern "C" fn kmerminhash_free(ptr: *mut KmerMinHash) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr);
    }
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

ffi_fn! {
unsafe fn kmerminhash_get_mins(ptr: *mut KmerMinHash) -> Result<*const u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let output = mh.mins.clone();
    let ptr = output.as_ptr();
    mem::forget(output);
    Ok(ptr)
}
}

ffi_fn! {
unsafe fn kmerminhash_get_min_idx(ptr: *mut KmerMinHash, idx: u64) -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    Ok(mh.mins[idx as usize])
}
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
pub extern "C" fn kmerminhash_mins_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins.push(val)
}

ffi_fn! {
unsafe fn kmerminhash_get_abund_idx(ptr: *mut KmerMinHash, idx: u64) -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
      Ok(abunds[idx as usize])
    } else {
      Ok(0)
    }
}
}

#[no_mangle]
pub extern "C" fn kmerminhash_get_abunds_size(ptr: *mut KmerMinHash) -> usize {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
        abunds.len()
    } else {
        0
    }
}

#[no_mangle]
pub extern "C" fn kmerminhash_abunds_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
        abunds.push(val)
    }
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

ffi_fn! {
unsafe fn kmerminhash_merge(ptr: *mut KmerMinHash, other: *const KmerMinHash) -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.merge(other_mh)?;
    Ok(())
}
}

ffi_fn! {
unsafe fn kmerminhash_add_from(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.add_from(other_mh)
}
}

ffi_fn! {
unsafe fn kmerminhash_count_common(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.count_common(other_mh)
}
}

ffi_fn! {
unsafe fn kmerminhash_intersection(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    if let Ok((common, size)) = mh.intersection(other_mh) {
        return Ok(size);
    }
    Ok(0)
}
}

ffi_fn! {
unsafe fn kmerminhash_compare(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<f64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.compare(other_mh)
}
}
