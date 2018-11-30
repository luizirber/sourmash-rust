all: test

check: build test bench

build:
	cargo build

bench:
	cargo bench

test:
	cargo test

include/sourmash.h: src/lib.rs src/ffi.rs src/errors.rs
	RUST_BACKTRACE=1 cbindgen --clean -c cbindgen.toml -o $@

.phony: test
