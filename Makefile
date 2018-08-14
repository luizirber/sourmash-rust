all: test

test:
	cargo test

target/sourmash.h: src/lib.rs src/ffi.rs src/errors.rs
	RUST_BACKTRACE=1 cbindgen --clean -c cbindgen.toml -o $@

.phony: test
