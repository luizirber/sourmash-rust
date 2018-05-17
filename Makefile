all: test

test:
	cargo test

target/sourmash.h: src/lib.rs src/ffi.rs src/errors.rs
	RUST_BACKTRACE=1 cbindgen --clean -c cbindgen.toml -o $@

wasm:
#	wasm-pack init --mode no-installs
	wasm-pack init

.phony: test
