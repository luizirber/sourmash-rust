all: test

test:
	cargo test

target/sourmash.h: src/lib.rs src/ffi.rs src/errors.rs
	$(eval tempdir := $(shell mktemp -d))
	RUST_BACKTRACE=1 CARGO_EXPAND_TARGET_DIR=${tempdir} cbindgen -c cbindgen.toml -o $@
	-rm -rf ${tempdir}

.phony: test
