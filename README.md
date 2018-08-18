<meta charset="utf-8"/>

# sourmash-rust

**sourmash core functionality implemented in Rust**

[![Build Status](https://travis-ci.org/luizirber/sourmash-rust.svg?branch=master)](https://travis-ci.org/luizirber/sourmash-rust)
[![Build status](https://ci.appveyor.com/api/projects/status/pk0druxv8eafq0r1?svg=true)](https://ci.appveyor.com/project/luizirber/sourmash-rust)
[![](http://meritbadge.herokuapp.com/sourmash)](https://crates.io/crates/sourmash)
[![](https://img.shields.io/crates/d/sourmash.svg)](https://crates.io/crates/sourmash)
[![API Documentation on docs.rs](https://docs.rs/sourmash/badge.svg)](https://docs.rs/sourmash)

[sourmash][sourmash] is a command-line tool and Python library for
computing MinHash sketches from DNA sequences, comparing them to each other, and
plotting the results. This allows you to estimate sequence similarity between
even very large data sets quickly and accurately. The core data structure is
implemented in C++.

[sourmash]: https://github.com/dib-lab/sourmash

There is a [PR in sourmash][sourmash_pr] to replace the C++ core with this
implementation (tests passing, yay!).

[sourmash_pr]: https://github.com/dib-lab/sourmash/pull/424

Another goal is to compile this code to webassembly and use it in the browser.
There is a [NPM package][package] already available, based on this [PR][wasm].
For an example usage of the NPM package check [wort-dnd][wort-dnd].

[package]: https://www.npmjs.com/package/sourmash
[wasm]: https://github.com/luizirber/sourmash-rust/pull/3
[wort-dnd]: https://github.com/luizirber/wort-dnd

For more details, check [luizirber/2018-python-rust][poster] for a poster
presented at [GCCBOSC][bosc] and [SciPy][scipy] 2018.

[poster]: https://github.com/luizirber/2018-python-rust
[bosc]: https://gccbosc2018.sched.com/
[scipy]: https://scipy2018.scipy.org/ehome/index.php?eventid=299527&tabid=712461&cid=2233543&sessionid=21618890&sessionchoice=1&

## License

This project is licensed under a [BSD 3-Clause License](LICENSE).
