# Third-party notices

This notice distinguishes one incorporated CC-BY-4.0 source adaptation and
two incorporated MIT source adaptations from scientific papers and external
data. These third-party materials do not change the Apache-2.0 license of
other repository-authored code.

## Incorporated source: Finsberg, Sundnes, and van den Brink cardiac benchmark

Affected file and license:
`examples/cardiac_benchmark/activation.py` — CC BY 4.0.

The activation and pressure histories are adapted from these files in the
exact Zenodo deposit credited to Henrik Finsberg, Joakim Sundnes (listed as
"sundnes" in the record), and Jonas van den Brink:

- [`activation_model.py`](https://github.com/finsberg/cardiac_benchmark/blob/325d17d850c2e2032abb85a4191a5795d3008ab7/src/cardiac_benchmark/activation_model.py)
- [`pressure_model.py`](https://github.com/finsberg/cardiac_benchmark/blob/325d17d850c2e2032abb85a4191a5795d3008ab7/src/cardiac_benchmark/pressure_model.py)
- `finsberg/cardiac_benchmark` v1.0.0, peeled commit
  `325d17d850c2e2032abb85a4191a5795d3008ab7`
- archived release: [Zenodo 10.5281/zenodo.10875818](https://doi.org/10.5281/zenodo.10875818)
- exact archive: `finsberg/cardiac_benchmark-v1.0.0.zip`, 345,885 bytes,
  Zenodo MD5 `be92da5dbc1fd26d424bf88ef7db13b4`, independently computed
  SHA-256 `df6e5f03e644cb055ba3649f901030ba1d18840ff2ea94f6c71fd12bded28185`
- archived `activation_model.py` SHA-256
  `aff5299d14ff63d61b091a581e720fe5585ac2d1b31076c0f8039da844156a63`
- archived `pressure_model.py` SHA-256
  `3ac374f97c2a122b527c59f6e4fa9bebf379449086b3ce155f26caba3b2aae43`

The Zenodo record labels the v1.0.0 software deposit CC BY 4.0. This project
uses that conservative license basis; the complete terms are retained at
`LICENSES/CC-BY-4.0.txt`. The archived README also says "MIT," while the
`LICENSE` referenced by the archived `pyproject.toml` is absent. This project
does not rely on that incomplete MIT record and does not supply a downstream-
created MIT copyright notice for it.

Modification indication: CoupFE-Cardiac combined the two upstream modules,
routed parameters through its benchmark registry, vectorized scalar math,
added input/integration/completion/shape/finiteness checks, and renamed the
public functions to `tau_of_t` and `p_of_t`. See the file header for the same
attribution and modification notice.

## Incorporated source: cardiac benchmark toolkit

Affected files and license:

- `examples/cardiac_benchmark/fiber_crosscheck.py` — MIT.
- `examples/cardiac_benchmark/structural_directions.py` — MIT.

The NumPy cross-check and runtime physical-coordinate reconstruction are
adapted from
[`ellipsoid_fiber_generation.py`](https://github.com/Reidmen/cardiac_benchmark_toolkit/blob/e8d47553cfc83eb274eba3e177de33148e7f441c/src/cardiac_benchmark_toolkit/ellipsoid_fiber_generation.py)
in `Reidmen/cardiac_benchmark_toolkit`, commit
`e8d47553cfc83eb274eba3e177de33148e7f441c`. The pinned tree contains the
[MIT license](https://github.com/Reidmen/cardiac_benchmark_toolkit/blob/e8d47553cfc83eb274eba3e177de33148e7f441c/LICENSE).

The exact upstream MIT file is retained at
`LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt`.

## External paper and dataset

The benchmark article and the separately downloaded Zenodo comparison dataset
are both credited in `NOTICE` and `examples/REFERENCES.md`. The 23.2 GB archive
and raw team pickle files are not vendored. Checked-in reports and SVG figures
do redistribute transformed CC-BY-4.0 material derived from the dataset:
interpolated all-team mean, standard deviation, minimum/maximum envelope,
selected Simula curves, and numerical comparison summaries. CoupFE-Cardiac
selected hash-pinned publisher files, mapped them to stated common time grids,
and computed those statistics and graphics; these are modifications, not raw
publisher files. The source is R. A. Arostica Barrera and Cristobal Bertoglio,
"A software benchmark for cardiac elastodynamics," Zenodo,
https://doi.org/10.5281/zenodo.14260459, licensed CC BY 4.0. The full license
text is in `LICENSES/CC-BY-4.0.txt`. These terms do not relicense
repository-authored source code.
