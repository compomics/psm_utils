# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2022-11-12

### Added

- `PSM`: Add `ion_mobility` field
- `PSMList`: Allow slicing with bool arrays (e.g., `psm_df[psm_df["retention_time"] < 2000]`)
- `rename_modifications`: Add support for fixed modifications
- Add example files
- Online: Add support for GZipped files
- Online: Add support for logarithmic score (e.g. e-values)
- Docs: Extend contributing with example contributions
- Docs: Add notes to `PSM.get_usi()` method
- Docs: Extend quickstart on PSMList
- Docs: Add "psm_utils tags" for file formats, as used in high-level read/write/convert functions
- Docs: Peptide Record: add notes on unsupported modification types; add example for C-terminal modification
- Docs: More clearly document conversion to DataFrame
- Docs: Add bioconda install instructions
- Docs: Add citation for preprint
- Tests: Added tests for PSMList `set_ranks` and `get_rank1_psms` methods

### Changed

- `PSMList`: Refactor `set_ranks` and `get_rank1_psms` methods
- Update `.vscode/settings.json`
- Typing: Replace Union with OR operator `|`
- Online: Use percentiles instead of randomly sampling for PP plot
- Docs: Force TOC-tree max depth
- Tests: Expand unit tests in general

### Fixed

- `PSMList`: Truncate __repr__ to first five entries only, avoiding crashing notebook output
- `Peptidoform`: Minor typing fix
- `add_fixed_modifications`: Allow input as dict as well as list of tuples
- `io`: Fix issue where the `NamedTemporaryFile` for `_supports_write_psm` was seen as invalid Percolator file
- `io.convert`: pass ` progressbar` argument to class, not `write_file`
- `io.mzid`: Add more supported MS-GF score names, make SpecEValue default
- `io.peptide_record`: `spec_id` is now a required column (`spectrum_id` is also required in PSM)
- `io.peptide_record`: Fix parsing of C-terminal modifications from proforma to peprec
- `io.percolator`: Fix Percolator peptide notation writing (fixes #18)
- `io.tsv`: Fix issue where `TSVReader` would not use string type for metadata
- `io.xtandem`: Fix issue where optional arguments were not accepted by `XTandemReader`
- `io.xtandem`: Do not split spectrum title on space
- `io.xtandem`: Fix issue where optional arguments were not accepted by `XTandemReader`
- Online: Fix pi-0 diagonal calculation
- Remove obsolete to do comments in code

## [0.1.0] - 2022-10-14

### Added

- Initial version
