# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **SpecHLA now runs on the HPC profiles.** `params.spechla_path` defaulted to
  `/app/SpecHLA`, a path inside the container image, while the `roihu` and `puhti`
  profiles set `container = null` for the SPECHLA processes so they run natively on
  the host. The command was therefore
  `bash /app/SpecHLA/script/whole/SpecHLA.sh ...`, which exits 127. Confirmed on
  Roihu 2026-09-18. `conf/roihu_params.yaml` and `conf/puhti_params.yaml` now carry
  `spechla_path: /projappl/project_2008084/SpecHLAx` and
  `use_local_spechla: true` -- the values the working `slurm_giab_*.sh` scripts
  passed on the command line.

  Verified end to end on Roihu against HG00096: SpecHLA produced a complete
  `hla.result.txt` typing all six genes, and the class-I calls
  (A\*01:01 / A\*29:02, B\*08:01 / B\*44:03, C\*07:01 / C\*16:01) match the 1000
  Genomes reference typing. This is the first call SpecHLA has produced in this
  framework.

- **SpecHLA no longer manufactures an empty result.** `modules/spechla.nf` wrote
  `# No results generated` into `${sample}_spechla.txt` and exited 0 when no
  `hla.result.txt` was found -- the 196-byte header-only file seen in every pilot
  stage, indistinguishable downstream from a caller that genuinely typed nothing.
  It now fails and reports where it looked. Both processes also check that
  `SpecHLA.sh` exists before spending twenty minutes on BAM sorting and FASTQ
  extraction.

- **The test suite runs, and CI with it.** `requirements.txt` could never be added
  because `.gitignore` carried a blanket `*.txt`, so CI's
  `pip install -r requirements.txt` failed before any test ran. The subprocess
  tests and `bin/run_1000g_benchmark.py` spawned `python3` rather than the
  interpreter running them; test fixtures hardcoded a CSC scratch path for files
  inside the repository; and relative config paths were re-anchored when the runner
  serialised a config into a temporary directory. 19 failures to 0.

### Fixed — Roihu

- **The pipeline can be started on Roihu at all.** Every `slurm_*.sh` said
  `module load nextflow`, a Puhti idiom that cannot work there:
  `/etc/profile.d/zz-csc-env.sh`, which defines `module`, returns immediately when
  `$PS1` is unset -- so in every SLURM batch script -- unless
  `CSC_ENV_INIT_NON_INTERACTIVE=yes` is exported first; and the bare name
  `nextflow` has no default version and sits behind a `bio-apps/v202603`
  prerequisite. `bin/roihu_env.sh` holds the working sequence (nextflow 25.10.2,
  openjdk 17.0.11, samtools 1.21) and `slurm_panelhla_roihu.sh` is the launcher.
  Nextflow was never missing from Roihu; nothing here knew how to ask for it.

- **HLA-HD could not run on Roihu.** `conf/roihu_params.yaml` pointed `hlahd_db` at
  a host directory that was never migrated and does not exist, so `hlahd.sh` was
  handed a missing `-f <db>/freq_data`. `main.nf` only checks the parameter is
  non-empty, so the wrong path passed the guard. The database ships inside
  `hlahd.sif`; the parameter now names it (`/app/hlahd.1.4.0`).

- **POLYSOLVER produced nothing on chr-prefixed BAMs.** Its own
  `shell_call_hla_type` hardcodes the hg38 HLA regions as `6:29941260-...` while
  GRCh38_full_analysis_set BAMs name the contig `chr6`, so every query matched
  nothing. The module now detects the naming from the BAM header and rewrites the
  queries, and records the exit status instead of discarding it with `|| true`.

- **Nine more placeholder-result writers removed**, across polysolver, kourami,
  locityper and immuannot, on top of the five fixed earlier. Each wrote a
  header-only file and exited 0 when its tool produced nothing -- one kourami
  branch called `exit 0` outright. No module manufactures an empty result now.

- **CI parses the Nextflow workflow.** It previously ran `pytest` and nothing else,
  so a malformed process block could reach `main` with the tests green. Pinned to
  Nextflow 25.10.2, the version Roihu provides.

- **Nextflow 26.x is incompatible** and the manifest now says so
  (`>=23.04.0, <26.0.0`). Its strict config parser rejects function definitions in
  `nextflow.config`, and `check_max()` is one. Lifting the bound means migrating to
  `process.resourceLimits`.

- **SpecHLA never reached the consensus.** It was the only caller with no parse
  step, so `modules/spechla.nf` copied its native wide table
  (`Sample | HLA_A_1 | HLA_A_2 | ...`) straight through, while
  `modules/aggregation.nf` reads only long format. Aggregation took "Sample" as a
  gene on one line and the sample id on the next, matched neither against
  `mv_genes`, and dropped both -- so SpecHLA ran, called correctly, and contributed
  nothing to the vote on every run. `bin/parse_spechla_results.py` converts wide to
  long using the same column convention `extract_wide_gene_columns()` in
  `bin/hla_benchmark.py` uses, and the native output is kept as
  `${sample}_spechla_raw.txt`. No published result is affected: the manuscript's
  numbers come from `hla_benchmark.py`, which parses SpecHLA correctly.

- **The consensus step called nothing.** `conf/tool_weights_{wgs,wes,rna}.json` are
  flat `{tool: weight}` maps, but `lookup_weight()` understood only a nested runtime
  schema and a legacy `raw_accuracy` one and returned `0.0` for everything else. On
  a sample five callers had typed correctly, every gene came back
  `no_call / no_nonzero_weight`. The keys also spelled tools `OptiType` and
  `HLA-HD` while `aggregated_calls.tsv` says `optitype` and `hlahd`. Both fixed, and
  `tests/test_majority_voting_weights.py` asserts every tool in every shipped weight
  file resolves non-zero.

- **seq2HLA's p-value filter was inverted.** Its "Confidence" column is a p-value --
  its README says the tool reports "a p-value for each call" -- but the parser
  skipped loci whose values were *below* 0.1, so it discarded seq2HLA's strongest
  calls (A, B, C at p < 0.05) and kept its weakest (DQA1 at p 0.135/0.217). The
  comparison flips and the constant is renamed `P_VALUE_THRESHOLD`.

- **The consensus step never ran.** `enable_majority_voting` defaults to false, so
  runs stopped at per-tool calls with no `aggregated_calls.tsv`. The Roihu launcher
  now requests it; the default is unchanged for anyone wanting caller outputs alone.

- `bin/roihu_env.sh` guards `set -u` while sourcing the CSC environment, which
  dereferences `$PS1` unguarded, and `slurm_panelhla_roihu.sh` takes the repository
  location from `PANELHLA_HOME` -- SLURM copies a batch script to its spool
  directory, so `$BASH_SOURCE` cannot find it.

### Known

- `conf/tool_weights_wgs.json` has no POLYSOLVER entry, so POLYSOLVER scores zero
  in a calibrated WGS vote. That is a missing calibrated measurement, not a code
  defect. `tool_weights_wes.json` -- the modality the paper reports -- does include
  it.
- The pipeline default is `weighting = calibrated`, while the project's primary
  method is the equal-weight plurality vote (`--weighting equal`).

- `.gitattributes` pins LF for scripts, configs and workflow files. The repository
  is edited from Windows and deployed to Linux.

### Added
- **Algorithm-family and correlation-aware challenger gates.**
  `champion_challenger.override_policy.gate_mode` selects `tool_count` (default,
  unchanged behaviour), `family_count`, `correlation` or `learned`.
  `conf/tool_families.yaml` groups callers into alignment, assembly-graph and k-mer
  families; a caller outside the map forms its own singleton family, so the gate
  never merges voters silently. `conf/champion_challenger_family_gate.yaml` enables
  it as an overlay. Reported in the manuscript as the FamilyGatedCC comparator.

### Changed
- Project renamed to **PanelHLA** throughout the documentation, and the README now
  leads with the equal-weight plurality vote -- the primary method since the
  2026-09-08 amendment -- rather than Champion-Challenger, which is a reported
  ablation.

## [3.1.0] - 2026-04-19

### Added
- `--weight-alpha` / `--weight-beta` CLI flags on `hla_benchmark.py` and `run_1000g_benchmark.py`
  to override the default 0.7/0.3 reliability-confidence split at runtime
- `--figures-dir` and `--tables-dir` are now required args in `generate_html_report.py`;
  hard-coded scratch paths and username removed for portability
- Weight formula sensitivity analysis: three variants (1.0/0.0, 0.5/0.5, 0.0/1.0) evaluated
  via `slurm_weight_sensitivity.sh`; results in `analysis/weight_sensitivity/`
- Wilson 95% CIs on all `overall_correct_call_rate` summary columns in benchmark output tables
- `analysis/weight_sensitivity/sensitivity_comparison.tsv`: cross-modality summary of
  WeightedConsensus performance across all weight variants

### Changed
- `bin/generate_html_report.py`: `--figures-dir` and `--tables-dir` are now required (no defaults)

## [3.0.0] - 2026-04-18

### Added
- 8-tool HLA typing ensemble: OptiType, ArcasHLA, SpecHLA, HLA-HD, Kourami, T1K, POLYSOLVER, Seq2HLA
- `read_confidence_v2` confidence weighting: `final_weight = 0.7 × base_reliability + 0.3 × effective_confidence_score`
- Guardrail system per tool × modality: `applied`, `poor_calibration`, `no_confidence` statuses with shrinkage factors
- Tri-modal benchmark (WGS n=99, WES n=51, RNA-seq n=50) with 60/20/20 training/validation/holdout split
- 1000 Genomes NYGC 30× CRAM integration: 49-sample wave-2 WGS batch downloaded from EBI and typed
- Publication-quality figure generation (`bin/generate_figures_v2.py`): 8 figures PNG 300dpi + PDF
- HTML benchmark report v4 (`analysis/generate_report_v4.py`): standalone 3MB report embedding all figures
- IMGT/HLA v3.59.0 pinned for reproducibility
- Abstention tradeoff analysis (`abstention_tradeoff.tsv`): callable rate vs accuracy at configurable thresholds
- Discordance taxonomy (`discordance_summary.tsv`): DNA/RNA discordance, expression bias, technical conflict tags
- SLURM array scripts with `%10` throttle for disk-safe batch processing

### Key benchmark results (holdout)
- WES (n=12): WeightedConsensus 91.7%, MajorityVote 94.4%
- RNA-seq (n=7): WeightedConsensus/MajorityVote 100%, T1K 90.5%
- WGS (n=20): WeightedConsensus 26.7%, callable rate 56.7% at threshold 0.45
- WGS best single tool: OptiType 43.3%

### Changed
- Benchmark scope: `partial_realdata_benchmark` (A/B/C loci; DRB1/DQB1 excluded pending truth source)
- Truth source: Gourraud et al. 2014 1000 Genomes (HLA-A, -B, -C only)
- Benchmark split expanded from 30/10/10 to 60/20/20 to accommodate 100-sample tri-modal cohort

---

## [2.0.0] - 2025-11-21

### Added
- Complete multi-tool HLA typing pipeline with Nextflow DSL2
- Support for three HLA typing tools:
  - OptiType (DNA/RNA-seq)
  - ArcasHLA (RNA-seq optimized, works with DNA)
  - SpecHLA (Exome/WGS with variant calling)
- Majority voting / consensus calling module
- Flexible input support (BAM, CRAM, FASTQ)
- HLA region extraction for space-efficient processing
- Chromosome 6 extraction optimization for SpecHLA
- CSC Puhti HPC profile with SLURM support
- Docker and Singularity container support
- Comprehensive error handling and logging
- Automatic BAM index creation
- Chromosome naming convention detection (chr6 vs 6)
- Exon-only mode for SpecHLA (critical for exome data)
- SLURM submission script with proper parameter handling
- Comprehensive documentation (README, QUICKSTART guide)
- Pipeline execution reports (timeline, trace, DAG)

### Features
- **Smart Input Handling**: Automatically detects and processes different input formats
- **Optimized for Exome Data**: Special exon-only mode for SpecHLA
- **Space Efficient**: Optional HLA region extraction to reduce disk usage
- **Resume Capability**: Nextflow's built-in resume functionality
- **Parallel Processing**: Efficiently processes multiple samples simultaneously
- **Comprehensive Logging**: Detailed logs for each tool and process
- **Resource Management**: Configurable CPU, memory, and time limits
- **Container Ready**: Pre-configured for Docker and Singularity

### Fixed
- SLURM script parameter handling (parameters now properly passed to Nextflow)
- Container paths for CSC Puhti pre-downloaded SIF files
- SpecHLA chromosome naming issues (chr6 vs 6)
- BAM indexing in various scenarios
- Error handling for missing or empty files

### Configuration
- CSC Puhti profile with pre-downloaded containers
- SLURM executor configuration
- Resource labels for different process types
- Automatic retry on failure
- Comprehensive parameter validation

### Documentation
- Complete README with usage examples
- Quick Start guide for new users
- Detailed parameter descriptions
- Troubleshooting section
- CSC Puhti specific instructions
- Performance benchmarks

### Known Issues
- SpecHLA requires exon-only mode (`--spechla_exon_only 1`) for exome data
- Some tools may produce empty results for low-coverage samples
- Chromosome 6 naming must be detected correctly from BAM headers

### Compatibility
- Nextflow: ≥23.04.0
- Docker/Singularity required
- Tested on CSC Puhti (SLURM)
- Compatible with hg19 and hg38 reference builds

## [1.0.0] - 2025-11-XX

### Initial Development
- Basic pipeline structure
- Single tool support
- Initial testing on 1000 Genomes data

---

## Planned for Future Releases

### [2.1.0] - Planned
- [ ] Additional HLA typing tools (HLA-HD, HLA*LA, xHLA)
- [ ] Advanced consensus algorithms
- [ ] Quality score integration
- [ ] Phasing information
- [ ] Interactive HTML reports
- [ ] Database storage support

### [2.2.0] - Planned
- [ ] Cloud deployment (AWS, Google Cloud)
- [ ] GUI interface
- [ ] Real-time monitoring dashboard
- [ ] Automated quality control
- [ ] Population frequency analysis

### [3.0.0] - Planned
- [ ] Machine learning-based consensus
- [ ] Novel allele detection
- [ ] Integration with clinical databases
- [ ] Multi-sample joint calling
