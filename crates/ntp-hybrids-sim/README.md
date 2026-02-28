# ntp-hybrids-sim

`ntp-hybrids-sim` is a Rust CLI for reproducing and exploring the conceptual hybrid nuclear thermal propulsion trade space from the paper _Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable, Routine Human Mars Missions_.

The crate focuses on three concepts from the paper:

- `EANTR`: electrostatic augmentation of a nuclear-thermal exhaust stream through partial ionization and downstream grid acceleration.
- `ARENTP`: acoustic resonance enhancement inside fuel-element channels to improve heat transfer and raise effective propulsion performance.
- `SF-NTR`: full-path supercritical propellant operation with higher chamber pressure and supercritical heat-transfer gains.

## Paper

Source paper:

- DOI: <https://doi.org/10.5281/zenodo.18809245>
- Repository: <https://github.com/infinityabundance/EANTR-ARENTP-SFNTR>

Suggested citation:

```text
Riaan de Beer (2026), "Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable,
Routine Human Mars Missions: The Electrostatic-Augmented Nuclear Thermal Rocket (EANTR),
Acoustic Resonance-Enhanced Nuclear Thermal Propulsion (ARENTP), and Supercritical-Fluid
Nuclear Thermal Rocket (SF-NTR)," Version 1.0, February 28, 2026. DOI: 10.5281/zenodo.18809245
```

## What This Crate Does

The crate implements the simplified 0-D and 1-D style numerical workflow described in the paper's software section.

For a selected hybrid and mission delta-v, it:

- Computes mass ratio with the Tsiolkovsky rocket equation for variable specific impulse.
- Compares baseline NTP propellant mass fraction to the hybrid comparison bands discussed in the paper.
- Uses hydrogen `sqrt(T/M)` scaling where the model needs temperature-derived Isp.
- Applies paper-aligned mechanism parameters for each hybrid instead of generic placeholder ranges.
- Simulates variable mission modes for trans-Mars injection, heliocentric cruise, and Mars-orbit capture.
- Runs exactly `360` Monte Carlo samples per hybrid to expose conceptual uncertainty in Isp, propellant mass fraction, and transit time.
- Writes all outputs into a timestamped folder so repeated runs never overwrite prior results.

## Why This Crate Exists

The paper is intentionally conceptual. Its main value is in comparing propulsion mechanisms and mission-level implications, not in claiming high-fidelity flight analysis.

This crate exists to turn those conceptual claims into a reproducible workflow:

- so the paper's hybrid assumptions can be inspected rather than only read narratively,
- so baseline-vs-hybrid tradeoffs can be rerun quickly for different mission delta-v values,
- so uncertainty can be visualized through repeatable Monte Carlo output,
- and so readers can use the resulting CSVs in Rust, Python, spreadsheets, or Google Colab without rebuilding the model themselves.

In practice, the CLI is the reproducible batch interface and the companion notebook is the interactive visualization layer.

## Modeling Scope

The crate is paper-aligned, but deliberately simplified.

It includes:

- Tsiolkovsky mass-ratio calculations.
- Hybrid-specific mechanism parameterizations for `EANTR`, `ARENTP`, and `SF-NTR`.
- Variable-mode phase profiles.
- Paper-band summary comparisons.
- CSV outputs for downstream plotting and review.

It does not include:

- detailed CFD,
- neutronics,
- structural or thermal stress analysis,
- radiation shielding analysis,
- high-fidelity trajectory optimization,
- or engineering qualification of the proposed concepts.

The results are conceptual estimates for comparison and visualization, not flight-certified performance predictions.

## How The Rust Crate Works

The crate operates as a command-line batch workflow with a fixed sequence:

1. Parse inputs with `clap`.
   - `--hybrid` selects `EANTR`, `ARENTP`, `SF-NTR`, or `all`.
   - `--dv` sets the mission delta-v in m/s.
   - `--runs` is fixed to the paper's `360` Monte Carlo runs per hybrid.
2. Create a dated output directory under `output-ntp-hybrids-sim/YYYY-MM-DD_HH-MM-SS`.
3. Build the concept-specific mechanism parameters for the selected hybrid.
4. Simulate a variable-mode mission profile across:
   - trans-Mars injection,
   - heliocentric cruise,
   - and Mars-orbit capture.
5. Evaluate the simplified propulsion model using:
   - the Tsiolkovsky rocket equation,
   - hydrogen `sqrt(T/M)` scaling where appropriate,
   - paper-aligned effective-Isp calculations for each concept,
   - and phase-level transit-time estimates.
6. Run an exactly `360`-sample Monte Carlo analysis that perturbs the main conceptual mechanism variables.
7. Write the resulting CSV files for later plotting, inspection, or downstream analysis.

The concept logic differs by hybrid:

- `EANTR` uses ionization fraction, electrostatic grid voltage, stage count, and acceleration efficiency to model augmented cruise behavior.
- `ARENTP` uses acoustic frequency and heat-transfer gain to model resonance-enhanced performance.
- `SF-NTR` uses chamber pressure, methane blend fraction, and supercritical heat-transfer gain to model high-pressure supercritical operation.

The crate is therefore not just producing one headline number. It is generating multiple views of the same conceptual trade study:

- a nominal mission summary,
- an Isp sensitivity sweep,
- a hybrid-specific mechanism sweep,
- a Monte Carlo uncertainty study,
- and a phase-by-phase thrust profile.

## How The Colab Notebook Works

The companion notebook uses the Rust crate directly rather than reimplementing the simulation in Python.

Its workflow is:

1. Ensure the repository exists in the Colab runtime, cloning it automatically if needed.
2. Install Rust if `cargo` is not already available.
3. Install Python plotting dependencies such as `pandas`, `matplotlib`, and `ipywidgets`.
4. Build the Rust crate in the notebook environment.
5. Expose widget controls so the user can choose the hybrid and mission delta-v.
6. Run the same CLI command the crate expects.
7. Load the CSV files produced by the crate.
8. Render visual outputs and compact summary tables for the key paper-facing metrics.
9. Save plots into the same dated output directory as the CSV files.

This matters because the notebook is purely a visualization and execution shell around the crate. The source of truth for the actual analysis remains the Rust code.

## Running The CLI

From this crate directory:

```bash
cargo run -- --hybrid EANTR --dv 5500 --runs 360
```

Supported hybrid selections:

- `EANTR`
- `ARENTP`
- `SF-NTR`
- `all`

The CLI creates a new output folder for every run:

```text
output-ntp-hybrids-sim/YYYY-MM-DD_HH-MM-SS
```

## Output Files

Each run writes:

- `summary.csv`
- `isp_sweep_<hybrid>.csv`
- `mechanism_sweep_<hybrid>.csv`
- `monte_carlo_<hybrid>.csv`
- `thrust_profile_<hybrid>.csv`

These files are designed to answer slightly different questions:

- `summary.csv` reports the main paper-facing mission figures and reference bands.
- `isp_sweep_<hybrid>.csv` shows how raw rocket-equation propellant fraction changes with Isp.
- `mechanism_sweep_<hybrid>.csv` shows how the primary concept driver affects effective Isp, PMF, and transit time.
- `monte_carlo_<hybrid>.csv` captures the 360-run uncertainty study for the selected concept.
- `thrust_profile_<hybrid>.csv` records the phase-by-phase variable-mode mission behavior.

## Companion Notebook

The repository also includes a Google Colab notebook for building the crate, running the CLI, and plotting the generated CSV outputs:

- Notebook path: `crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb`
- Colab link: <https://colab.research.google.com/github/infinityabundance/EANTR-ARENTP-SFNTR/blob/main/crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb>

The notebook exists for readers who want a guided, visual workflow rather than a terminal-first one.

## License

Licensed under Apache License 2.0.
