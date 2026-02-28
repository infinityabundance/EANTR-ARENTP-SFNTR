# ntp-hybrids-sim

`ntp-hybrids-sim` is a Rust CLI for exploring the conceptual hybrid nuclear thermal propulsion trade space described in the paper _Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable, Routine Human Mars Missions_.

It models three concepts:

- `EANTR` for partial ionization plus downstream electrostatic acceleration.
- `ARENTP` for acoustic resonance-enhanced heat transfer in fuel-element channels.
- `SF-NTR` for full-path supercritical propellant flow and high-pressure operation.

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

## What The Crate Does

The crate implements the same simplified 0-D and 1-D style analysis scope described in the paper's software section.

It provides:

- A Tsiolkovsky rocket-equation model with variable specific impulse.
- Propellant mass fraction comparisons against a 65% baseline NTR case and the paper's hybrid comparison bands.
- Hydrogen `sqrt(T/M)` Isp scaling.
- Paper-aligned mechanism models for `EANTR`, `ARENTP`, and `SF-NTR`.
- Variable-mode phase profiles for trans-Mars injection, cruise, and Mars-orbit capture.
- Monte Carlo uncertainty studies with exactly `360` runs per hybrid.
- CSV outputs written to timestamped folders under `output-ntp-hybrids-sim`.

## Install And Run

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

## Companion Notebook

The repository also includes a Google Colab notebook for building the crate, running the CLI, and plotting the generated CSV outputs:

- Notebook path: `crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb`
- Colab link: <https://colab.research.google.com/github/infinityabundance/EANTR-ARENTP-SFNTR/blob/main/crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb>

## License

Licensed under Apache License 2.0.

## Scope And Limitations

This crate is intentionally conceptual. It does not attempt detailed CFD, neutronics, thermal transients, structural analysis, radiation shielding analysis, or full mission trajectory optimization.

The numbers produced by the tool are paper-aligned estimates meant for reproducible comparison and visualization, not flight-ready performance predictions.
