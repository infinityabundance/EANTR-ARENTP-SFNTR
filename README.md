# EANTR-ARENTP-SFNTR
Conceptual hybrid nuclear thermal propulsion (NTP) architectures for affordable, routine human Mars missions: EANTR, ARENTP, and SF-NTR

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/infinityabundance/EANTR-ARENTP-SFNTR/blob/main/crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb)

## Paper

Source paper:

- DOI: <https://doi.org/10.5281/zenodo.18809245>
- Title: _Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable, Routine Human Mars Missions_

Suggested citation:

```text
Riaan de Beer (2026), "Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable,
Routine Human Mars Missions: The Electrostatic-Augmented Nuclear Thermal Rocket (EANTR),
Acoustic Resonance-Enhanced Nuclear Thermal Propulsion (ARENTP), and Supercritical-Fluid
Nuclear Thermal Rocket (SF-NTR)," Version 1.0, February 28, 2026. DOI: 10.5281/zenodo.18809245
```

## What This Repository Does

This repository packages the paper together with a reproducible software workflow for exploring its conceptual propulsion trade study.

At the repository level, it does three things:

- documents the proposed hybrid NTP concepts and the paper they come from,
- provides a Rust simulation crate that reproduces the paper's simplified numerical comparisons,
- and provides a Google Colab workflow for readers who want to run the model and inspect outputs visually instead of only from the command line.

The software focuses on the three propulsion concepts introduced in the paper:

- `EANTR`: electrostatic augmentation of nuclear-thermal exhaust through partial ionization and downstream grid acceleration,
- `ARENTP`: acoustic resonance-enhanced heat transfer inside fuel-element channels,
- `SF-NTR`: full-path supercritical propellant operation with higher chamber pressure and improved thermal-fluid behavior.

The repository is intentionally designed so a reader can move from paper to executable model without needing to reconstruct the assumptions from scratch.

## Why This Repository Exists

The paper is conceptual by design. Its purpose is to propose and compare hybrid nuclear thermal propulsion mechanisms, not to present a flight-qualified engineering program.

This repository exists to make that conceptual work inspectable and reproducible:

- so the paper's assumptions are not locked inside prose alone,
- so the mission-level implications can be rerun for different delta-v values,
- so outputs like effective Isp, propellant mass fraction, and transit time can be explored with the same simplified workflow the paper describes,
- and so readers, reviewers, and collaborators can test the trade space directly in Rust or in Google Colab.

In practical terms, the repository bridges the gap between the paper's narrative claims and a usable computational artifact.

## Modeling Scope

The repository's simulation workflow is paper-aligned, but simplified.

It includes:

- Tsiolkovsky mass-ratio calculations,
- paper-aligned mechanism parameterizations for `EANTR`, `ARENTP`, and `SF-NTR`,
- variable-mode mission phase profiles,
- summary comparisons against baseline and paper-reference bands,
- Monte Carlo uncertainty analysis,
- and CSV plus plot-friendly outputs.

It does not include:

- detailed CFD,
- detailed neutronics,
- structural qualification,
- full mission trajectory optimization,
- radiation shielding design,
- or any claim of flight-readiness.

The outputs are conceptual estimates meant to support inspection, comparison, and visualization of the paper's ideas.

## How The Rust Crate Works

The Rust crate lives in `crates/ntp-hybrids-sim` and acts as the reproducible batch-analysis engine for the repository.

At a code level, it works in a fixed sequence:

1. It parses CLI inputs with `clap`, primarily the selected hybrid, mission delta-v, and Monte Carlo run count.
2. It builds a dated output directory under `output-ntp-hybrids-sim/YYYY-MM-DD_HH-MM-SS` so each run is isolated.
3. It selects the requested propulsion concept (`EANTR`, `ARENTP`, `SF-NTR`, or all three).
4. For each concept, it constructs paper-aligned mechanism parameters:
   - `EANTR`: ionization fraction, electrostatic grid voltage, stage count, and acceleration efficiency.
   - `ARENTP`: acoustic frequency and heat-transfer gain.
   - `SF-NTR`: chamber pressure, methane blend fraction, and supercritical heat-transfer gain.
5. It evaluates a variable-mode mission profile with phase-specific behavior for trans-Mars injection, heliocentric cruise, and Mars-orbit capture.
6. It computes concept outputs using the simplified models described in the paper:
   - Tsiolkovsky mass ratio,
   - hydrogen `sqrt(T/M)` Isp scaling where applicable,
   - paper-aligned effective-Isp models for each hybrid,
   - mission propellant mass fraction,
   - and phase-level transit-time estimates.
7. It performs an exactly `360`-run Monte Carlo study per hybrid to sample the main conceptual uncertainty drivers.
8. It writes multiple CSV artifacts so the results can be inspected outside the terminal.

Those CSV outputs are intentionally split by purpose:

- `summary.csv`: paper-facing mission summary values and comparison bands.
- `isp_sweep_<hybrid>.csv`: raw rocket-equation sensitivity to Isp.
- `mechanism_sweep_<hybrid>.csv`: sensitivity to the primary hybrid-specific mechanism parameter.
- `monte_carlo_<hybrid>.csv`: uncertainty-study output for the 360 sampled runs.
- `thrust_profile_<hybrid>.csv`: the phase-by-phase mission profile used for the selected concept.

The Rust crate is the authoritative computation layer. Everything else in the repository builds on those generated outputs.

## How The Colab Notebook Works

The Colab notebook lives at `crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb` and acts as the interactive visualization layer on top of the Rust crate.

Its workflow is:

1. Bootstrap the environment:
   - detect whether the repository already exists under `/content`,
   - clone the repository from GitHub if it does not,
   - install Rust if `cargo` is missing,
   - and install the Python plotting dependencies used for visualization.
2. Build the Rust crate directly in Colab.
3. Let the user choose a hybrid and delta-v value from notebook widgets.
4. Run the same CLI command the crate expects, so the notebook is not using a separate shadow implementation.
5. Load the generated CSV outputs from the dated output folder.
6. Render the paper-facing plots:
   - mass fraction versus Isp,
   - mechanism-driver versus effective Isp,
   - Monte Carlo propellant mass fraction distribution,
   - transit time histogram,
   - and a phase-by-phase mission profile.
7. Display compact summary tables instead of raw wide dataframe dumps, so the key mission and paper-comparison values remain readable in Colab.
8. Save the plots back into the same dated output folder as the CSVs.

This notebook exists so someone can inspect the model quickly, visually, and reproducibly without needing a local Rust development setup.

## Repository Layout

The main simulation assets live here:

- CLI crate: `crates/ntp-hybrids-sim`
- Colab notebook: `crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb`
- Output root: `output-ntp-hybrids-sim`

The Rust crate is the reproducible batch interface. It generates timestamped CSV outputs for the paper-aligned analysis.

The Colab notebook is the interactive interface. It can clone the repository, install Rust, run the crate, and visualize the generated outputs directly in the browser.

The notebook button above opens that Colab workflow directly.
