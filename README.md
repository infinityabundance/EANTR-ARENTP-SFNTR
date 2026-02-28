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

## Repository Layout

The main simulation assets live here:

- CLI crate: `crates/ntp-hybrids-sim`
- Colab notebook: `crates/ntp-hybrids-sim/ntp-hybrids-sim-visualization.ipynb`
- Output root: `output-ntp-hybrids-sim`

The Rust crate is the reproducible batch interface. It generates timestamped CSV outputs for the paper-aligned analysis.

The Colab notebook is the interactive interface. It can clone the repository, install Rust, run the crate, and visualize the generated outputs directly in the browser.

The notebook button above opens that Colab workflow directly.
