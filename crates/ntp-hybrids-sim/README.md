# ntp-hybrids-sim

`ntp-hybrids-sim` is a conceptual Rust simulator for the paper _Novel Hybrid Nuclear Thermal Propulsion Concepts for Enabling Affordable, Routine Human Mars Missions_.

## What It Does

The crate models three hybrid nuclear thermal propulsion concepts:

- `EANTR`
- `ARENTP`
- `SF-NTR`

It provides:

- A Tsiolkovsky rocket-equation model with variable specific impulse.
- Propellant mass fraction comparisons against a 65% baseline NTR case and a 35-45% hybrid target band.
- A hydrogen temperature-to-Isp scaling model using `sqrt(T/M)`.
- Paper-aligned mechanism models for:
  - `EANTR`: partial ionization plus electrostatic acceleration.
  - `ARENTP`: acoustic heat-transfer enhancement and resonance-driven thrust boost.
  - `SF-NTR`: supercritical pressure tuning and chamber-pressure uplift.
- Variable-mode thrust profiles for high-thrust departure, cruise augmentation, thrust boosting, and chamber-pressure tuning.
- Monte Carlo uncertainty studies with exactly 360 runs per hybrid.
- CSV summary tables written into timestamped folders under `output-ntp-hybrids-sim`.

## Why It Exists

The project is intended to make the paper's hybrid-NTP trade space easier to inspect, rerun, and visualize.

Instead of leaving the concepts at the narrative level, the crate gives a concrete workflow for:

- comparing baseline and hybrid propellant requirements,
- exploring how Isp and thermal state influence performance,
- estimating transit-time sensitivity across propulsion modes,
- and generating reproducible output files for further analysis in Python or Colab.

## Running It

From this crate directory:

```bash
cargo run -- --hybrid EANTR --dv 5500 --runs 360
```

Supported hybrid selections:

- `EANTR`
- `ARENTP`
- `SF-NTR`
- `all`

The CLI creates a new output folder each run using:

```text
output-ntp-hybrids-sim/YYYY-MM-DD_HH-MM-SS
```

That keeps CSVs and generated artifacts from being overwritten.

## Outputs

Each run writes:

- `summary.csv`
- `isp_sweep_<hybrid>.csv`
- `mechanism_sweep_<hybrid>.csv`
- `monte_carlo_<hybrid>.csv`
- `thrust_profile_<hybrid>.csv`

The companion notebook, `ntp-hybrids-sim-visualization.ipynb`, can build the crate in Google Colab, run the same CLI, and generate saved plots from those CSV outputs.
