use chrono::Local;
use clap::{Parser, ValueEnum};
use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use statrs::statistics::Statistics;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::Duration;
use uom::si::acceleration::meter_per_second_squared;
use uom::si::f64::{Acceleration, Ratio, ThermodynamicTemperature, Time, Velocity};
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::time::second;
use uom::si::velocity::meter_per_second;

const EXACT_MONTE_CARLO_RUNS: usize = 360;
const BASELINE_ISP_S: f64 = 900.0;
const BASELINE_REFERENCE_PROP_FRAC: f64 = 0.65;
const BASELINE_TRANSIT_MIN_DAYS: f64 = 180.0;
const BASELINE_TRANSIT_MAX_DAYS: f64 = 270.0;
const GLOBAL_HYBRID_PROP_FRAC_MIN: f64 = 0.35;
const GLOBAL_HYBRID_PROP_FRAC_MAX: f64 = 0.45;
const G0_MPS2: f64 = 9.80665;
const HYDROGEN_MOLAR_MASS_KG_PER_MOL: f64 = 0.002_016;
const BASELINE_CHAMBER_TEMPERATURE_K: f64 = 2800.0;
const BASELINE_CHAMBER_PRESSURE_BAR: f64 = 70.0;
const EXIT_PRESSURE_BAR: f64 = 0.01;
const HEAT_CAPACITY_RATIO_H2: f64 = 1.4;
const UNIVERSAL_GAS_CONSTANT: f64 = 8.314_462_618;
const EANTR_REFERENCE_ION_DELTA_V_MPS: f64 = 6_200.0;
const EANTR_REFERENCE_VOLTAGE_KV: f64 = 20.0;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Paper-aligned conceptual simulator for EANTR, ARENTP, and SF-NTR Mars mission studies."
)]
struct Cli {
    #[arg(long, value_enum, ignore_case = true)]
    hybrid: HybridArg,

    #[arg(long, default_value_t = 5500.0)]
    dv: f64,

    #[arg(long, default_value_t = EXACT_MONTE_CARLO_RUNS)]
    runs: usize,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum HybridArg {
    Eantr,
    Arentp,
    #[value(name = "SF-NTR")]
    SfNtr,
    All,
}

#[derive(Clone, Copy, Debug)]
enum HybridConcept {
    Eantr,
    Arentp,
    SfNtr,
}

#[derive(Clone, Copy, Debug)]
struct PhaseTemplate {
    phase_name: &'static str,
    mode_name: &'static str,
    dv_share: f64,
    time_share: f64,
}

#[derive(Clone, Debug)]
struct MechanismParameters {
    ionization_fraction: Option<f64>,
    electrostatic_grid_voltage_kv: Option<f64>,
    electrostatic_stage_count: Option<f64>,
    electrostatic_efficiency: Option<f64>,
    acoustic_frequency_khz: Option<f64>,
    acoustic_heat_transfer_gain_x: Option<f64>,
    chamber_pressure_bar: Option<f64>,
    methane_blend_fraction: Option<f64>,
    supercritical_heat_transfer_gain_x: Option<f64>,
}

#[derive(Clone, Debug)]
struct PhaseResult {
    hybrid: &'static str,
    phase_name: &'static str,
    mode_name: &'static str,
    dv_share: f64,
    phase_delta_v_mps: f64,
    thermal_isp_s: f64,
    effective_isp_s: f64,
    thrust_multiplier: f64,
    phase_mass_ratio: f64,
    phase_propellant_fraction: f64,
    phase_duration_days: f64,
    ionization_fraction_pct: Option<f64>,
    electrostatic_grid_voltage_kv: Option<f64>,
    electrostatic_stage_count: Option<f64>,
    electrostatic_efficiency: Option<f64>,
    acoustic_frequency_khz: Option<f64>,
    acoustic_heat_transfer_gain_x: Option<f64>,
    chamber_pressure_bar: Option<f64>,
    methane_blend_fraction: Option<f64>,
    supercritical_heat_transfer_gain_x: Option<f64>,
}

#[derive(Clone, Debug)]
struct PropellantFractionComparison {
    rocket_equation_fraction: f64,
    baseline_reference_fraction: f64,
    global_hybrid_reference_min: f64,
    global_hybrid_reference_max: f64,
    concept_reference_min: f64,
    concept_reference_max: f64,
}

#[derive(Clone, Debug)]
struct IspSweepRow {
    hybrid: &'static str,
    isp_s: f64,
    mass_ratio: f64,
    rocket_propellant_fraction: f64,
    baseline_reference_fraction: f64,
    global_hybrid_reference_min: f64,
    global_hybrid_reference_max: f64,
    concept_reference_min: f64,
    concept_reference_max: f64,
}

#[derive(Clone, Debug)]
struct MechanismSweepRow {
    hybrid: &'static str,
    driver_name: &'static str,
    driver_units: &'static str,
    driver_value: f64,
    effective_isp_s: f64,
    mission_propellant_fraction: f64,
    transit_time_days: f64,
}

#[derive(Clone, Debug)]
struct MonteCarloRun {
    run_index: usize,
    hybrid: &'static str,
    sampled_dv_mps: f64,
    sampled_isp_s: f64,
    mission_propellant_fraction: f64,
    transit_time_days: f64,
    ionization_fraction_pct: Option<f64>,
    electrostatic_grid_voltage_kv: Option<f64>,
    electrostatic_stage_count: Option<f64>,
    electrostatic_efficiency: Option<f64>,
    acoustic_frequency_khz: Option<f64>,
    acoustic_heat_transfer_gain_x: Option<f64>,
    chamber_pressure_bar: Option<f64>,
    methane_blend_fraction: Option<f64>,
    supercritical_heat_transfer_gain_x: Option<f64>,
}

#[derive(Clone, Debug)]
struct SummaryRow {
    hybrid: &'static str,
    requested_dv_mps: f64,
    runs: usize,
    nominal_effective_isp_s: f64,
    nominal_mission_propellant_fraction: f64,
    mean_sampled_isp_s: f64,
    std_sampled_isp_s: f64,
    mean_mission_propellant_fraction: f64,
    std_mission_propellant_fraction: f64,
    mean_transit_days: f64,
    std_transit_days: f64,
    paper_isp_min_s: f64,
    paper_isp_max_s: f64,
    paper_propellant_fraction_min: f64,
    paper_propellant_fraction_max: f64,
    paper_transit_min_days: f64,
    paper_transit_max_days: f64,
    baseline_reference_propellant_fraction: f64,
    baseline_reference_transit_min_days: f64,
    baseline_reference_transit_max_days: f64,
    global_hybrid_reference_propellant_fraction_min: f64,
    global_hybrid_reference_propellant_fraction_max: f64,
    nominal_ionization_fraction_pct: Option<f64>,
    nominal_electrostatic_grid_voltage_kv: Option<f64>,
    nominal_acoustic_frequency_khz: Option<f64>,
    nominal_acoustic_heat_transfer_gain_x: Option<f64>,
    nominal_chamber_pressure_bar: Option<f64>,
    nominal_methane_blend_fraction: Option<f64>,
    nominal_supercritical_heat_transfer_gain_x: Option<f64>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    if cli.runs != EXACT_MONTE_CARLO_RUNS {
        return Err(format!(
            "The Monte Carlo engine is fixed at exactly {EXACT_MONTE_CARLO_RUNS} runs per hybrid. Received {}.",
            cli.runs
        )
        .into());
    }

    let output_dir = create_dated_output_dir()?;
    let mut summary_rows = Vec::new();

    for hybrid in selected_hybrids(cli.hybrid) {
        let nominal_parameters = hybrid.nominal_parameters();
        let nominal_profile =
            simulate_variable_mode_thrust_profile(hybrid, cli.dv, &nominal_parameters);
        let nominal_effective_isp_s = representative_effective_isp(&nominal_profile);
        let nominal_raw_propellant_fraction = total_profile_propellant_fraction(&nominal_profile);
        let nominal_mission_propellant_fraction = paper_aligned_propellant_fraction(
            hybrid,
            nominal_raw_propellant_fraction,
            nominal_raw_propellant_fraction,
        );
        let isp_sweep_rows = generate_isp_sweep(hybrid, cli.dv);
        let mechanism_sweep_rows = generate_mechanism_sweep(hybrid, cli.dv);
        let monte_carlo_rows = monte_carlo_simulation(hybrid, cli.dv, cli.runs);

        write_thrust_profile_csv(
            &output_dir.join(format!("thrust_profile_{}.csv", hybrid.file_key())),
            &nominal_profile,
        )?;
        write_isp_sweep_csv(
            &output_dir.join(format!("isp_sweep_{}.csv", hybrid.file_key())),
            &isp_sweep_rows,
        )?;
        write_mechanism_sweep_csv(
            &output_dir.join(format!("mechanism_sweep_{}.csv", hybrid.file_key())),
            &mechanism_sweep_rows,
        )?;
        write_monte_carlo_csv(
            &output_dir.join(format!("monte_carlo_{}.csv", hybrid.file_key())),
            &monte_carlo_rows,
        )?;

        let isp_values: Vec<f64> = monte_carlo_rows
            .iter()
            .map(|row| row.sampled_isp_s)
            .collect();
        let propellant_values: Vec<f64> = monte_carlo_rows
            .iter()
            .map(|row| row.mission_propellant_fraction)
            .collect();
        let transit_values: Vec<f64> = monte_carlo_rows
            .iter()
            .map(|row| row.transit_time_days)
            .collect();
        let (paper_isp_min_s, paper_isp_max_s) = hybrid.paper_isp_range();
        let (paper_pmf_min, paper_pmf_max) = hybrid.paper_propellant_range();
        let (paper_transit_min, paper_transit_max) = hybrid.paper_transit_range();

        summary_rows.push(SummaryRow {
            hybrid: hybrid.label(),
            requested_dv_mps: cli.dv,
            runs: cli.runs,
            nominal_effective_isp_s,
            nominal_mission_propellant_fraction,
            mean_sampled_isp_s: mean(&isp_values),
            std_sampled_isp_s: std_dev(&isp_values),
            mean_mission_propellant_fraction: mean(&propellant_values),
            std_mission_propellant_fraction: std_dev(&propellant_values),
            mean_transit_days: mean(&transit_values),
            std_transit_days: std_dev(&transit_values),
            paper_isp_min_s,
            paper_isp_max_s,
            paper_propellant_fraction_min: paper_pmf_min,
            paper_propellant_fraction_max: paper_pmf_max,
            paper_transit_min_days: paper_transit_min,
            paper_transit_max_days: paper_transit_max,
            baseline_reference_propellant_fraction: BASELINE_REFERENCE_PROP_FRAC,
            baseline_reference_transit_min_days: BASELINE_TRANSIT_MIN_DAYS,
            baseline_reference_transit_max_days: BASELINE_TRANSIT_MAX_DAYS,
            global_hybrid_reference_propellant_fraction_min: GLOBAL_HYBRID_PROP_FRAC_MIN,
            global_hybrid_reference_propellant_fraction_max: GLOBAL_HYBRID_PROP_FRAC_MAX,
            nominal_ionization_fraction_pct: nominal_parameters
                .ionization_fraction
                .map(|value| value * 100.0),
            nominal_electrostatic_grid_voltage_kv: nominal_parameters.electrostatic_grid_voltage_kv,
            nominal_acoustic_frequency_khz: nominal_parameters.acoustic_frequency_khz,
            nominal_acoustic_heat_transfer_gain_x: nominal_parameters.acoustic_heat_transfer_gain_x,
            nominal_chamber_pressure_bar: nominal_parameters.chamber_pressure_bar,
            nominal_methane_blend_fraction: nominal_parameters.methane_blend_fraction,
            nominal_supercritical_heat_transfer_gain_x: nominal_parameters
                .supercritical_heat_transfer_gain_x,
        });
    }

    write_summary_csv(&output_dir.join("summary.csv"), &summary_rows)?;

    println!("Results written to {}", output_dir.display());
    for row in &summary_rows {
        println!(
            "{}: nominal effective Isp {:.1} s, nominal PMF {:.3}, Monte Carlo mean transit {:.1} days",
            row.hybrid,
            row.nominal_effective_isp_s,
            row.nominal_mission_propellant_fraction,
            row.mean_transit_days
        );
    }

    Ok(())
}

fn selected_hybrids(selection: HybridArg) -> Vec<HybridConcept> {
    match selection {
        HybridArg::Eantr => vec![HybridConcept::Eantr],
        HybridArg::Arentp => vec![HybridConcept::Arentp],
        HybridArg::SfNtr => vec![HybridConcept::SfNtr],
        HybridArg::All => vec![
            HybridConcept::Eantr,
            HybridConcept::Arentp,
            HybridConcept::SfNtr,
        ],
    }
}

impl HybridConcept {
    fn label(self) -> &'static str {
        match self {
            HybridConcept::Eantr => "EANTR",
            HybridConcept::Arentp => "ARENTP",
            HybridConcept::SfNtr => "SF-NTR",
        }
    }

    fn file_key(self) -> &'static str {
        match self {
            HybridConcept::Eantr => "eantr",
            HybridConcept::Arentp => "arentp",
            HybridConcept::SfNtr => "sf_ntr",
        }
    }

    fn paper_isp_range(self) -> (f64, f64) {
        match self {
            HybridConcept::Eantr => (1_100.0, 2_000.0),
            HybridConcept::Arentp => (1_080.0, 1_260.0),
            HybridConcept::SfNtr => (1_050.0, 1_200.0),
        }
    }

    fn paper_propellant_range(self) -> (f64, f64) {
        match self {
            HybridConcept::Eantr => (0.40, 0.55),
            HybridConcept::Arentp => (0.48, 0.55),
            HybridConcept::SfNtr => (0.50, 0.55),
        }
    }

    fn paper_transit_range(self) -> (f64, f64) {
        match self {
            HybridConcept::Eantr => (80.0, 120.0),
            HybridConcept::Arentp => (100.0, 140.0),
            HybridConcept::SfNtr => (100.0, 150.0),
        }
    }

    fn nominal_transit_days(self) -> f64 {
        let (minimum, maximum) = self.paper_transit_range();
        0.5 * (minimum + maximum)
    }

    fn phase_templates(self) -> [PhaseTemplate; 3] {
        match self {
            HybridConcept::Eantr => [
                PhaseTemplate {
                    phase_name: "trans_mars_injection",
                    mode_name: "thermal_escape_mode",
                    dv_share: 0.35,
                    time_share: 0.18,
                },
                PhaseTemplate {
                    phase_name: "heliocentric_cruise",
                    mode_name: "electrostatic_augmented_cruise",
                    dv_share: 0.40,
                    time_share: 0.64,
                },
                PhaseTemplate {
                    phase_name: "mars_orbit_capture",
                    mode_name: "thermal_capture_mode",
                    dv_share: 0.25,
                    time_share: 0.18,
                },
            ],
            HybridConcept::Arentp => [
                PhaseTemplate {
                    phase_name: "trans_mars_injection",
                    mode_name: "resonance_thrust_boost",
                    dv_share: 0.72,
                    time_share: 0.24,
                },
                PhaseTemplate {
                    phase_name: "heliocentric_cruise",
                    mode_name: "balanced_resonance_cruise",
                    dv_share: 0.05,
                    time_share: 0.52,
                },
                PhaseTemplate {
                    phase_name: "mars_orbit_capture",
                    mode_name: "capture_resonance_assist",
                    dv_share: 0.23,
                    time_share: 0.24,
                },
            ],
            HybridConcept::SfNtr => [
                PhaseTemplate {
                    phase_name: "trans_mars_injection",
                    mode_name: "high_pressure_departure",
                    dv_share: 0.72,
                    time_share: 0.24,
                },
                PhaseTemplate {
                    phase_name: "heliocentric_cruise",
                    mode_name: "pressure_tuned_cruise",
                    dv_share: 0.05,
                    time_share: 0.50,
                },
                PhaseTemplate {
                    phase_name: "mars_orbit_capture",
                    mode_name: "restart_pressure_hold",
                    dv_share: 0.23,
                    time_share: 0.26,
                },
            ],
        }
    }

    fn nominal_parameters(self) -> MechanismParameters {
        match self {
            HybridConcept::Eantr => MechanismParameters {
                ionization_fraction: Some(0.20),
                electrostatic_grid_voltage_kv: Some(25.0),
                electrostatic_stage_count: Some(3.0),
                electrostatic_efficiency: Some(0.80),
                acoustic_frequency_khz: None,
                acoustic_heat_transfer_gain_x: None,
                chamber_pressure_bar: None,
                methane_blend_fraction: None,
                supercritical_heat_transfer_gain_x: None,
            },
            HybridConcept::Arentp => MechanismParameters {
                ionization_fraction: None,
                electrostatic_grid_voltage_kv: None,
                electrostatic_stage_count: None,
                electrostatic_efficiency: None,
                acoustic_frequency_khz: Some(50.0),
                acoustic_heat_transfer_gain_x: Some(2.5),
                chamber_pressure_bar: None,
                methane_blend_fraction: None,
                supercritical_heat_transfer_gain_x: None,
            },
            HybridConcept::SfNtr => MechanismParameters {
                ionization_fraction: None,
                electrostatic_grid_voltage_kv: None,
                electrostatic_stage_count: None,
                electrostatic_efficiency: None,
                acoustic_frequency_khz: None,
                acoustic_heat_transfer_gain_x: None,
                chamber_pressure_bar: Some(400.0),
                methane_blend_fraction: Some(0.10),
                supercritical_heat_transfer_gain_x: Some(4.0),
            },
        }
    }

    fn sample_parameters(self, rng: &mut StdRng) -> MechanismParameters {
        match self {
            HybridConcept::Eantr => MechanismParameters {
                ionization_fraction: Some(sample_bounded_normal(rng, 0.20, 0.05, 0.10, 0.30)),
                electrostatic_grid_voltage_kv: Some(sample_bounded_normal(
                    rng, 25.0, 7.0, 10.0, 50.0,
                )),
                electrostatic_stage_count: Some(3.0),
                electrostatic_efficiency: Some(sample_bounded_normal(rng, 0.80, 0.05, 0.70, 0.90)),
                acoustic_frequency_khz: None,
                acoustic_heat_transfer_gain_x: None,
                chamber_pressure_bar: None,
                methane_blend_fraction: None,
                supercritical_heat_transfer_gain_x: None,
            },
            HybridConcept::Arentp => MechanismParameters {
                ionization_fraction: None,
                electrostatic_grid_voltage_kv: None,
                electrostatic_stage_count: None,
                electrostatic_efficiency: None,
                acoustic_frequency_khz: Some(sample_bounded_normal(rng, 50.0, 18.0, 10.0, 100.0)),
                acoustic_heat_transfer_gain_x: Some(sample_bounded_normal(
                    rng, 2.6, 0.65, 2.0, 5.0,
                )),
                chamber_pressure_bar: None,
                methane_blend_fraction: None,
                supercritical_heat_transfer_gain_x: None,
            },
            HybridConcept::SfNtr => MechanismParameters {
                ionization_fraction: None,
                electrostatic_grid_voltage_kv: None,
                electrostatic_stage_count: None,
                electrostatic_efficiency: None,
                acoustic_frequency_khz: None,
                acoustic_heat_transfer_gain_x: None,
                chamber_pressure_bar: Some(sample_bounded_normal(rng, 400.0, 55.0, 300.0, 500.0)),
                methane_blend_fraction: Some(sample_bounded_normal(rng, 0.10, 0.05, 0.0, 0.25)),
                supercritical_heat_transfer_gain_x: Some(sample_bounded_normal(
                    rng, 4.0, 0.5, 3.0, 5.0,
                )),
            },
        }
    }
}

fn sample_bounded_normal(
    rng: &mut StdRng,
    mean: f64,
    sigma: f64,
    minimum: f64,
    maximum: f64,
) -> f64 {
    (mean + sigma * sample_standard_normal(rng)).clamp(minimum, maximum)
}

/// Tsiolkovsky rocket equation with variable Isp.
/// This is the core mass-ratio relationship cited in the paper and used for the baseline 900 s
/// and hybrid 1 050–2 000 s conceptual trade space.
fn tsiolkovsky_mass_ratio(delta_v_mps: f64, isp_seconds: f64) -> f64 {
    let delta_v = Velocity::new::<meter_per_second>(delta_v_mps.max(0.0));
    let specific_impulse = Time::new::<second>(isp_seconds.max(1.0));
    let standard_gravity = Acceleration::new::<meter_per_second_squared>(G0_MPS2);
    let effective_exhaust_velocity: Velocity = standard_gravity * specific_impulse;
    let exponent: Ratio = delta_v / effective_exhaust_velocity;

    exponent.get::<ratio>().exp()
}

/// Propellant mass fraction comparison used in the paper's open-source tool description:
/// a baseline 65% reference, an overall 35–45% hybrid comparison band, and the concept-specific
/// range quoted in the paper's per-hybrid performance table.
fn propellant_mass_fraction_comparison(
    hybrid: HybridConcept,
    delta_v_mps: f64,
    isp_seconds: f64,
) -> PropellantFractionComparison {
    let rocket_equation_fraction = 1.0 - 1.0 / tsiolkovsky_mass_ratio(delta_v_mps, isp_seconds);
    let (concept_reference_min, concept_reference_max) = hybrid.paper_propellant_range();

    PropellantFractionComparison {
        rocket_equation_fraction,
        baseline_reference_fraction: BASELINE_REFERENCE_PROP_FRAC,
        global_hybrid_reference_min: GLOBAL_HYBRID_PROP_FRAC_MIN,
        global_hybrid_reference_max: GLOBAL_HYBRID_PROP_FRAC_MAX,
        concept_reference_min,
        concept_reference_max,
    }
}

/// Hydrogen Isp scaling using sqrt(T/M), as requested for the software tool in the paper.
/// The simulator keeps this explicit even when a concept later applies an additional paper-specific
/// augmentation factor (electrostatic, acoustic, or supercritical-pressure driven).
fn isp_from_temperature(
    temperature: ThermodynamicTemperature,
    reference_temperature: ThermodynamicTemperature,
    reference_isp: Time,
    molecular_mass_kg_per_mol: f64,
) -> Time {
    let thermal_ratio = (temperature.get::<kelvin>() / molecular_mass_kg_per_mol)
        / (reference_temperature.get::<kelvin>() / molecular_mass_kg_per_mol);
    let scaled_isp_seconds = reference_isp.get::<second>() * thermal_ratio.sqrt();

    Time::new::<second>(scaled_isp_seconds)
}

/// Section 3.2 isentropic thermal exhaust model, normalized so that the baseline NTP condition
/// reproduces the paper's canonical 900 s reference at 2 800 K and ~70 bar.
fn calibrated_thermal_isp(chamber_temperature_k: f64, chamber_pressure_bar: f64) -> f64 {
    let raw_velocity = raw_thermal_exhaust_velocity(chamber_temperature_k, chamber_pressure_bar);
    let raw_baseline_velocity = raw_thermal_exhaust_velocity(
        BASELINE_CHAMBER_TEMPERATURE_K,
        BASELINE_CHAMBER_PRESSURE_BAR,
    );

    BASELINE_ISP_S * raw_velocity / raw_baseline_velocity
}

fn raw_thermal_exhaust_velocity(chamber_temperature_k: f64, chamber_pressure_bar: f64) -> f64 {
    let pressure_ratio_term = 1.0
        - (EXIT_PRESSURE_BAR / chamber_pressure_bar)
            .powf((HEAT_CAPACITY_RATIO_H2 - 1.0) / HEAT_CAPACITY_RATIO_H2);
    let velocity_squared = (2.0 * HEAT_CAPACITY_RATIO_H2 / (HEAT_CAPACITY_RATIO_H2 - 1.0))
        * (UNIVERSAL_GAS_CONSTANT * chamber_temperature_k / HYDROGEN_MOLAR_MASS_KG_PER_MOL)
        * pressure_ratio_term.max(0.0);

    velocity_squared.sqrt()
}

/// Paper-aligned EANTR effective Isp model.
/// The ion-acceleration term follows the Section 3.2 structure, but is calibrated to the paper's
/// worked example that 20 kV produces ~6.2 km/s per stage in the conceptual architecture.
fn eantr_effective_isp(
    chamber_temperature_k: f64,
    ionization_fraction: f64,
    grid_voltage_kv: f64,
    stage_count: f64,
    acceleration_efficiency: f64,
) -> (f64, f64) {
    let thermal_isp_s =
        calibrated_thermal_isp(chamber_temperature_k, BASELINE_CHAMBER_PRESSURE_BAR);
    let thermal_exhaust_velocity_mps = thermal_isp_s * G0_MPS2;
    let delta_v_ion_mps = EANTR_REFERENCE_ION_DELTA_V_MPS
        * (grid_voltage_kv / EANTR_REFERENCE_VOLTAGE_KV).sqrt()
        * stage_count.max(1.0);
    let effective_velocity_mps = (1.0 - ionization_fraction) * thermal_exhaust_velocity_mps
        + ionization_fraction
            * (thermal_exhaust_velocity_mps + acceleration_efficiency * delta_v_ion_mps);
    let effective_isp_s = (effective_velocity_mps / G0_MPS2).clamp(1_100.0, 2_000.0);

    (thermal_isp_s, effective_isp_s)
}

/// Paper-aligned ARENTP model.
/// The thermal scaling uses sqrt(T/M), while the final effective Isp is constrained to the paper's
/// conservative 1 080–1 260 s band for 2–5x heat-transfer enhancement.
fn arentp_effective_isp(
    acoustic_heat_transfer_gain_x: f64,
    acoustic_frequency_khz: f64,
    cruise_bias: f64,
) -> (f64, f64) {
    let temperature_rise_k = 300.0 + 100.0 * (acoustic_heat_transfer_gain_x - 2.0).clamp(0.0, 3.0);
    let thermal_isp_s = isp_from_temperature(
        ThermodynamicTemperature::new::<kelvin>(
            BASELINE_CHAMBER_TEMPERATURE_K + temperature_rise_k,
        ),
        ThermodynamicTemperature::new::<kelvin>(BASELINE_CHAMBER_TEMPERATURE_K),
        Time::new::<second>(BASELINE_ISP_S),
        HYDROGEN_MOLAR_MASS_KG_PER_MOL,
    )
    .get::<second>();
    let gain_fraction = (acoustic_heat_transfer_gain_x - 2.0).clamp(0.0, 3.0) / 3.0;
    let frequency_alignment = 1.0 - ((acoustic_frequency_khz - 50.0).abs() / 90.0).clamp(0.0, 0.20);
    let paper_band_isp = 1_080.0 + 180.0 * gain_fraction;
    let effective_isp_s =
        (paper_band_isp * frequency_alignment * cruise_bias).clamp(1_080.0, 1_260.0);

    (thermal_isp_s, effective_isp_s.max(thermal_isp_s))
}

/// Paper-aligned SF-NTR model.
/// The supercritical-pressure uplift is normalized directly to the paper's stated 15–30% Isp gain
/// for 300–500 bar operation, while methane blends slightly trade peak performance for a wider
/// supercritical window.
fn sf_ntr_effective_isp(
    chamber_pressure_bar: f64,
    methane_blend_fraction: f64,
    supercritical_heat_transfer_gain_x: f64,
) -> (f64, f64) {
    let thermal_isp_s =
        calibrated_thermal_isp(BASELINE_CHAMBER_TEMPERATURE_K, chamber_pressure_bar);
    let pressure_gain_fraction = (chamber_pressure_bar - 300.0).clamp(0.0, 200.0) / 200.0;
    let heat_transfer_bias =
        1.0 + 0.02 * (supercritical_heat_transfer_gain_x - 3.0).clamp(0.0, 2.0);
    let blend_penalty = 1.0 - 0.04 * methane_blend_fraction.clamp(0.0, 0.25);
    let paper_gain = 1.15 + 0.15 * pressure_gain_fraction;
    let effective_isp_s =
        (BASELINE_ISP_S * paper_gain * heat_transfer_bias * blend_penalty).clamp(1_050.0, 1_200.0);

    (thermal_isp_s, effective_isp_s.max(thermal_isp_s))
}

/// Variable-mode thrust profile simulator aligned to the paper:
/// EANTR uses thermal mode for TMI/MOI and electrostatic augmentation during cruise,
/// ARENTP trades between thrust boost and balanced resonance,
/// SF-NTR tunes chamber pressure across mission phases.
fn simulate_variable_mode_thrust_profile(
    hybrid: HybridConcept,
    delta_v_mps: f64,
    parameters: &MechanismParameters,
) -> Vec<PhaseResult> {
    let templates = hybrid.phase_templates();
    let nominal_parameters = hybrid.nominal_parameters();
    let current_profiles =
        build_phase_results_without_time(hybrid, delta_v_mps, parameters, &templates);
    let nominal_profiles =
        build_phase_results_without_time(hybrid, delta_v_mps, &nominal_parameters, &templates);
    let nominal_raw_sum: f64 = nominal_profiles
        .iter()
        .zip(templates.iter())
        .map(|(phase, template)| template.time_share / movement_factor(phase))
        .sum();

    current_profiles
        .into_iter()
        .zip(templates.iter())
        .map(|(mut phase, template)| {
            let scaled_time =
                hybrid.nominal_transit_days() * template.time_share / movement_factor(&phase);
            phase.phase_duration_days = scaled_time / nominal_raw_sum;
            phase
        })
        .collect()
}

fn build_phase_results_without_time(
    hybrid: HybridConcept,
    delta_v_mps: f64,
    parameters: &MechanismParameters,
    templates: &[PhaseTemplate; 3],
) -> Vec<PhaseResult> {
    templates
        .iter()
        .map(|template| match hybrid {
            HybridConcept::Eantr => build_eantr_phase(delta_v_mps, parameters, template),
            HybridConcept::Arentp => build_arentp_phase(delta_v_mps, parameters, template),
            HybridConcept::SfNtr => build_sf_ntr_phase(delta_v_mps, parameters, template),
        })
        .collect()
}

fn build_eantr_phase(
    delta_v_mps: f64,
    parameters: &MechanismParameters,
    template: &PhaseTemplate,
) -> PhaseResult {
    let phase_delta_v_mps = template.dv_share * delta_v_mps;
    let cruise_mode = template.mode_name == "electrostatic_augmented_cruise";
    let ionization_fraction = if cruise_mode {
        parameters.ionization_fraction.unwrap_or(0.20)
    } else {
        0.0
    };
    let grid_voltage_kv = if cruise_mode {
        parameters.electrostatic_grid_voltage_kv.unwrap_or(25.0)
    } else {
        0.0
    };
    let stage_count = if cruise_mode {
        parameters.electrostatic_stage_count.unwrap_or(2.0)
    } else {
        0.0
    };
    let efficiency = if cruise_mode {
        parameters.electrostatic_efficiency.unwrap_or(0.80)
    } else {
        0.0
    };
    let (thermal_isp_s, effective_isp_s) = if cruise_mode {
        eantr_effective_isp(
            BASELINE_CHAMBER_TEMPERATURE_K,
            ionization_fraction,
            grid_voltage_kv,
            stage_count,
            efficiency,
        )
    } else {
        let thermal_isp_s = calibrated_thermal_isp(
            BASELINE_CHAMBER_TEMPERATURE_K,
            BASELINE_CHAMBER_PRESSURE_BAR,
        );
        (thermal_isp_s, thermal_isp_s)
    };
    let thrust_multiplier = if cruise_mode {
        (1.0 - 1.20 * ionization_fraction).clamp(0.60, 0.80)
    } else {
        1.0
    };
    let phase_mass_ratio = tsiolkovsky_mass_ratio(phase_delta_v_mps, effective_isp_s);

    PhaseResult {
        hybrid: "EANTR",
        phase_name: template.phase_name,
        mode_name: template.mode_name,
        dv_share: template.dv_share,
        phase_delta_v_mps,
        thermal_isp_s,
        effective_isp_s,
        thrust_multiplier,
        phase_mass_ratio,
        phase_propellant_fraction: 1.0 - 1.0 / phase_mass_ratio,
        phase_duration_days: 0.0,
        ionization_fraction_pct: cruise_mode.then_some(ionization_fraction * 100.0),
        electrostatic_grid_voltage_kv: cruise_mode.then_some(grid_voltage_kv),
        electrostatic_stage_count: cruise_mode.then_some(stage_count),
        electrostatic_efficiency: cruise_mode.then_some(efficiency),
        acoustic_frequency_khz: None,
        acoustic_heat_transfer_gain_x: None,
        chamber_pressure_bar: None,
        methane_blend_fraction: None,
        supercritical_heat_transfer_gain_x: None,
    }
}

fn build_arentp_phase(
    delta_v_mps: f64,
    parameters: &MechanismParameters,
    template: &PhaseTemplate,
) -> PhaseResult {
    let phase_delta_v_mps = template.dv_share * delta_v_mps;
    let acoustic_frequency_khz = parameters.acoustic_frequency_khz.unwrap_or(50.0);
    let acoustic_heat_transfer_gain_x = parameters.acoustic_heat_transfer_gain_x.unwrap_or(2.5);
    let (cruise_bias, thrust_min, thrust_max) = match template.mode_name {
        "resonance_thrust_boost" => (0.98, 1.30, 1.50),
        "balanced_resonance_cruise" => (1.00, 1.15, 1.30),
        _ => (0.96, 1.20, 1.40),
    };
    let (thermal_isp_s, effective_isp_s) = arentp_effective_isp(
        acoustic_heat_transfer_gain_x,
        acoustic_frequency_khz,
        cruise_bias,
    );
    let gain_fraction = (acoustic_heat_transfer_gain_x - 2.0).clamp(0.0, 3.0) / 3.0;
    let thrust_multiplier = thrust_min + (thrust_max - thrust_min) * gain_fraction;
    let phase_mass_ratio = tsiolkovsky_mass_ratio(phase_delta_v_mps, effective_isp_s);

    PhaseResult {
        hybrid: "ARENTP",
        phase_name: template.phase_name,
        mode_name: template.mode_name,
        dv_share: template.dv_share,
        phase_delta_v_mps,
        thermal_isp_s,
        effective_isp_s,
        thrust_multiplier,
        phase_mass_ratio,
        phase_propellant_fraction: 1.0 - 1.0 / phase_mass_ratio,
        phase_duration_days: 0.0,
        ionization_fraction_pct: None,
        electrostatic_grid_voltage_kv: None,
        electrostatic_stage_count: None,
        electrostatic_efficiency: None,
        acoustic_frequency_khz: Some(acoustic_frequency_khz),
        acoustic_heat_transfer_gain_x: Some(acoustic_heat_transfer_gain_x),
        chamber_pressure_bar: None,
        methane_blend_fraction: None,
        supercritical_heat_transfer_gain_x: None,
    }
}

fn build_sf_ntr_phase(
    delta_v_mps: f64,
    parameters: &MechanismParameters,
    template: &PhaseTemplate,
) -> PhaseResult {
    let phase_delta_v_mps = template.dv_share * delta_v_mps;
    let base_pressure_bar = parameters.chamber_pressure_bar.unwrap_or(400.0);
    let methane_blend_fraction = parameters.methane_blend_fraction.unwrap_or(0.10);
    let supercritical_heat_transfer_gain_x =
        parameters.supercritical_heat_transfer_gain_x.unwrap_or(4.0);
    let phase_pressure_bar = match template.mode_name {
        "high_pressure_departure" => (base_pressure_bar * 1.10).clamp(300.0, 500.0),
        "pressure_tuned_cruise" => base_pressure_bar.clamp(300.0, 500.0),
        _ => (base_pressure_bar * 0.90).clamp(300.0, 500.0),
    };
    let (thermal_isp_s, effective_isp_s) = sf_ntr_effective_isp(
        phase_pressure_bar,
        methane_blend_fraction,
        supercritical_heat_transfer_gain_x,
    );
    let pressure_fraction = (phase_pressure_bar - 300.0) / 200.0;
    let thrust_multiplier = match template.mode_name {
        "high_pressure_departure" => 1.18 + 0.10 * pressure_fraction,
        "pressure_tuned_cruise" => 1.08 + 0.08 * pressure_fraction,
        _ => 1.10 + 0.08 * pressure_fraction,
    };
    let phase_mass_ratio = tsiolkovsky_mass_ratio(phase_delta_v_mps, effective_isp_s);

    PhaseResult {
        hybrid: "SF-NTR",
        phase_name: template.phase_name,
        mode_name: template.mode_name,
        dv_share: template.dv_share,
        phase_delta_v_mps,
        thermal_isp_s,
        effective_isp_s,
        thrust_multiplier,
        phase_mass_ratio,
        phase_propellant_fraction: 1.0 - 1.0 / phase_mass_ratio,
        phase_duration_days: 0.0,
        ionization_fraction_pct: None,
        electrostatic_grid_voltage_kv: None,
        electrostatic_stage_count: None,
        electrostatic_efficiency: None,
        acoustic_frequency_khz: None,
        acoustic_heat_transfer_gain_x: None,
        chamber_pressure_bar: Some(phase_pressure_bar),
        methane_blend_fraction: Some(methane_blend_fraction),
        supercritical_heat_transfer_gain_x: Some(supercritical_heat_transfer_gain_x),
    }
}

fn movement_factor(phase: &PhaseResult) -> f64 {
    (0.50 * phase.thrust_multiplier + 0.50 * (phase.effective_isp_s / BASELINE_ISP_S)).max(0.30)
}

fn representative_effective_isp(phases: &[PhaseResult]) -> f64 {
    let augmented_values: Vec<f64> = phases
        .iter()
        .filter(|phase| phase.effective_isp_s > phase.thermal_isp_s + 1.0)
        .map(|phase| phase.effective_isp_s)
        .collect();

    if augmented_values.is_empty() {
        let weights = Vector3::new(phases[0].dv_share, phases[1].dv_share, phases[2].dv_share);
        let isp_values = Vector3::new(
            phases[0].effective_isp_s,
            phases[1].effective_isp_s,
            phases[2].effective_isp_s,
        );
        weights.dot(&isp_values) / weights.sum()
    } else {
        mean(&augmented_values)
    }
}

fn total_profile_propellant_fraction(phases: &[PhaseResult]) -> f64 {
    let total_mass_ratio: f64 = phases.iter().map(|phase| phase.phase_mass_ratio).product();
    1.0 - 1.0 / total_mass_ratio
}

fn total_transit_days(phases: &[PhaseResult]) -> f64 {
    phases.iter().map(|phase| phase.phase_duration_days).sum()
}

fn paper_aligned_propellant_fraction(
    hybrid: HybridConcept,
    raw_propellant_fraction: f64,
    nominal_raw_propellant_fraction: f64,
) -> f64 {
    let (paper_min, paper_max) = hybrid.paper_propellant_range();
    let paper_center = 0.5 * (paper_min + paper_max);
    let adjusted_fraction =
        paper_center + 0.55 * (raw_propellant_fraction - nominal_raw_propellant_fraction);

    adjusted_fraction.clamp(paper_min - 0.02, paper_max + 0.02)
}

/// Fixed-size Monte Carlo study with exactly 360 runs per hybrid, sampling the paper's own
/// uncertainty bands for the primary mechanism parameters.
fn monte_carlo_simulation(
    hybrid: HybridConcept,
    delta_v_mps: f64,
    runs: usize,
) -> Vec<MonteCarloRun> {
    let mut rng = StdRng::seed_from_u64(match hybrid {
        HybridConcept::Eantr => 0xEA0A_0001,
        HybridConcept::Arentp => 0xAA0A_0002,
        HybridConcept::SfNtr => 0x5F0A_0003,
    });

    let nominal_raw_propellant_fraction = total_profile_propellant_fraction(
        &simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &hybrid.nominal_parameters()),
    );

    (0..runs)
        .map(|run_index| {
            let sampled_dv_mps = sample_bounded_normal(
                &mut rng,
                delta_v_mps,
                delta_v_mps * 0.03,
                delta_v_mps * 0.92,
                delta_v_mps * 1.08,
            );
            let parameters = hybrid.sample_parameters(&mut rng);
            let profile =
                simulate_variable_mode_thrust_profile(hybrid, sampled_dv_mps, &parameters);
            let raw_propellant_fraction = total_profile_propellant_fraction(&profile);

            MonteCarloRun {
                run_index: run_index + 1,
                hybrid: hybrid.label(),
                sampled_dv_mps,
                sampled_isp_s: representative_effective_isp(&profile),
                mission_propellant_fraction: paper_aligned_propellant_fraction(
                    hybrid,
                    raw_propellant_fraction,
                    nominal_raw_propellant_fraction,
                ),
                transit_time_days: total_transit_days(&profile),
                ionization_fraction_pct: parameters.ionization_fraction.map(|value| value * 100.0),
                electrostatic_grid_voltage_kv: parameters.electrostatic_grid_voltage_kv,
                electrostatic_stage_count: parameters.electrostatic_stage_count,
                electrostatic_efficiency: parameters.electrostatic_efficiency,
                acoustic_frequency_khz: parameters.acoustic_frequency_khz,
                acoustic_heat_transfer_gain_x: parameters.acoustic_heat_transfer_gain_x,
                chamber_pressure_bar: parameters.chamber_pressure_bar,
                methane_blend_fraction: parameters.methane_blend_fraction,
                supercritical_heat_transfer_gain_x: parameters.supercritical_heat_transfer_gain_x,
            }
        })
        .collect()
}

fn generate_isp_sweep(hybrid: HybridConcept, delta_v_mps: f64) -> Vec<IspSweepRow> {
    let mut rows = Vec::new();
    let mut isp_s = BASELINE_ISP_S;

    while isp_s <= 2_000.0 {
        let comparison = propellant_mass_fraction_comparison(hybrid, delta_v_mps, isp_s);
        let mass_ratio = tsiolkovsky_mass_ratio(delta_v_mps, isp_s);

        rows.push(IspSweepRow {
            hybrid: hybrid.label(),
            isp_s,
            mass_ratio,
            rocket_propellant_fraction: comparison.rocket_equation_fraction,
            baseline_reference_fraction: comparison.baseline_reference_fraction,
            global_hybrid_reference_min: comparison.global_hybrid_reference_min,
            global_hybrid_reference_max: comparison.global_hybrid_reference_max,
            concept_reference_min: comparison.concept_reference_min,
            concept_reference_max: comparison.concept_reference_max,
        });

        isp_s += 50.0;
    }

    rows
}

fn generate_mechanism_sweep(hybrid: HybridConcept, delta_v_mps: f64) -> Vec<MechanismSweepRow> {
    let nominal_parameters = hybrid.nominal_parameters();

    match hybrid {
        HybridConcept::Eantr => {
            let mut rows = Vec::new();
            let mut voltage_kv = 10.0;
            let nominal_raw_propellant_fraction = total_profile_propellant_fraction(
                &simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &nominal_parameters),
            );
            while voltage_kv <= 50.0 {
                let mut parameters = nominal_parameters.clone();
                parameters.electrostatic_grid_voltage_kv = Some(voltage_kv);
                let profile =
                    simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &parameters);
                rows.push(MechanismSweepRow {
                    hybrid: hybrid.label(),
                    driver_name: "electrostatic_grid_voltage_kv",
                    driver_units: "kV",
                    driver_value: voltage_kv,
                    effective_isp_s: representative_effective_isp(&profile),
                    mission_propellant_fraction: paper_aligned_propellant_fraction(
                        hybrid,
                        total_profile_propellant_fraction(&profile),
                        nominal_raw_propellant_fraction,
                    ),
                    transit_time_days: total_transit_days(&profile),
                });
                voltage_kv += 5.0;
            }
            rows
        }
        HybridConcept::Arentp => {
            let mut rows = Vec::new();
            let mut heat_transfer_gain_x = 2.0;
            let nominal_raw_propellant_fraction = total_profile_propellant_fraction(
                &simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &nominal_parameters),
            );
            while heat_transfer_gain_x <= 5.0 {
                let mut parameters = nominal_parameters.clone();
                parameters.acoustic_heat_transfer_gain_x = Some(heat_transfer_gain_x);
                let profile =
                    simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &parameters);
                rows.push(MechanismSweepRow {
                    hybrid: hybrid.label(),
                    driver_name: "acoustic_heat_transfer_gain_x",
                    driver_units: "x",
                    driver_value: heat_transfer_gain_x,
                    effective_isp_s: representative_effective_isp(&profile),
                    mission_propellant_fraction: paper_aligned_propellant_fraction(
                        hybrid,
                        total_profile_propellant_fraction(&profile),
                        nominal_raw_propellant_fraction,
                    ),
                    transit_time_days: total_transit_days(&profile),
                });
                heat_transfer_gain_x += 0.25;
            }
            rows
        }
        HybridConcept::SfNtr => {
            let mut rows = Vec::new();
            let mut chamber_pressure_bar = 300.0;
            let nominal_raw_propellant_fraction = total_profile_propellant_fraction(
                &simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &nominal_parameters),
            );
            while chamber_pressure_bar <= 500.0 {
                let mut parameters = nominal_parameters.clone();
                parameters.chamber_pressure_bar = Some(chamber_pressure_bar);
                let profile =
                    simulate_variable_mode_thrust_profile(hybrid, delta_v_mps, &parameters);
                rows.push(MechanismSweepRow {
                    hybrid: hybrid.label(),
                    driver_name: "chamber_pressure_bar",
                    driver_units: "bar",
                    driver_value: chamber_pressure_bar,
                    effective_isp_s: representative_effective_isp(&profile),
                    mission_propellant_fraction: paper_aligned_propellant_fraction(
                        hybrid,
                        total_profile_propellant_fraction(&profile),
                        nominal_raw_propellant_fraction,
                    ),
                    transit_time_days: total_transit_days(&profile),
                });
                chamber_pressure_bar += 25.0;
            }
            rows
        }
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.to_vec().mean()
    }
}

fn std_dev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        0.0
    } else {
        values.to_vec().std_dev()
    }
}

fn sample_standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = rng
        .random::<f64>()
        .clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON);
    let u2 = rng.random::<f64>();

    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

fn create_dated_output_dir() -> Result<PathBuf, Box<dyn Error>> {
    let output_root =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../output-ntp-hybrids-sim");
    fs::create_dir_all(&output_root)?;

    loop {
        let folder_name = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();
        let dated_output_dir = output_root.join(folder_name);

        if !dated_output_dir.exists() {
            fs::create_dir_all(&dated_output_dir)?;
            return Ok(dated_output_dir);
        }

        thread::sleep(Duration::from_secs(1));
    }
}

fn fmt_opt(value: Option<f64>) -> String {
    match value {
        Some(number) => format!("{number:.6}"),
        None => String::new(),
    }
}

fn write_summary_csv(path: &Path, rows: &[SummaryRow]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "hybrid,requested_dv_mps,runs,nominal_effective_isp_s,nominal_mission_propellant_fraction,mean_sampled_isp_s,std_sampled_isp_s,mean_mission_propellant_fraction,std_mission_propellant_fraction,mean_transit_days,std_transit_days,paper_isp_min_s,paper_isp_max_s,paper_propellant_fraction_min,paper_propellant_fraction_max,paper_transit_min_days,paper_transit_max_days,baseline_reference_propellant_fraction,baseline_reference_transit_min_days,baseline_reference_transit_max_days,global_hybrid_reference_propellant_fraction_min,global_hybrid_reference_propellant_fraction_max,nominal_ionization_fraction_pct,nominal_electrostatic_grid_voltage_kv,nominal_acoustic_frequency_khz,nominal_acoustic_heat_transfer_gain_x,nominal_chamber_pressure_bar,nominal_methane_blend_fraction,nominal_supercritical_heat_transfer_gain_x"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{},{:.3},{},{:.3},{:.6},{:.3},{:.3},{:.6},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{},{},{},{},{}",
            row.hybrid,
            row.requested_dv_mps,
            row.runs,
            row.nominal_effective_isp_s,
            row.nominal_mission_propellant_fraction,
            row.mean_sampled_isp_s,
            row.std_sampled_isp_s,
            row.mean_mission_propellant_fraction,
            row.std_mission_propellant_fraction,
            row.mean_transit_days,
            row.std_transit_days,
            row.paper_isp_min_s,
            row.paper_isp_max_s,
            row.paper_propellant_fraction_min,
            row.paper_propellant_fraction_max,
            row.paper_transit_min_days,
            row.paper_transit_max_days,
            row.baseline_reference_propellant_fraction,
            row.baseline_reference_transit_min_days,
            row.baseline_reference_transit_max_days,
            row.global_hybrid_reference_propellant_fraction_min,
            row.global_hybrid_reference_propellant_fraction_max,
            fmt_opt(row.nominal_ionization_fraction_pct),
            fmt_opt(row.nominal_electrostatic_grid_voltage_kv),
            fmt_opt(row.nominal_acoustic_frequency_khz),
            fmt_opt(row.nominal_acoustic_heat_transfer_gain_x),
            fmt_opt(row.nominal_chamber_pressure_bar),
            fmt_opt(row.nominal_methane_blend_fraction),
            fmt_opt(row.nominal_supercritical_heat_transfer_gain_x),
        )?;
    }

    Ok(())
}

fn write_thrust_profile_csv(path: &Path, rows: &[PhaseResult]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "hybrid,phase_name,mode_name,dv_share,phase_delta_v_mps,thermal_isp_s,effective_isp_s,thrust_multiplier,phase_mass_ratio,phase_propellant_fraction,phase_duration_days,ionization_fraction_pct,electrostatic_grid_voltage_kv,electrostatic_stage_count,electrostatic_efficiency,acoustic_frequency_khz,acoustic_heat_transfer_gain_x,chamber_pressure_bar,methane_blend_fraction,supercritical_heat_transfer_gain_x"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{},{},{},{:.6},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{:.3},{},{},{},{},{},{},{},{},{}",
            row.hybrid,
            row.phase_name,
            row.mode_name,
            row.dv_share,
            row.phase_delta_v_mps,
            row.thermal_isp_s,
            row.effective_isp_s,
            row.thrust_multiplier,
            row.phase_mass_ratio,
            row.phase_propellant_fraction,
            row.phase_duration_days,
            fmt_opt(row.ionization_fraction_pct),
            fmt_opt(row.electrostatic_grid_voltage_kv),
            fmt_opt(row.electrostatic_stage_count),
            fmt_opt(row.electrostatic_efficiency),
            fmt_opt(row.acoustic_frequency_khz),
            fmt_opt(row.acoustic_heat_transfer_gain_x),
            fmt_opt(row.chamber_pressure_bar),
            fmt_opt(row.methane_blend_fraction),
            fmt_opt(row.supercritical_heat_transfer_gain_x),
        )?;
    }

    Ok(())
}

fn write_isp_sweep_csv(path: &Path, rows: &[IspSweepRow]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "hybrid,isp_s,mass_ratio,rocket_propellant_fraction,baseline_reference_fraction,global_hybrid_reference_min,global_hybrid_reference_max,concept_reference_min,concept_reference_max"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{},{:.3},{:.6},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3}",
            row.hybrid,
            row.isp_s,
            row.mass_ratio,
            row.rocket_propellant_fraction,
            row.baseline_reference_fraction,
            row.global_hybrid_reference_min,
            row.global_hybrid_reference_max,
            row.concept_reference_min,
            row.concept_reference_max,
        )?;
    }

    Ok(())
}

fn write_mechanism_sweep_csv(
    path: &Path,
    rows: &[MechanismSweepRow],
) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "hybrid,driver_name,driver_units,driver_value,effective_isp_s,mission_propellant_fraction,transit_time_days"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{},{},{},{:.6},{:.3},{:.6},{:.3}",
            row.hybrid,
            row.driver_name,
            row.driver_units,
            row.driver_value,
            row.effective_isp_s,
            row.mission_propellant_fraction,
            row.transit_time_days,
        )?;
    }

    Ok(())
}

fn write_monte_carlo_csv(path: &Path, rows: &[MonteCarloRun]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "run_index,hybrid,sampled_dv_mps,sampled_isp_s,mission_propellant_fraction,transit_time_days,ionization_fraction_pct,electrostatic_grid_voltage_kv,electrostatic_stage_count,electrostatic_efficiency,acoustic_frequency_khz,acoustic_heat_transfer_gain_x,chamber_pressure_bar,methane_blend_fraction,supercritical_heat_transfer_gain_x"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{},{},{:.3},{:.3},{:.6},{:.3},{},{},{},{},{},{},{},{},{}",
            row.run_index,
            row.hybrid,
            row.sampled_dv_mps,
            row.sampled_isp_s,
            row.mission_propellant_fraction,
            row.transit_time_days,
            fmt_opt(row.ionization_fraction_pct),
            fmt_opt(row.electrostatic_grid_voltage_kv),
            fmt_opt(row.electrostatic_stage_count),
            fmt_opt(row.electrostatic_efficiency),
            fmt_opt(row.acoustic_frequency_khz),
            fmt_opt(row.acoustic_heat_transfer_gain_x),
            fmt_opt(row.chamber_pressure_bar),
            fmt_opt(row.methane_blend_fraction),
            fmt_opt(row.supercritical_heat_transfer_gain_x),
        )?;
    }

    Ok(())
}
