import argparse
import os
import time
from pathlib import Path

import numpy as np
from Utils import (
    add_angle_accuracy,
    build_grating_output_dir,
    build_scatterer_output_dir,
    ensure_dir,
    find_farfield_peak,
    has_scatterer_results,
    load_scatterer_results,
    plot_candidate_summary,
    save_candidate_report,
    save_farfield_data_single_wavelength,
    save_json,
    save_lsf_file,
    score_scatterer_wavelengths,
    wavelength_tag,
)


DEFAULT_SCATTERER_FSP = "./LC_simulation.fsp"
DEFAULT_GRATING_FSP = "./Grating_simulation.fsp"
SCATTERER_SCHEME_CHOICES = ("SiN on Top", "SiN on Bottom", "TiO2 on Top", "TiO2 on Bottom")

OUTPUT_ROOT = r"./output"

DEFAULT_PERIOD_LIST = [0.34e-6, 0.36e-6, 0.38e-6]


def period_radius_pairs():
    """Return the default period sweep with radius scaled from r=0.12 um at p=0.36 um."""
    period_list = DEFAULT_PERIOD_LIST
    r_scatter_list = [0.12 / 0.36 * p for p in period_list]
    return list(zip(period_list, r_scatter_list))


def regime_schemes(selected_scheme):
    """Return one requested scheme or all four supported regimes by default."""
    if selected_scheme is not None:
        return (selected_scheme,)
    return SCATTERER_SCHEME_CHOICES


def scatterer_namespace_config(args):
    return {
        "lsf_fsp_file": str(args.scatterer_fsp),
        "field_monitor_name": "monitor",
        "trans_monitor_name": "monitor",
        "hide_lumerical": bool(args.hide_lumerical),
        "period": float(args.period),
        "t_Al": float(args.t_al),
        "t_spacer": float(args.t_spacer),
        "t_LC": float(args.t_lc),
        "t_ITO": float(args.t_ito),
        "t_glass": float(args.t_glass),
        "r_scatter": float(args.r_scatter),
        "t_scatter": float(args.t_scatter),
        "scatterer_scheme": args.scatterer_scheme,
        "lambda_start": float(args.lambda_start_nm * 1e-9),
        "lambda_stop": float(args.lambda_stop_nm * 1e-9),
        "index_values": np.arange(args.index_start, args.index_stop + 1e-12, args.index_step).tolist(),
        "num_wave": int(args.num_wave),
        "wave_nm": np.linspace(args.lambda_start_nm, args.lambda_stop_nm, args.num_wave).tolist(),
        "max_retry": int(args.max_retry),
    }


def get_scatterer_phase_profile(args):
    import Sweep_param_scatter

    if args.scatterer_output:
        data_folder = Path(args.scatterer_output)
    else:
        data_folder = build_scatterer_output_dir(
            args.scatterer_output_base,
            scatterer_namespace_config(args),
        )
    args.scatterer_output = str(data_folder)

    if args.force_scatterer_sweep or not has_scatterer_results(data_folder):
        print("Running scatterer sweep to generate phase/transmission profile.")
        Sweep_param_scatter.run_sweep(
            lsf_fsp_file=args.scatterer_fsp,
            output_root=data_folder,
            output_base=args.scatterer_output_base,
            hide_lumerical=args.hide_lumerical,
            period=args.period,
            t_Al=args.t_al,
            t_spacer=args.t_spacer,
            t_LC=args.t_lc,
            t_ITO=args.t_ito,
            t_glass=args.t_glass,
            r_scatter=args.r_scatter,
            t_scatter=args.t_scatter,
            lambda_start=args.lambda_start_nm * 1e-9,
            lambda_stop=args.lambda_stop_nm * 1e-9,
            index_values=np.arange(args.index_start, args.index_stop + 1e-12, args.index_step),
            num_wave=args.num_wave,
            max_retry=args.max_retry,
            only_create_model=False,
            scatterer_scheme=args.scatterer_scheme,
        )
    else:
        print(f"Using existing scatterer results: {data_folder}")

    results = load_scatterer_results(data_folder)
    config_scheme = results["config"].get("scatterer_scheme")
    if config_scheme != args.scatterer_scheme:
        raise ValueError(
            "Scatterer result scheme mismatch: "
            f"folder has '{config_scheme}', but this run requested '{args.scatterer_scheme}'."
        )

    return results


def design_and_simulate_grating(candidate, args, output_dir):
    import lumapi
    import Grating_fdtd

    wavelength_nm = candidate["wavelength_nm"]
    target_lambda_m = wavelength_nm * 1e-9
    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    Grating_fdtd.SCATTERER_SCHEME = args.scatterer_scheme
    grating_data = Grating_fdtd.load_scatterer_config_and_data(
        Path(args.scatterer_output),
        target_wavelength_nm=wavelength_nm,
    )

    wl_tag = wavelength_tag(wavelength_nm, prefix="", decimals=1)
    lsf_path = output_dir / f"grating_{wl_tag}.lsf"
    save_lsf_file(grating_data["lsf_content"], lsf_path)

    fdtd = lumapi.FDTD(hide=args.hide_lumerical)
    try:
        if os.path.exists(args.grating_fsp):
            fdtd.load(args.grating_fsp)

        fdtd.eval(grating_data["lsf_content"])
        fdtd.eval(f"createmodel({args.grating_cells},{args.steering_angle_deg});")

        try:
            fdtd.getnamed(args.farfield_monitor, "name")
        except Exception as exc:
            debug_fsp = output_dir / "debug_monitor_missing_before_run.fsp"
            fdtd.save(str(debug_fsp))
            raise RuntimeError(
                f"Monitor object '{args.farfield_monitor}' was not created by createmodel(). "
                f"Saved debug file: {debug_fsp}"
            ) from exc

        # fdtd.save(str(output_dir / "debug_before_run.fsp"))
        fdtd.eval("run;")
        time.sleep(1.0)
        # fdtd.save(str(output_dir / "debug_after_run.fsp"))

        farfield_dir = output_dir / "farfield"
        save_farfield_data_single_wavelength(
            fdtd=fdtd,
            monitor_name=args.farfield_monitor,
            target_lambda_m=target_lambda_m,
            output_dir=farfield_dir,
            prefix=f"grating_{wl_tag}",
            make_plot=True,
            save_angle_axes=True,
        )

        peak = find_farfield_peak(fdtd, args.farfield_monitor, target_lambda_m)
        peak = add_angle_accuracy(peak, args.steering_angle_deg)
        peak["designed_wavelength_nm"] = float(wavelength_nm)
        peak["candidate_score"] = candidate["score"]
        peak["candidate_mean_transmission"] = candidate["mean_transmission"]
        peak["candidate_coverage_fraction_lut"] = candidate["coverage_fraction_lut"]
        peak["candidate_uniformity_score_lut"] = candidate["uniformity_score_lut"]
        peak["grating_cells"] = int(args.grating_cells)
        peak["scatterer_scheme"] = args.scatterer_scheme

        save_json(output_dir / "farfield_peak_result.json", peak)

        return peak
    finally:
        try:
            fdtd.close()
        except Exception:
            pass


def candidates_for_requested_wavelengths(requested_wavelengths_nm, scored_candidates):
    """Create simulation candidates from user-requested wavelengths."""
    scored_wavelengths = np.array([candidate["wavelength_nm"] for candidate in scored_candidates])

    requested_candidates = []
    for requested_wavelength_nm in requested_wavelengths_nm:
        nearest_idx = int(np.argmin(np.abs(scored_wavelengths - requested_wavelength_nm)))
        nearest = dict(scored_candidates[nearest_idx])
        nearest["requested_wavelength_nm"] = float(requested_wavelength_nm)
        nearest["wavelength_nm"] = float(scored_wavelengths[nearest_idx])
        nearest["selection_mode"] = "user_requested_nearest_available"
        nearest["request_error_nm"] = float(nearest["wavelength_nm"] - requested_wavelength_nm)
        requested_candidates.append(nearest)

    return requested_candidates


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run/load a scatterer phase sweep, identify wavelengths with good 0-2pi "
            "phase coverage and transmission, then design grating simulations and "
            "check far-field steering peaks."
        )
    )

    parser.add_argument("--scatterer-output", default=None, help="Explicit scatterer data folder. If omitted, use the automatic scatterer namespace.")
    parser.add_argument("--scatterer-output-base", default=os.path.join(OUTPUT_ROOT, "scatterer_datas"))
    parser.add_argument("--grating-output-base", default=os.path.join(OUTPUT_ROOT, "grating_datas"))
    parser.add_argument("--pipeline-output", default=os.path.join(OUTPUT_ROOT, "pipeline_reports"))
    parser.add_argument("--scatterer-fsp", default=DEFAULT_SCATTERER_FSP)
    parser.add_argument("--grating-fsp", default=DEFAULT_GRATING_FSP)
    parser.add_argument(
        "--scatterer-scheme",
        default=None,
        choices=SCATTERER_SCHEME_CHOICES,
        help="Run only one scatterer regime. If omitted, sweep all four regimes.",
    )
    parser.add_argument("--force-scatterer-sweep", action="store_true")
    parser.add_argument(
        "--show-lumerical",
        dest="hide_lumerical",
        action="store_false",
        help="Show the Lumerical GUI. By default, simulations run hidden.",
    )
    parser.set_defaults(hide_lumerical=False)

    parser.add_argument("--period", type=float, default=0.36e-6)
    parser.add_argument("--t-al", type=float, default=0.30e-6)
    parser.add_argument("--t-spacer", type=float, default=0.17e-6)
    parser.add_argument("--t-lc", type=float, default=0.50e-6)
    parser.add_argument("--t-ito", type=float, default=0.05e-6)
    parser.add_argument("--t-glass", type=float, default=3.0e-6)
    parser.add_argument("--r-scatter", type=float, default=0.12e-6)
    parser.add_argument("--t-scatter", type=float, default=0.20e-6)

    parser.add_argument("--lambda-start-nm", type=float, default=400.0)
    parser.add_argument("--lambda-stop-nm", type=float, default=700.0)
    parser.add_argument("--num-wave", type=int, default=100)
    parser.add_argument("--index-start", type=float, default=1.55)
    parser.add_argument("--index-stop", type=float, default=1.75)
    parser.add_argument("--index-step", type=float, default=0.01)
    parser.add_argument("--max-retry", type=int, default=6)

    parser.add_argument("--min-mean-transmission", type=float, default=0.5)
    parser.add_argument("--candidate-count", type=int, default=5)
    parser.add_argument("--simulate-top-n", type=int, default=1)
    parser.add_argument(
        "--test-wavelength-nm",
        type=float,
        nargs="+",
        default=None,
        help=(
            "One or more wavelengths to simulate directly. Values are snapped to "
            "the nearest wavelength available in the scatterer sweep. If omitted, "
            "the top scored wavelengths are simulated."
        ),
    )
    parser.add_argument("--grating-cells", type=int, default=10)
    parser.add_argument("--steering-angle-deg", type=float, default=8.0)
    parser.add_argument("--farfield-monitor", default="R_monitor")

    return parser


def run_single_regime_period_radius(args, scheme, period, r_scatter):
    args.scatterer_scheme = scheme
    args.period = period
    args.r_scatter = r_scatter
    args.scatterer_output = None

    output_root = Path(args.pipeline_output)
    ensure_dir(output_root)

    print("\n" + "=" * 72)
    print(
        f"Running pipeline for scheme='{scheme}', "
        f"period={period:.4e} m, r_scatter={r_scatter:.4e} m"
    )

    scatterer_results = get_scatterer_phase_profile(args)

    run_report_dir = output_root / Path(args.scatterer_output).name
    ensure_dir(run_report_dir)
    save_json(run_report_dir / "pipeline_config.json", vars(args))

    candidates = score_scatterer_wavelengths(
        scatterer_results,
        min_mean_transmission=args.min_mean_transmission,
    )
    save_candidate_report(candidates, run_report_dir)
    plot_candidate_summary(candidates, run_report_dir)

    print("\nTop wavelength candidates:")
    for rank, candidate in enumerate(candidates[: args.candidate_count], start=1):
        print(
            f"{rank:2d}. {candidate['wavelength_nm']:.2f} nm | "
            f"score={candidate['score']:.3f}, "
            f"coverage={candidate['coverage_fraction_lut']:.3f}, "
            f"uniformity={candidate['uniformity_score_lut']:.3f}, "
            f"mean T={candidate['mean_transmission']:.3f}, "
            f"LUT points={candidate['usable_lut_points']}/{candidate['total_index_points']}"
        )

    if args.test_wavelength_nm:
        simulation_candidates = candidates_for_requested_wavelengths(args.test_wavelength_nm, candidates)
        print("\nUser-requested wavelength candidates:")
        for candidate in simulation_candidates:
            print(
                f"requested={candidate['requested_wavelength_nm']:.2f} nm, "
                f"nearest={candidate['wavelength_nm']:.2f} nm, "
                f"error={candidate['request_error_nm']:.3f} nm, "
                f"score={candidate['score']:.3f}"
            )
    else:
        simulation_candidates = candidates[: args.simulate_top_n]

    peak_results = []
    for candidate in simulation_candidates:
        grating_config = {
            **scatterer_results["config"],
            "wavelength_nm": float(candidate["wavelength_nm"]),
            "grating_cells": int(args.grating_cells),
            "steering_angle_deg": float(args.steering_angle_deg),
            "farfield_monitor": args.farfield_monitor,
        }
        sim_dir = build_grating_output_dir(args.grating_output_base, grating_config)
        print(f"\nDesigning and simulating grating for {candidate['wavelength_nm']:.2f} nm")
        peak_results.append(design_and_simulate_grating(candidate, args, sim_dir))

    save_json(run_report_dir / "grating_farfield_peak_results.json", peak_results)

    if peak_results:
        print("\nFar-field peak checks:")
        for result in peak_results:
            print(
                f"{result['designed_wavelength_nm']:.2f} nm | "
                f"peak theta_x={result['peak_theta_x_deg']:.2f} deg, "
                f"theta_y={result['peak_theta_y_deg']:.2f} deg, "
                f"best error={result['best_abs_error_deg']:.2f} deg "
                f"({result['best_sign_match']} sign)"
            )

    print(f"\nPipeline report saved to: {run_report_dir.resolve()}")
    return {
        "period": float(period),
        "r_scatter": float(r_scatter),
        "scatterer_scheme": scheme,
        "scatterer_output": args.scatterer_output,
        "report_dir": str(run_report_dir),
        "top_candidates": candidates[: args.candidate_count],
        "simulation_candidates": simulation_candidates,
        "peak_results": peak_results,
    }


def main():
    args = build_arg_parser().parse_args()
    output_root = ensure_dir(args.pipeline_output)

    sweep_results = []
    for scheme in regime_schemes(args.scatterer_scheme):
        for period, r_scatter in period_radius_pairs():
            sweep_results.append(run_single_regime_period_radius(args, scheme, period, r_scatter))

    save_json(output_root / "sweep_summary.json", sweep_results)
    print(f"\nAll pipeline sweep reports saved under: {output_root.resolve()}")


if __name__ == "__main__":
    main()
