import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import lumapi
from Utils import (
    ensure_dir,
    extract_complex_field_from_monitor,
    average_phase_over_xy,
    try_get_transmission,
    plot_summary_figures,
)
from Scatterer_lsf_gen import scatterer_lsf_gen


# =========================================================
# 1) 用户配置
# =========================================================

# Lumerical function 脚本

# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./LC_simulation.fsp"

# 输出根目录
OUTPUT_ROOT = Path("automated_sweep_results")

# monitor names
FIELD_MONITOR_NAME = "monitor"       # 用来取 E
TRANS_MONITOR_NAME = "monitor"     # 用来取 transmission / reflectance

# Scatterer scheme: "SiN on Top" or "SiN on Bottom" or "TiO2 on Top" or "TiO2 on Bottom"
SCATTERER_SCHEME = "TiO2 on Top"
SCATTERER_SCHEME_CHOICES = ("SiN on Top", "SiN on Bottom", "TiO2 on Top", "TiO2 on Bottom")


# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = True

# =========================================================
# 2) 建模参数
# =========================================================

period = 0.36e-6
t_Al = 0.30e-6
t_spacer = 0.17e-6
t_LC = 0.50e-6
t_ITO = 0.05e-6
t_glass = 3.0e-6
r_scatter = 0.12e-6
t_scatter = 0.20e-6
lambda_start = 0.40e-6
lambda_stop = 0.70e-6

# 折射率扫描
index_values = np.arange(1.55, 1.75 + 1e-12, 0.01)

# 波长轴
num_wave = 100
wave_nm = np.linspace(400, 700, num_wave)

# 最大重试次数
MAX_RETRY = 6


def run_sweep(
    lsf_fsp_file=LSF_FSP_FILE,
    output_root=OUTPUT_ROOT,
    field_monitor_name=FIELD_MONITOR_NAME,
    trans_monitor_name=TRANS_MONITOR_NAME,
    hide_lumerical=HIDE_LUMERICAL,
    period=period,
    t_Al=t_Al,
    t_spacer=t_spacer,
    t_LC=t_LC,
    t_ITO=t_ITO,
    t_glass=t_glass,
    r_scatter=r_scatter,
    t_scatter=t_scatter,
    lambda_start=lambda_start,
    lambda_stop=lambda_stop,
    index_values=None,
    num_wave=num_wave,
    max_retry=MAX_RETRY,
    only_create_model=False,
    scatterer_scheme=SCATTERER_SCHEME,
):
    """Run the LC index sweep with configurable parameters."""
    output_root = Path(output_root)
    lsf_fsp_file = Path(lsf_fsp_file)

    if index_values is None:
        index_values = np.arange(1.55, 1.75 + 1e-12, 0.01)
    index_values = np.asarray(index_values, dtype=float)

    wave_nm = np.linspace(400, 700, num_wave)
    n_index = len(index_values)
    n_wave = len(wave_nm)

    ensure_dir(output_root)
    config = {
        "lsf_fsp_file": str(lsf_fsp_file),
        "output_root": str(output_root),
        "field_monitor_name": field_monitor_name,
        "trans_monitor_name": trans_monitor_name,
        "hide_lumerical": bool(hide_lumerical),
        "period": float(period),
        "t_Al": float(t_Al),
        "t_spacer": float(t_spacer),
        "t_LC": float(t_LC),
        "t_ITO": float(t_ITO),
        "t_glass": float(t_glass),
        "r_scatter": float(r_scatter),
        "t_scatter": float(t_scatter),
        "scatterer_scheme": scatterer_scheme,
        "lambda_start": float(lambda_start),
        "lambda_stop": float(lambda_stop),
        "index_values": index_values.tolist(),
        "num_wave": int(num_wave),
        "wave_nm": wave_nm.tolist(),
        "max_retry": int(max_retry),
    }
    with open(output_root / "config.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    fig_dir = output_root / "fig_phase_trans"
    ensure_dir(fig_dir)

    phase_matrix = np.zeros((n_index, n_wave))
    trans_matrix = np.zeros((n_index, n_wave))

    fdtd = lumapi.FDTD(hide=hide_lumerical)

    try:
        if not os.path.isfile(lsf_fsp_file):
            raise FileNotFoundError(f"FSP file NOT found: {lsf_fsp_file}")
        fdtd.load(str(lsf_fsp_file))

        lsf_code = scatterer_lsf_gen(scatterer_scheme)

        if len(lsf_code.strip()) == 0:
            raise ValueError(f"LSF code has error: {lsf_code}")

        print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
        fdtd.eval(lsf_code)

        fdtd.switchtolayout()
        fdtd.eval(f'setglobalmonitor("frequency points",{num_wave});')
        print(f"[DEBUG] Set global monitor frequency points to {num_wave}")
        

        for i, lc_index in enumerate(index_values):
            print(f"[{i+1}/{n_index}] Running LC index = {lc_index:.3f}")
            success = False

            for retry in range(max_retry):
                print(f"  Attempt {retry+1}/{max_retry}")
                try:
                    fdtd.eval("switchtolayout;")
                    fdtd.eval("deleteall;")

                    build_cmd = (
                        f"build_unit_cell_model("
                        f"{period},"
                        f"{t_Al},"
                        f"{t_spacer},"
                        f"{t_LC},"
                        f"{t_ITO},"
                        f"{t_glass},"
                        f"{r_scatter},"
                        f"{t_scatter},"
                        f"{lambda_start},"
                        f"{lambda_stop}"
                        f");"
                    )
                    fdtd.eval(build_cmd)

                    fdtd.select("liquid_crystal")
                    fdtd.set("index", float(lc_index))
                    
                    if only_create_model:
                        print("Only creating model, skipping simulation.")
                        time.sleep(2000.0)
                        return

                    time.sleep(1.0)
                    fdtd.run()
                    time.sleep(1.0)

                    Ex, wavelength_from_monitor = extract_complex_field_from_monitor(
                        fdtd, field_monitor_name
                    )
                    phase, _, _ = average_phase_over_xy(Ex)

                    trans = try_get_transmission(fdtd, trans_monitor_name)

                    if len(phase) != n_wave:
                        raise RuntimeError(
                            f"Phase length {len(phase)} != expected {n_wave}"
                        )
                    if len(trans) != n_wave:
                        raise RuntimeError(
                            f"Transmission length {len(trans)} != expected {n_wave}"
                        )

                    phase_matrix[i, :] = phase
                    trans_matrix[i, :] = trans

                    success = True
                    break
                except Exception as e:
                    print(f"    Attempt {retry+1} failed: {e}")
                    time.sleep(1.0)

            if not success:
                raise RuntimeError(
                    f"Failed to complete simulation for LC index {lc_index:.3f} after {max_retry} attempts."
                )

        np.save(output_root / "index_values.npy", index_values)
        np.save(output_root / "wave_nm.npy", wave_nm)
        np.save(output_root / "phase_matrix.npy", phase_matrix)
        np.save(output_root / "trans_matrix.npy", trans_matrix)

        plot_summary_figures(
            outdir=fig_dir,
            phase_matrix_raw=phase_matrix,
            trans_matrix=trans_matrix,
            index_values=index_values,
            wave_nm=wave_nm,
        )

        print("Done.")
        print(f"Results saved to: {output_root.resolve()}")

        return {
            "index_values": index_values,
            "wave_nm": wave_nm,
            "phase_matrix": phase_matrix,
            "trans_matrix": trans_matrix,
        }
    finally:
        try:
            fdtd.close()
        except Exception:
            pass


# =========================================================
# 4) 主流程
# =========================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Run a Lumerical LC index sweep.")
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT), help="Output directory")
    parser.add_argument("--num-wave", type=int, default=num_wave, help="Number of wavelength points")
    parser.add_argument("--max-retry", type=int, default=MAX_RETRY, help="Retry attempts per index")
    parser.add_argument("--index-start", type=float, default=1.55, help="Start LC index")
    parser.add_argument("--index-stop", type=float, default=1.75, help="Stop LC index")
    parser.add_argument("--index-step", type=float, default=0.01, help="Step size for LC index")
    parser.add_argument("--hide-gui", action="store_false", help="Hide Lumerical GUI")
    parser.add_argument(
        "--scatterer-scheme",
        default=SCATTERER_SCHEME,
        choices=SCATTERER_SCHEME_CHOICES,
        help="Scatterer material and z-placement scheme",
    )
    args = parser.parse_args()
    
    # r_scatter_list = [0.08e-6, 0.10e-6, 0.12e-6, 0.14e-6, 0.16e-6]
    # output_root_list = [f"{OUTPUT_ROOT}_r{int(r*1e9)}nm" for r in r_scatter_list]

    # for i, output_root in enumerate(output_root_list):
    #     run_sweep(
    #         r_scatter=r_scatter_list[i],
    #         output_root=output_root,
    #         hide_lumerical=args.hide_gui,
    #         num_wave=args.num_wave,
    #         max_retry=args.max_retry,
    #         index_values=np.arange(args.index_start, args.index_stop + 1e-12, args.index_step),
    #     )
    
    # period_list = [0.34e-6, 0.36e-6, 0.38e-6]
    # output_root_list = [f"{OUTPUT_ROOT}_p{int(p*1e9)}nm" for p in period_list]
    
    # for i, output_root in enumerate(output_root_list):
    #     run_sweep(
    #         period=period_list[i],
    #         output_root=output_root,
    #         hide_lumerical=args.hide_gui,
    #         num_wave=args.num_wave,
    #         max_retry=args.max_retry,
    #         index_values=np.arange(args.index_start, args.index_stop + 1e-12, args.index_step),
    #     )
    
    # period_list = [0.25e-6, 0.26e-6, 0.27e-6, 0.28e-6, 0.29e-6, 0.30e-6, 0.31e-6, 0.32e-6, 0.33e-6, 0.34e-6, 0.35e-6, 0.36e-6, 0.37e-6, 0.38e-6]
    # period_list = [0.322e-6, 0.324e-6, 0.326e-6, 0.328e-6, 0.33e-6]
    period_list = [0.34e-6]
    r_scatter_list = [0.12/0.36*p for p in period_list]
    output_root_list = [f"./automated_result_glass_0.3/{OUTPUT_ROOT}_r{int(r*1e9)}nm_p{int(p*1e9)}nm" for r, p in zip(r_scatter_list, period_list)]

    
    for i, output_root in enumerate(output_root_list):
        print(f"Would run sweep with period={period_list[i]:.3e} m and r_scatter={r_scatter_list[i]:.3e} m, output_root='{output_root}'")

        run_sweep(
            period=period_list[i],
            r_scatter=r_scatter_list[i],
            output_root=output_root,
            hide_lumerical=False,
            num_wave=args.num_wave,
            max_retry=args.max_retry,
            index_values=np.arange(args.index_start, args.index_stop + 1e-12, args.index_step),
            only_create_model=True,
            scatterer_scheme=args.scatterer_scheme,
        )
        
if __name__ == "__main__":
    main()
