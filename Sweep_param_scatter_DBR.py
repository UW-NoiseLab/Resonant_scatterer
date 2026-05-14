import os
import time
from pathlib import Path

import numpy as np
import lumapi
from Utils import (
    build_scatterer_output_dir,
    ensure_dir,
    extract_complex_field_from_monitor,
    extract_monitor_wavelengths,
    average_phase_over_xy,
    try_get_transmission,
    plot_summary_figures,
    hash_from_namespace,
    save_json,
    stable_config_hash,
)
from Scatterer_lsf_gen import scatter_DBR_lsf_gen


# =========================================================
# 1) 用户配置
# =========================================================

# Lumerical function 脚本

# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./LC_simulation_DBR.fsp"

# Scatterer sweeps default to output/scatterer_datas/{scheme}_p{period}_r{radius}_{hash}.
OUTPUT_BASE = Path("output_single_scatterer/scatterer_datas_DBR")

# monitor names
FIELD_MONITOR_NAME = "monitor"       # 用来取 E
TRANS_MONITOR_NAME = "monitor"     # 用来取 transmission / reflectance

# DBR high-index material: "TiO2" or "SiN"
DBR_MATERIAL_CHOICES = ("TiO2", "SiN")


# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = True
DBR_DESIGN_WAVELENGTH = 0.620e-6
DBR_MATERIAL = "TiO2"
ONLY_CREATE_MODEL = False

t_Al = 0.30e-6
t_spacer = 0.17e-6
t_LC = 0.50e-6
t_scatter = 0.20e-6


period = 0.36e-6
t_ITO = 0.2e-6
t_glass = 0.1e-6
r_scatter = 0.11e-6


lambda_start = 0.400e-6
lambda_stop = 0.70e-6

# 折射率扫描
index_start = 1.522
index_stop = 1.813
index_step = 0.01

# 波长轴
num_wave = 151

# 最大重试次数
MAX_RETRY = 6


def run_sweep(
    lsf_fsp_file=LSF_FSP_FILE,
    output_root=None,
    output_base=OUTPUT_BASE,
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
    dbr_design_wavelength=DBR_DESIGN_WAVELENGTH,
    dbr_material=DBR_MATERIAL,
):
    """Run the LC index sweep with configurable parameters."""
    lsf_fsp_file = Path(lsf_fsp_file)
    dbr_design_wavelength = float(dbr_design_wavelength)
    if dbr_material not in DBR_MATERIAL_CHOICES:
        valid_materials = ", ".join(DBR_MATERIAL_CHOICES)
        raise ValueError(f"Unknown DBR material '{dbr_material}'. Valid materials: {valid_materials}")

    if index_values is None:
        index_values = np.arange(1.55, 1.75 + 1e-12, 0.01)
    index_values = np.asarray(index_values, dtype=float)

    wave_nm = np.linspace(lambda_start * 1e9, lambda_stop * 1e9, num_wave)  # Convert to nm for output
    n_index = len(index_values)
    n_wave = len(wave_nm)
    dbr_scheme_name = f"DBR {dbr_material} {dbr_design_wavelength * 1e9:.0f}nm"

    namespace_config = {
        "lsf_fsp_file": str(lsf_fsp_file),
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
        "scatterer_scheme": dbr_scheme_name,
        "dbr_design_wavelength": float(dbr_design_wavelength),
        "dbr_material": dbr_material,
        "lambda_start": float(lambda_start),
        "lambda_stop": float(lambda_stop),
        "index_values": index_values.tolist(),
        "num_wave": int(num_wave),
        "wave_nm": wave_nm.tolist(),
        "max_retry": int(max_retry),
    }
    if output_root is None:
        output_root = build_scatterer_output_dir(output_base, namespace_config)
    output_root = Path(output_root)
    config_hash = hash_from_namespace(output_root) or stable_config_hash(namespace_config)

    ensure_dir(output_root)
    config = {
        **namespace_config,
        "output_root": str(output_root),
        "output_namespace": str(output_root),
        "config_hash": config_hash,
    }
    save_json(output_root / "config.json", config)

    fig_dir = output_root / "fig_phase_trans"
    ensure_dir(fig_dir)

    phase_matrix = np.zeros((n_index, n_wave))
    trans_matrix = np.zeros((n_index, n_wave))

    fdtd = lumapi.FDTD(hide=hide_lumerical)

    try:
        if not os.path.isfile(lsf_fsp_file):
            raise FileNotFoundError(f"FSP file NOT found: {lsf_fsp_file}")
        fdtd.load(str(lsf_fsp_file))

        lsf_code = scatter_DBR_lsf_gen(dbr_design_wavelength, dbr_material)

        if len(lsf_code.strip()) == 0:
            raise ValueError(f"LSF code has error: {lsf_code}")

        print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
        fdtd.eval(lsf_code)

        fdtd.switchtolayout()
        fdtd.eval(f"setglobalmonitor('use wavelength spacing',1);")
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

                    # time.sleep(1.0)
                    fdtd.run()
                    time.sleep(1.0)

                    Ex, wavelength_from_monitor = extract_complex_field_from_monitor(
                        fdtd, field_monitor_name
                    )
                    
                    tol = 1e-10  # 可以调

                    if np.allclose(wavelength_from_monitor, wave_nm*1e-9, atol=tol, rtol=0):
                        print("[DEBUG] Wavelengths from monitor match expected wave_nm.")

                    elif np.allclose(wavelength_from_monitor[::-1], wave_nm*1e-9, atol=tol, rtol=0):
                        print("[DEBUG] Wavelengths match in reverse order. Reversing data.")
                        Ex = Ex[..., ::-1]
                        wavelength_from_monitor = wavelength_from_monitor[::-1]

                    else:
                        raise ValueError("[WARNING] Wavelength mismatch!")
                    
                    
                    phase, _, _ = average_phase_over_xy(Ex)

                    trans = try_get_transmission(fdtd, trans_monitor_name)
                    trans_wavelength = extract_monitor_wavelengths(fdtd, trans_monitor_name)

                    if np.allclose(trans_wavelength, wave_nm * 1e-9, atol=tol, rtol=0):
                        print("[DEBUG] Reflection wavelengths match expected wave_nm.")
                    elif np.allclose(trans_wavelength[::-1], wave_nm * 1e-9, atol=tol, rtol=0):
                        print("[DEBUG] Reflection wavelengths match in reverse order. Reversing data.")
                        trans = trans[::-1]
                        trans_wavelength = trans_wavelength[::-1]
                    else:
                        raise ValueError("[WARNING] Reflection wavelength mismatch!")

                    if len(phase) != n_wave:
                        raise RuntimeError(
                            f"Phase length {len(phase)} != expected {n_wave}"
                        )
                    if len(trans) != n_wave:
                        raise RuntimeError(
                            f"Transmission length {len(trans)} != expected {n_wave}"
                        )
                    if not np.allclose(wavelength_from_monitor, trans_wavelength, atol=tol, rtol=0):
                        raise RuntimeError(
                            "Field monitor wavelengths do not match reflection wavelengths after ordering."
                        )

                    phase_matrix[i, :] = phase
                    trans_matrix[i, :] = trans

                    success = True
                    break
                except Exception as e:
                    print(f"    Attempt {retry+1} failed: {e}")
                    print(lsf_code)
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
    parser.add_argument(
        "--output-root",
        default=None,
        help="Explicit output directory. If omitted, use output/scatterer_datas/{scheme}_p{period}_r{radius}_{hash}.",
    )
    parser.add_argument("--output-base", default=str(OUTPUT_BASE), help="Base directory for automatic scatterer namespaces")
    parser.add_argument("--num-wave", type=int, default=num_wave, help="Number of wavelength points")
    parser.add_argument("--max-retry", type=int, default=MAX_RETRY, help="Retry attempts per index")
    parser.add_argument("--index-start", type=float, default=index_start, help="Start LC index")
    parser.add_argument("--index-stop", type=float, default=index_stop, help="Stop LC index")
    parser.add_argument("--index-step", type=float, default=index_step, help="Step size for LC index")
    parser.add_argument("--dbr-wavelength", type=float, default=DBR_DESIGN_WAVELENGTH, help="DBR design wavelength in meters")
    parser.add_argument("--dbr-material", choices=DBR_MATERIAL_CHOICES, default=DBR_MATERIAL, help="High-index DBR material")

    args = parser.parse_args()
    
    run_sweep(
        output_root=Path(args.output_root) if args.output_root else None,
        output_base=Path(args.output_base),
        hide_lumerical=HIDE_LUMERICAL,
        num_wave=args.num_wave,
        max_retry=args.max_retry,
        index_values=np.arange(args.index_start, args.index_stop + 1e-12, args.index_step),
        only_create_model=ONLY_CREATE_MODEL,
        dbr_design_wavelength=args.dbr_wavelength,
        dbr_material=args.dbr_material,
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
    )

def overwrite_default_params(params: dict):
    """Overwrite default parameters with provided ones."""
    global DBR_DESIGN_WAVELENGTH, DBR_MATERIAL, period, t_Al, t_spacer, t_LC, t_ITO, t_glass, r_scatter, t_scatter, lambda_start, lambda_stop, index_values, num_wave
    DBR_DESIGN_WAVELENGTH = params.get("dbr_design_wavelength", DBR_DESIGN_WAVELENGTH)
    DBR_MATERIAL = params.get("dbr_material", DBR_MATERIAL)
    period = params.get("period", period)
    t_Al = params.get("t_Al", t_Al)
    t_spacer = params.get("t_spacer", t_spacer)
    t_LC = params.get("t_LC", t_LC)
    t_ITO = params.get("t_ITO", t_ITO)
    t_glass = params.get("t_glass", t_glass)
    r_scatter = params.get("r_scatter", r_scatter)
    t_scatter = params.get("t_scatter", t_scatter)
    lambda_start = params.get("lambda_start", lambda_start)
    lambda_stop = params.get("lambda_stop", lambda_stop)
    
    
if __name__ == "__main__":
    from archived_scatterer_params import SiN_b_565nm, SiN_t_490nm, ACSNano_TiO2_665nm, \
    TiO2_t_Zhihao, SiN_b_Zhihao, Meta_TiO2_opt_a, Meta_TiO2_opt_b
    overwrite_default_params(Meta_TiO2_opt_a)
    main()
