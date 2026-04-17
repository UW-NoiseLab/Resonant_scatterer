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
from Grating_lsf_gen import grating_lsf_gen

# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./Grating_simulation.fsp"

SCATTERER_SCHEME = "TiO2 on Top"

# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = False

grating_wavelength_nm = 533.0

scatterer_data_folder = Path("./automated_result_glass_3/automated_sweep_results_r107nm_p322nm")

# 输出根目录
OUTPUT_ROOT = Path(os.path.join("./automated_grating_results", f"scatterer_wavelength_{int(grating_wavelength_nm)}nm"))

def load_scatterer_config_and_data(data_folder, target_wavelength_nm=533.0):
    """Load configuration and lookup tables from scatterer results at target wavelength."""

    config_file = data_folder / "config.json"
    with open(config_file, "r") as f:
        config = json.load(f)

    config_scheme = config.get("scatterer_scheme")
    if config_scheme is None:
        raise ValueError(
            f"Scatterer config {config_file} does not contain 'scatterer_scheme'. "
            "Regenerate the scatterer LUT with the updated Sweep_param_scatter.py "
            "so Grating_fdtd.py can verify the geometry scheme."
        )
    if config_scheme != SCATTERER_SCHEME:
        raise ValueError(
            "Scatterer scheme mismatch: "
            f"LUT folder was generated with '{config_scheme}', "
            f"but Grating_fdtd.py is configured for '{SCATTERER_SCHEME}'."
        )
    print(f"Verified scatterer scheme: {config_scheme}")

    # Load numpy arrays
    wave_nm = np.load(data_folder / "wave_nm.npy")
    index_values = np.load(data_folder / "index_values.npy")
    phase_matrix = np.load(data_folder / "phase_matrix.npy")
    trans_matrix = np.load(data_folder / "trans_matrix.npy")

    # Find the index closest to target wavelength
    wl_idx = np.argmin(np.abs(wave_nm - target_wavelength_nm))
    actual_wl_nm = wave_nm[wl_idx]

    print(f"Target wavelength: {target_wavelength_nm} nm")
    print(f"Actual wavelength used: {actual_wl_nm} nm")

    # Extract lookup tables at target wavelength
    # phase_matrix shape: (num_index, num_wavelength)
    # trans_matrix shape: (num_index, num_wavelength)
    lc_phase_data = phase_matrix[:, wl_idx]
    lc_phase_data = (lc_phase_data - lc_phase_data[0]) % (2 * np.pi)  # Normalize to [0, 2pi)

    lc_trans_data = trans_matrix[:, wl_idx]
    lc_index_data = index_values.copy()

    # Keep only strictly increasing phase points
    keep_mask = np.zeros(len(lc_phase_data), dtype=bool)
    keep_mask[0] = True

    current_max = lc_phase_data[0]
    for i in range(1, len(lc_phase_data)):
        if lc_phase_data[i] > current_max:
            keep_mask[i] = True
            current_max = lc_phase_data[i]

    num_removed = len(lc_phase_data) - np.sum(keep_mask)

    lc_phase_data = lc_phase_data[keep_mask]
    lc_trans_data = lc_trans_data[keep_mask]
    lc_index_data = lc_index_data[keep_mask]

    print(f"Removed {num_removed} non-monotonic LUT points.")
    print(f"Remaining LUT points: {len(lc_phase_data)}")

    # Convert to Python lists for LSF generation
    lc_phase_data = lc_phase_data.tolist()
    lc_trans_data = lc_trans_data.tolist()
    lc_index_data = lc_index_data.tolist()
    
    # Generate LSF content using the creategratinglsf function
    lsf_content = grating_lsf_gen(
        period=config["period"],
        t_Al=config["t_Al"],
        t_spacer=config["t_spacer"],
        t_LC=config["t_LC"],
        t_ITO=config["t_ITO"],
        t_glass=config["t_glass"],
        r_scatter=config["r_scatter"],
        t_scatter=config["t_scatter"],
        lambda_design=target_wavelength_nm * 1e-9,
        lambda_span=100e-9,  # 100 nm span for the LSF function
        LC_index_data=lc_index_data,
        LC_phase_data=lc_phase_data,
        scheme=SCATTERER_SCHEME
    )
    
    return {
        "config": config,
        "wave_nm": wave_nm,
        "index_values": index_values,
        "phase_matrix": phase_matrix,
        "trans_matrix": trans_matrix,
        "target_wavelength_nm": actual_wl_nm,
        "lc_index_data": lc_index_data,
        "lc_phase_data": lc_phase_data,
        "lsf_content": lsf_content,
    }


def save_lsf_file(lsf_content, output_file):
    """Save generated LSF content to file."""
    ensure_dir(Path(output_file).parent)
    with open(output_file, 'w') as f:
        f.write(lsf_content)
    print(f"LSF file saved to: {output_file}")

def save_farfield_data_single_wavelength(
    fdtd,
    monitor_name,
    target_lambda_m,
    output_dir,
    prefix="farfield",
    make_plot=True,
    save_angle_axes=False,
):
    """
    Save far-field data from a monitor at the wavelength closest to target_lambda_m.

    Save only E2 + metadata by default, to avoid too many duplicated files.
    Plot axes are converted to angle (deg).
    """
    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    c0 = 299792458.0

    # frequency vector
    f_vec = np.array(fdtd.getdata(monitor_name, "f")).squeeze()
    if f_vec.size == 0:
        raise ValueError(f"No frequency data found in monitor: {monitor_name}")

    # convert to wavelength
    lambda_vec = c0 / f_vec

    # find closest wavelength index
    f_index_py = int(np.argmin(np.abs(lambda_vec - target_lambda_m)))
    actual_lambda_m = float(lambda_vec[f_index_py])
    f_index_lsf = f_index_py + 1  # Lumerical is 1-based

    print(f"Target wavelength: {target_lambda_m*1e9:.3f} nm")
    print(f"Actual wavelength used: {actual_lambda_m*1e9:.3f} nm")
    print(f"Frequency index used in Lumerical: {f_index_lsf}")

    # far field
    E2 = np.array(fdtd.farfield3d(monitor_name, f_index_lsf))
    ux = np.array(fdtd.farfieldux(monitor_name, f_index_lsf)).squeeze()
    uy = np.array(fdtd.farfielduy(monitor_name, f_index_lsf)).squeeze()
    E2 = np.squeeze(E2)

    # convert to angle axis (degree)
    theta_x_deg = np.degrees(np.arcsin(np.clip(ux, -1.0, 1.0)))
    theta_y_deg = np.degrees(np.arcsin(np.clip(uy, -1.0, 1.0)))

    wl_tag = f"{actual_lambda_m * 1e9:.1f}nm".replace(".", "p")

    # save only E2 by default
    np.save(output_dir / f"{prefix}_E2_{wl_tag}.npy", E2)

    # optionally save angle axes once if needed
    if save_angle_axes:
        np.save(output_dir / f"{prefix}_theta_x_deg.npy", theta_x_deg)
        np.save(output_dir / f"{prefix}_theta_y_deg.npy", theta_y_deg)

    meta = {
        "monitor_name": monitor_name,
        "target_lambda_m": float(target_lambda_m),
        "actual_lambda_m": actual_lambda_m,
        "target_lambda_nm": float(target_lambda_m * 1e9),
        "actual_lambda_nm": float(actual_lambda_m * 1e9),
        "f_index_python": int(f_index_py),
        "f_index_lumerical": int(f_index_lsf),
        "E2_shape": list(E2.shape),
        "theta_x_range_deg": [float(theta_x_deg.min()), float(theta_x_deg.max())],
        "theta_y_range_deg": [float(theta_y_deg.min()), float(theta_y_deg.max())],
        "theta_x_size": int(theta_x_deg.size),
        "theta_y_size": int(theta_y_deg.size),
    }

    with open(output_dir / f"{prefix}_meta_{wl_tag}.json", "w") as f:
        json.dump(meta, f, indent=2)

    if make_plot:
        plt.figure(figsize=(6, 5))

        # 用 pcolormesh 比 imshow 更适合非线性角度轴
        THETA_X, THETA_Y = np.meshgrid(theta_x_deg, theta_y_deg)
        plt.pcolormesh(THETA_X, THETA_Y, E2, shading="auto")

        plt.xlabel("Angle X (deg)")
        plt.ylabel("Angle Y (deg)")
        plt.title(f"|E|^2 at {actual_lambda_m*1e9:.1f} nm")
        plt.colorbar(label="|E|^2")
        plt.tight_layout()
        plt.savefig(output_dir / f"{prefix}_{wl_tag}.png", dpi=300)
        plt.close()

    return {
        "target_lambda_m": float(target_lambda_m),
        "actual_lambda_m": actual_lambda_m,
        "actual_lambda_nm": actual_lambda_m * 1e9,
        "f_index_python": f_index_py,
        "f_index_lumerical": f_index_lsf,
        "E2_shape": list(E2.shape),
    }


def save_farfield_data_multiple_wavelengths(
    fdtd,
    monitor_name,
    target_wavelengths_m,
    output_dir,
    prefix="farfield",
    make_plot=True,
):
    """
    Save far-field data for multiple target wavelengths.
    Axes in plots are angle (deg).
    Avoid duplicated saving if multiple target wavelengths map to the same actual frequency index.
    """
    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    c0 = 299792458.0
    f_vec = np.array(fdtd.getdata(monitor_name, "f")).squeeze()
    if f_vec.size == 0:
        raise ValueError(f"No frequency data found in monitor: {monitor_name}")

    lambda_vec = c0 / f_vec

    results = []
    used_f_indices = set()
    saved_angle_axes = False

    for target_lambda_m in target_wavelengths_m:
        print("\n" + "=" * 60)
        print(f"Processing target wavelength: {target_lambda_m * 1e9:.3f} nm")

        f_index_py = int(np.argmin(np.abs(lambda_vec - target_lambda_m)))
        actual_lambda_m = float(lambda_vec[f_index_py])

        if f_index_py in used_f_indices:
            print(
                f"Skipped: target {target_lambda_m*1e9:.3f} nm maps to already-saved "
                f"actual wavelength {actual_lambda_m*1e9:.3f} nm"
            )
            continue

        used_f_indices.add(f_index_py)

        result = save_farfield_data_single_wavelength(
            fdtd=fdtd,
            monitor_name=monitor_name,
            target_lambda_m=target_lambda_m,
            output_dir=output_dir,
            prefix=prefix,
            make_plot=make_plot,
            save_angle_axes=not saved_angle_axes,
        )
        saved_angle_axes = True
        results.append(result)

    with open(output_dir / f"{prefix}_summary.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nAll far-field data saved.")
    return results

# Load scatterer data and generate LSF file at 533nm
if __name__ == "__main__":
    ensure_dir(OUTPUT_ROOT)
    
    result = load_scatterer_config_and_data(scatterer_data_folder, target_wavelength_nm=grating_wavelength_nm)
    
    
    # Save the generated LSF content
    # save_lsf_file(result["lsf_content"], LSF_FSP_FILE.replace(".fsp", ".lsf"))
    
    print("\nGenerated LSF parameters:")
    print(f"  Period: {result['config']['period']:.4e} m")
    print(f"  Scatterer scheme: {result['config']['scatterer_scheme']}")
    print(f"  Wavelength: {result['target_wavelength_nm']} nm")
    print(f"  LC Index values: {len(result['lc_index_data'])} points")
    print(f"  LC Phase data: {len(result['lc_phase_data'])} points")
    
    plt.figure(figsize=(6, 4))
    plt.plot(result["lc_index_data"], np.array(result["lc_phase_data"]), marker='o')
    plt.xlabel("LC Refractive Index")
    plt.ylabel("Phase Delay (rad)")
    plt.title(f"LC Phase vs Index at {result['target_wavelength_nm']} nm")
    plt.grid(True)
    plt.savefig(os.path.join(OUTPUT_ROOT, f"lc_phase_vs_index_{int(result['target_wavelength_nm'])}nm.png"), dpi=300)
    print(f"Lookup table plot saved to: {os.path.abspath(os.path.join(OUTPUT_ROOT, f'lc_phase_vs_index_{int(result["target_wavelength_nm"])}nm.png'))}")
    
    fdtd = lumapi.FDTD(hide=HIDE_LUMERICAL)
    
    if os.path.exists(LSF_FSP_FILE):
        fdtd.load(LSF_FSP_FILE)
    
    # fdtd.eval(result["lsf_content"])
    
    # fdtd.eval("createmodel(10,10);")

    lsf_code = result["lsf_content"]

    if len(lsf_code.strip()) == 0:
        raise ValueError(f"LSF code is empty: {lsf_code}")

    print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
    fdtd.eval(lsf_code)
    
    fdtd.eval("createmodel(15,10);")
    
    time.sleep(1000.0)  # 等待模型创建完成

    
    # fdtd.eval("run;")
    # # fdtd.farfield3d()
    # target_wavelengths_nm = np.linspace(grating_wavelength_nm - 20, grating_wavelength_nm + 20, 5)  # Example: 5 wavelengths from 500nm to 600nm
    # target_wavelengths_m = [wl * 1e-9 for wl in target_wavelengths_nm]

    # farfield_results = save_farfield_data_multiple_wavelengths(
    #     fdtd=fdtd,
    #     monitor_name="R_1.7",
    #     target_wavelengths_m=target_wavelengths_m,
    #     output_dir=OUTPUT_ROOT / "farfield_data",
    #     prefix="grating_farfield",
    #     make_plot=True,
    # )
    
    # fdtd.save(LSF_FSP_FILE)
