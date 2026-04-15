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


# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./Grating_simulation.fsp"


# 输出根目录
OUTPUT_ROOT = Path("./automated_grating_results")


# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = False

grating_wavelength_nm = 533.0

scatterer_data_folder = Path("./automated_result_glass_3/automated_sweep_results_r107nm_p322nm")


def load_scatterer_config_and_data(data_folder, target_wavelength_nm=533.0):
    """Load configuration and lookup tables from scatterer results at target wavelength."""
    from Grating_lsf_gen import grating_lsf_gen
    
    config_file = data_folder / "config.json"
    with open(config_file, 'r') as f:
        config = json.load(f)
    
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
    lc_phase_data = phase_matrix[:, wl_idx].tolist()
    lc_trans_data = trans_matrix[:, wl_idx].tolist()
    lc_index_data = index_values.tolist()
    
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


# Load scatterer data and generate LSF file at 533nm
if __name__ == "__main__":
    ensure_dir(OUTPUT_ROOT)
    
    result = load_scatterer_config_and_data(scatterer_data_folder, target_wavelength_nm=grating_wavelength_nm)
    
    
    # Save the generated LSF content
    save_lsf_file(result["lsf_content"], LSF_FSP_FILE.replace(".fsp", ".lsf"))
    
    print("\nGenerated LSF parameters:")
    print(f"  Period: {result['config']['period']:.4e} m")
    print(f"  Wavelength: {result['target_wavelength_nm']} nm")
    print(f"  LC Index values: {len(result['lc_index_data'])} points")
    print(f"  LC Phase data: {len(result['lc_phase_data'])} points")
    
    plt.figure(figsize=(6, 4))
    plt.plot(result["lc_index_data"], (np.array(result["lc_phase_data"]) - np.array(result["lc_phase_data"][0])) % (2*np.pi), marker='o')
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
    
    # print(result["lsf_content"])
    
    # fdtd.eval("createmodel(10,10);")
    with open(LSF_FSP_FILE.replace(".fsp", ".lsf"), "r", encoding="utf-8") as f:
        lsf_code = f.read()

    if len(lsf_code.strip()) == 0:
        raise ValueError(f"LSF file is empty: {LSF_FSP_FILE.replace(".fsp", ".lsf")}")

    print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
    fdtd.eval(lsf_code)
    
    fdtd.eval("createmodel(10,10);")
    
    fdtd.save(LSF_FSP_FILE)


