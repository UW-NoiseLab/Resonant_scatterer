import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import lumapi
from Utils import (
    build_grating_output_dir,
    ensure_dir,
    extract_lut_at_wavelength,
    load_scatterer_results,
    save_farfield_data_multiple_wavelengths,
    save_farfield_data_single_wavelength,
    save_lsf_file,
    validate_scatterer_scheme,
)
from Grating_lsf_gen import grating_DBR_lsf_gen

# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./Grating_simulation_DBR.fsp"



# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = True

DBR_DESIGN_WAVELENGTH = 0.620e-6
DBR_MATERIAL = "TiO2"
DBR_MATERIAL_CHOICES = ("TiO2", "SiN")
DBR_PAIRS = 6

scatterer_data_folder = Path("./output_single_scatterer/scatterer_datas_DBR/Pairs_2_DBR_TiO2_620nm_p350nm_r120nm_2155c068b4")
grating_wavelength_nm = 660.0

SOURCE_POLARIZATION = "x"  # "x" or "y", x is the direction of gradient structure, y is perpendicular to the gradient structure

# 输出根目录
if SOURCE_POLARIZATION == "x":
    OUTPUT_BASE = Path("output_single_grating/grating_datas_DBR_farfield_index3")
else:
    OUTPUT_BASE = Path("output_single_grating/grating_datas_DBR_y_polarization")


def dbr_scatterer_scheme(material, wavelength, dbr_pairs):
    if material not in DBR_MATERIAL_CHOICES:
        valid_materials = ", ".join(DBR_MATERIAL_CHOICES)
        raise ValueError(f"Unknown DBR material '{material}'. Valid materials: {valid_materials}")
    dbr_pairs = int(dbr_pairs)
    if dbr_pairs <= 0:
        raise ValueError(f"DBR pair count must be positive. Got: {dbr_pairs}")
    
    return f"Pairs {dbr_pairs} DBR {material} {float(wavelength) * 1e9:.0f}nm"


def legacy_dbr_scatterer_scheme(material, wavelength):
    return f"DBR {material} {float(wavelength) * 1e9:.0f}nm"


def load_scatterer_config_and_data(data_folder, target_wavelength_nm):
    """Load configuration and lookup tables from scatterer results at target wavelength."""

    data_folder = Path(data_folder)
    scatterer_results = load_scatterer_results(data_folder)
    config = scatterer_results["config"]
    dbr_material = config.get("dbr_material", DBR_MATERIAL)
    dbr_design_wavelength = float(config.get("dbr_design_wavelength", DBR_DESIGN_WAVELENGTH))
    dbr_pairs = int(config.get("dbr_pairs", DBR_PAIRS))
    if dbr_material not in DBR_MATERIAL_CHOICES:
        valid_materials = ", ".join(DBR_MATERIAL_CHOICES)
        raise ValueError(f"Unknown DBR material '{dbr_material}'. Valid materials: {valid_materials}")
    if dbr_pairs <= 0:
        raise ValueError(f"DBR pair count must be positive. Got: {dbr_pairs}")

    expected_schemes = {dbr_scatterer_scheme(dbr_material, dbr_design_wavelength, dbr_pairs)}
    if "dbr_pairs" not in config:
        expected_schemes.add(legacy_dbr_scatterer_scheme(dbr_material, dbr_design_wavelength))
    if config.get("scatterer_scheme") not in expected_schemes:
        validate_scatterer_scheme(
            config,
            dbr_scatterer_scheme(dbr_material, dbr_design_wavelength, dbr_pairs),
            data_folder / "config.json",
        )
    print(f"Verified scatterer scheme: {config['scatterer_scheme']}")
    print(f"DBR material: {dbr_material}")
    print(f"DBR design wavelength: {dbr_design_wavelength * 1e9:.1f} nm")
    print(f"DBR pairs: {dbr_pairs}")

    lut = extract_lut_at_wavelength(scatterer_results, target_wavelength_nm)
    print(f"Target wavelength: {target_wavelength_nm} nm")
    print(f"Actual wavelength used: {lut['actual_wavelength_nm']} nm")
    print(f"Removed {lut['num_removed']} non-monotonic LUT points.")
    print(f"Remaining LUT points: {len(lut['lc_phase_data'])}")
    
    # Generate LSF content using the creategratinglsf function
    lsf_content = grating_DBR_lsf_gen(
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
        LC_index_data=lut["lc_index_data"],
        LC_phase_data=lut["lc_phase_data"],
        wavelength=dbr_design_wavelength,
        material=dbr_material,
        dbr_pairs=dbr_pairs,
    )
    
    return {
        "config": config,
        "wave_nm": scatterer_results["wave_nm"],
        "index_values": scatterer_results["index_values"],
        "phase_matrix": scatterer_results["phase_matrix"],
        "trans_matrix": scatterer_results["trans_matrix"],
        "target_wavelength_nm": lut["actual_wavelength_nm"],
        "dbr_design_wavelength": dbr_design_wavelength,
        "dbr_material": dbr_material,
        "dbr_pairs": dbr_pairs,
        "lc_index_data": lut["lc_index_data"],
        "lc_phase_data": lut["lc_phase_data"],
        "lc_trans_data": lut["lc_trans_data"],
        "lsf_content": lsf_content,
    }

if __name__ == "__main__":
    result = load_scatterer_config_and_data(scatterer_data_folder, target_wavelength_nm=grating_wavelength_nm)
    grating_cells = 15
    steering_angle_deg = 5
    for steering_angle_deg in range(0, 26, 1):  # Example: 0, 5, 10 degrees
        output_config = {
            **result["config"],
            "wavelength_nm": float(result["target_wavelength_nm"]),
            "grating_cells": grating_cells,
            "steering_angle_deg": steering_angle_deg,
            "farfield_monitor": "R_monitor",
            "dbr_design_wavelength": result["dbr_design_wavelength"],
            "dbr_material": result["dbr_material"],
            "dbr_pairs": result["dbr_pairs"],
        }
        
        output_root = ensure_dir(build_grating_output_dir(OUTPUT_BASE, output_config))
        output_root = Path(output_root) / f"steering_{steering_angle_deg}deg_grating_cells_{grating_cells}"
        output_root = ensure_dir(output_root)
        
        if output_root.exists():
            if (output_root / "farfield_data").exists():
                if not any((output_root / "farfield_data").iterdir()):
                    print(f"Empty output directory: {output_root}. Running simulation...")
                else:
                    print(f"Output directory already contains data: {output_root}. Skipping simulation.")
                    continue
        
        # Save the generated LSF content
        # save_lsf_file(result["lsf_content"], LSF_FSP_FILE.replace(".fsp", ".lsf"))
        
        print("\nGenerated LSF parameters:")
        print(f"  Period: {result['config']['period']:.4e} m")
        print(f"  Scatterer scheme: {result['config']['scatterer_scheme']}")
        print(f"  DBR material: {result['dbr_material']}")
        print(f"  DBR design wavelength: {result['dbr_design_wavelength'] * 1e9:.1f} nm")
        print(f"  DBR pairs: {result['dbr_pairs']}")
        print(f"  Wavelength: {result['target_wavelength_nm']} nm")
        print(f"  LC Index values: {len(result['lc_index_data'])} points")
        print(f"  LC Phase data: {len(result['lc_phase_data'])} points")
        
        plt.figure(figsize=(6, 4))
        plt.plot(result["lc_index_data"], np.array(result["lc_phase_data"]), marker='o')
        plt.xlabel("LC Refractive Index")
        plt.ylabel("Phase Delay (rad)")
        plt.title(f"LC Phase vs Index at {result['target_wavelength_nm']} nm")
        plt.grid(True)
        lut_plot_path = output_root / f"lc_phase_vs_index_{int(result['target_wavelength_nm'])}nm.png"
        plt.savefig(lut_plot_path, dpi=300)
        print(f"Lookup table plot saved to: {os.path.abspath(lut_plot_path)}")
        
        fdtd = lumapi.FDTD(hide=HIDE_LUMERICAL)
        
        if os.path.exists(LSF_FSP_FILE):
            fdtd.load(LSF_FSP_FILE)

        lsf_code = result["lsf_content"]

        if len(lsf_code.strip()) == 0:
            raise ValueError(f"LSF code is empty: {lsf_code}")

        print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
        fdtd.eval(lsf_code)
        
        fdtd.eval(f"createmodel({grating_cells},{steering_angle_deg});")
        
        if SOURCE_POLARIZATION == "x":
            pass
        elif SOURCE_POLARIZATION == "y":
            fdtd.eval("setnamed('source', 'polarization angle', 90);")
        
        fdtd.eval("run;")

        target_wavelengths_nm = np.linspace(grating_wavelength_nm - 20, grating_wavelength_nm + 20, 5)  # Example: 5 wavelengths from 500nm to 600nm
        target_wavelengths_m = [wl * 1e-9 for wl in target_wavelengths_nm]

        farfield_results = save_farfield_data_multiple_wavelengths(
            fdtd=fdtd,
            monitor_name="R_monitor",
            target_wavelengths_m=target_wavelengths_m,
            output_dir=output_root / "farfield_data",
            prefix="grating_farfield",
            make_plot=True,
        )
