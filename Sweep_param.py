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


# =========================================================
# 1) 用户配置
# =========================================================

# Lumerical function 脚本
LSF_FUNCTION_FILE = r"./create_scatterer.lsf"

# 预先建好的 fsp 模型
LSF_FSP_FILE = r"./LC_simulation.fsp"

# 输出根目录
OUTPUT_ROOT = Path("automated_sweep_results")

# monitor names
FIELD_MONITOR_NAME = "monitor"       # 用来取 E
TRANS_MONITOR_NAME = "T_monitor"     # 用来取 transmission / reflectance

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
t_glass = 0.1e-6
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



# =========================================================
# 4) 主流程
# =========================================================

def main():
    ensure_dir(OUTPUT_ROOT)
    fig_dir = OUTPUT_ROOT / "fig_phase_trans"
    ensure_dir(fig_dir)

    n_index = len(index_values)
    n_wave = len(wave_nm)

    phase_matrix = np.zeros((n_index, n_wave))
    trans_matrix = np.zeros((n_index, n_wave))

    fdtd = lumapi.FDTD(hide=HIDE_LUMERICAL)

    try:
        # 读 fsp
        if not os.path.isfile(LSF_FSP_FILE):
            raise FileNotFoundError(f"FSP file NOT found: {LSF_FSP_FILE}")
        fdtd.load(str(LSF_FSP_FILE))

        # 读 lsf function
        if not os.path.isfile(LSF_FUNCTION_FILE):
            raise FileNotFoundError(f"LSF file NOT found: {LSF_FUNCTION_FILE}")

        with open(LSF_FUNCTION_FILE, "r", encoding="utf-8") as f:
            lsf_code = f.read()

        if len(lsf_code.strip()) == 0:
            raise ValueError(f"LSF file is empty: {LSF_FUNCTION_FILE}")

        print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
        fdtd.eval(lsf_code)

        fdtd.switchtolayout()
        fdtd.eval(f'setglobalmonitor("frequency points",{num_wave});')
        print(f"[DEBUG] Set global monitor frequency points to {num_wave}")

        for i, lc_index in enumerate(index_values):
            print(f"[{i+1}/{n_index}] Running LC index = {lc_index:.3f}")

            success = False

            for retry in range(MAX_RETRY):
                print(f"  Attempt {retry+1}/{MAX_RETRY}")

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

                    time.sleep(1.0)
                    fdtd.run()
                    time.sleep(1.0)

                    # 读场
                    Ex, wavelength_from_monitor = extract_complex_field_from_monitor(
                        fdtd, FIELD_MONITOR_NAME
                    )
                    phase, _, _ = average_phase_over_xy(Ex)

                    # 读 transmission / reflectance
                    trans = try_get_transmission(fdtd, TRANS_MONITOR_NAME)

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
                    f"Failed to complete simulation for LC index {lc_index:.3f} after {MAX_RETRY} attempts."
                )

        # 只保存全局 npy
        np.save(OUTPUT_ROOT / "index_values.npy", index_values)
        np.save(OUTPUT_ROOT / "wave_nm.npy", wave_nm)
        np.save(OUTPUT_ROOT / "phase_matrix.npy", phase_matrix)
        np.save(OUTPUT_ROOT / "trans_matrix.npy", trans_matrix)

        # 画图
        plot_summary_figures(
            outdir=fig_dir,
            phase_matrix_raw=phase_matrix,
            trans_matrix=trans_matrix,
            index_values=index_values,
            wave_nm=wave_nm
        )

        print("Done.")
        print(f"Results saved to: {OUTPUT_ROOT.resolve()}")

    finally:
        try:
            fdtd.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()