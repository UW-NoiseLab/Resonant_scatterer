import os
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import time
# =========================================================
# 1) Lumerical lumapi path
#    这里改成你自己电脑上的 Lumerical API 路径
# =========================================================
# Windows 例子：
# sys.path.append(r"C:\Program Files\Lumerical\v242\api\python")

# Linux 例子：
# sys.path.append("/opt/lumerical/v242/api/python")

import lumapi

# =========================================================
# 2) 用户配置
# =========================================================

# 你的 Lumerical function 脚本
LSF_FUNCTION_FILE = r"./create_scatterer.lsf"

# Lumericial fsp 模型文件（可选，如果你在 LSF 里已经建模了，就不需要了）
LSF_FSP_FILE = r"./LC_simulation.fsp"

# 输出根目录
OUTPUT_ROOT = Path("automated_sweep_results")

# monitor names
FIELD_MONITOR_NAME = "monitor"     # 用来取 E
TRANS_MONITOR_NAME = "T_monitor"    # 用来取 transmission

# 是否隐藏 Lumerical GUI
HIDE_LUMERICAL = True

# =========================================================
# 3) 建模参数
# =========================================================

period = 0.36e-6
t_Al = 0.30e-6
t_spacer = 0.17e-6
t_LC = 0.50e-6
t_ITO = 0.05e-6
t_glass = 2.0e-6
r_scatter = 0.12e-6
t_scatter = 0.20e-6
lambda_start = 0.40e-6
lambda_stop = 0.70e-6

# 你现在想扫的 LC index 范围
index_values = np.arange(1.55, 1.75, 0.01)

# wavelength axis for plotting / bookkeeping
num_wave = 151
wave_nm = np.linspace(400, 700, num_wave)
wave_m = wave_nm * 1e-9


# =========================================================
# 4) 工具函数
# =========================================================

def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def save_array_csv(path: Path, array: np.ndarray, header: str = ""):
    np.savetxt(path, array, delimiter=",", header=header, comments="")


def extract_complex_field_from_monitor(fdtd, monitor_name: str):
    """
    从 field monitor 读取 E，并返回 Ex, wavelength
    假设结果里存在 Ex 和 lambda / f 信息
    """
    E = fdtd.getresult(monitor_name, "E")

    # if "Ex" not in E:
    #     print(f"Monitor '{monitor_name}' result {list(E.keys())}")

    Ex = np.array(E["E"])[...,0].squeeze()  # 取 Ex 分量，去掉单维度

    # 常见情况下 wavelength 可能叫 lambda
    if "lambda" in E:
        wavelength = np.array(E["lambda"]).flatten()
    elif "f" in E:
        c = 299792458.0
        freq = np.array(E["f"]).flatten()
        wavelength = c / freq
    else:
        wavelength = None

    return Ex, wavelength, E


def average_phase_over_xy(Ex: np.ndarray):
    """
    对每个 wavelength，把 Ex 在横向平面平均后取相位。
    兼容常见维度：
    - (Nx, Ny, Nlambda)
    - (Nx, Ny, 1, Nlambda)
    - (Nx, Ny, Nz, Nlambda) 但通常 2D monitor 不会这样
    """
    arr = np.squeeze(Ex)

    if arr.ndim == 3:
        # (Nx, Ny, Nlambda)
        complex_mean = arr.mean(axis=(0, 1))
    elif arr.ndim == 2:
        # 只有一个波长
        complex_mean = np.array([arr.mean()])
    else:
        raise ValueError(f"Unexpected Ex shape after squeeze: {arr.shape}")

    phase = np.angle(complex_mean)
    amplitude = np.abs(complex_mean)
    return phase, amplitude, complex_mean


def try_get_transmission(fdtd, monitor_name: str):
    """
    优先直接用 transmission(monitor_name)；
    如果不行，再尝试从 getresult 里找 T。
    """
    try:
        t = fdtd.transmission(monitor_name)
        t = np.array(t).flatten()
        return t
    except Exception:
        pass

    try:
        res = fdtd.getresult(monitor_name, "T")
        if "T" in res:
            return np.array(res["T"]).flatten()
    except Exception:
        pass

    raise RuntimeError(f"Failed to get transmission from monitor '{monitor_name}'.")


def plot_heatmap(data, x_vals, y_vals, xlabel, ylabel, title, save_path: Path):
    ''' 绘制热图 
    data: 2D array, shape (len(y_vals), len(x_vals))
    '''
    
    plt.figure(figsize=(8, 5))
    plt.imshow(
        data,
        aspect="auto",
        origin="lower",
        extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
    )
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def plot_selected_curves(wave_nm, data, index_values, selected_indices, ylabel, title, save_path: Path):
    plt.figure(figsize=(8, 5))
    for idx in selected_indices:
        plt.plot(wave_nm, data[idx], label=f"n={index_values[idx]:.2f}")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


# =========================================================
# 5) 主流程
# =========================================================

def main():
    ensure_dir(OUTPUT_ROOT)
    ensure_dir(OUTPUT_ROOT / "figures")
    ensure_dir(OUTPUT_ROOT / "raw_per_index")

    # 初始化结果矩阵
    n_index = len(index_values)
    n_wave = len(wave_nm)

    phase_matrix = np.zeros((n_index, n_wave))
    trans_matrix = np.zeros((n_index, n_wave))
    amp_matrix = np.zeros((n_index, n_wave))

    # 打开 Lumerical FDTD
    fdtd = lumapi.FDTD(hide=HIDE_LUMERICAL)
    fdtd.load(str(LSF_FSP_FILE))  # 如果你有预先建好的 fsp 模型，可以直接加载
    # time.sleep(5)  # 等待模型加载完成


    # 读入你定义函数的 lsf 文件
    if not os.path.isfile(LSF_FUNCTION_FILE):
        raise FileNotFoundError(f"LSF file NOT found: {LSF_FUNCTION_FILE}")
    with open(LSF_FUNCTION_FILE, "r", encoding="utf-8") as f:
        lsf_code = f.read()
        
    if len(lsf_code.strip()) == 0:
        raise ValueError(f"LSF file is empty: {LSF_FUNCTION_FILE}")

    print(f"[DEBUG] Successfully read LSF file ({len(lsf_code)} characters)")
    fdtd.eval(lsf_code)
    fdtd.switchtolayout()
    
    fdtd.eval(f"""setglobalmonitor("frequency points",{num_wave});""")
    print(f"[DEBUG] Set global monitor frequency points to {num_wave}")
    

    # 可选：保存初始文件
    # fdtd.save(str((OUTPUT_ROOT / "initial_model.fsp").resolve()))
    
    MAX_RETRY = 6  # 最多重试次数

    for i, lc_index in enumerate(index_values):
        print(f"[{i+1}/{n_index}] Running LC index = {lc_index:.3f}")
    
        SUCCESS = False
        for retry in range(MAX_RETRY):
            fdtd.eval("switchtolayout;")  # 确保在布局模式下建模
            fdtd.eval("deleteall;")  # 每次循环先清空布局，确保干净的建模环境
            bd_model = f"""build_unit_cell_model({period}, {t_Al}, {t_spacer}, {t_LC}, {t_ITO}, {t_glass}, {r_scatter}, {t_scatter}, {lambda_start}, {lambda_stop});"""

            fdtd.eval(bd_model)
            fdtd.select("liquid_crystal")
            fdtd.set("index", float(lc_index))
            time.sleep(1.5)  # 等建模命令执行完
            print(f"  Attempt {retry+1} of {MAX_RETRY}...")
            try:
                fdtd.run()
                time.sleep(1.5)
                # -------- 读取场 --------
                Ex, wavelength_from_monitor, _ = extract_complex_field_from_monitor(fdtd, FIELD_MONITOR_NAME)
                phase, amp, complex_mean = average_phase_over_xy(Ex)

                # -------- 读取 transmission --------
                trans = try_get_transmission(fdtd, TRANS_MONITOR_NAME)

                SUCCESS = True
                break
            except Exception as e:
                print(f"    Attempt {retry+1} failed: {e}")
                if retry == MAX_RETRY - 1:
                    raise
                time.sleep(1)

        if not SUCCESS:
            raise RuntimeError(
                f"Failed to complete simulation for LC index {lc_index:.3f} after {MAX_RETRY} attempts."
            )

        # -------- 读取 transmission --------
        trans = try_get_transmission(fdtd, TRANS_MONITOR_NAME)

        # 长度校验
        if len(phase) != n_wave:
            raise RuntimeError(
                f"Phase length {len(phase)} != expected n_wave {n_wave}. "
                f"Check monitor spectral points."
            )
        if len(trans) != n_wave:
            raise RuntimeError(
                f"Transmission length {len(trans)} != expected n_wave {n_wave}. "
                f"Check source/monitor spectral points."
            )

        phase_matrix[i, :] = phase
        trans_matrix[i, :] = trans
        amp_matrix[i, :] = amp

        # 每个 index 单独存一份
        per_index_dir = OUTPUT_ROOT / "raw_per_index" / f"index_{lc_index:.3f}"
        ensure_dir(per_index_dir)

        sio.savemat(
            per_index_dir / "result.mat",
            {
                "lc_index": lc_index,
                "wave_nm": wave_nm,
                "phase": phase,
                "transmission": trans,
                "amplitude": amp,
                "complex_mean_real": np.real(complex_mean),
                "complex_mean_imag": np.imag(complex_mean),
            },
        )

        save_array_csv(
            per_index_dir / "curves.csv",
            np.column_stack([wave_nm, phase, trans, amp]),
            header="wavelength_nm,phase_rad,transmission,field_mean_amplitude",
        )

        # 可选：每次都保存一个 fsp
        # fdtd.save(str((per_index_dir / "model.fsp").resolve()))

    # =====================================================
    # 全局保存
    # =====================================================
    sio.savemat(
        OUTPUT_ROOT / "sweep_all.mat",
        {
            "index_values": index_values,
            "wave_nm": wave_nm,
            "phase_matrix": phase_matrix,
            "trans_matrix": trans_matrix,
            "amp_matrix": amp_matrix,
        },
    )

    np.save(OUTPUT_ROOT / "index_values.npy", index_values)
    np.save(OUTPUT_ROOT / "wave_nm.npy", wave_nm)
    np.save(OUTPUT_ROOT / "phase_matrix.npy", phase_matrix)
    np.save(OUTPUT_ROOT / "trans_matrix.npy", trans_matrix)
    np.save(OUTPUT_ROOT / "amp_matrix.npy", amp_matrix)

    save_array_csv(
        OUTPUT_ROOT / "phase_matrix.csv",
        phase_matrix,
        header="phase matrix, rows=index sweep, cols=wavelength",
    )
    save_array_csv(
        OUTPUT_ROOT / "trans_matrix.csv",
        trans_matrix,
        header="transmission matrix, rows=index sweep, cols=wavelength",
    )
    save_array_csv(
        OUTPUT_ROOT / "amp_matrix.csv",
        amp_matrix,
        header="field amplitude matrix, rows=index sweep, cols=wavelength",
    )

    # =====================================================
    # 画图
    # =====================================================
    fig_dir = OUTPUT_ROOT / "figures"

    plot_heatmap(
        phase_matrix,
        x_vals=wave_nm,
        y_vals=index_values,
        xlabel="Wavelength (nm)",
        ylabel="LC index",
        title="Phase vs Wavelength and LC Index",
        save_path=fig_dir / "phase_heatmap.png",
    )

    plot_heatmap(
        trans_matrix,
        x_vals=wave_nm,
        y_vals=index_values,
        xlabel="Wavelength (nm)",
        ylabel="LC index",
        title="Transmission vs Wavelength and LC Index",
        save_path=fig_dir / "transmission_heatmap.png",
    )

    plot_heatmap(
        amp_matrix,
        x_vals=wave_nm,
        y_vals=index_values,
        xlabel="Wavelength (nm)",
        ylabel="LC index",
        title="Mean Field Amplitude vs Wavelength and LC Index",
        save_path=fig_dir / "field_amplitude_heatmap.png",
    )

    # 选几条曲线画出来
    selected_indices = np.linspace(0, len(index_values) - 1, 5, dtype=int)

    plot_selected_curves(
        wave_nm,
        phase_matrix,
        index_values,
        selected_indices,
        ylabel="Phase (rad)",
        title="Phase Curves for Selected LC Indices",
        save_path=fig_dir / "phase_selected_curves.png",
    )

    plot_selected_curves(
        wave_nm,
        trans_matrix,
        index_values,
        selected_indices,
        ylabel="Transmission",
        title="Transmission Curves for Selected LC Indices",
        save_path=fig_dir / "transmission_selected_curves.png",
    )

    print("Done.")
    print(f"Results saved to: {OUTPUT_ROOT.resolve()}")

    
    fdtd.close()
     


if __name__ == "__main__":
    main()