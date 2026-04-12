import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
# =========================================================
# 3) 工具函数
# =========================================================

def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def wrap_to_2pi(x: np.ndarray) -> np.ndarray:
    return np.mod(x, 2 * np.pi)


def extract_complex_field_from_monitor(fdtd, monitor_name: str):
    """
    从 field monitor 读取 E，并取 Ex 分量。
    假设 E["E"] 最后一维前三个分量对应 Ex/Ey/Ez。
    """
    # E = fdtd.getresult(monitor_name, "E")
    # field = np.array(E["E"])

    # if field.size == 0:
    #     raise RuntimeError(f"Empty field data from monitor '{monitor_name}'.")

    # 取 Ex 分量
    Ex = fdtd.getdata(monitor_name, "Ex").squeeze()
    c = 299792458.0
    freq = np.array(fdtd.getdata(monitor_name,'f')).flatten()
    wavelength = c / freq
    return Ex, wavelength


def average_phase_over_xy(Ex: np.ndarray):
    """
    对每个 wavelength，把 Ex 在横向平面平均后取相位。
    支持 squeeze 后:
    - (Nx, Ny, Nlambda)
    - (Nx, Ny) 只有一个波长
    """
    arr = np.squeeze(Ex)

    if arr.ndim == 3:
        complex_mean = arr.mean(axis=(0, 1))
    elif arr.ndim == 2:
        complex_mean = np.array([arr.mean()])
    else:
        raise ValueError(f"Unexpected Ex shape after squeeze: {arr.shape}")

    phase = np.angle(complex_mean)
    amplitude = np.abs(complex_mean)
    return phase, amplitude, complex_mean


def try_get_transmission(fdtd, monitor_name: str):
    """
    优先直接用 transmission(monitor_name)，失败再尝试 getresult(..., "T")
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


def plot_summary_figures(
    outdir: Path,
    phase_matrix_raw: np.ndarray,
    trans_matrix: np.ndarray,
    index_values: np.ndarray,
    wave_nm: np.ndarray
):
    """
    生成你 MATLAB 版本对应的图：
    1) 每个波长一张，左轴 phase，右轴 reflectance
    2) phase map
    3) reflectance map
    """
    ensure_dir(outdir)

    # 相位归一化：以第一行作为参考
    phase_matrix = phase_matrix_raw - phase_matrix_raw[0:1, :]
    phase_matrix = wrap_to_2pi(phase_matrix)

    # -----------------------------------------------------
    # 每个 wavelength 一张双 y 轴图
    # -----------------------------------------------------
    for i in range(len(wave_nm)):
        fig, ax1 = plt.subplots(figsize=(6, 4.5))

        ax1.plot(index_values, phase_matrix[:, i], linewidth=1.5)
        ax1.set_xlabel("Refractive Index", fontsize=12)
        ax1.set_ylabel("Phase (rad)", fontsize=12)
        ax1.set_ylim(0, 2 * np.pi)
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(index_values, trans_matrix[:, i], linestyle="--", linewidth=1.5)
        ax2.set_ylabel("Reflectance", fontsize=12)
        ax2.set_ylim(0, 1)

        ax1.set_title(f"$\\lambda$ = {wave_nm[i]:.1f} nm", fontsize=13)
        ax1.tick_params(labelsize=11)
        ax2.tick_params(labelsize=11)

        fig.tight_layout()
        fig.savefig(outdir / f"phase_trans_{int(round(wave_nm[i]))}nm.png", dpi=300)
        plt.close(fig)

    # -----------------------------------------------------
    # Phase map
    # -----------------------------------------------------
    fig = plt.figure(figsize=(7, 5))
    plt.imshow(
        phase_matrix,
        aspect="auto",
        origin="lower",
        extent=[wave_nm[0], wave_nm[-1], index_values[0], index_values[-1]]
    )
    plt.xlabel("Wavelength (nm)", fontsize=12)
    plt.ylabel("Refractive Index", fontsize=12)
    plt.title("Phase Map", fontsize=13)
    cbar = plt.colorbar()
    cbar.set_label("Phase (rad)", fontsize=12)
    plt.clim(0, 2 * np.pi)
    plt.tight_layout()
    plt.savefig(outdir / "phase_map.png", dpi=300)
    plt.close(fig)

    # -----------------------------------------------------
    # Reflectance map
    # -----------------------------------------------------
    fig = plt.figure(figsize=(7, 5))
    plt.imshow(
        trans_matrix,
        aspect="auto",
        origin="lower",
        extent=[wave_nm[0], wave_nm[-1], index_values[0], index_values[-1]]
    )
    plt.xlabel("Wavelength (nm)", fontsize=12)
    plt.ylabel("Refractive Index", fontsize=12)
    plt.title("Reflectance Map", fontsize=13)
    cbar = plt.colorbar()
    cbar.set_label("Reflectance", fontsize=12)
    plt.clim(0, 1)
    plt.tight_layout()
    plt.savefig(outdir / "reflectance_map.png", dpi=300)
    plt.close(fig)

