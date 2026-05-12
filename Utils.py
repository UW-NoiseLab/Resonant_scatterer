import hashlib
import json
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


C0 = 299792458.0
TWO_PI = 2 * np.pi
SCATTERER_RESULT_FILES = (
    "config.json",
    "index_values.npy",
    "wave_nm.npy",
    "phase_matrix.npy",
    "trans_matrix.npy",
)


# =========================================================
# Filesystem and namespace helpers
# =========================================================

def ensure_dir(path: Path):
    """Create a directory tree if needed and keep call sites concise."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _json_ready(value):
    """Convert numpy/path objects into stable JSON-serializable values."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return [_json_ready(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_ready(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def stable_config_hash(config, length=10):
    """Return a short deterministic hash for a config dictionary."""
    payload = json.dumps(_json_ready(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:length]


def sanitize_token(value):
    """Make a compact filesystem-safe token from a human-readable string."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")
    return token or "unnamed"


def nm_tag(value_m, prefix, decimals=0):
    """Format a length in meters as a compact nanometer token, such as p340nm."""
    value_nm = float(value_m) * 1e9
    if decimals == 0:
        return f"{prefix}{int(round(value_nm))}nm"
    formatted = f"{value_nm:.{decimals}f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"{prefix}{formatted}nm"


def wavelength_tag(wavelength_nm, prefix="wavelength", decimals=1):
    """Format a wavelength already in nm as a compact token."""
    formatted = f"{float(wavelength_nm):.{decimals}f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"{prefix}{formatted}nm"


def build_scatterer_output_dir(base_dir, config):
    """
    Build the canonical scatterer output namespace:
    output/scatterer_datas/{scheme}_p{period}_r{radius}_{hash}.
    """
    scheme = sanitize_token(config["scatterer_scheme"])
    period = nm_tag(config["period"], "p")
    radius = nm_tag(config["r_scatter"], "r")
    config_hash = stable_config_hash(config)
    return Path(base_dir) / f"{scheme}_{period}_{radius}_{config_hash}"


def hash_from_namespace(path):
    """Extract the trailing hash token from one of our managed output folders."""
    if path is None:
        return None
    name = Path(path).name
    parts = name.rsplit("_", 1)
    if len(parts) != 2:
        return None
    suffix = parts[1]
    if re.fullmatch(r"[0-9a-fA-F]+", suffix):
        return suffix
    return None


def build_grating_output_dir(base_dir, config):
    """
    Build the canonical grating output namespace:
    output/grating_datas/{scheme}_p{period}_r{radius}_wavelength{wl}_{hash}.

    The hash intentionally reuses the scatterer sweep hash when config_hash is
    present, so grating folders can be traced back to the exact scatterer LUT
    that produced them. Wavelength/grating settings stay readable in the folder
    name instead of changing the shared hash.
    """
    scheme = sanitize_token(config["scatterer_scheme"])
    period = nm_tag(config["period"], "p")
    radius = nm_tag(config["r_scatter"], "r")
    wavelength = wavelength_tag(config["wavelength_nm"])
    config_hash = (
        config.get("scatterer_config_hash")
        or hash_from_namespace(config.get("output_namespace"))
        or hash_from_namespace(config.get("output_root"))
        or config.get("config_hash")
    )
    if config_hash is None:
        config_hash = stable_config_hash(config)
    return Path(base_dir) / f"{scheme}_{period}_{radius}_{wavelength}_{config_hash}"


def save_json(path, data):
    """Write JSON with the same formatting everywhere in the project."""
    path = Path(path)
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_json_ready(data), f, indent=2)


def save_lsf_file(lsf_content, output_file):
    """Save generated Lumerical script text and create the parent directory."""
    output_file = Path(output_file)
    ensure_dir(output_file.parent)
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(lsf_content)
    print(f"LSF file saved to: {output_file}")


def wrap_to_2pi(x: np.ndarray) -> np.ndarray:
    """Map phase values to the [0, 2pi) interval."""
    return np.mod(x, TWO_PI)


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
    freq = np.array(fdtd.getdata(monitor_name,'f')).flatten()
    wavelength = C0 / freq
    return Ex, wavelength


def extract_monitor_wavelengths(fdtd, monitor_name: str):
    """Return the wavelength axis of a monitor in meters."""
    freq = np.array(fdtd.getdata(monitor_name, "f")).flatten()
    if freq.size == 0:
        raise RuntimeError(f"Empty frequency data from monitor '{monitor_name}'.")
    return C0 / freq


def average_phase_over_xy(Ex: np.ndarray):
    """
    对每个 wavelength，把 Ex 在横向平面平均后取相位。
    支持 squeeze 后:
    - (Nx, Ny, Nlambda)
    """
    arr = np.squeeze(Ex)

    if arr.ndim == 3:
        complex_mean = arr.mean(axis=(0, 1))
    else:
        raise ValueError(f"Unexpected Ex shape after squeeze: {arr.shape}")

    phase = np.angle(complex_mean)
    amplitude = np.abs(complex_mean)
    # print(f"Average amplitude over xy: {amplitude}")
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


# =========================================================
# Scatterer result and wavelength-selection helpers
# =========================================================

def load_scatterer_results(data_folder):
    """Load a completed scatterer sweep folder and validate required files."""
    data_folder = Path(data_folder)
    missing = [
        str(data_folder / filename)
        for filename in SCATTERER_RESULT_FILES
        if not (data_folder / filename).exists()
    ]
    if missing:
        raise FileNotFoundError("Missing scatterer result files:\n" + "\n".join(missing))

    with open(data_folder / "config.json", "r", encoding="utf-8") as f:
        config = json.load(f)

    return {
        "config": config,
        "index_values": np.load(data_folder / "index_values.npy"),
        "wave_nm": np.load(data_folder / "wave_nm.npy"),
        "phase_matrix": np.load(data_folder / "phase_matrix.npy"),
        "trans_matrix": np.load(data_folder / "trans_matrix.npy"),
    }


def has_scatterer_results(data_folder):
    """Return True only when a scatterer sweep folder is complete enough to reuse."""
    try:
        load_scatterer_results(data_folder)
        return True
    except FileNotFoundError:
        return False


def validate_scatterer_scheme(config, expected_scheme, config_file=None):
    """Fail early if a LUT folder was generated with a different scatterer scheme."""
    config_scheme = config.get("scatterer_scheme")
    location = f" {config_file}" if config_file is not None else ""
    if config_scheme is None:
        raise ValueError(
            f"Scatterer config{location} does not contain 'scatterer_scheme'. "
            "Regenerate the scatterer LUT so grating code can verify the geometry scheme."
        )
    if config_scheme != expected_scheme:
        raise ValueError(
            "Scatterer scheme mismatch: "
            f"LUT folder was generated with '{config_scheme}', "
            f"but this run is configured for '{expected_scheme}'."
        )


def monotonic_lut_mask(phase_values):
    """Keep the strictly increasing phase samples required by Lumerical interp()."""
    keep_mask = np.zeros(len(phase_values), dtype=bool)
    keep_mask[0] = True
    current_max = phase_values[0]
    for i in range(1, len(phase_values)):
        if phase_values[i] > current_max:
            keep_mask[i] = True
            current_max = phase_values[i]
    return keep_mask


def extract_lut_at_wavelength(scatterer_results, target_wavelength_nm):
    """Extract monotonic LC-index/phase LUT data nearest to a target wavelength."""
    wave_nm = scatterer_results["wave_nm"]
    index_values = scatterer_results["index_values"]
    phase_matrix = scatterer_results["phase_matrix"]
    trans_matrix = scatterer_results["trans_matrix"]

    wl_idx = int(np.argmin(np.abs(wave_nm - target_wavelength_nm)))
    actual_wl_nm = float(wave_nm[wl_idx])

    lc_phase_data = wrap_to_2pi(phase_matrix[:, wl_idx] - phase_matrix[0, wl_idx])
    lc_trans_data = trans_matrix[:, wl_idx]
    lc_index_data = index_values.copy()

    keep_mask = monotonic_lut_mask(lc_phase_data)
    return {
        "target_wavelength_nm": float(target_wavelength_nm),
        "actual_wavelength_nm": actual_wl_nm,
        "wavelength_index": wl_idx,
        "num_removed": int(len(lc_phase_data) - np.sum(keep_mask)),
        "lc_phase_data": lc_phase_data[keep_mask].tolist(),
        "lc_trans_data": lc_trans_data[keep_mask].tolist(),
        "lc_index_data": lc_index_data[keep_mask].tolist(),
        "keep_mask": keep_mask,
    }


def phase_uniformity_metrics(phase_values):
    """Measure circular 0-to-2pi phase coverage and spacing uniformity."""
    phase_values = wrap_to_2pi(np.asarray(phase_values, dtype=float))
    sorted_phase = np.sort(phase_values)
    circular_gaps = np.diff(np.r_[sorted_phase, sorted_phase[0] + TWO_PI])

    ideal_gap = TWO_PI / len(sorted_phase)
    gap_cv = float(np.std(circular_gaps) / ideal_gap) if ideal_gap > 0 else np.inf
    uniformity = float(np.clip(1.0 - gap_cv, 0.0, 1.0))
    coverage = float((TWO_PI - np.max(circular_gaps)) / TWO_PI)

    return {
        "coverage_fraction": coverage,
        "largest_gap_rad": float(np.max(circular_gaps)),
        "gap_cv": gap_cv,
        "uniformity_score": uniformity,
    }


def score_scatterer_wavelengths(scatterer_results, min_mean_transmission):
    """Rank wavelengths by phase coverage, phase uniformity, and transmission."""
    wave_nm = scatterer_results["wave_nm"]
    index_values = scatterer_results["index_values"]
    phase_matrix = scatterer_results["phase_matrix"]
    trans_matrix = scatterer_results["trans_matrix"]

    rows = []
    for wl_idx, wavelength_nm in enumerate(wave_nm):
        phase = wrap_to_2pi(phase_matrix[:, wl_idx] - phase_matrix[0, wl_idx])
        trans = np.asarray(trans_matrix[:, wl_idx], dtype=float)
        trans_power = np.abs(trans)
        keep_mask = monotonic_lut_mask(phase)
        kept_phase = phase[keep_mask]

        all_metrics = phase_uniformity_metrics(phase)
        kept_metrics = phase_uniformity_metrics(kept_phase)
        mean_trans = float(np.mean(trans_power))
        min_trans = float(np.min(trans_power))
        raw_mean_trans = float(np.mean(trans))
        raw_min_trans = float(np.min(trans))

        trans_score = float(np.clip(mean_trans, 0.0, 1.0))
        if mean_trans < min_mean_transmission:
            trans_score *= mean_trans / max(min_mean_transmission, 1e-12)

        score = (
            0.45 * kept_metrics["coverage_fraction"]
            + 0.25 * kept_metrics["uniformity_score"]
            + 0.20 * trans_score
            + 0.10 * (len(kept_phase) / len(index_values))
        )

        rows.append(
            {
                "wavelength_nm": float(wavelength_nm),
                "score": float(score),
                "mean_transmission": mean_trans,
                "min_transmission": min_trans,
                "raw_mean_transmission": raw_mean_trans,
                "raw_min_transmission": raw_min_trans,
                "coverage_fraction_all_points": all_metrics["coverage_fraction"],
                "coverage_fraction_lut": kept_metrics["coverage_fraction"],
                "uniformity_score_all_points": all_metrics["uniformity_score"],
                "uniformity_score_lut": kept_metrics["uniformity_score"],
                "largest_gap_lut_rad": kept_metrics["largest_gap_rad"],
                "gap_cv_lut": kept_metrics["gap_cv"],
                "usable_lut_points": int(np.sum(keep_mask)),
                "total_index_points": int(len(index_values)),
            }
        )

    rows.sort(key=lambda item: item["score"], reverse=True)
    return rows


def save_candidate_report(candidates, output_dir):
    """Save wavelength candidate ranking as JSON and spreadsheet-friendly CSV."""
    output_dir = ensure_dir(output_dir)
    save_json(output_dir / "wavelength_candidates.json", candidates)

    header = [
        "rank",
        "wavelength_nm",
        "score",
        "mean_transmission",
        "min_transmission",
        "coverage_fraction_lut",
        "uniformity_score_lut",
        "largest_gap_lut_rad",
        "usable_lut_points",
        "total_index_points",
    ]
    with open(output_dir / "wavelength_candidates.csv", "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for rank, row in enumerate(candidates, start=1):
            values = [rank] + [row[key] for key in header[1:]]
            f.write(",".join(str(value) for value in values) + "\n")


def plot_candidate_summary(candidates, output_dir):
    """Plot wavelength score, phase coverage, and transmission on shared x-axis."""
    output_dir = ensure_dir(output_dir)
    ordered = sorted(candidates, key=lambda item: item["wavelength_nm"])
    wave_nm = np.array([item["wavelength_nm"] for item in ordered])

    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    axes[0].plot(wave_nm, [item["score"] for item in ordered], marker="o")
    axes[0].set_ylabel("Score")
    axes[0].grid(True)

    axes[1].plot(wave_nm, [item["coverage_fraction_lut"] for item in ordered], marker="o")
    axes[1].set_ylabel("LUT coverage")
    axes[1].set_ylim(0, 1.02)
    axes[1].grid(True)

    axes[2].plot(wave_nm, [item["mean_transmission"] for item in ordered], marker="o")
    axes[2].set_ylabel("Mean T")
    axes[2].set_xlabel("Wavelength (nm)")
    axes[2].set_ylim(0, 1.02)
    axes[2].grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "wavelength_candidate_scores.png", dpi=300)
    plt.close(fig)


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


# =========================================================
# Far-field helpers
# =========================================================

def _parse_steering_angle_deg(folder_name):
    """Extract the designed steering angle from names like steering_12deg_..."""
    match = re.search(r"steering_(-?\d+(?:p\d+|\.\d+)?)deg", folder_name)
    if match is None:
        return None
    return float(match.group(1).replace("p", "."))


def _parse_farfield_profile_wavelength_nm(csv_path):
    """Extract the wavelength from grating_farfield_620p2nm_x_profile.csv."""
    match = re.search(r"grating_farfield_(\d+(?:p\d+|\.\d+)?)nm_x_profile\.csv$", csv_path.name)
    if match is None:
        return None
    return float(match.group(1).replace("p", "."))


def _closest_farfield_x_profile_csv(farfield_dir, target_wavelength_nm):
    """Return the x-profile CSV in farfield_dir nearest to target_wavelength_nm."""
    candidates = []
    for csv_path in Path(farfield_dir).glob("grating_farfield_*nm_x_profile.csv"):
        wavelength_nm = _parse_farfield_profile_wavelength_nm(csv_path)
        if wavelength_nm is not None:
            candidates.append((abs(wavelength_nm - target_wavelength_nm), wavelength_nm, csv_path))

    if not candidates:
        return None, None

    _, wavelength_nm, csv_path = min(candidates, key=lambda item: item[0])
    return csv_path, wavelength_nm


def plot_grating_farfield_steering_map(
    grating_dir,
    wavelength_nm,
    glass_refractive_index=1.5,
    air_refractive_index=1.0,
    output_file=None,
    save_matrix=True,
):
    """
    Stack each steering subfolder's nearest wavelength x-profile into a 2D map.

    The CSV theta_x_deg values are treated as glass-side angles. The plotted
    x-axis is converted to the corresponding air angle by Snell's law so it can
    be compared directly with the designed beam steering angle.
    """
    grating_dir = Path(grating_dir)
    if not grating_dir.exists():
        raise FileNotFoundError(f"Grating directory not found: {grating_dir}")

    rows = []
    for subdir in sorted(item for item in grating_dir.iterdir() if item.is_dir()):
        steering_angle_deg = _parse_steering_angle_deg(subdir.name)
        if steering_angle_deg is None:
            continue

        farfield_dir = subdir / "farfield_data"
        if not farfield_dir.is_dir():
            continue

        csv_path, actual_wavelength_nm = _closest_farfield_x_profile_csv(
            farfield_dir, float(wavelength_nm)
        )
        if csv_path is None:
            continue

        data = np.genfromtxt(csv_path, delimiter=",", names=True)
        theta_x_deg = np.asarray(data["theta_x_deg"], dtype=float).ravel()
        intensity = np.asarray(data["intensity"], dtype=float).ravel()
        if theta_x_deg.size != intensity.size:
            raise ValueError(f"Malformed x-profile CSV: {csv_path}")

        rows.append(
            {
                "steering_angle_deg": steering_angle_deg,
                "actual_wavelength_nm": float(actual_wavelength_nm),
                "csv_path": csv_path,
                "theta_x_deg": theta_x_deg,
                "intensity": intensity,
            }
        )

    if not rows:
        raise FileNotFoundError(
            f"No grating_farfield_*nm_x_profile.csv files found under {grating_dir}"
        )

    rows.sort(key=lambda item: item["steering_angle_deg"])
    theta_x_glass_deg = rows[0]["theta_x_deg"]
    for row in rows[1:]:
        if row["theta_x_deg"].shape != theta_x_glass_deg.shape or not np.allclose(
            row["theta_x_deg"], theta_x_glass_deg
        ):
            raise ValueError(
                "Cannot concatenate x-profiles because theta_x_deg axes differ. "
                f"First mismatch: {row['csv_path']}"
            )

    steering_angles_deg = np.array([row["steering_angle_deg"] for row in rows])
    actual_wavelengths_nm = np.array([row["actual_wavelength_nm"] for row in rows])
    intensity_map = np.vstack([row["intensity"] for row in rows])
    sin_theta_air = (
        float(glass_refractive_index)
        / float(air_refractive_index)
        * np.sin(np.deg2rad(theta_x_glass_deg))
    )
    air_valid_mask = np.abs(sin_theta_air) <= 1.0
    theta_x_air_deg = np.rad2deg(np.arcsin(sin_theta_air[air_valid_mask]))
    intensity_map_air = intensity_map[:, air_valid_mask]

    if output_file is None:
        wl_tag = wavelength_tag(wavelength_nm, prefix="", decimals=1)
        output_file = grating_dir / f"farfield_steering_map_{wl_tag}.png"
    output_file = Path(output_file)
    ensure_dir(output_file.parent)

    plt.figure(figsize=(8, 5))
    plt.pcolormesh(theta_x_air_deg, steering_angles_deg, intensity_map_air, shading="auto")
    plt.xlabel("theta_x_air_deg")
    plt.ylabel("Designed beam steering angle (deg)")
    plt.title(
        f"Far-field X profiles near {float(wavelength_nm):.1f} nm "
        f"(actual {np.mean(actual_wavelengths_nm):.1f} nm, n_glass={float(glass_refractive_index):g})"
    )
    plt.colorbar(label="Integrated intensity")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    if save_matrix:
        np.save(output_file.with_suffix(".npy"), intensity_map_air)
        np.save(output_file.with_name(output_file.stem + "_theta_x_air_deg.npy"), theta_x_air_deg)
        np.save(output_file.with_name(output_file.stem + "_theta_x_glass_deg.npy"), theta_x_glass_deg)
        np.save(
            output_file.with_name(output_file.stem + "_steering_angles_deg.npy"),
            steering_angles_deg,
        )

    print(f"Saved far-field steering map to: {output_file}")
    return {
        "output_file": str(output_file),
        "theta_x_air_deg": theta_x_air_deg,
        "theta_x_glass_deg": theta_x_glass_deg,
        "steering_angles_deg": steering_angles_deg,
        "actual_wavelengths_nm": actual_wavelengths_nm,
        "intensity_map": intensity_map_air,
        "glass_refractive_index": float(glass_refractive_index),
        "air_refractive_index": float(air_refractive_index),
        "source_files": [str(row["csv_path"]) for row in rows],
    }


def _monitor_wavelength_index(fdtd, monitor_name, target_lambda_m):
    """Find the monitor frequency index nearest to target wavelength."""

    f_vec = np.array(fdtd.getdata(monitor_name, "f")).squeeze()
    
    if f_vec.size == 0:
        raise ValueError(f"No frequency data found in monitor: {monitor_name}")

    lambda_vec = C0 / f_vec
    f_index_py = int(np.argmin(np.abs(lambda_vec - target_lambda_m)))
    actual_lambda_m = float(lambda_vec[f_index_py])
    return f_index_py, f_index_py + 1, actual_lambda_m


def _farfield_angle_data(fdtd, monitor_name, f_index_lsf):
    """Read far-field intensity and convert ux/uy axes to angles in degrees."""
    e2 = np.squeeze(np.array(fdtd.farfield3d(monitor_name, f_index_lsf)))
    ux = np.squeeze(np.array(fdtd.farfieldux(monitor_name, f_index_lsf)))
    uy = np.squeeze(np.array(fdtd.farfielduy(monitor_name, f_index_lsf)))
    theta_x_deg = np.degrees(np.arcsin(np.clip(ux, -1.0, 1.0)))
    theta_y_deg = np.degrees(np.arcsin(np.clip(uy, -1.0, 1.0)))
    return e2, theta_x_deg, theta_y_deg


def farfield_x_profile(e2, theta_x_deg, theta_y_deg, reduction="sum"):
    """Collapse a 2D far-field map into Angle X versus intensity."""
    if e2.shape == (theta_x_deg.size, theta_y_deg.size):
        # Rows are theta_x, columns are theta_y.
        axis_y = 1
    elif e2.shape == (theta_y_deg.size, theta_x_deg.size):
        # Rows are theta_y, columns are theta_x.
        axis_y = 0
    else:
        raise ValueError(
            f"Unexpected far-field shape {e2.shape}; "
            f"theta_x has {theta_x_deg.size}, theta_y has {theta_y_deg.size}."
        )

    if reduction == "sum":
        intensity_x = np.nansum(e2, axis=axis_y)
    elif reduction == "max":
        intensity_x = np.nanmax(e2, axis=axis_y)
    else:
        raise ValueError(f"Unknown far-field profile reduction: {reduction}")

    return np.asarray(theta_x_deg), np.asarray(intensity_x)


def save_farfield_x_profile(e2, theta_x_deg, theta_y_deg, output_dir, prefix):
    """Save a 1D Angle X profile plot and arrays for locating steering maxima."""
    output_dir = ensure_dir(output_dir)
    theta_x_deg, intensity_x = farfield_x_profile(e2, theta_x_deg, theta_y_deg)

    np.save(output_dir / f"{prefix}_theta_x_profile_deg.npy", theta_x_deg)
    np.save(output_dir / f"{prefix}_intensity_x_profile.npy", intensity_x)
    np.savetxt(
        output_dir / f"{prefix}_x_profile.csv",
        np.column_stack([theta_x_deg, intensity_x]),
        delimiter=",",
        header="theta_x_deg,intensity",
        comments="",
    )

    peak_idx = int(np.nanargmax(intensity_x))
    peak_info = {
        "peak_theta_x_deg": float(theta_x_deg[peak_idx]),
        "peak_intensity_x": float(intensity_x[peak_idx]),
        "profile_reduction": "sum_over_theta_y",
    }
    save_json(output_dir / f"{prefix}_x_profile_peak.json", peak_info)

    plot_mask = (theta_x_deg >= -25.0) & (theta_x_deg <= 25.0)
    plot_theta_x = theta_x_deg[plot_mask]
    plot_intensity_x = intensity_x[plot_mask]

    plt.figure(figsize=(7, 4))
    plt.plot(plot_theta_x, plot_intensity_x, linewidth=1.5)
    plt.axvline(theta_x_deg[peak_idx], linestyle="--", linewidth=1.0)
    plt.xlabel("Angle X (deg)")
    plt.ylabel("Integrated intensity")
    plt.title(f"Far-field X profile, peak at {theta_x_deg[peak_idx]:.2f} deg")
    plt.xlim(-25, 25)
    plt.xticks(np.arange(-25, 25.1, 2.5))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / f"{prefix}_x_profile.png", dpi=300)
    plt.close()

    return peak_info


def save_farfield_data_single_wavelength(
    fdtd,
    monitor_name,
    target_lambda_m,
    output_dir,
    prefix="farfield",
    make_plot=True,
    save_angle_axes=False,
):
    """Save far-field E2 data and metadata for the monitor wavelength nearest target."""
    output_dir = ensure_dir(output_dir)
    f_index_py, f_index_lsf, actual_lambda_m = _monitor_wavelength_index(
        fdtd, monitor_name, target_lambda_m
    )
    e2, theta_x_deg, theta_y_deg = _farfield_angle_data(fdtd, monitor_name, f_index_lsf)

    print(f"Target wavelength: {target_lambda_m*1e9:.3f} nm")
    print(f"Actual wavelength used: {actual_lambda_m*1e9:.3f} nm")
    print(f"Frequency index used in Lumerical: {f_index_lsf}")

    wl_tag = wavelength_tag(actual_lambda_m * 1e9, prefix="", decimals=1)
    np.save(output_dir / f"{prefix}_E2_{wl_tag}.npy", e2)
    x_profile_peak = save_farfield_x_profile(
        e2=e2,
        theta_x_deg=theta_x_deg,
        theta_y_deg=theta_y_deg,
        output_dir=output_dir,
        prefix=f"{prefix}_{wl_tag}",
    )

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
        "E2_shape": list(e2.shape),
        "theta_x_range_deg": [float(theta_x_deg.min()), float(theta_x_deg.max())],
        "theta_y_range_deg": [float(theta_y_deg.min()), float(theta_y_deg.max())],
        "theta_x_size": int(theta_x_deg.size),
        "theta_y_size": int(theta_y_deg.size),
        "x_profile_peak": x_profile_peak,
    }
    save_json(output_dir / f"{prefix}_meta_{wl_tag}.json", meta)

    if make_plot:
        plt.figure(figsize=(6, 5))
        theta_x_grid, theta_y_grid = np.meshgrid(theta_x_deg, theta_y_deg)
        if e2.shape == (theta_x_deg.size, theta_y_deg.size):
            plot_e2 = e2.T
        elif e2.shape == (theta_y_deg.size, theta_x_deg.size):
            plot_e2 = e2
        else:
            raise ValueError(
                f"Unexpected far-field shape {e2.shape}; "
                f"theta_x has {theta_x_deg.size}, theta_y has {theta_y_deg.size}."
            )
        plt.pcolormesh(theta_x_grid, theta_y_grid, plot_e2, shading="auto")
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
        "E2_shape": list(e2.shape),
    }


def save_farfield_data_multiple_wavelengths(
    fdtd,
    monitor_name,
    target_wavelengths_m,
    output_dir,
    prefix="farfield",
    make_plot=True,
):
    """Save far-field data for several wavelengths without duplicate frequency indices."""
    output_dir = ensure_dir(output_dir)
    results = []
    used_f_indices = set()
    saved_angle_axes = False

    for target_lambda_m in target_wavelengths_m:
        print("\n" + "=" * 60)
        print(f"Processing target wavelength: {target_lambda_m * 1e9:.3f} nm")

        f_index_py, _, actual_lambda_m = _monitor_wavelength_index(
            fdtd, monitor_name, target_lambda_m
        )
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

    save_json(output_dir / f"{prefix}_summary.json", results)
    print("\nAll far-field data saved.")
    return results


def find_farfield_peak(fdtd, monitor_name, target_lambda_m):
    """Return the strongest far-field pixel and its angular coordinates."""
    _, f_index_lsf, actual_lambda_m = _monitor_wavelength_index(fdtd, monitor_name, target_lambda_m)
    e2, theta_x_deg, theta_y_deg = _farfield_angle_data(fdtd, monitor_name, f_index_lsf)

    peak_index = np.unravel_index(int(np.nanargmax(e2)), e2.shape)
    if e2.shape == (theta_x_deg.size, theta_y_deg.size):
        peak_theta_x = float(theta_x_deg[peak_index[0]])
        peak_theta_y = float(theta_y_deg[peak_index[1]])
    elif e2.shape == (theta_y_deg.size, theta_x_deg.size):
        peak_theta_x = float(theta_x_deg[peak_index[1]])
        peak_theta_y = float(theta_y_deg[peak_index[0]])
    else:
        raise ValueError(
            f"Unexpected far-field shape {e2.shape}; "
            f"theta_x has {theta_x_deg.size}, theta_y has {theta_y_deg.size}."
        )

    return {
        "target_lambda_m": float(target_lambda_m),
        "actual_lambda_m": actual_lambda_m,
        "actual_lambda_nm": actual_lambda_m * 1e9,
        "f_index_lumerical": int(f_index_lsf),
        "peak_theta_x_deg": peak_theta_x,
        "peak_theta_y_deg": peak_theta_y,
        "peak_intensity": float(np.nanmax(e2)),
        "total_intensity": float(np.nansum(e2)),
        "peak_fraction_of_total": float(np.nanmax(e2) / max(np.nansum(e2), 1e-30)),
        "e2_shape": list(e2.shape),
    }


def add_angle_accuracy(peak_result, target_theta_deg):
    """Annotate a far-field peak with same-sign and sign-flipped angle errors."""
    peak_theta_x = peak_result["peak_theta_x_deg"]
    same_sign_error = abs(peak_theta_x - target_theta_deg)
    opposite_sign_error = abs(peak_theta_x + target_theta_deg)

    peak_result["target_theta_deg"] = float(target_theta_deg)
    peak_result["same_sign_error_deg"] = float(same_sign_error)
    peak_result["opposite_sign_error_deg"] = float(opposite_sign_error)
    peak_result["best_abs_error_deg"] = float(min(same_sign_error, opposite_sign_error))
    peak_result["best_sign_match"] = "same" if same_sign_error <= opposite_sign_error else "opposite"
    return peak_result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot a stacked far-field x-profile map across grating steering folders."
    )
    parser.add_argument(
        "grating_dir",
        nargs="?",
        default=(
            "output_single_grating/grating_datas/"
            "TiO2_on_Top_p350nm_r135nm_wavelength620nm_da79e378a7"
        ),
        help="Directory containing steering_* subfolders.",
    )
    parser.add_argument(
        "wavelength_nm",
        nargs="?",
        type=float,
        default=620.0,
        help="Target wavelength in nm.",
    )
    parser.add_argument(
        "glass_refractive_index",
        nargs="?",
        type=float,
        default=1.5,
        help="Refractive index of glass used to convert theta_x_deg to air angle.",
    )
    parser.add_argument(
        "--output-file",
        default=None,
        help="Optional output PNG path. Defaults inside grating_dir.",
    )
    args = parser.parse_args()

    plot_grating_farfield_steering_map(
        grating_dir=args.grating_dir,
        wavelength_nm=args.wavelength_nm,
        glass_refractive_index=args.glass_refractive_index,
        output_file=args.output_file,
    )
