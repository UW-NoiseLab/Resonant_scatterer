template = """
############################################
# Unit-cell metasurface model for Lumerical
# Geometry + FDTD + monitors
############################################
function build_unit_cell_model( 
    period,
    t_Al,
    t_spacer,
    t_LC,
    t_ITO,
    t_glass,
    r_scatter,
    t_scatter,
    lambda_start,
    lambda_stop
)
{
############################################
# Hyperparameters
############################################
deleteall;

# Substrate
substrate_z_min = -4.0e-6;
substrate_z_max = -0.15e-6;

# Al layer
Al_z_center = 0;
Al_z_min = Al_z_center - t_Al/2;
Al_z_max = Al_z_center + t_Al/2;

# Spacer (SiO2), directly above Al
spacer_z_min = Al_z_max;
spacer_z_max = spacer_z_min + t_spacer;
spacer_z_center = 0.5*(spacer_z_min + spacer_z_max);

# Liquid crystal, directly above spacer
LC_z_min = spacer_z_max;
LC_z_max = LC_z_min + t_LC;
LC_z_center = 0.5*(LC_z_min + LC_z_max);

# ITO, directly above LC
ITO_z_min = LC_z_max;
ITO_z_max = ITO_z_min + t_ITO;
ITO_z_center = 0.5*(ITO_z_min + ITO_z_max);

# Glass layer above ITO
glass_z_min = ITO_z_max;
glass_z_max = glass_z_min + t_glass;
glass_z_center = 0.5*(glass_z_min + glass_z_max);

# __SCATTER_LABEL__ cylindrical scatterer, __SCATTER_POSITION_COMMENT__
__SCATTER_Z_POSITION__
scatter_z_center = 0.5*(scatter_z_min + scatter_z_max);

# FDTD margins
z_margin_bottom = 0.85e-6;
z_margin_top = 0.98e-6;

fdtd_z_min = Al_z_min - 0.3e-6;
fdtd_z_max = glass_z_min + 0.8e-6;

############################################
# Materials note
# Please change material database names if needed
############################################

mat_substrate = "SiO2 (Glass) - Palik";
mat_Al        = "Al (Aluminium) - Palik";
mat_spacer    = "SiO2 (Glass) - Palik";
mat_glass     = "SiO2 (Glass) - Palik";
mat_scatter   = "__SCATTER_MATERIAL__";

############################################
# Geometry
############################################

# Substrate
addrect;
set("name", "substrate");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z min", substrate_z_min);
set("z max", substrate_z_max);
set("material", mat_substrate);

# Al layer
addrect;
set("name", "Al_layer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_Al);
set("z", Al_z_center);
set("material", mat_Al);

# SiO2 spacer
addrect;
set("name", "spacer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_spacer);
set("z", spacer_z_center);
set("material", mat_spacer);

# Liquid crystal block (object-defined dielectric, n = 1.75)
addrect;
set("name", "liquid_crystal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_LC);
set("z", LC_z_center);
set("material", "<Object defined dielectric>");
set("index", 1.75);

# ITO layer
addrect;
set("name", "ITO_layer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_ITO);
set("z", ITO_z_center);
set("material", "ITO");

# Glass superstrate
addrect;
set("name", "glass_top");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z min", glass_z_min);
set("z max", glass_z_max);
set("material", mat_glass);

# Cylindrical __SCATTER_LABEL__ scatterer, __SCATTER_POSITION_COMMENT__
addcircle;
set("name", "scatterer");
set("x", 0);
set("y", 0);
set("radius", r_scatter);
set("z min", scatter_z_min);
set("z max", scatter_z_max);
set("material", mat_scatter);

############################################
# FDTD region
############################################

addfdtd;
set("dimension", "3D");

set("x", 0);
set("y", 0);
set("z", 0.5*(fdtd_z_min + fdtd_z_max));

set("x span", period);
set("y span", period);
set("z span", fdtd_z_max - fdtd_z_min);

set("allow symmetry on all boundaries", 1);
# Boundary conditions for unit cell
set("x min bc", "Anti-Symmetric");
set("x max bc", "Anti-Symmetric");
set("y min bc", "Symmetric");
set("y max bc", "Symmetric");
set("z min bc", "PML");
set("z max bc", "PML");

# Simulation control
set("mesh accuracy", 3);
set("simulation time", 1000e-15);

############################################
# Plane wave source
# Incidence from top to bottom
############################################

addplane;
set("name", "source");
set("injection axis", "z");
set("direction", "Backward");    # from +z toward -z
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", glass_z_min + 0.28e-6);
set("wavelength start", lambda_start);
set("wavelength stop", lambda_stop);

############################################
# Power monitors for R/T
############################################

# Reflection monitor (above structure, below source)
#adddftmonitor;
#set("name", "R_monitor");
#set("monitor type", "2D Z-normal");
#set("x", 0);
#set("y", 0);
#set("x span", period);
#set("y span", period);
#set("z", glass_z_min + 0.10e-6);

# Transmission monitor (below structure)
adddftmonitor;
set("name", "T_monitor");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", substrate_z_max - 0.10e-6);

############################################
# Index monitors
############################################

# Horizontal cross-section index monitor (XY plane)
addindex;
set("name", "index_xy");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", scatter_z_center);

# Vertical cross-section index monitor (XZ plane)
addindex;
set("name", "index_xz");
set("monitor type", "2D Y-normal");
set("x", 0);
set("y", 0);
set("z", 0.5*(fdtd_z_min + fdtd_z_max));
set("x span", period);
set("z span", fdtd_z_max - fdtd_z_min);

############################################
# Field monitors (DFT monitors for E-field)
############################################

# Top monitor
adddftmonitor;
set("name", "monitor");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", ITO_z_max + 0.48e-6);
#set("override global monitor settings", 1);
#set("use source limits", 1);


# Side View monitor
adddftmonitor;
set("name", "x_normal_monitor");
set("monitor type", "2D x-normal");

set("y", 0);
set("x", 0);
set("y span", period);
set("z min", Al_z_min - 0.15e-6);
set("z max", glass_z_min + 0.18e-6);

# Horizontal cross-section field monitor (same position as index_xy)
#adddftmonitor;
#set("name", "field_xy");
#set("monitor type", "2D Z-normal");
#set("x", 0);
#set("y", 0);
#set("x span", period);
#set("y span", period);
#set("z", scatter_z_center);
#set("override global monitor settings", 1);
#set("use source limits", 1);

# Vertical cross-section field monitor (same position as index_xz)
#adddftmonitor;
#set("name", "field_xz");
#set("monitor type", "2D Y-normal");
#set("x", 0);
#set("y", 0);
#set("z", 0.5*(fdtd_z_min + fdtd_z_max));
#set("x span", period);
#set("z span", fdtd_z_max - fdtd_z_min);
#set("override global monitor settings", 1);
#set("use source limits", 1);

############################################
# Optional mesh override around scatterer / thin layers
############################################

#addmesh;
#set("name", "mesh_fine");
#set("x", 0);
#set("y", 0);
#set("z", 0.5*(scatter_z_min + ITO_z_max));
#set("x span", period);
#set("y span", period);
#set("z span", (ITO_z_max - scatter_z_min) +Sy 0.10e-6);
#set("dx", 5e-9);
#set("dy", 5e-9);
#set("dz", 5e-9);
}
############################################
# Save
############################################

# Uncomment if you want to save automatically:
# save("unit_cell_metasurface.fsp");

############################################
# Notes
############################################
# 1. If material names do not match your local database,
#    replace the strings above with your exact material names.
# 2. The scatterer scheme is __SCATTER_SCHEME__.
# 3. The source is normal incidence from top.
# 4. Periodic BC is used in x/y for single unit-cell simulation.
# 5. You can directly run the script, then run the simulation.
############################################
"""

from typing import Literal, Union

_SCATTERER_SCHEMES = {
    "SiN on Top": {
        "label": "Si3N4",
        "material": "Si3N4 (Silicon Nitride) - Phillip",
        "position_comment": "directly below ITO",
        "z_position": """scatter_z_max = ITO_z_min;
scatter_z_min = scatter_z_max - t_scatter;""",
    },
    "SiN on Bottom": {
        "label": "Si3N4",
        "material": "Si3N4 (Silicon Nitride) - Phillip",
        "position_comment": "directly above spacer SiO2",
        "z_position": """scatter_z_min = spacer_z_max;
scatter_z_max = scatter_z_min + t_scatter;""",
    },
    "TiO2 on Top": {
        "label": "TiO2",
        "material": "TiO2 (Titanium Dioxide) - Siefke",
        "position_comment": "directly below ITO",
        "z_position": """scatter_z_max = ITO_z_min;
scatter_z_min = scatter_z_max - t_scatter;""",
    },
    "TiO2 on Bottom": {
        "label": "TiO2",
        "material": "TiO2 (Titanium Dioxide) - Siefke",
        "position_comment": "directly above spacer SiO2",
        "z_position": """scatter_z_min = spacer_z_max;
scatter_z_max = scatter_z_min + t_scatter;""",
    },
}


_DBR_MATERIALS = {
    "SiN": {
        "label": "Si3N4",
        "material": "Si3N4 (Silicon Nitride) - Phillip",
        "index": 2.0,
    },
    "TiO2": {
        "label": "TiO2",
        "material": "TiO2 (Titanium Dioxide) - Siefke",
        "index": 2.4,
    },
}


_DBR_TEMPLATE = """
############################################
# Unit-cell LC + DBR model for Lumerical
# Geometry + FDTD + monitors
############################################
function build_unit_cell_model(
    period,
    t_Al,
    t_spacer,
    t_LC,
    t_ITO,
    t_glass,
    r_scatter,
    t_scatter,
    lambda_start,
    lambda_stop
)
{
############################################
# Hyperparameters
############################################
deleteall;

# Substrate
substrate_z_min = -4.0e-6;
substrate_z_max = -0.15e-6;

# Al layer
Al_z_center = 0;
Al_z_min = Al_z_center - t_Al/2;
Al_z_max = Al_z_center + t_Al/2;

# Spacer (SiO2), directly above Al
spacer_z_min = Al_z_max;
spacer_z_max = spacer_z_min + t_spacer;
spacer_z_center = 0.5*(spacer_z_min + spacer_z_max);

# Liquid crystal, directly above spacer
LC_z_min = spacer_z_max;
LC_z_max = LC_z_min + t_LC;
LC_z_center = 0.5*(LC_z_min + LC_z_max);

# ITO, directly above LC
ITO_z_min = LC_z_max;
ITO_z_max = ITO_z_min + t_ITO;
ITO_z_center = 0.5*(ITO_z_min + ITO_z_max);

# Quarter-wave DBR above ITO.
# High-index layer is __DBR_LABEL__; low-index layer is SiO2.
dbr_wavelength = __DBR_WAVELENGTH__;
dbr_pairs = __DBR_PAIRS__;
n_dbr_high = __DBR_INDEX__;
n_dbr_low = 1.45;
dbr_high_thickness = dbr_wavelength/(4*n_dbr_high);
dbr_low_thickness = dbr_wavelength/(4*n_dbr_low);
dbr_z_min = ITO_z_max;
dbr_z_max = dbr_z_min + dbr_pairs*(dbr_high_thickness + dbr_low_thickness);
dbr_z_center = 0.5*(dbr_z_min + dbr_z_max);

# Glass layer above DBR
glass_z_min = dbr_z_max;
glass_z_max = glass_z_min + t_glass;
glass_z_center = 0.5*(glass_z_min + glass_z_max);

# FDTD margins
z_margin_bottom = 0.85e-6;
z_margin_top = 0.98e-6;

fdtd_z_min = Al_z_min - 0.3e-6;
fdtd_z_max = glass_z_min + 0.8e-6;

############################################
# Materials note
# Please change material database names if needed
############################################

mat_substrate = "SiO2 (Glass) - Palik";
mat_Al        = "Al (Aluminium) - Palik";
mat_spacer    = "SiO2 (Glass) - Palik";
mat_glass     = "SiO2 (Glass) - Palik";
mat_dbr_high  = "__DBR_MATERIAL__";
mat_dbr_low   = "SiO2 (Glass) - Palik";

############################################
# Geometry
############################################

# Substrate
addrect;
set("name", "substrate");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z min", substrate_z_min);
set("z max", substrate_z_max);
set("material", mat_substrate);

# Al layer
addrect;
set("name", "Al_layer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_Al);
set("z", Al_z_center);
set("material", mat_Al);

# SiO2 spacer
addrect;
set("name", "spacer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_spacer);
set("z", spacer_z_center);
set("material", mat_spacer);

# Liquid crystal block (object-defined dielectric, n = 1.75)
addrect;
set("name", "liquid_crystal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_LC);
set("z", LC_z_center);
set("material", "<Object defined dielectric>");
set("index", 1.75);

# ITO layer
addrect;
set("name", "ITO_layer");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z span", t_ITO);
set("z", ITO_z_center);
set("material", "ITO");

# DBR stack above ITO: (__DBR_LABEL__/SiO2)^dbr_pairs
current_z = dbr_z_min;
for (i=1:dbr_pairs) {
    addrect;
    set("name", "__DBR_LABEL___DBR_high_"+num2str(i));
    set("x", 0);
    set("y", 0);
    set("x span", period);
    set("y span", period);
    set("z min", current_z);
    set("z max", current_z + dbr_high_thickness);
    set("material", mat_dbr_high);
    current_z = current_z + dbr_high_thickness;

    addrect;
    set("name", "SiO2_DBR_low_"+num2str(i));
    set("x", 0);
    set("y", 0);
    set("x span", period);
    set("y span", period);
    set("z min", current_z);
    set("z max", current_z + dbr_low_thickness);
    set("material", mat_dbr_low);
    current_z = current_z + dbr_low_thickness;
}

# Glass superstrate
addrect;
set("name", "glass_top");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z min", glass_z_min);
set("z max", glass_z_max);
set("material", mat_glass);

############################################
# FDTD region
############################################

addfdtd;
set("dimension", "3D");

set("x", 0);
set("y", 0);
set("z", 0.5*(fdtd_z_min + fdtd_z_max));

set("x span", period);
set("y span", period);
set("z span", fdtd_z_max - fdtd_z_min);

set("allow symmetry on all boundaries", 1);
set("x min bc", "Anti-Symmetric");
set("x max bc", "Anti-Symmetric");
set("y min bc", "Symmetric");
set("y max bc", "Symmetric");
set("z min bc", "PML");
set("z max bc", "PML");

set("mesh accuracy", 3);
set("simulation time", 1000e-15);

############################################
# Plane wave source
# Incidence from top to bottom
############################################

addplane;
set("name", "source");
set("injection axis", "z");
set("direction", "Backward");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", glass_z_min + 0.28e-6);
set("wavelength start", lambda_start);
set("wavelength stop", lambda_stop);

############################################
# Power monitors for R/T
############################################

adddftmonitor;
set("name", "T_monitor");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", substrate_z_max - 0.10e-6);

############################################
# Index monitors
############################################

addindex;
set("name", "index_xy");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", dbr_z_center);

addindex;
set("name", "index_xz");
set("monitor type", "2D Y-normal");
set("x", 0);
set("y", 0);
set("z", 0.5*(fdtd_z_min + fdtd_z_max));
set("x span", period);
set("z span", fdtd_z_max - fdtd_z_min);

############################################
# Field monitors (DFT monitors for E-field)
############################################

adddftmonitor;
set("name", "monitor");
set("monitor type", "2D Z-normal");
set("x", 0);
set("y", 0);
set("x span", period);
set("y span", period);
set("z", dbr_z_max + 0.48e-6);

adddftmonitor;
set("name", "x_normal_monitor");
set("monitor type", "2D x-normal");
set("y", 0);
set("x", 0);
set("y span", period);
set("z min", Al_z_min - 0.15e-6);
set("z max", glass_z_min + 0.18e-6);
}

############################################
# Notes
############################################
# 1. This DBR model intentionally contains no cylindrical scatterer in the LC.
# 2. The DBR material is __DBR_LABEL__; the DBR is designed for __DBR_WAVELENGTH__ m.
# 3. The source is normal incidence from top.
# 4. Periodic BC is used in x/y for single unit-cell simulation.
############################################
"""


def scatterer_lsf_gen(
    scheme: Union[Literal["SiN on Top", "SiN on Bottom", "TiO2 on Top", "TiO2 on Bottom"], None] = None,
):
    if scheme is None:
        scheme = "SiN on Top"

    if scheme not in _SCATTERER_SCHEMES:
        valid_schemes = ", ".join(_SCATTERER_SCHEMES)
        raise ValueError(f"Unknown scatterer scheme '{scheme}'. Valid schemes: {valid_schemes}")

    scheme_config = _SCATTERER_SCHEMES[scheme]

    return (
        template
        .replace("__SCATTER_SCHEME__", scheme)
        .replace("__SCATTER_LABEL__", scheme_config["label"])
        .replace("__SCATTER_MATERIAL__", scheme_config["material"])
        .replace("__SCATTER_POSITION_COMMENT__", scheme_config["position_comment"])
        .replace("__SCATTER_Z_POSITION__", scheme_config["z_position"])
    )
    

def scatter_DBR_lsf_gen(
    wavelength: float,
    material: Literal["TiO2", "SiN"] = "TiO2",
    dbr_pairs: int = 6,
):
    """Generate an LC unit-cell script with a quarter-wave DBR above ITO.

    Parameters
    ----------
    wavelength:
        DBR design wavelength in meters.
    material:
        High-index DBR material. Must be "TiO2" or "SiN".
    """
    if material not in _DBR_MATERIALS:
        valid_materials = ", ".join(_DBR_MATERIALS)
        raise ValueError(f"Unknown DBR material '{material}'. Valid materials: {valid_materials}")

    wavelength = float(wavelength)
    if wavelength <= 0:
        raise ValueError(f"DBR wavelength must be positive. Got: {wavelength}")
    dbr_pairs = int(dbr_pairs)
    if dbr_pairs <= 0:
        raise ValueError(f"DBR pair count must be positive. Got: {dbr_pairs}")

    dbr_config = _DBR_MATERIALS[material]

    return (
        _DBR_TEMPLATE
        .replace("__DBR_LABEL__", dbr_config["label"])
        .replace("__DBR_MATERIAL__", dbr_config["material"])
        .replace("__DBR_INDEX__", f'{dbr_config["index"]:.16g}')
        .replace("__DBR_WAVELENGTH__", f"{wavelength:.16g}")
        .replace("__DBR_PAIRS__", str(dbr_pairs))
    )


def scatterer_DBR_lsf_gen(
    wavelength: float,
    material: Literal["TiO2", "SiN"] = "TiO2",
    dbr_pairs: int = 6,
):
    return scatter_DBR_lsf_gen(wavelength, material, dbr_pairs)
