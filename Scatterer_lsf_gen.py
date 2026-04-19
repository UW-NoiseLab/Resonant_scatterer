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
