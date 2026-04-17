import numpy as np
from typing import Literal, Union

_SCATTERER_SCHEMES = {
    "SiN on Top": {
        "label": "SiN",
        "material": "Si3N4 (Silicon Nitride) - Phillip",
        "patch_name": "SiN_patch_",
        "group_name": "SiN_patches",
        "position_comment": "directly below ITO",
        "z_position": """scatter_z_max = ITO_z_min;
scatter_z_min = scatter_z_max - __T_SCATTER__;""",
    },
    "SiN on Bottom": {
        "label": "SiN",
        "material": "Si3N4 (Silicon Nitride) - Phillip",
        "patch_name": "SiN_patch_",
        "group_name": "SiN_patches",
        "position_comment": "directly above spacer SiO2",
        "z_position": """scatter_z_min = spacer_z_max;
scatter_z_max = scatter_z_min + __T_SCATTER__;""",
    },
    "TiO2 on Top": {
        "label": "TiO2",
        "material": "TiO2 (Titanium Dioxide) - Siefke",
        "patch_name": "TiO2_patch_",
        "group_name": "TiO2_patches",
        "position_comment": "directly below ITO",
        "z_position": """scatter_z_max = ITO_z_min;
scatter_z_min = scatter_z_max - __T_SCATTER__;""",
    },
    "TiO2 on Bottom": {
        "label": "TiO2",
        "material": "TiO2 (Titanium Dioxide) - Siefke",
        "patch_name": "TiO2_patch_",
        "group_name": "TiO2_patches",
        "position_comment": "directly above spacer SiO2",
        "z_position": """scatter_z_min = spacer_z_max;
scatter_z_max = scatter_z_min + __T_SCATTER__;""",
    },
}

def grating_lsf_gen(
    period,
    t_Al,
    t_spacer,
    t_LC,
    t_ITO,
    t_glass,
    r_scatter,
    t_scatter,
    lambda_design,
    lambda_span,
    
    LC_index_data,
    LC_phase_data,
    scheme: Union[Literal["SiN on Top", "SiN on Bottom", "TiO2 on Top", "TiO2 on Bottom"], None] = None,
):
    '''
    Generate LSF content for the grating structure based on the provided parameters and lookup tables.
    The unit of all length parameters is meter, and the unit of wavelength parameters is meter as well.
    LC_index_data and LC_phase_data should be lists or 1D numpy arrays of the same length, 
    representing the lookup table of LC refractive index vs phase delay at the design wavelength.
    '''
    template = """############################################
# PARAMETERS
############################################
function createmodel(M, theta){
    
#theta is in degree

theta=theta/180.0*pi;

M = M;
N = 1;

period = __PERIOD__;

wavelength=__LAMBDA_DESIGN__;
wl_span=__LAMBDA_SPAN__;

############################################
# DOMAIN SIZE
############################################

Lx = M*period;
Ly = N*period;

############################################
# Hyperparameters
############################################
switchtolayout;
deleteall;

# Substrate
substrate_z_min = -4.0e-6;
substrate_z_max = -0.15e-6;

# Al layer
Al_z_center = 0;
Al_z_min = Al_z_center - 0.5*__T_AL__;
Al_z_max = Al_z_center + 0.5*__T_AL__;

# Spacer (SiO2), directly above Al
spacer_z_min = Al_z_max;
spacer_z_max = spacer_z_min + __T_SPACER__;
spacer_z_center = 0.5*(spacer_z_min + spacer_z_max);

# Liquid crystal, directly above spacer
LC_z_min = spacer_z_max;
LC_z_max = LC_z_min + __T_LC__;
LC_z_center = 0.5*(LC_z_min + LC_z_max);

# ITO, directly above LC
ITO_z_min = LC_z_max;
ITO_z_max = ITO_z_min + __T_ITO__;
ITO_z_center = 0.5*(ITO_z_min + ITO_z_max);

# Glass layer above ITO
glass_z_min = ITO_z_max;
glass_z_max = glass_z_min + __T_GLASS__;
glass_z_center = 0.5*(glass_z_min + glass_z_max);

# __SCATTER_LABEL__ cylindrical scatterer, __SCATTER_POSITION_COMMENT__
__SCATTER_Z_POSITION__
scatter_z_center = 0.5*(scatter_z_min + scatter_z_max);

# FDTD margins
z_margin_bottom = 0.3e-6;
z_margin_top = 0.3e-6;

fdtd_z_min = Al_z_min - 0.3e-6;
fdtd_z_max = glass_z_min + 2.0e-6;

############################################
# SUBSTRATE
############################################

patch_size = __PERIOD__;  # This is the pixelated Al patch size, not the actual Al patch size. We can set it to be the same as period to make sure the Al patch is a solid block without holes. The actual Al patch size is determined by the LC cell size, which is also set to be the same as period.


addrect;
set("name","substrate");
set("material","SiO2 (Glass) - Palik");

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z min",substrate_z_min);
set("z max",substrate_z_max);


############################################
# ITO
############################################

addrect;
set("name","ITO");
set("material","ITO");

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z min",ITO_z_min);
set("z max",ITO_z_max);

############################################
# TOP GLASS
############################################

addrect;
set("name","glass");
set("material","SiO2 (Glass) - Palik");

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z min",glass_z_min);
set("z max",glass_z_max);


############################################
# LIC n
############################################

LC_index_data = __LC_INDEX_DATA__;

LC_phase_data = __LC_PHASE_DATA__;

# Deciding the phase gradient with theta

k0 = 2*pi/wavelength;

# phase increment between adjacent cells
dphi = k0 * period * sin(theta);

# construct phase distribution
LC_phase = matrix(M);

for(i=1:M){
    LC_phase(i) = mod((i-1)*dphi, 2*pi);
}

LC_index = interp(LC_index_data,LC_phase_data,LC_phase);

?LC_index;


############################################
# UNIT CELLS
############################################

for(i=1:M){

    xpos = -Lx/2 + (i-0.5)*period;
    
    ####################################
    # Liquid crystal cell
    ####################################

    addrect;

    set("name","LC_cell_"+num2str(i));

    set("material","<Object defined dielectric>");
    set("index",LC_index(i));

    set("x",xpos);
    set("y",0);

    set("x span",period);
    set("y span",Ly);

    set("z min",LC_z_min);
    set("z max",LC_z_max);
    

    ####################################
    # Aluminum patch
    ####################################

    addrect;
    set("name","Al_patch_"+num2str(i));

    set("material","Al (Aluminium) - Palik");

    set("x",xpos);
    set("y",0);

    set("x span",patch_size);
    set("y span",patch_size);

    set("z min",Al_z_min);
    set("z max",Al_z_max);
    #addtogroup("Al_patches");

    ####################################
    # SiO2 spacer
    ####################################

    addrect;
    set("name","Oxide_patch_"+num2str(i));

    set("material","SiO2 (Glass) - Palik");

    set("x",xpos);
    set("y",0);

    set("x span",patch_size);
    set("y span",patch_size);

    set("z min",spacer_z_min);
    set("z max",spacer_z_max);
    #addtogroup("Oxide_patches");

    ####################################
    # __SCATTER_LABEL__ metasurface cylinder
    ####################################

    addcircle;
    set("name","__SCATTER_PATCH_NAME__"+num2str(i));

    set("material","__SCATTER_MATERIAL__");

    set("x",xpos);
    set("y",0);

    set("radius",__R_SCATTER__);

    set("z min",scatter_z_min);
    set("z max",scatter_z_max);
    #addtogroup("__SCATTER_GROUP_NAME__");

}

############################################
# SOURCE
############################################

addplane;

set("name","source");

set("injection axis","z");
set("direction","Backward");

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z",glass_z_min + 0.5e-6);

set("wavelength start",wavelength-wl_span/2);
set("wavelength stop",wavelength+wl_span/2);

set("angle theta",0);
set("angle phi",0);

############################################
# REFLECTION MONITOR
############################################

adddftmonitor;

set("name","R_monitor");

set("monitor type","2D Z-normal");

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z", glass_z_min + 1.0e-6);

set("override global monitor settings",1);
set("frequency points",200);

############################################
# FDTD REGION
############################################

addfdtd;

set("x",0);
set("y",0);      # added

set("x span",Lx);
set("y span",Ly);

set("z min",fdtd_z_min);
set("z max",fdtd_z_max);

set("x min bc","Bloch");
set("x max bc","Bloch");

set("y min bc","Periodic");
set("y max bc","Periodic");

set("z min bc","PML");
set("z max bc","PML");

set("set based on source angle",1);


############################################
# MONITORS AT y = 0
############################################

####################################
# DFT field monitor
####################################
adddftmonitor;

set("name","field_monitor_y0");

set("monitor type","2D Y-normal");

set("x",0);
set("y",0);
set("z min",fdtd_z_min);

set("z max",fdtd_z_max);  

set("x span",Lx);

set("override global monitor settings",1);
set("frequency points",200);


####################################
# Index monitor
####################################
addindex;

set("name","index_monitor_y0");

set("monitor type","2D Y-normal");

set("x",0);
set("y",0);
set("z min",fdtd_z_min);
set("z max",fdtd_z_max);  
set("x span",Lx);

}
"""

    if scheme is None:
        scheme = "SiN on Top"

    if scheme not in _SCATTERER_SCHEMES:
        valid_schemes = ", ".join(_SCATTERER_SCHEMES)
        raise ValueError(f"Unknown scatterer scheme '{scheme}'. Valid schemes: {valid_schemes}")

    scheme_config = _SCATTERER_SCHEMES[scheme]

    assert len(LC_index_data) == len(LC_phase_data), "LC_index_data and LC_phase_data must have the same length."
    assert np.all(np.diff(np.array(LC_phase_data)) > 0), f"LC_phase_data must be strictly increasing. Got: {LC_phase_data}"

    LC_index_data = [float(x) for x in LC_index_data]
    LC_phase_data = [float(x) for x in LC_phase_data]
    # index_data_str = ",\n    ".join(f"{x:.16g}" for x in LC_index_data)
    # phase_data_str = ",\n    ".join(f"{x:.16g}" for x in LC_phase_data)
    
    index_data_str = "[\n    " + ",\n    ".join(f"{x:.16g}" for x in LC_index_data) + "\n]"

    phase_data_str = "[\n    " + ",\n    ".join(f"{x:.16g}" for x in LC_phase_data) + "\n]"
    
    template = template.replace("__SCATTER_Z_POSITION__", scheme_config["z_position"]) \
            .replace("__PERIOD__", f"{period:.16g}") \
            .replace("__T_AL__", f"{t_Al:.16g}") \
            .replace("__T_SPACER__", f"{t_spacer:.16g}") \
            .replace("__T_LC__", f"{t_LC:.16g}") \
            .replace("__T_ITO__", f"{t_ITO:.16g}") \
            .replace("__T_GLASS__", f"{t_glass:.16g}") \
            .replace("__R_SCATTER__", f"{r_scatter:.16g}") \
            .replace("__T_SCATTER__", f"{t_scatter:.16g}") \
            .replace("__LAMBDA_DESIGN__", f"{lambda_design:.16g}") \
            .replace("__LAMBDA_SPAN__", f"{lambda_span:.16g}") \
            .replace("__LC_INDEX_DATA__", index_data_str) \
            .replace("__LC_PHASE_DATA__", phase_data_str) \
            .replace("__SCATTER_LABEL__", scheme_config["label"]) \
            .replace("__SCATTER_MATERIAL__", scheme_config["material"]) \
            .replace("__SCATTER_PATCH_NAME__", scheme_config["patch_name"]) \
            .replace("__SCATTER_GROUP_NAME__", scheme_config["group_name"]) \
            .replace("__SCATTER_POSITION_COMMENT__", scheme_config["position_comment"]) \

    return template
