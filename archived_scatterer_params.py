SiN_b_565nm = {
    'scatterer_scheme': "SiN on Bottom",
    'period': 0.36e-6,
    't_Al': 0.30e-6,
    't_spacer': 0.17e-6,
    't_LC': 0.50e-6,
    't_ITO': 0.2e-6,
    't_glass': 0.1e-6,
    'r_scatter': 0.11e-6,
    't_scatter': 0.20e-6,
} # 565nm 时 t_ITO 从 2e-6 调整到 0.2e-6，其他参数保持不变


SiN_t_490nm = {
    'scatterer_scheme': "SiN on Top",
    'period': 0.36e-6,
    't_Al': 0.30e-6,
    't_spacer': 0.17e-6,
    't_LC': 0.50e-6,
    't_ITO': 0.05e-6,
    't_glass': 0.1e-6,
    'r_scatter': 0.12e-6,
    't_scatter': 0.20e-6,
} # 490nm 时 t_ITO 从 0.2e-6 调整到 0.05e-6，其他参数保持不变


ACSNano_TiO2_665nm = {
    'scatterer_scheme': "TiO2 on Bottom",
    'period': 0.36e-6,
    't_Al': 0.15e-6,
    't_spacer': 0.17e-6,
    't_LC': 0.70e-6,
    't_ITO': 0.023e-6,
    't_glass': 0.1e-6,
    'r_scatter': 0.13e-6,
    't_scatter': 0.20e-6,
    'lambda_start': 550e-9,
    'lambda_stop': 700e-9,
} 
# Below is from ACS Nano: https://pubs.acs.org/doi/10.1021/acsnano.3c04071
# a planar bottom metallic (Al) reflector (hAl = 150 nm), 
# a silicon dioxide (SiO2) spacer (hSpacer = 170 nm), 
# an array of titanium dioxide (TiO2) disk resonators (hTiO2) = 200 nm, 
# dTiO2 = 260 nm, PTiO2 = 360 nm), an LC cell surrounding the metasurface (hLC = 500 nm), 
# and a top glass superstrate coated by a thin indium tin oxide (ITO, hITO = 23 nm) layer in a bottom-up configuration, 
# where h, d, and P stand for height, diameter, and periodicity, respectively. In our finite-difference-time-domain (FDTD, Lumerical) simulation model, 
# when the incident light propagates through the LC cell along the z-axis (k̅||ẑ) and is polarized along the same direction as the in-plane LC directors (x-axis in Figure 1a, i.e., Eı = E0 x̂), 
# it “sees” the extraordinary refractive index (ne) of LC. The LC directors can rotate in xz-plane, subtending an angle θLC with the z axis (i.e., θLC = 90°along x-axis and θLC = 0° along z-axis) that can be controlled with an applied bias. 
# At their vertical orientation, the LC directors are said to be switched completely, 
# and light “sees” an ordinary refractive index (n0) while passing through the LC cell. 
# We used commercially available LC with ne = 1.5221 and n0 = 1.522 (i.e., birefringence of Δn = 0.291). 


TiO2_t_Zhihao = {
    'scatterer_scheme': "TiO2 on Top",
    'period': 0.36e-6,
    't_Al': 0.30e-6,
    't_spacer': 0.17e-6,
    't_LC': 0.50e-6,
    't_ITO': 0.023e-6,
    't_glass': 0.1e-6,
    'r_scatter': 0.13e-6,
    't_scatter': 0.22e-6,
    'lambda_start': 550e-9,
    'lambda_stop': 700e-9,
}


SiN_b_Zhihao = {
    'scatterer_scheme': "SiN on Bottom",
    'period': 0.36e-6,
    't_Al': 0.30e-6,
    't_spacer': 0.17e-6,
    't_LC': 0.50e-6,
    't_ITO': 0.023e-6,
    't_glass': 0.1e-6,
    'r_scatter': 0.13e-6,
    't_scatter': 0.22e-6,
    'lambda_start': 550e-9,
    'lambda_stop': 700e-9,
}