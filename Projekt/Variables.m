E_smd = 105e9; %Young's modulus [Pa]
E_pcb = 105e9;
E_sol = 50e9;

nu_smd = 0.118; %Poisson's ratio
nu_pcb = 0.136;
nu_sol = 0.36;

k_smd = 0.29; %thermal conductivity [W/mK]
k_pcb = 1.059;
k_sol = 66.8;

rho_smd = 1850; %mass density [kg/m^3]
rho_pcb = 1850;
rho_sol = 7625;

c_smd = 950; %specfic heat conductivity [J/(kgK]
c_pcb = 950;
c_sol = 210;

alpha_smd = 1.2e-5; %termal expansion coefficient  [1/K]
alpha_pcb = 2e-5;
alpha_sol = 1.2e-5;

T_0 = 30; %initial temperature [C]
alpha_c = 40; %convection coefficient [W/(m^2K)]
q_el = 9e3; %heat flow [W/m^3]
T_inf = 20 %surrounding temperature [C]