E_smd = 105e9; %Young's modulus [Pa]
E_pcb = 105e9;
E_sol = 50e9;

nu_smd = 0.118; %Poisson's ratio
nu_pcb = 0.136;
nu_sol = 0.36;

k_smd = 0.29; %thermal conductivity [W/mK]
k_pcb = 1.059;
k_sol = 66.8;

ro_smd = 1850; %mass density [kg/m^3]
ro_pcb = 1850;
ro_sol = 7265;

c_smd = 950; %specfic heat conductivity [J/(kgK]
c_pcb = 950;
c_sol = 210;

alpha_smd = 1.2e-5; %termal expansion coefficient  [1/K]
alpha_pcb = 2e-5;
alpha_sol = 1.2e-5;

%Create the three possible D's for the different materials
D_smd = hooke(2, E_smd, nu_smd);
D_smd = red(D_smd, 3);
D_pcb = hooke(2, E_pcb, nu_pcb);
D_pcb = red(D_pcb, 3);
D_sol = hooke(2, E_sol, nu_sol);
D_sol = red(D_sol, 3);

T_0 = 30; %initial temperature [C]
alpha_c = 40; %convection coefficient [W/(m^2K)]
q_el = 9e3; %heat flow [W/m^2]
T_inf = 20; %surrounding temperature [C]