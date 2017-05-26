function Stress(Temp_stat)
constants;

Topology;


%We have twice as many dof - correct edof och ndof
edof_old = edof;
edof = [edof(:, 1:2) edof(:,2)+ndof, edof(:, 3), edof(:,3)+ndof, ...
    edof(:,4) edof(:,4)+ndof];
ndof = 2*ndof;

%Skapa tom K, f och f0
K = sparse(ndof,ndof);
f0 = sparse(ndof, 1);

 [E,nu,~, ~, ~, alpha] = decideElementproperties(elements, nelem);
 
%Calculate D depending on material, and assembla Ke och fe0
for i = 1:nelem
 
 De = hooke(2, E(i), nu(i));
 De = red(De, 3);
 
 %elementets x-koordinater i noderna    
 x = Ex(i,:)';
 %elementets y-koordinater i noderna
 y = Ey(i,:)';
 
 %Vi vill veta vilka noder som ar kopplade...
 nodes = edof_old(i,2:4);
 %för attt kunna plocka ut rätt temperatur ur a_stat
 T = Temp_stat(nodes);
 
 %Räkna ut Ke
 Ke = plante(Ex(i,:), Ey(i,:), [2 1], De);
 
 %Räkna ut f0e
 T_avg = sum(T)/3;
 a_tmp_constant=0.5*alpha(i)*E(i)/(1-2*nu(i))*(T_avg-T_0);
 
 e0 = (1 + nu(i)) * alpha(i) * (T_avg-T_0) * [1; 1; 0];
 
 Be_110 = [y(2)-y(3); x(3)- x(2); y(3)-y(1); x(1)-x(3); y(1)-y(2); x(2)-x(1)];
 f0e = a_tmp_constant * Be_110;
 [K, f0] = assem(edof(i, :), K, Ke, f0, f0e);
 
end

%Find nodes on boundary with displacement conditions
ux0_edge = find(coord(:,1) > 1e-3 - 1e-6);

ux0_middle = find(coord(:,1) <  1e-6);

ux0_bottom = find(coord(:,2) <  1e-6);
uy0_bottom = ux0_bottom + ndof/2;

bc_dof = unique([ux0_edge; ux0_middle; uy0_bottom]);

%Defines boundary condition
bc = [bc_dof zeros(length(bc_dof), 1)];


%Find displacements
u = solveq(K,f0,bc);

%von Mises effective stress, elementwise
Seff_el = zeros(nelem,1);

for i = 1:nelem
    
    De = hooke(2, E(i), nu(i));
    De = red(De, 3);

    e_nodes = edof(i, 2:7);

    %calculate elementwise stress and strains
    [es, et] = plants(Ex(i,:), Ey(i,:), [2 1], De, u(e_nodes)' );
    
    %avarage temperature increase in an element
    enod2 = edof_old(i, 2:4);
    T_nodes = Temp_stat(enod2);
    dT_avg = (T_nodes(1)+T_nodes(2)+T_nodes(3))/3 - T_0;
    
    De_e0 = alpha(i)*E(i)*dT_avg/(1-2*nu(i)) * [1; 1; 0];
    es = es-De_e0';
    
    sxx=es(1);
    syy=es(2);
    
    
    G = 0.5*E(i)/(1+nu(i));
    tau_xy = G*et(3); 
    szz = nu(i)*(sxx+syy)-alpha(i)*E(i)*dT_avg; %13.42
    
    Seff_el(i) = sqrt(sxx^2+syy^2+szz^2-sxx*syy-sxx*szz-syy*szz+3*(tau_xy)^2); 
end

%get nodewise Seffs 
Seff_nod = zeros(ndof/2, 1);

for i=1: size(coord,1)
    [c0,c1] = find(edof_old(:,2:4)==i);
    Seff_nod(i) = sum(Seff_el(c0))/size(c0,1);
end

%PRINTA!!!
figure(1)
hold on
Ed = extract(edof,u);
eldraw2(Ex,Ey,[1,4,1]);
eldraw2(-Ex,Ey,[1,4,1]);

Ed_neg = [-Ed(:,1),Ed(:,2), -Ed(:,3), Ed(:,4), -Ed(:,5), Ed(:,6)];
[sfac] = eldisp2(Ex,Ey,Ed, [2 1 1]);
         eldisp2(-Ex,Ey, Ed_neg, [2 1 1]);
    xlabel('x [m]');
    ylabel('y [m]');
scale = 'auto';
draw_temps(Ex, Ey, edof_old, Seff_nod, 7, scale);
    xlabel('x [m]');
    ylabel('y [m]');
    h = colorbar;
    xlabel(h, 'Stress [Pa]')