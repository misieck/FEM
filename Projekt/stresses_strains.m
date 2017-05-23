clear all;
project;

%Define stationary temperature
Temp_stat = a_stat;

%vi har fler (dubbelt s� m�nga) frihetsgrader - redigera edof och ndof
%Kan detta g�rasp� snyggare s�tt?

edof_old = edof;
edof = [edof(:, 1:2) edof(:,2)+ndof, edof(:, 3), edof(:,3)+ndof, ...
    edof(:,4) edof(:,4)+ndof];
ndof = 2*ndof;

%Skapa tom K, f och f0
K = sparse(ndof,ndof);
f = sparse(ndof, 1);
f0 = sparse(ndof, 1);

% %kolla s� elementnummer st�mmer med r�tt D
% %Tillskriv r�tt D beroende p� material samt assembla Ke och fe0
for i = 1:nelem
 [E,nu,k,ro,c,alpha,D] = decideElementproperties(elements, i);
 
 %elementets x-koordinater i noderna    
 x = Ex(i,:)';
 %elementets y-koordinater i noderna
 y = Ey(i,:)';
 
 %Vi vill veta vilka noder som �r kopplade...
 nodes = edof_old(i,2:4);
 %för attt kunna plocka ut rätt temperatur ur a_stat
 T = Temp_stat(nodes);
 
 %Räkna ut Ke
 Ke = plante(Ex(i,:), Ey(i,:), [2 1], D);
 
 %Räkna ut f0e, lite OSÄKER PÅ DETTA UTTRYCK!!!!
 T_avg = sum(T)/3;
 constant=0.5*alpha*E/(1-2*nu)*(T_avg-T_0);
 Be110 = [y(2)-y(3); x(3)- x(2); y(3)-y(1); x(1)-x(3); y(1)-y(2); x(2)-x(1)];
 f0e = constant * Be110;
 [K, f0] = assem(edof(i, :), K, Ke, f0, f0e);
 
end

%Find nodes on boundary with displacement conditions
ux0_edge = find(coord(:,1) > 1e-3 - 1e-6);
uy0_edge = ux0_edge + ndof/2;

ux0_middle = find(coord(:,1) <  1e-6);

ux0_bottom = find(coord(:,2)<0 + 1e-6);
uy0_bottom = ux0_bottom + ndof/2;

bc_dof = unique([ux0_edge; uy0_bottom]);

%Defines boundary condition
bc = [bc_dof zeros(length(bc_dof), 1)];
f = f + f0;


%Find displacements
a = solveq(K,f,bc);


%stress
Seff_el = zeros(nelem,1);

for i = 1:nelem
    [E,nu,k,ro,c,alpha,D] = decideElementproperties(elements, i);
     
    enod = edof(i, 2:7);
    enod2 = edof_old(i, 2:4);
    Tnodes = Temp_stat(enod2);
    [es, et] = plants(Ex(i,:), Ey(i,:), [2 1], D, a(enod)' );
    dTavg = (Tnodes(1)+Tnodes(2)+Tnodes(3))/3 - 30;
    De0 = alpha*E*dTavg/(1-nu) * [1; 1; 0];
    es = es-De0';
    
    sxx=es(1);
    syy=es(2);
    G = 0.5*E/(1-nu);
    tao_xy = G*et(3); 
    
    szz = nu*(sxx-syy)-alpha*E*dTavg;
    
    Seff_temp=sqrt(sxx^2+syy^2+szz^2-sxx*syy-sxx*szz-syy*szz+3*tao_xy^2);
    Seff_el(i) = Seff_temp;
end

Seff_nod = zeros(ndof/2, 1);

for i=1: size(coord,1)
    [c0,c1] = find(edof_old(:,2:4)==i);
    Seff_nod(i) = sum(Seff_el(c0))/size(c0,1);
end


%Ed1=extract(Edof1,a);
%Ed2=extract(Edof2,a);
%[sfac]=scalfact2(Ex1,Ey1,Ed1,0.2);

%eldisp2(Ex1,Ey1,Ed1,[2 1 1],sfac);
%eldisp2(Ex2,Ey2,Ed2,[2 1 1],sfac);


%PRINTA!!!
figure(6)
eldraw2(Ex,Ey,[1,4,1]);
Ed = extract(edof,a);
[sfac] = eldisp2(Ex,Ey,Ed, [2 1 1]);
eldisp2(Ex,Ey,Ed,[2 1 1],sfac);

draw_temps(Ex, Ey, edof_old, Seff_nod, 7, scale);