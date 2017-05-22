clear all;
project;

%Define stationary temperature
Temp_stat = stationary_temps;

%vi har fler (dubbelt s� m�nga) frihetsgrader - redigera edof och ndof
%Kan detta g�rasp� snyggare s�tt?
add = ndof*ones(nelem,1);
edof_old = edof;
edof = [edof(:, 1:2) edof(:,2)+add, edof(:, 3), edof(:,3)+add, edof(:,4) edof(:,4)+add];
ndof = 2*ndof;

%Skapa tom K, f och f0
K = sparse(ndof,ndof);
f = sparse(ndof, 1);
f0 = sparse(ndof, 1);


%Skapa de tre m�jliga D som finns f�r de olika elementen
D_smd = hooke(2, E_smd, nu_smd);
D_smd = red(D_smd, 3);
D_pcb = hooke(2, E_pcb, nu_pcb);
D_pcb = red(D_pcb, 3);
D_sol = hooke(2, E_sol, nu_sol);
D_sol = red(D_sol, 3);


%kolla s� elementnummer st�mmer med r�tt D
%Tillskriv r�tt D beroende p� material samt assembla Ke och fe0
for i = 1:nelem
     switch elements(4, i)
        case 1
            D = D_pcb;
            nu = nu_pcb;
            E = E_pcb;
            alpha = alpha_pcb;
            
        case 2
            D = D_smd;
            nu = nu_smd;
            E = E_smd;
            alpha = alpha_smd;
            
        case 3
            D = D_sol;
            nu = nu_sol;
            E = E_sol;
            alpha = alpha_sol;
     end
 
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

%Hitta noder med "randvillkor"
%Kolla så att dessa punkter stämmer, vad menar man med bara ux och xh uy?
ux0_edge = find(coord(:,1) > 1e-3 - 1e-6);
uy0_edge = ux0_edge + ndof/2;

ux0_middle = find(coord(:,1) >  - 1e-6);

ux0_bottom = find(coord(:,2)<0 + 1e-6);
uy0_bottom = ux0_bottom + ndof/2;

bc_dof = unique([ux0_edge; uy0_bottom; ux0_middle]);

%Tänker jag rätt i att vi har randvillkor t=0 där förskjutningen ej är
%given? JA!

%definera randvillkor
bc = [bc_dof zeros(length(bc_dof), 1)];
f = f0;


%R�kna ut f�rskjutning
a = solveq(K,f,bc);


%stress
D_smd = hooke(2, E_smd, nu_smd);
D_pcb = hooke(2, E_pcb, nu_pcb);
D_sol = hooke(2, E_sol, nu_sol);

Seff_el = zeros(nelem,1);

for i = 1:nelem
         switch elements(4, i)
        case 1
            D = D_smd;
            
        case 2
            D = D_smd;
            
        case 3
            D = D_sol;
         end
    
    enod = edof(i, 2:7);    
    [es, et] = plants(Ex(i,:), Ey(i,:), [2 1], D, a(enod)' );
    sxx=es(1);
    syy=es(2);
    szz= es(3);
    sxy=es(4);
    Seff_temp=sqrt(sxx^2+syy^2+szz^2-sxx*syy-sxx*szz-syy*szz+3*sxy^2);
    Seff_el(i) = Seff_temp;
end

Seff_nod = zeros(ndof/2, 1);

for i=1: size(coord,1)
    [c0,c1] = find(edof(:,2:4)==i);
    Seff_nod(i,1) = sum(Seff_el(c0))/size(c0,1);
end


%Ed1=extract(Edof1,a);
%Ed2=extract(Edof2,a);
%[sfac]=scalfact2(Ex1,Ey1,Ed1,0.2);

%eldisp2(Ex1,Ey1,Ed1,[2 1 1],sfac);
%eldisp2(Ex2,Ey2,Ed2,[2 1 1],sfac);


%PRINTA!!!
Ed = extract(edof,a);
[sfac] = eldisp2(Ex,Ey,Ed, [2 1 1])
