clear all;
project;


%Här ska den stationära lösningen för temperaturen komma!!!!!!!!!!!!!!!
%Importera från project
Temp_stat = Tsnap(:,3);

%vi har fler (dubbelt så många) frihetsgrader - redigera edof och ndof
%Kan detta göraspå snyggare sätt?
add = ndof*ones(nelem,1);
edof_old = edof;
edof = [edof(:, 1:2) edof(:,2)+add, edof(:, 3), edof(:,3)+add, edof(:,4) edof(:,4)+add];
ndof = 2*ndof;

%Skapa tom K, f och f0
K = sparse(ndof,ndof);
f = sparse(ndof, 1);
f0 = sparse(ndof, 1);


%Skapa de tre möjliga D som finns för de olika elementen
D_smd = hooke(2, E_smd, nu_smd);
D_smd = red(D_smd, 3);
D_pcb = hooke(2, E_pcb, nu_pcb);
D_pcb = red(D_pcb, 3);
D_sol = hooke(2, E_sol, nu_sol);
D_sol = red(D_sol, 3);


%kolla så elementnummer stämmer med rätt D
%Tillskriv rätt D beroende på material samt assembla Ke och fe0
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
 
 %Vi vill veta vilka noder som är kopplade...
 nodes = edof_old(i,2:4);
 %för attt kunna plocka ut rätt temperatur ur a_stat
 T = [Temp_stat(nodes(1)); Temp_stat(nodes(2)); Temp_stat(nodes(3))];
 
 %Räkna ut Ke
 Ke = plante(Ex(i,:), Ey(i,:), [2 1], D);
 
 %Räkna ut f0e MYCKET OSÄKER PÅ DETTA UTTRYCK!!!!
 Ttot = sum(T);
 constant=0.5*alpha*E/(1-2*nu)*(1/3*Ttot-T_0);
 Be110 = [y(2)-y(3); x(3)- x(2); y(3)-y(1); x(1)-x(3); y(1)-y(2); x(2)-x(1)];
 f0e = constant * Be110;
 [K, f0] = assem(edof(i, :), K, Ke, f0, f0e);
 
end

%Hitta noder med "randvillkor"
%Kolla så att dessa punkter stämmer, vad menar man med bara ux och xh uy?
ux0 = find(coord(:,1) > 1e-3 - 1e-6);
ux00 = ux0 + ndof/2*ones(length(ux0),1);
uy0 = find(coord(:,2)<0 + 1e-6);
uy00 = uy0 + ndof/2*ones(length(uy0),1);
bc_dof = [ux0; ux00; uy0; uy00];

%Tänker jag rätt i att vi har randvillkor t=0 där förskjutningen ej är
%given?

%definera randvillkor
bc = [bc_dof zeros(length(bc_dof), 1)];
f = f0;


%Räkna ut förskjutning
a = solveq(K,f, bc);

%DENNA ÄR INTE KLAR ÄN!
%Räkna ut stresses and strains
t = [];
s = [];
vonMises = [];
for i = 1:nelem
         switch elements(4, i)
        case 1
            D = D_pcb;
            
        case 2
            D = D_smd;
            
        case 3
            D = D_sol;
         end
         
    enod = edof(i, 2:7);    
    [es, et] = plants(Ex(i,:), Ey(i,:), [2,1], D, enod);
    s = [s;es];
    t = [t;et];
    t2 = et.^2;
    t2 = sum(t2');
    vonMises=[vonMises;sqrt(t2)];
end

%PRINTA!!!
aEd = extract(edof,a);
[sfac] = eldisp2(Ex,Ey,Ed);
