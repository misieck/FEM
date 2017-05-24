
clear all;
constants;
coord=[0 0;
        1 0;
        0 1;
        2 1;
        3 0;
        4 0;
        4 1;
        3 2;
        1 2;
        0 3;
        2 3; 
        4 3];
 coord = coord/10000;
 edof=[1 1 2 3;
       2 2 4 3;
       3 2 5 4;
       4 5 7 4;
       5 6 7 5
       6 3 4 9;
       7 4 7 8;
       8 3 9 10;
       9 4 8 9;
       10 7 12 8;
       11 9 11 10;
       12 9 8 11;
       13 8 12 11]
   
%coord=points' ;
ndof=max(max(edof(:, 2:4)));
nelem = length(edof(:,1));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
eldraw2(Ex,Ey,[1,4,1]);

   
dt = 0.1;
Time = 6;

K = sparse(ndof,ndof);
K_c = sparse(ndof,ndof);
f = sparse(ndof,1);
f_b = f;
C = sparse(ndof, ndof);
edge_length = abs(coord(10,1) - coord(11,1));
k = 0;
c_ro = 0;

f_b(10) = q_el * edge_length/2;
f_b(11) = 2 * q_el * edge_length/2;
f_b(12) = q_el * edge_length/2;

for elem = 1:nelem
    k = k_pcb;
    c_ro = c_pcb * ro_pcb;
    D= [k 0;
        0 k];
    Ke = flw2te(Ex(elem,:), Ey(elem, :),1, D);
    K_old = K;
    K = assem(edof(elem,:),K, Ke);
    Ce = plantml(Ex(elem,:), Ey(elem,:), c_ro);
    C = assem(edof(elem,:),C, Ce); 
    if max(max(K)) > 300
        
        disp 'hej'
    end
end
K_1 = K;

d0 = T_0 * ones(ndof,1);
ip = [dt, Time, 0.4, [3, 5, [0.1 0.3 2.8], [1, 2, 6, 4, 11] ]  ];
Tsnap = step1(K, C, d0, ip, f_b, []);