

clear all
constants;

stuff = load('mesh-lo.mat');
edges = stuff.edges;

%correct scale in stuff.points 
points = stuff.points/1000;

elements = stuff.triangles;
nelem=length(elements(1,:));
edof(:,1)=1:nelem ;
edof(:,2:4)=elements(1:3,:)' ;
coord=points' ;
ndof=max(max(elements(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
eldraw2(Ex,Ey,[1,4,1]);

dt = 0.1;
Time = 1;

K = sparse(ndof,ndof);
f = sparse(ndof,1);
f_b = f;
C = sparse(ndof, ndof);

k = 0;
c_ro = 0;


q_edges = find (edges(5,:)==3 & points (1, edges(1,:)) < 0.0002 + 1e-6 & points (1, edges(2,:)) < 0.0002 + 1e-6  );

conv_edges = setdiff (find ( edges(5,:)==3 ), q_edges);
conv_edges = union (conv_edges, find ( edges(5,:)==4 ));



for elem = 1:nelem
    switch elements(4, elem)
        case 1
            k = k_pcb;
            c_ro = c_pcb * ro_pcb;
        case 2
            k = k_smd;
            c_ro = c_smd * ro_smd;
        case 3
            k = k_sol;
            c_ro = c_sol * ro_sol;
    end
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

for edge = q_edges
    nodes = edges(1:2, edge);
    edge_length = distance( points (:, edges(1, edge)), points (:, edges(2, edge)) );
    f_b_e = sparse(ndof, 1);
    f_b_e(nodes) = -q_el * edge_length/2;
    f_b = f_b + f_b_e;
end

for edge = conv_edges
    nodes = edges(1:2, edge);
    edge_length = distance( points (:, edges(1, edge)), points (:, edges(2, edge)) );
    f_b_e = sparse(ndof, 1);
    f_b_e(nodes) = T_inf * alpha_c * edge_length/2;
    f_b = f_b + f_b_e;
end

for edge = conv_edges
    nodes = edges(1:2, edge);
    edge_length = distance( points (:, edges(1, edge)), points (:, edges(2, edge)) );
    K_c_e = sparse(ndof, ndof);
    K_c_e (nodes(1), nodes(2)) = alpha_c * edge_length/6;
    K_c_e (nodes(2), nodes(1)) = K_c_e(nodes(1), nodes(2));
    K_c_e (nodes(1), nodes(1)) = alpha_c * edge_length/3;
    K_c_e (nodes(2), nodes(2)) = K_c_e (nodes(1), nodes(1)) ;
    K = K + K_c_e;
end
diag(K)


d0 = T_0 * ones(ndof,1);
ip = [dt, Time, 0, [3, 3, [0, 0.5, 1], [58, 53, 74] ]  ];
Tsnap = step1(K, C, d0, ip, f, []);


