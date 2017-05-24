clear all;
close all;

constants;

stuff = load('mesh-lo.mat');
edges = stuff.edges;

%correct scale in stuff.points 
points = stuff.points/1000;
elements = stuff.triangles;
nelem=length(elements(1,:));
edof = zeros(nelem,4);
edof(:,1)=1:nelem ;
edof(:,2:4)=elements(1:3,:)' ;
coord=points' ;
ndof=max(max(elements(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

[~,~,k,ro,c,~] = decideElementproperties(elements, nelem);

%eldraw2(Ex,Ey,[1,4,1]);

dt = 0.2;
Time = 240;

K = sparse(ndof,ndof);
K_c = sparse(ndof,ndof);
f = sparse(ndof,1);
f_b = f;
C = sparse(ndof, ndof);


q_edges = find (edges(5,:)==3 & points (1, edges(1,:)) < 0.0002 + 1e-6 & points (1, edges(2,:)) < 0.0002 + 1e-6  );

conv_edges = setdiff (find ( edges(5,:)==3 ), q_edges);
conv_edges = union (conv_edges, find ( edges(5,:)==4 ));



for elem = 1:nelem
    
    D = [k(elem) 0;
        0 k(elem)];
    
    c_ro = c(elem)*ro(elem);
    
    Ke = flw2te(Ex(elem,:), Ey(elem, :),1, D);
    K = assem(edof(elem,:),K, Ke);
    
    Ce = plantml(Ex(elem,:), Ey(elem,:), c_ro);
    C = assem(edof(elem,:),C, Ce);
end


for edge = q_edges
    nodes = edges(1:2, edge);
    edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );
    f_be = q_el * edge_length/2*ones(2,1);
    f_b(nodes) = f_b(nodes)+f_be;
end


for edge = conv_edges
    nodes = edges(1:2, edge);
     edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );

        
        f_be = T_inf * alpha_c * edge_length/2*ones(2,1);
        K_ce = edge_length*alpha_c/6 * [2, 1; 1, 2];
       [K_c, f_b] = assem([1, nodes(1), nodes(2)], K_c, K_ce, f_b, f_be);
 
end

K = K + K_c;

pbound = unique (reshape(edges(1:2, conv_edges), 1,[]));
pbound = [pbound' T_0*ones(length(pbound),1)];

d0 = T_0 * ones(ndof,1);
ip = [dt, Time, 1, [4, ndof, [1 10 60 240], 1:ndof ]  ];

[Tsnap, D, V] = step1(K, C, d0, ip, f_b, []);

scale = 'auto';
for i = 1:3
    draw_temps(Ex, Ey, edof, Tsnap(:,i), i+1,scale);
end

a_stat = solveq(K, f_b);
draw_temps(Ex, Ey, edof, a_stat, 5, scale);
