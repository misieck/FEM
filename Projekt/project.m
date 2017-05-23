%clear all
constants;

stuff = load('mesh.mat');
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

%close all;
eldraw2(Ex,Ey,[1,4,1]);

dt = 0.2;
Time = 2400;

K = sparse(ndof,ndof);
K_c = sparse(ndof,ndof);
f = sparse(ndof,1);
f_b = f;
C = sparse(ndof, ndof);

k = 0;
c_ro = 0;


q_edges = find (edges(5,:)==3 & points (1, edges(1,:)) < 0.0002 + 1e-6 & points (1, edges(2,:)) < 0.0002 + 1e-6  );

conv_edges = setdiff (find ( edges(5,:)==3 ), q_edges);
conv_edges = union (conv_edges, find ( edges(5,:)==4 ));



for elem = 1:nelem
    [E,nu,k,ro,c,alpha,De] = decideElementproperties(elements, elem);
    D= [k 0;
        0 k];
    Ke = flw2te(Ex(elem,:), Ey(elem, :),1, D);
    K_old = K;
    K = assem(edof(elem,:),K, Ke);
    Ce = plantml(Ex(elem,:), Ey(elem,:), c);
    C = assem(edof(elem,:),C, Ce);
end


for edge = q_edges
    nodes = edges(1:2, edge);
    edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );
    f_b_e = q_el * edge_length/2*ones(2,1);
    f_b(nodes) = f_b(nodes)+f_b_e;
end

% for edge = conv_edges
%     nodes = edges(1:2, edge);
%     edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );
%     f_b_e = sparse(ndof, 1);
%     f_b_e(nodes) = T_inf * alpha_c * edge_length/2;
%     f_b = f_b + f_b_e;
% end

for edge = conv_edges
    nodes = edges(1:2, edge);
     edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );
%     K_c_e = sparse(ndof, ndof);
%     K_c_e (nodes(1), nodes(2)) = alpha_c * edge_length/6;
%     K_c_e (nodes(2), nodes(1)) = K_c_e(nodes(1), nodes(2));
%     K_c_e (nodes(1), nodes(1)) = alpha_c * edge_length/3;
%     K_c_e (nodes(2), nodes(2)) = K_c_e (nodes(1), nodes(1)) ;
%     K_c = K_c + K_c_e;
        
        f_be = T_inf * alpha_c * edge_length/2*ones(2,1);
        K_ce = edge_length*alpha_c/6 * [2, 1; 1, 2];
        [K_c, f_b] = assem([1, nodes(1), nodes(2)], K_c, K_ce, f_b, f_be);
 
end

K = K + K_c;

pbound = unique (reshape(edges(1:2, conv_edges), 1,[]));
pbound = [pbound' T_0*ones(length(pbound),1)];
time_history_of_the_load = f_b;

d0 = T_0 * ones(ndof,1);
ip = [dt, Time, 1, [4, ndof, [1 10 60 240], [1:ndof] ]  ];
[Tsnap] = step1(K, C, d0, ip, f_b, []);


%stationary_index = find_stationary(V, 0.0001, 50);
%stationary_temps = D(:, stationary_index);

scale = 'auto';
for i = 1:3
    draw_temps(Ex, Ey, edof, Tsnap(:,i), i+1, scale);
end

%draw_temps(Ex, Ey, edof, stationary_temps, 5, scale);
a_stat = solveq(K, f_b);
draw_temps(Ex, Ey, edof, a_stat, 5, scale);

% 
% % Starting temperature
% a0 = T_0*ones(ndof,1);
% % Solve equation system over time
% dt = 1; % Time step
% Kt = C/dt + K;
% for i = 1:250 % Run for 10 s
%     ft = f_b + C*a0/dt;
%     a=solveq(Kt,ft);
%    
%     a0 = a;
% end


