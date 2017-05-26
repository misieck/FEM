function [a_stat] = Temperature
constants;

Topology;

[~,~,k,ro,c,~] = decideElementproperties(elements, nelem);

K = sparse(ndof,ndof);
K_c = sparse(ndof,ndof);
f_b = sparse(ndof,1);
C = sparse(ndof, ndof);

%Defines K_e and C_e for each element and assembles them to K and C
for elem = 1:nelem
    
    D = [k(elem) 0;
        0 k(elem)];
    
    c_ro = c(elem)*ro(elem);
    
    Ke = flw2te(Ex(elem,:), Ey(elem, :),1, D);
    K = assem(edof(elem,:),K, Ke);
    
    Ce = plantml(Ex(elem,:), Ey(elem,:), c_ro);
    C = assem(edof(elem,:),C, Ce);
end

%Finds nodal points and corresponding element boundaries where the flux is q_el, i. e the nodes on boundary no. 3 where x is less or equal to 0.2
q_edges = find (edges(5,:)==3 & points (1, edges(1,:)) < 0.0002 + 1e-6 ...
    & points (1, edges(2,:)) < 0.0002 + 1e-6  );

%Assembles the element load vector f_be with the global load vector f_b for the elements with an boundary where the flux is q_el
for edge = q_edges
    nodes = edges(1:2, edge);
    edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );
    f_be = q_el * edge_length/2*ones(2,1);
    f_b(nodes) = f_b(nodes)+f_be;
end

%Finds nodal points and corresonding element boundaries where convection is present, i.e  the remainder of boundary no. 3 and boundary no. 5
conv_edges = setdiff (find ( edges(5,:)==3 ), q_edges);
conv_edges = union (conv_edges, find ( edges(5,:)==4 ));

%Assembles element K_ce and f_be with the global K_c and f_b for the elements with a boundary where convection is present
for edge = conv_edges
    nodes = edges(1:2, edge);
     edge_length = distance_proj( points(:, edges(1, edge)), points (:, edges(2, edge)) );

        
        f_be = T_inf * alpha_c * edge_length/2*ones(2,1);
        K_ce = edge_length*alpha_c/6 * [2, 1; 1, 2];
       [K_c, f_b] = assem([1, nodes(1), nodes(2)], K_c, K_ce, f_b, f_be);
 
end

%Defines new stiffness matrix K
K = K + K_c;

%Finds the stationary solution
a_stat = solveq(K, f_b);


%Solves the transient heat problem by using the CALFEM-function step1
dt = 0.2;
Time = 240;


d0 = T_0 * ones(ndof,1);
ip = [dt, Time, 1, [4, ndof, [1 10 60 240], 1:ndof ]  ];

[Tsnap] = step1(K, C, d0, ip, f_b, []);


%Draws the temperature at for different points in time
scale = 'auto';
for i = 1:4
    hold on
    draw_temps(Ex, Ey, edof, Tsnap(:,i), i+1,scale);
    xlabel('x [m]');
    ylabel('y [m]');
    h = colorbar;
    xlabel(h, 'Temperature [C]')

end

%Draws the stationary solution
hold on;
draw_temps(Ex, Ey, edof, a_stat, i+2, scale);
    xlabel('x [m]');
    ylabel('y [m]');
    h = colorbar;
    xlabel(h, 'Temperature [C]')
