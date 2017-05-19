

clear all
constants;
boundary_geometry_lo;
stuff = load('mesh-lo.mat');
edges = stuff.edges;
points = stuff.points;
triangles = stuff.triangles;
nelm=length(triangles(1,:));
edof(:,1)=1:nelm ;
edof(:,2:4)=triangles(1:3,:)' ;
coord=points' ;
ndof=max(max(triangles(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);
eldraw2(Ex,Ey,[1,4,1])   


K = zeros(ndof);
f = zeros(ndof,1);

for elem = 1:nelem
   
    
    Ke = spring1e();
    K = assem
    
end
