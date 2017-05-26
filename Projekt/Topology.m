
%loads information from mesh generated from pdetool
stuff = load('mesh.mat');

%vector that identifies how nodal point on edges are connected
edges = stuff.edges;

%Our mesh has been generated in [mm] and needs to be converted to [m]
points = stuff.points/1000;

%Containd information about elements nodal point and which subdomain it
%belongs to
elements = stuff.triangles;

%Number of elements
nelem=length(elements(1,:));

%Creates topologymatrix, i.e which nodes the element has
edof = zeros(nelem,4);
edof(:,1)=1:nelem ;
edof(:,2:4)=elements(1:3,:)' ;

%Defines coordinates for nodal points
coord=points';

%ndof = number of nodal points
ndof=max(max(elements(1:3,:)));

%Ex cotains the x-coordinates for each element, and Ey contains the
%y-coordinates
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);