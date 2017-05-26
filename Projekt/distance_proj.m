%Calculates the length of between to points with coordinated (a_x, a_y) and
%(b_x, b_y)
function d=distance_proj(a, b)

x1 = a(1);
y1 = a(2);

x2 = b(1);
y2 = b(2);

d = sqrt((x2 -x1)^2  + (y2 - y1)^2 );
