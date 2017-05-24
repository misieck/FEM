function draw_temps(Ex, Ey, edof, Temps, figurenr,scale)

figure(figurenr);
hold on;
ed=extract(edof,Temps);
%colormap('jet')
fill(Ex',Ey',ed');
fill(-Ex', Ey', ed');
colorbar;
caxis(scale);

