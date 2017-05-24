function draw_temps(Ex, Ey, edof, Temps, figurenr)

figure(figurenr);
hold on;
scale = 'auto';
ed=extract(edof,Temps);
colormap('jet')
fill(Ex',Ey',ed');
fill(-Ex', Ey', ed');
colorbar;
caxis(scale);

