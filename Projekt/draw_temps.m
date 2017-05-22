function draw_temps(Ex, Ey, edof, Temps, figurenr, colormapscale)
figure(figurenr);
ed=extract(edof,Temps);
%colormap(hot)
fill(Ex',Ey',ed')
colorbar;
caxis(colormapscale);