function draw_temps(Ex, Ey, edof, Temps, figurenr, colormapscale)

figure(figurenr);
hold on;
ed=extract(edof,Temps);
%colormap('jet')
fill(Ex',Ey',ed');
fill(-Ex', Ey', ed');
colorbar;
caxis(colormapscale);

% figure(1)
% hold on
% ed = extract(edof,astat);
% fill(Ex', Ey', ed');
% fill(-Ex', Ey', ed');