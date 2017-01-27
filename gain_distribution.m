load('Ducati_jam2_LHCP.mat')
close all;
figure(1)
surf(azimvals-180, elevvals, gainvals','EdgeColor', 'None', 'facecolor', 'interp');
view(2);
colormap parula
% caxis auto
caxis([min(min(gainvals)), max(max(gainvals))])
colorbar
h = colorbar;
ylabel(h, 'Transmitter Gain [-]', 'Fontsize', 16)
axis([-180, 180, 0 90])
xlabel('Azimuth \theta [deg]', 'Fontsize', 16);
ylabel('Elevation \phi [deg]', 'Fontsize', 16);
fig = gcf;
fig.PaperUnits = 'normalized';
fig.PaperPosition = [0 0 1 0.4];

iptsetpref('ImshowBorder','tight');
print('-r300', '-opengl', '-depsc', '../report/figures/gain_map')