% Plot colored streamlins for flow past a sphere
% of radius a located centred at c in the domain

clear all; 
close all;

a = 0.2; % Sphere radius
c = [1.0 0.5]; % Sphere centre

% Setup grids
dx = 0.01;
x = 0:dx:2;
y = 0:dx:1;
[X,Y] = meshgrid(x,y);

% Compute streamfunction
R = sqrt((X-c(1)).^2 + (Y-c(2)).^2);
theta = atan((Y-c(2))./(X-c(1)));
Streamfunction = R.^2 .* (1-a^3./R.^3).*sin(theta).^2;
Streamfunction(R <= a) = 0;

% Specify contour spacing
maxPsi = max(max(abs(Streamfunction)));
v = -maxPsi:maxPsi/20:maxPsi;

% Setup figure and colormap
h = figure();
h.Position = [200 200 800 400];
ax = axes;
colormap(ax, autumn);

% Make colored streamfunction plot
cbar = streamfunctioncolor(X,Y,Streamfunction,v);

% Ensure proper aspect ratio
daspect([1 1 1]);

% Add some labels and decorations
box on;
title('Flow past a sphere, $\psi=r^2$sin$(\theta)^2 [1-(a/r)^3] $');
cbar.Label.String = 'Scaled |U|';

% Save image
savename = [mfilename('fullpath'), '.png'];
fprintf('Writing plotfile %s \n', savename);
print(h,savename,'-dpng','-r200')
