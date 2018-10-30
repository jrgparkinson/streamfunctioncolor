function cVel = streamfunctioncolor(X,Y,U,V,Streamfunction,v)
%STREAMFUNCTIONCOLOR  Color streamlines according to magnitude of the velocity
%   Inspired by STREAMCOLOR, by Bertrand Dano
%   (https://uk.mathworks.com/matlabcentral/fileexchange/24049-streamcolor)
%   The difference between this function and STREAMCOLOR is that here the 
%   streamlines are evenly spaced in streamfunction space, which can be
%   useful in some circumstances.
%
%    cbar = STREAMCOLOR(X,Y,U,V,Streamfunction,v)
%    cbar = STREAMCOLOR(X,Y,Streamfunction,v) 
%             Computes U,V from Streamfunction and X,Y. Assumes dx=dy.
% 
%    Assumes X,Y are meshgrids
%            U,V are velocities defined on the meshgrid
%            Streamfunction is the Stokes Streamfunction defined on the meshgrid
%            v contains the contour levels, as used by contourc
%
%    If you do not have the streamfunction you can calculate it from U,V data by solving
%    \nabla^2 \psi = \nabla \times (U,V), where \psi is the streamfunction.
%    There are plenty of options on the MATLAB exchange for taking the curl of a vector field
%    and solving poisson's equation.
%
%    Example:
%     [X,Y] = meshgrid(0:0.01:1, 0:0.01:1);
%     Streamfunction = cos(pi*X).*sin(Y*pi);
%     v = -1:0.05:1;
%     ax = axes; colormap(ax, autumn);
%     cbar = streamfunctioncolor(X,Y,Streamfunction,v);
%
%    See also STREAMCOLOR, STREAKARROW, COUNTOURC
%
%   Jamie Parkinson 29-10-2018
%   Copyright 1984-2018 The MathWorks, Inc.

if nargin < 5
    % Assuming inputs of form STREAMCOLOR(X,Y,Streamfunction,v)
    Streamfunction = U; v = V;
    
	dx = X(1,2) - X(1,1); dy = Y(2,1) - Y(1,1);
	[U,V] = gradient(Streamfunction, dx, dy);
	
	% Correct the sign, not that this matters
	U = -U;
end

% Line width for contours
lw = 0.5;

% Find the colormap we're using
cmap = colormap(gca);

% Need vectors of x and y for contourc
x = X(1, :);
y = Y(:, 1);

% Get contour data
C = contourc(x,y,Streamfunction, v);

% Make colorbar
cVel = colorbar('Location', 'eastoutside');

% Compute speed field, normalised between 0 and 1, in order to decide how
% to color the lines
magU = sqrt(U.^2 + V.^2);
maxU = max(max(abs(magU)));
relUField = magU/maxU;

% Compile vertices into a cell array
contours = {};
thisContour = [];

% C is a ContourMatrix (see
% https://uk.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.contour-properties.html#budgut_-ContourMatrix)
% Each pair is either an (x,y) co-ordinate, or specifies the contour level
% and number of vertices in the contour

% Number of vertices in the first contour
numVertices = C(2,1);
vertexCount = 1;

for i=2:length(C)
    c = C(:, i); 
    
    if vertexCount > numVertices
        % We've reached the end of a contour - save it and start the next
        % one
        contours{end+1} = thisContour;
        thisContour = [];
        vertexCount = 1;
        numVertices = c(2);
        
    else
        % Add vertex to current contour
        thisContour(end+1, :) = c;
        vertexCount = vertexCount + 1;
    end
    
end

% Plot contours in line segments, colored by the size of U 
hold on;

for i=1:length(contours)
    
    % c contains an array of xy pairs
    c = contours{i};
    
    numVertices = size(c, 1);
    
    for j=2:numVertices
        relU = interp2(X, Y, relUField, c(j, 1), c(j, 2));
        relU = min(relU, 1); % Ensure relU is <= 1
        color = cmap(ceil(relU*length(cmap)), :); % Get color from colormap
        plot([c(j, 1), c(j-1, 1)], [c(j, 2), c(j-1, 2)], 'Color', color, 'LineWidth', lw);
    end
end

hold off;


end


function [U,V] = computeVelocity(X,Y,Streamfunction)

end
