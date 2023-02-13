%% Chebyshev Curve Repulsion
%
%% Curve Interpolation
% For testing purposes, create an octagon with the ends not connected.
xpts = zeros(8,1);
ypts = zeros(8,1);
zpts = zeros(8,1);
for i = 1:8
    theta = 2 * pi / 8 * i;
    xpts(i) = cos(theta);
    ypts(i) = sin(theta);
    zpts(i) = 0;
end

% Creating interpolant chebfun
x = interpolateAtChebpts(xpts);
y = interpolateAtChebpts(ypts);
z = interpolateAtChebpts(zpts);