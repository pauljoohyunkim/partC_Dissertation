%% Chebyshev Curve Repulsion
%
%% Curve Interpolation
% For testing purposes, create an octagon with the ends not connected.
pointnum=16
xpts = zeros(pointnum,1);
ypts = zeros(pointnum,1);
zpts = zeros(pointnum,1);
for i = 1:pointnum
    theta = 2 * pi / pointnum * i;
    xpts(i) = cos(theta);
    ypts(i) = sin(theta);
    zpts(i) = 0;
end

% Creating interpolant chebfun
x = interpolateAtChebpts(xpts);
y = interpolateAtChebpts(ypts);
z = interpolateAtChebpts(zpts);

%%
% Energy Computation
tangentPointEnergy(x, y, z, pointnum, 2, 4)

%%
% Gradient Computation
tangentPointEnergyGradient(x, y, z, pointnum, 0.00001, 2, 4)

%% Gradient Descent
stepsize=0.0005;
resolution=90;
perturbation=0.0001;
epsilon = 0.0001;
alpha=2;
beta=4;
J = length(x);
for t = 1:10000
    coeffsx = chebcoeffs(x);
    coeffsy = chebcoeffs(y);
    coeffsz = chebcoeffs(z);
    coeffs = [coeffsx; coeffsy; coeffsz];
    coeffs = coeffs - stepsize * tangentPointEnergyGradient(x, y, z, resolution, perturbation, alpha, beta);
    %tangIentPointEnergy(x,y,z,resolution, alpha, beta)
    %Deflation
    i = 0;
    while norm(coeffs(:,end-i)) < epsilon
        i = i + 1;
    end
    J = J - i;

    x = chebfun(coeffs(1:J), 'coeffs');
    y = chebfun(coeffs(J+1:2*J), 'coeffs');
    z = chebfun(coeffs(2*J + 1:3*J), 'coeffs');
    plot(x,y);
    timeIndex = sprintf("t=%d",t);
    title(timeIndex)
    drawnow;
end