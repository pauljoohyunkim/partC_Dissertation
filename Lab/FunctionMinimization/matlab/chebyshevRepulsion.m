%% Chebyshev Curve Repulsion
%
%% Curve Interpolation
% For testing purposes, create an octagon with the ends not connected.
pointnum=16;
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
%% Initial Curve: Elliptic Curve
t = chebfun('t');
x = 0.000001*t^3 + (1.1*t)^2 - 1;
y = (1.1*t)^3 - (1.1*t);
z = 0.000001*t^3 + t;
%%
% Energy Computation
tangentPointEnergy(x, y, z, pointnum, 2, 4)

%%
% Gradient Computation
tangentPointEnergyGradient(x, y, z, pointnum, 0.00001, 2, 4)

%% Gradient Descent
stepsize=0.0005;
resolution=40;
perturbation=0.01;
epsilon = 0.0001;
alpha=2;
beta=4;
viewsize=1;
linewidth=4;
J = length(x);

snapshotInterval = 100;
snapshotNumber = 0;
for t = 1:301
    coeffsx = chebcoeffs(x);
    coeffsy = chebcoeffs(y);
    coeffsz = chebcoeffs(z);
    coeffs = [coeffsx; coeffsy; coeffsz];
    coeffs = coeffs - stepsize * tangentPointEnergyGradient(x, y, z, resolution, perturbation, alpha, beta);
    %tangIentPointEnergy(x,y,z,resolution, alpha, beta)
    

    coeffsx = coeffs(1:J);
    coeffsy = coeffs(J+1:2*J);
    coeffsz = coeffs(2*J + 1:3*J);
    %Deflation
    i = 0;
    while norm(coeffs(:,end-i)) < epsilon
        i = i + 1;
    end
    J = J - i;
    coeffsx = coeffsx(1:J);
    coeffsy = coeffsy(1:J);
    coeffsz = coeffsz(1:J);

    x = chebfun(coeffs(1:J), 'coeffs');
    y = chebfun(coeffs(J+1:2*J), 'coeffs');
    z = chebfun(coeffs(2*J + 1:3*J), 'coeffs');
    plot3(x,y,z,'LineWidth',linewidth);
    xlim([-viewsize,viewsize])
    ylim([-viewsize,viewsize])
    zlim([-viewsize,viewsize])
    
    timeIndex = sprintf("t=%d",t);
    %title(timeIndex)
    %view([-53+90 80])
    drawnow;
    if mod(t, snapshotInterval) == 1
        exportgraphics(gcf, sprintf("%d.png",snapshotNumber),"Resolution",300, "BackgroundColor","none");
        snapshotNumber = snapshotNumber + 1;
    end
end