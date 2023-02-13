function gradient = tangentPointEnergyGradient(x, y, z, resolution, perturbation, alpha, beta)
%TANGENTPOINTENERGYGRADIENT Computes gradient by central difference
    J = length(x);
    coeffs = [chebcoeffs(x); chebcoeffs(y); chebcoeffs(z)];
    gradient = zeros(numel(coeffs), 1);
    for i=1:numel(coeffs)
        % Store quantity before perturbation
        temp = coeffs(i);
        % Perturb
        coeffs(i) = coeffs(i) + perturbation;
        xcoeffs = coeffs(1:J);
        ycoeffs = coeffs(J+1:2*J);
        zcoeffs = coeffs(2*J + 1: 3*J);
        xf = chebfun(xcoeffs, 'coeffs');
        yf = chebfun(ycoeffs, 'coeffs');
        zf = chebfun(zcoeffs, 'coeffs');
        energyP = tangentPointEnergy(xf, yf, zf, resolution, alpha, beta);

        % Restore
        coeffs(i) = temp;

        % Perturb the other way around
        coeffs(i) = coeffs(i) - perturbation;
        xcoeffs = coeffs(1:J);
        ycoeffs = coeffs(J+1:2*J);
        zcoeffs = coeffs(2*J + 1: 3*J);
        xf = chebfun(xcoeffs, 'coeffs');
        yf = chebfun(ycoeffs, 'coeffs');
        zf = chebfun(zcoeffs, 'coeffs');
        energyN = tangentPointEnergy(xf, yf, zf, resolution, alpha, beta);

        % Restore
        coeffs(i) = temp;
        
        gradient(i) = (energyP - energyN) / 2;
    end
end

