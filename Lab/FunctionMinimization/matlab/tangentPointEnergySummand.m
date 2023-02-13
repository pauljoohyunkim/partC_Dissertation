function energyfrag = tangentPointEnergySummand(xi, xip1, xj, xjp1, alpha, beta)
%TANGENTPOINTENERGYSUMMAND Computes the summand used for quadrature for
%tangent point energy
    energyfrag = 0;
    xI = xip1 - xi;
    xJ = xjp1 - xj;
    xINorm = norm(xI);
    xJNorm = norm(xJ);
    TI = xI / xINorm;
    energyfrag = energyfrag + (norm(cross(TI, xi - xj))^alpha) / (norm(xi - xj)^beta);
    energyfrag = energyfrag + (norm(cross(TI, xi - xjp1))^alpha) / (norm(xi - xjp1)^beta);
    energyfrag = energyfrag + (norm(cross(TI, xip1 - xj))^alpha) / (norm(xip1 - xj)^beta);
    energyfrag = energyfrag + (norm(cross(TI, xip1 - xjp1))^alpha) / (norm(xip1 - xjp1)^beta);
    energyfrag = energyfrag / 4 * xINorm * xJNorm;
end

