function energy = tangentPointEnergy(x, y, z, resolution, alpha, beta)
%TANGENTPOINTENERGY Computes the tangent point energy. (x,y,z are chebfuns)
    energy = 0;
    %t = linspace(-1, 1, resolution);
    t = chebpts(resolution);
    xpoints = x(t);
    ypoints = y(t);
    zpoints = z(t);
    parfor i = 2:resolution-1
        for j = 2:resolution-1
            if abs(i - j) > 1 && abs(i - j + resolution) > 1 && abs(i - j - resolution) > 1
                xi = [xpoints(i); ypoints(i); zpoints(i)];
                xip1 = [xpoints(i+1); ypoints(i+1); zpoints(i+1)];
                xj = [xpoints(j); ypoints(j); zpoints(j)];
                xjp1 = [xpoints(j+1); ypoints(j+1); zpoints(j+1)];
                energy = energy + tangentPointEnergySummand(xi, xip1, xj, xjp1, alpha, beta);
            end
        end
    end
end

