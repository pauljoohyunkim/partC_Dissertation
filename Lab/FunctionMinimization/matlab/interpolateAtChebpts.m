function interpolantF = interpolateAtChebpts(yVals)
%INTERPOLATEATCHEBPTS Outputs chebfun from interpolation at chebpts of
%yVals
    N = numel(yVals);
    interpolantF = interp1(chebpts(N), yVals, domain(-1,1));
end

