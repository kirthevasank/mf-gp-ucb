function [mff, sff] = getCurrinExpFn(numFidels, numDims)
% numFidels isn't really needed now. Keeping it since I might want to change it in 
% the future.

  f2 = @(x) (1 - exp(-0.5./x(:,2))) .* ...
        (2300* x(:,1).^3 + 1900* x(:,1).^2 + 2092*x(:,1) + 60) ./ ...
        (100* x(:,1).^3 + 500* x(:,1).^2 + 4*x(:,1) + 20);
  f1 = @(x) -f2(x);

  numFidels = 2;
  funcHs = {f1; f2};
  hfMaxPt = [];
  hfMaxVal = 13.7986850;
  bounds = repmat([0 1], 2, 1);
  costs = (10.^(0:(numFidels-1)))';

  mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
  sff = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
                   hfMaxPt, hfMaxVal);

end

