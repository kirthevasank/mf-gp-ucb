function [mff, sff] = getMFHart3Fn(numFidels, numDims)

  alpha0 = [1.0, 1.2, 3.0, 3.2]';
  delta = [0.01 -0.01 -0.1 0.1]';

  funcHs = cell(numFidels, 1);
  costs = round(10.^(0:(numFidels-1)))';
  bounds = repmat([0, 1], 3, 1);

  for i = 1:numFidels
    alphai = alpha0 + delta * (numFidels-i);
    funcHs{i} = @(x) hart3Wrap(x, alphai);
  end

  % Highest fidelity maximum.
  hfMaxPt = [0.114614, 0.555649, 0.852547];
  hfMaxVal = funcHs{numFidels}(hfMaxPt);

  % A multifidelity function and a single fidelity function 
  mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
  sff = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
                   hfMaxPt, hfMaxVal);

end

function vals = hart3Wrap(X, alphai)
  numPts = size(X, 1);
  vals = zeros(numPts, 1);
  for i = 1:numPts
    vals(i) = hart3( X(i,:), alphai);
  end
end


function [y] = hart3(xx, alpha)

A = [3.0, 10, 30;
     0.1, 10, 35;
     3.0, 10, 30;
     0.1, 10, 35];
P = 10^(-4) * [3689, 1170, 2673;
               4699, 4387, 7470;
               1091, 8732, 5547;
               381, 5743, 8828];

outer = 0;
for ii = 1:4
  inner = 0;
  for jj = 1:3
    xj = xx(jj);
    Aij = A(ii, jj);
    Pij = P(ii, jj);
    inner = inner + Aij*(xj-Pij)^2;
  end
  new = alpha(ii) * exp(-inner);
  outer = outer + new;
end

y = outer;

end

