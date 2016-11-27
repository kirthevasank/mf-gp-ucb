function [mff, sff] = get1D2ModalMFF(numFidels)
% Returns a function of 2 modes in 1D

  FUNC_LB = -700;
  bounds = [0, 1];
  gaussVar = 0.005;
  centres = [0.25; 0.7];

  probs = [0.15; 0.85];
  delta = 10*[0.01; -0.01];
  cDelta = 10*[-0.01; 0.01];
  bDelta = -1.0;

%   % NO difference
%   delta = 0*[0.01; -0.01];
%   cDelta = 0*[-0.01; 0.01];
%   bDelta = -0;

  funcs = cell(numFidels, 1);
  for f = 1:numFidels
    fProbs = probs + delta * (numFidels - f);
    fCentres = centres + delta * (numFidels-f);
    fBias = bDelta * (numFidels -f);
    funcs{f} = @(t) max( - 700, fBias + log( ...
        fProbs(1) * mvnpdf(t, fCentres(1), gaussVar) + ...
        fProbs(2) * mvnpdf(t, fCentres(2), gaussVar) ) );
  end

  [~,maxCentreIdx] = max(probs);
  hfMaxPt = centres(maxCentreIdx, :);
  hfMaxVal = funcs{numFidels}(hfMaxPt);

  costs = (10.^(0:(numFidels-1)))';

  % Determine eta
  N = 10000;
  th = linspace(0,1,N)';
  F = zeros(N, numFidels);
  for f = 1:numFidels
    F(:,f) = funcs{f}(th); 
  end
  eta = max(abs(F(:,1) - F(:,2)));
  for f = 2:(numFidels-1)
    eta = max( eta, max( abs(F(:,f) - F(:,f+1)) ) );
  end

  mff = mfFunction(funcs, bounds, costs, eta, hfMaxPt, hfMaxVal);
  sff = mfFunction({funcs{numFidels}}, bounds, costs(numFidels), eta, ...
                    hfMaxPt, hfMaxVal);

  if hfMaxVal < FUNC_LB
    error('Decrease Lower Bound. MaxVal = %.4f\n', hfMaxVal);
  end

end

