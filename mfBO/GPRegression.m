function [teMean, teK, teStd, funcH] = ...
  GPRegression(Xte, Xtr, Ytr, bw, scale, meanFunc, noiseVar)
% Outputs the posterior mean (nTe x 1), standard deviations (nTe x 1) and covariance
% matrix (nTe x nTe) of the test data Xte. In addition, returns a function Handle
% for the GPs.

  numTrData = size(Xtr, 1);

  if numTrData == 0
    Ytr = zeros(0, 1);
  end

  Ktrtr = sqExpKernel(bw, scale, Xtr, Xtr) + diag(noiseVar * ones(numTrData, 1));
  Y_ = Ytr - meanFunc(Xtr);
  L = stableCholesky(Ktrtr);
  alpha = L' \ (L \ Y_);

  % obtain the function handle
  funcH = @(X) GPComputeOutputs(X, Xtr, L, alpha, bw, scale, meanFunc);
  % Compute outputs for the test data
  if ~isempty(Xte)
    [teMean, teK, teStd] = funcH(Xte);
  else
    teMean = []; teK = []; teStd = [];
  end

end


function [yMu, yK, yStd] = GPComputeOutputs(Xte, Xtr, L, alpha, bw, scale, meanFunc)

  meanXte = meanFunc(Xte);
  Ktetr = sqExpKernel(bw, scale, Xte, Xtr);
  Ktete = sqExpKernel(bw, scale, Xte, Xte);

  % Predictive Mean
  yMu = meanXte + Ktetr * alpha;
  % Predictive Variance
  V = L \ (Ktetr)';
  yK = Ktete - V'*V;
  yStd = sqrt(real(diag(yK)));

end

