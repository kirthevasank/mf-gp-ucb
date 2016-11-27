function [acq, uncerts] = acqMFGPUCB(x, funcHs, t, zetas)

  numDims = size(x, 2);
  numFidels = numel(funcHs);
  if ~isa(funcHs, 'cell') 
    funcHs = {funcHs};
  end

%   beta_t = 0.2 * numDims * log(2*numDims*numFidels*t);
  beta_t = 0.2 * numDims * log(2*numDims*t);
  uncerts = zeros(numFidels, 1);
  augZetas = [zetas; 0];

  indUCBs = zeros(numFidels, 1);
  for i = 1:numFidels
    [mu, ~, sigma] = funcHs{i}(x);
    uncerts(i) = sqrt(beta_t) * sigma;
    indUCBs(i) = mu + uncerts(i) + augZetas(i);
  end
  acq = min(indUCBs);

end

