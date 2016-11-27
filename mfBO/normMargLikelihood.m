function nlml = normMargLikelihood(bw, scale, X, Y, meanFunc, noiseVar)
% Computes the normalised Marginal Likelihood for a single GP
  numData = size(X, 1);
  Ky = sqExpKernel(bw, scale, X, X) + diag(noiseVar * ones(numData, 1));
  Y_ = Y - meanFunc(X);
  L = stableCholesky(Ky);
  alpha = L' \ (L \ Y_);
  nlml = -0.5 * Y_' * alpha - sum( log(diag(L)) ) -  log(2*pi)*numData/2;
end

