function Gamma = getUnNormParams(Theta, bounds)
  Alpha = bsxfun(@times, Theta, (bounds(:,2) - bounds(:,1))'); 
  Gamma = bsxfun(@plus, Alpha, bounds(:,1)');
end

