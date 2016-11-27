function Theta = mapToCube(Gamma, bounds)
  Beta = bsxfun(@minus, Gamma, bounds(:,1)');
  Theta = bsxfun(@rdivide, Beta, (bounds(:,2) - bounds(:,1))');
end

