function [L, diagPower] = stableCholesky(K)
% Sometimes nominally psd matrices are not psd due to numerical issues. By
% adding a small value to the diagonal we can make it psd. This is what this
% function does.
% Use this iff you know that K should be psd. We do not check for errors

  % First see if the usual Cholesky decomposition works
  try
    L = chol(K, 'lower');
    return;
  catch err
  end
  

  diagPower = min( ceil(log10(abs(min(diag(K)))))-1, -11);
  if ~(abs(diagPower) < inf)
    diagPower = -10;
  end

  % Now keep trying until Cholesky is successful
  success = false;
  K = K + (10^diagPower) * eye(size(K));
  while ~success
    try
      L = chol(K, 'lower');
      success = true;
    catch err
      if diagPower > 1e-3,
        fprintf('CHOL failed with diagPower = %d\n', diagPower);
        K, 
        error('Failed !');
      end
      diagPower = diagPower + 1; 
      K = K + (10^diagPower) * eye(size(K));
    end
  end

end

