function K = sqExpKernel(bw, scale, X, Y)
% Returns the Kernel Matrix for a Gaussian Kernel of bandwidth h.
% X is an nxd matrix. K is an nxn matrix.
% If Y is nonempty then returns the gaussian kernel for XxY

  if ~exist('Y', 'var') 
    Y = X;
  end
  D = distSquaredGP(X, Y);
  K = scale * exp( -D / (2*bw^2) );
end

