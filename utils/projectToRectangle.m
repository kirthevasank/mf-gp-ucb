function pt = projectBackToRectangle(pt, bounds)
% pt is a dx1 vector. Bounds is a dx2 matrix with the lower bounds on the first
% column and the upper bound on the second column.
% If bounds is not given, we will take [0,1]^d

  numDims = size(pt, 1);
  if ~exist('bounds', 'var')
    bounds = repmat([0,1], numDims, 1);
  end

  belowLBs = (pt - bounds(:,1)) < 0;
  pt(belowLBs) = bounds(belowLBs,1);

  aboveUBs = (pt - bounds(:,2)) > 0;
  pt(aboveUBs) = bounds(aboveUBs,2);

end
