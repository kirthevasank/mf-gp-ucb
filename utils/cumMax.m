function ret = cumMax(A)

  if size(A, 1) > 1, 
    A = A';
    isColVec = false;
  else,
    isColVec = true;
  end

  ret = max(bsxfun(@(~,y)([A(1:y),-inf(1,numel(A)-y)]'), A', 1:numel(A)));

  if ~isColVec, ret = ret';
  end

end

