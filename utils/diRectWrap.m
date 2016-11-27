function [maxVal, maxPt, history, queries, queryVals] = ...
  diRectWrap(func,bounds,opts,varargin)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

  Problem.f = @(t, varargin) -func(t');
  [ret_minval, final_xatmin, history, queries, queryVals] = ...
    diRect(Problem, bounds, opts, varargin);
  queryVals = -queryVals;
  history(:,3) = -history(:,3);

  [maxVal, maxIdx] = max(queryVals);
  maxPt = queries(maxIdx, :);

end

