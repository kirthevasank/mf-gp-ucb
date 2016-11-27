function [maxVal, maxPt, queries, vals, history] = sffDiRectWrap(sff, bounds, budget)

  func = @(x) sff.evalAtFidel(1, x);
  opts.maxevals = ceil(budget/sff.costs(1));
  opts.showits = true;

  initTime = cputime;
  [maxVal, maxPt, ~, queries, vals] = diRectWrap(func, bounds, opts);
  endTime = cputime - initTime;
  numQueries = size(queries, 1);
  evalTimes = linspace(0, endTime, numQueries + 1); evalTimes = evalTimes(2:end);

  history.totalNumQueries = numQueries;
  history.evalVals = vals;
  history.evalPts = queries;
  history.evalFidels = ones(numQueries, 1);
  history.hfMaxVal = maxVal;
  history.evalTimes = linspace(0, endTime, numQueries);

end

