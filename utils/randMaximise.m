function [maxVal, maxPt, queries, vals, history] =  ...
  randMaximise(sff, bounds, budget, params)

  numDims = size(bounds, 1);
  if ~isfield(params, 'budgetType'), params.budgetType = 'givenCost';
  end

  queries = zeros(0, numDims);
  vals = zeros(0, 1);
  spentBudget = 0;

  initTime = cputime;
  currCost = 0;
  evalTimes = zeros(0, 1);
  spentBudget = 0;

  while spentBudget < budget
    currQuery = rand(1, numDims) .* (bounds(:,2)-bounds(:,1))' + bounds(:,1)';
    currVal = sff.evalAtFidel(1, currQuery);
    currTime = cputime - initTime;
    currCost = currCost + sff.costs(end);

    queries = [queries; currQuery];
    vals = [vals; currVal];
    evalTimes = [evalTimes; currTime];
    
    if strcmp(params.budgetType, 'givenCost'), spentBudget = currCost;  
    elseif strcmp(params.budgetType, 'time'), spentBudget = currTime;
    else, error('unknown budget type');
    end

  end

  numQueries = size(queries, 1);
  [maxVal, maxIdx] = max(vals);
  maxPt = queries(maxIdx);

%   startTime = cputime;
%   vals = sff.evalAtFidel(1, queries);
%   endTime = cputime - startTime;
%   evalTimes = linspace(0, endTime, numQueries + 1); evalTimes = evalTimes(2:end);

  history.totalNumQueries = numQueries;
  history.evalVals = vals;
  history.evalPts = queries;
  history.evalFidels = ones(numQueries, 1);
  history.hfMaxVal = maxVal;
  history.evalTimes = evalTimes;

  fprintf('Random: #evals: %d, time: %0.3f, maxVal: %0.6f\n', ...
    numQueries, currTime, maxVal);
end

