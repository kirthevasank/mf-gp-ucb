function [hfMaxVals, cumRewards, cumAvgRewards, cumCosts]= ...
  getSimCumRewards(mfFunc, history, params)

  numFidels = mfFunc.numFidels;

  if ~exist('params', 'var'), params = struct();
  end
  if ~isfield(params, 'cumRewardHighestFidelity')
    params.cumRewardHighestFidelity = true;
  end

  if params.cumRewardHighestFidelity
    instRewards = mfFunc.evalAtFidel( numFidels, history.evalPts);
  else
    instRewards = history.evalVals;
  end

  % Cumulative costs
  if strcmp(params.budgetType, 'givenCost')
    evalCosts = mfFunc.costs(history.evalFidels);
    cumCosts = cumsum(evalCosts);
  elseif strcmp(params.budgetType, 'time')
    cumCosts = history.evalTimes;
  else
    error('Unknown Cost type!');
  end

  % First do the cumulative rewards
  weightedInstRewards = instRewards .* evalCosts;
  cumRewards = cumsum(weightedInstRewards);
  cumAvgRewards = cumRewards ./ cumCosts;

  % Instantaneous Rewards
  hfQueryVals = history.evalVals .* (history.evalFidels == numFidels);
  hfQueryVals(hfQueryVals == 0) = -inf;
  hfMaxVals = cumMax(hfQueryVals);

end

