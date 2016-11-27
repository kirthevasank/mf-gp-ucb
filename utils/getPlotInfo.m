function [MR, SE, eMR, eSE, expGridVals] = getPlotInfo(methRegrets, methCumCosts, ...
  gridPts, errBarIdxs, params)

  % prelims
  numGridPts = numel(gridPts);
  numExperiments = numel(methRegrets);
  % range for outliers
  if ~exist('params', 'var'), params = struct; end
  if ~isfield(params, 'outlierLowIdx'), params.outlierLowIdx = 1; end
  if ~isfield(params, 'outlierHighIdx'), params.outlierHighIdx = numExperiments; end
  inlierExpIdxs = params.outlierLowIdx:params.outlierHighIdx;
  if numExperiments == 1, inlierExpIdxs = 1;
  end
  numInlierExperiments = params.outlierHighIdx - params.outlierLowIdx + 1;

  expGridVals = zeros(numExperiments, numGridPts);
%   methCumCosts, methRegrets, size(gridPts),

  for expIter = 1:numExperiments
    numPts = size(methCumCosts{expIter, 1});
    addErr = 1e-9 * (1:numPts)';
    cumCosts = methCumCosts{expIter,1} + addErr;
    expGridVals(expIter, :) = interp1( ...
      cumCosts, methRegrets{expIter, 1}, gridPts);
  end

  % Sort them in case we need to remove outliers
  sortedExpGridVals = sort(expGridVals, 1);
  inlierGridVals = sortedExpGridVals(inlierExpIdxs, :);
  MR = nanmean(inlierGridVals, 1);
  if numExperiments == 1
    SE = nan * ones(1, numGridPts);
  else
    SE = nanstd(inlierGridVals) / sqrt(numInlierExperiments);
  end

  eMR = MR(errBarIdxs);
  eSE = SE(errBarIdxs);

end

