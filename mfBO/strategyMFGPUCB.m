function [nextPt, nextFidel, nextPtAcq] = ...
  strategyMFGPUCB(t, funcHs, bounds, params)

  % Preliminaries and preprocessing
  numFidels = numel(funcHs);
  numDims = size(bounds, 1);

  % First maximise the MF-GP-UCB acquisition function.
  acquisition = @(arg) acqMFGPUCB(arg, funcHs, t, params.zetas);
  [nextPtAcq, nextPt] = ...
    diRectWrap( acquisition, bounds, params.diRectParams);
  [nextPtAcq, uncerts] = acquisition(nextPt);

  % Now determine the fidelity
  nextFidel = 1;
  for i = 2:numFidels
    threshold = params.gammas(i-1);
%     threshold = params.gammas(i-1) * ...
%       (params.rangeY * params.costs(i-1)/params.costs(end) * 1/t)^(1/(numDims+2));
    if uncerts(i-1) < threshold, %& i <= (params.currMaxFidel+1), 
      nextFidel = i;
    else, break; end
  end

end

