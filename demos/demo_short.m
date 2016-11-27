% A script to test the MF-GP-UCB framework on synthetic functions.

clear all;
close all;
clc;

randMethods = {'MF-GP-UCB', 'GP-UCB', 'RAND'};
numRandMethods = numel(randMethods);

% PLOT_EVERY_ITER = true;
PLOT_EVERY_ITER = false;
plotColours = {'b', 'g', 'r'};
plotMarkers = {'o-', '+-', '*-'};

numFidels = nan; numDims = nan;
% Pick toy problem --------------------------------------------------------------
toyFuncHandle = @getMFHartFn; toyFuncDesc = 'Hartmann'; numFidels = 4; 
% toyFuncHandle = @getMFHart3Fn; toyFuncDesc = 'Hartmann'; numFidels = 3; 
% toyFuncHandle = @getCurrinExpFn; toyFuncDesc = 'CurrinExp';
% toyFuncHandle = @getBoreHoleFn; toyFuncDesc = 'BoreHole';
% -------------------------------------------------------------------------------

% Experiment parameters
[mff, sff] = toyFuncHandle(numFidels, numDims);
budget = 100 * mff.costs(end);
params.budgetType = 'givenCost';
numFidels = mff.numFidels;
numDims = mff.numDims;
toyFuncDesc = sprintf('%s-%dD', toyFuncDesc, numDims);

% Set problem parameters
bounds = mff.bounds;
numFidels = mff.numFidels;
numDims = mff.numDims;
costs = mff.costs;

% To store the results for the random methods
simRegrets = cell(numRandMethods);
cumCosts = cell(numRandMethods);


fprintf(['(%s, max: %0.4f) \n', ...
  '=======================================================================\n'], ...
  toyFuncDesc, mff.hfMaxVal);

 
for methIter = 1:numRandMethods

  fprintf('\n');

  switch randMethods{methIter}

    case 'MF-GP-UCB'
      params.acqStrategy = 'MF-GP-UCB';
      mfFuncObj = mff;
      [maxVal, maxPt, queries, vals, history] = ...
        mfBO( mfFuncObj, bounds, budget, params);
      mfFuncObj = mff;

    case 'GP-UCB'
      params.acqStrategy = 'GP-UCB';
      mfFuncObj = mff;
      [maxVal, maxPt, queries, vals, history] = ...
        mfBO(mfFuncObj, bounds, budget, params);

    case 'RAND'
      [maxVal, maxPt, queries, vals, history] = ...
        randMaximise(sff, bounds, budget, struct());
      mfFuncObj = sff;

    otherwise
      error('Unknown Method');


  end % End Switch methods --------------------------------------------------------

  [sR, ~, ~, cC, ~] = getSimCumRegrets(mfFuncObj, history, params);
  simRegrets{methIter} = sR;
  cumCosts{methIter} = cC;
  fprintf('Method: %s, Max Val: %.7f\n', randMethods{methIter}, maxVal);

end % End for loop ----------------------------------------------------------------


% Plot the results out -------------------------------------------------------------------
figure;
for methIter = 1:numRandMethods
  plot(cumCosts{methIter}, simRegrets{methIter}, ...
    plotMarkers{methIter}, 'Color', plotColours{methIter});
  hold on;
end
ylabel('$S_T$', 'Interpreter', 'Latex');
xlabel('Cost', 'Interpreter', 'Latex');
titleStr = sprintf('Simple Regret: Costs = %s', mat2str(mff.costs));
title(titleStr);
legend(randMethods);
% ---------------------------------------------------------------------------------

