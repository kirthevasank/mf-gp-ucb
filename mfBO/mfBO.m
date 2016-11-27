function [hfMaxVal, hfMaxPt, boQueries, boVals, history] = mfBO( ...
  mfFunc, bounds, budget, params)
% Performs Bayesian optimisation with multiple fidelity evaluations.
% Inputs:
%   mfFunc: Is an mfFunction object which contains the multi-fidelity functions. If
%     mfFunc is just a function, we will do vanilla BO.
%   bounds: The bounds over which we should optimise the function.
% budget: The budget we have for performing BO.
% params: A struct which (optionally) contains the following fields:
%   - acquisition: which acquisition function to use.
%   - gammas: The thresholds at which we switch from one fidelity to another.
%   - initPts, initVals: A numFidels cell set of initial points and function values.
%   - numInitPts: The number of initialisation points.
%   - diRectParams: Configuration for optimising acquisition function.
%   - bwRange: The range for learning the bandwidth.
%
% Outputs:
%   maxVal, maxPt: The maximum value and point of the highest fidelity evaluations.
%     If no highest fidelity evaluations were performed, it returns nans.
%   boQueries, boVals: numFidel x 1 structs containing the points and function values
%     of each fidelity.
%   history: A struct which has the following fields.
%     hfMaxVal: The maximum value among the highest fidelity evaluations.
%     totalNumQueries: the total number of queries including initialisation.
%     The following are totalNumQueries x 1 vectors.
%       numIters: The number of iterations we ran the method. 
%       evalFidels: The fidelity at which each query was carried out.
%       evalVals: The function value of those values.
%       evalPts: The points at which each query was evaluated.
%       cumRewards: The cumRewards at each time step

  % CONSTANTS 
  NUM_ITERS_PER_PARAM_RELEARN = 25;
  NUM_ITERS_PER_INFO_PRINTOUT = 10;
  NUM_THRESHOLD_EXCEEDS = inf;
  GAMMA_INC_COEFF = 5;
  ZETA_INC_COEFF = 2;
  NUM_INIT_ZETA_CHECK_TIMES = 0;

  if isa(mfFunc, 'function_handle')
    mfFunc = mfFunction({mfFunc}, bounds, 1);
  end

  % Prelims
  numDims = mfFunc.numDims;
  numFidels = mfFunc.numFidels;

  % Preprocessing and setting up hyper- parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initTime = cputime;
  params.initTime = initTime;
  params = mfboPreProcessParams(params, mfFunc, bounds);
  % Initialisation
  numAllInitFuncEvals = params.numAllInitFuncEvals;
  numIndFuncEvals = params.numIndInitFuncEvals;
  incCost = params.initCost;

  % Set the parameters of the initial GP.
  gpBw = params.initGPBw;
  gpScale = params.initGPScale;
  funcHs = params.initGPFuncHs;
  % Set the hyper parameters of the GP.
  gpHPs.bwRange = params.bwRange;
  gpHPs.scaleRange = params.scaleRange;
  gpHPs.meanFuncs = params.gpMeanFuncs;
  gpHPs.diRectOpts = params.gpDirectOpts;
  gpHPs.noiseVars = params.gpNoiseVars;
  gpHPs.multGPLearnStrategy = params.multGPLearnStrategy;

  % Set up initial history fields
  evalVals = params.initEvalVals;
  evalPts = params.initEvalPts;
  cumRewards = params.initCumRewards;
  currCumReward = cumRewards(end);
  evalFidels = params.initEvalFidels;
  evalTimes = params.initEvalTimes;
  allNextVals = [];

  % Now do Bayesian Optimisation
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initialise some variables.
  currMaxFidel = 1; % The maximum fidelity queried by BO up to now.
  threshExceededCounter = 0; % To keep track of whether we are over-exploiting.
  gammaThreshExcCounter = 0; % To keep track of whether we are stuck at low fidels.
  boIter = 0; % Keeps track of the current iteration of BO after initialisation.
  nextTime = evalTimes(end);
  boQueries = params.initPts;
  boVals = params.initVals;
  [hfMaxVal, hfMaxIdx] = max(boVals{numFidels});
  if isempty(hfMaxVal), hfMaxPt = -inf; hfMaxVal = -inf;
  else, hfMaxPt = boQueries{numFidels}(hfMaxIdx);
  end

  % Set the spent Budget based on the given costs or time.
  if strcmp(params.budgetType, 'givenCost'), spentBudget = incCost;
  elseif strcmp(params.budgetType, 'time'), spentBudget = evalTimes(end);
  else, error('unknown budget type');
  end

  fprintf(['Running Bayesian Optimisation: %s\n', ...
    '------------------------------------------------------\n'], params.acqStrategy);

  % The Main loop which does BO. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while (spentBudget < budget)

    % Prelims -----------------------------------------------------------
    boIter = boIter + 1;
    totalNumQueries = boIter + numAllInitFuncEvals;
    params.currMaxFidel = currMaxFidel;
    % Print details out
    if mod(boIter, NUM_ITERS_PER_INFO_PRINTOUT) == 0
      fprintf(['MF-BO #%03d, cost: %d/%d: Q:%s, hfMaxVal: %0.4f\n'], ...
        boIter, round(incCost), round(budget), mat2str(numIndFuncEvals'), hfMaxVal);
    end


    % Update the GP Hyper-parameters ------------------------------------------------
    % -------------------------------------------------------------------------------
    if (mod(boIter, NUM_ITERS_PER_PARAM_RELEARN) == 1 || ...
       threshExceededCounter >= NUM_THRESHOLD_EXCEEDS )
      
      if (threshExceededCounter >= NUM_THRESHOLD_EXCEEDS)
        fprintf('  ***  Re-learning GP hyper-params since thresholds exceeded !\n');
        gpHPs.bwRange(2) = max(1.1*gpHPs.bwRange(1), 0.9*gpHPs.bwRange(2));
        threshExceededCounter = 0;
      end

      % Learn the GP
      [funcHs, learnedHPs] = multipleGPRegressionML(boQueries, boVals, gpHPs);
      gpBw = learnedHPs.bw; 
      gpScale = learnedHPs.scale;

    else
      % Use the same hyper parameters for the GP but with updated observations.
      % Only update the GP of the last fidelities played.
      if ~isempty(allNextFidels)
        for i = 1:numel(allNextFidels)
          fidel = allNextFidels(i);
          [~, ~, ~, currFuncH] = GPRegression([], boQueries{fidel}, ...
            boVals{fidel}, gpBw(fidel), gpScale(fidel), ...
            gpHPs.meanFuncs{fidel}, gpHPs.noiseVars(fidel));
          funcHs{fidel} = currFuncH;
        end
      end
    end
    % -------------------------------------------------------------------------------

    % Determine the next point ------------------------------------------
    params.noiseVars = learnedHPs.noiseVars; % need this for MFSKO
    params.costs = mfFunc.costs;
    [nextPt, nextFidel, ~, nextPtStd] = getNextQuery(boIter, funcHs, bounds, ...
      boQueries, boVals, params);

    % Perturb Query if the selected point is already very close to an existing
    % point, perturb the query a little bit.
    if sqrt(sum(bsxfun(@minus, boQueries{nextFidel}, nextPt).^2, 2)) < ...
       1e-4 * max(gpBw)
      fprintf(' ******** Perturbing point as its too close to an existing point\n');
      nextPt = perturbPt(nextPt, bounds, max(gpBw));
    end

    % -------------------------------------------------------------------------------
    nextPtVal = mfFunc.evalAtFidel(nextFidel, nextPt); % EVALUATE THE FUNCTION
    % -------------------------------------------------------------------------------

    % Some Book-keeping -------------------------------------------------
    if (nextFidel == numFidels) & (nextPtVal > hfMaxVal)
      hfMaxVal = nextPtVal;
      hfMaxPt = nextPt;
    end
    nextTime = cputime - initTime;

    % Keep track of the maximum fidelity played so far.    
    if nextFidel > currMaxFidel, currMaxFidel = nextFidel;
    end
    % Check if nextPtStd is too small
    if nextPtStd < params.optPtStdThresholdCoeff * max(gpScale);
      threshExceededCounter = threshExceededCounter + 1;
    else, threshExceededCounter = 0; % reset the counter
    end

    % Keep track of gamma and zeta
    % ===============================================================================
    % Keep track of whether we are increasing the fidelity.
%     params.fixGamma, params.fixZeta,
    if numFidels > 1 & ...
      (strcmp(params.acqStrategy, 'MF-GP-UCB') | ...
        strcmp(params.acqStrategy, 'Illus-MF-GP-UCB' ) ) & ~params.fixGamma
      if (nextFidel == numFidels) | (nextFidel > currMaxFidel)
        gammaThreshExcCounter = 0;      
      else % if nextFidel <= currMaxFidel,
        gammaThreshExcCounter = gammaThreshExcCounter + 1;
        if gammaThreshExcCounter >= ...
                  mfFunc.costs(nextFidel+1)/mfFunc.costs(nextFidel)
          % increase gamma and reset the counter.
          params.gammas(nextFidel) = ...
            params.gammas(nextFidel) * GAMMA_INC_COEFF;   
          params.gammas = max(params.gammas, params.gammas(nextFidel));
%           fprintf('Increased gammas(%d) to %0.4f\n', ...
%             nextFidel, params.gammas(nextFidel));
          gammaThreshExcCounter = 0;
        end
      end
    end

    allNextPts = nextPt;
    allNextFidels = nextFidel;
    allNextVals = nextPtVal;
    allNextTimes = nextTime;
    % Keep track of whether zetas are correct
    if nextFidel > 1 & ...
      (strcmp(params.acqStrategy, 'MF-GP-UCB') | ...
        strcmp(params.acqStrategy, 'Illus-MF-GP-UCB' ) ) & ~params.fixZeta

      diffZetas = [-diff(params.zetas); params.zetas(end)];
      diffEst = abs(nextPtVal - funcHs{nextFidel-1}(nextPt));

      if (numIndFuncEvals(nextFidel) < NUM_INIT_ZETA_CHECK_TIMES) | ...
        (diffEst > 0.9 * diffZetas(nextFidel-1))

%         fprintf('Checking for zeta violation at %s.\n', mat2str(nextPt));
        % First query at the lower fidelity and check that this is the case
        lfVal = mfFunc.evalAtFidel(nextFidel-1, nextPt);
        trueDiff = nextPtVal - lfVal;
        allNextPts = [allNextPts; nextPt];
        allNextVals = [allNextVals; lfVal];
        allNextFidels = [allNextFidels; (nextFidel-1)];
        nextTime = cputime - initTime;
        allNextTimes = [allNextTimes; nextTime];
        
        % If indeed different you need to update eta
        if trueDiff > 0.9 * diffZetas(nextFidel-1)
          trueDiff = ceil(1e4 * trueDiff)/1e4;
          diffZetas(nextFidel-1) = ZETA_INC_COEFF*trueDiff;
          diffZetas = max(diffZetas, diffZetas(nextFidel-1));
          newZetas = flipud( cumsum(flipud(diffZetas)) );
          params.zetas = newZetas;
%           fprintf('Increased Zetas: %s\n', mat2str(params.zetas));
        end

      end
    end

    % gamma/ zeta ends here ===================================================

    for q = 1:size(allNextPts,1)
      f = allNextFidels(q);
      p = allNextPts(q, :);
      v = allNextVals(q);
      t = allNextTimes(q);

      boQueries{f} = [boQueries{f}; p];
      boVals{f} = [boVals{f}; v];
      currCumReward = currCumReward + mfFunc.costs(f) * v;

      % Store these for history
      evalVals = [evalVals; v];
      evalPts = [evalPts; p];
      cumRewards = [cumRewards; currCumReward];
      evalFidels = [evalFidels; f];
      evalTimes = [evalTimes; t];

      % Update the cost
      numIndFuncEvals(f) = numIndFuncEvals(f) + 1;
      incCost = incCost + mfFunc.costs(f);
    end

    % Set the spent Budget based on the given costs or time.
    if strcmp(params.budgetType, 'givenCost'), spentBudget = incCost;
    elseif strcmp(params.budgetType, 'time'), spentBudget = evalTimes(end);
    else, error('unknown budget type');
    end

  end
  % BO Loop ends here
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  history.totalNumQueries = totalNumQueries;
  history.evalVals = evalVals;
  history.evalPts = evalPts;
  history.cumRewards = cumRewards;
  history.evalFidels = evalFidels;
  history.evalTimes = evalTimes;
  history.hfMaxVal = hfMaxVal;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% A wrapper to maximise the relevant acquisition function and determine the fidelity.
function [nextPt, nextFidel, nextPtMean, nextPtStd, nextPtAcq] = ...
  getNextQuery(boIter, funcHs, bounds, boQueries, boVals, params)

  % Set DiRect Parameters
  numDims = size(bounds, 1);
  numFidels = numel(funcHs);
  if ~isfield(params, 'diRectParams')
    diRectParams.maxevals = ceil(7 * min(5,numDims)^2 * sqrt(min(boIter, 1000)));
    diRectParams.maxits = inf;
%     fprintf('t = %d, diREctEvals: %d\n', boIter, diRectParams.maxevals);
  end
  params.diRectParams = diRectParams;

  switch params.acqStrategy

    case 'MF-GP-UCB'
      [nextPt, nextFidel, nextPtAcq] = ...
        strategyMFGPUCB(boIter, funcHs, bounds, params);

    case 'GP-UCB'
      params.zetas = zeros(0,1);
      [nextPt, ~, nextPtAcq] = strategyMFGPUCB(...
        boIter, {funcHs{numFidels}}, bounds, params);
      nextFidel = numFidels;

    case 'GP-EI'
      currMax = max(boVals{numFidels});
      if isempty(currMax), currMax = 0; end
      [nextPt, nextPtAcq] = strategyGPEI(funcHs{numFidels}, bounds, ...
        currMax, params);
      nextFidel = numFidels; 

    case 'GP-PI'
      currMax = max(boVals{numFidels});
      if isempty(currMax), currMax = 0; end
      [nextPt, nextPtAcq] = strategyGPPI(funcHs{numFidels}, bounds, ...
        currMax, boIter, params);
      nextFidel = numFidels; 

    case 'Illus-MF-GP-UCB'
      [nextPt, nextFidel, nextPtAcq] = ...
        strategyIllusMFGPUCB(boIter, funcHs, bounds, boQueries, params);

    case 'MFSKO'
      currMax = max(boVals{numFidels});
      if isempty(currMax), currMax = 0; end
      for f = 1:numFidels
        params.numFidelEvals(f) = size(boQueries{f}, 1);
      end
      [nextPt, nextFidel, nextPtAcq] = ...
        strategyMFSKO(boIter, funcHs, bounds, currMax, params);

    otherwise
      error('Unknown Acquisition');
  end

  % Obtain the posterior mean and std of the next point
  [nextPtMean, ~, nextPtStd] = funcHs{nextFidel}(nextPt);

end


function x = perturbPt(x, bounds, scale)
  numDims = size(bounds, 1);
  perturbation =  scale * ...
    max(bounds(2) - bounds(1)) * (rand(1, numDims) - 0.5);
  x = x + perturbation;
  x = projectToRectangle(x', bounds)';
end

