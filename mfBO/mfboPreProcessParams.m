function params = mfboProcessParams(params, mfFunc, bounds)
% Check for parameters and if not use default values.
% Prelims
  numDims = size(bounds, 1);
  numFidels = mfFunc.numFidels;

%   INIT_STRATEGY = 'ComplexMFVersion';
  INIT_STRATEGY = 'SimpleVersion';

  % Initialisation
  if ~isfield(params, 'initPts') | isempty(params.initPts)

    if isfield(params, 'initBudget'), initBudget = params.initBudget;
    else, initBudget = 10 * mfFunc.costs(end);
    end

    % Create these structs for initialisation;
    initPts = repmat({zeros(0,numDims)}, numFidels, 1);
    initVals = repmat({zeros(0,1)}, numFidels, 1);
    initTimes = repmat({zeros(0,1)}, numFidels, 1);

    switch INIT_STRATEGY

      case 'ComplexMFVersion'

        if numFidels > 1 & strcmp(params.acqStrategy, 'MF-GP-UCB')
          numF2InitPts = min(10, floor(0.5*initBudget/mfFunc.costs(2)));
          remBudget = initBudget - numF2InitPts*mfFunc.costs(2);
          numF1InitPts = min( ...
            [floor(remBudget/mfFunc.costs(1)), numDims*numFidels*5, 200] );
          numF2InitPts = min(numF2InitPts, numF1InitPts);
        elseif numFidels > 1
          numF1InitPts = ...
            min([200, 5*numFidels*numDims, floor(0.5*initBudget/mfFunc.costs(1))]);
        else
          numF1InitPts = floor(initBudget/mfFunc.costs(end));
        end

        F1InitPts = boGetInitPts(bounds, numF1InitPts);
        F1InitTimes = zeros(numF1InitPts, 1);
        F1InitVals = zeros(numF1InitPts, 1);
        for i = 1:numF1InitPts
          F1InitVals(i) = mfFunc.evalAtFidel(1, F1InitPts(i,:));
          F1InitTimes(i) = cputime - params.initTime;
        end
        F1InitVals = mfFunc.evalAtFidel(1, F1InitPts);
        
        initPts{1} = F1InitPts;
        initVals{1} = F1InitVals;
        initTimes{1} = F1InitTimes;

        % Now for MF-GP-UCB
        if numFidels > 1 & strcmp(params.acqStrategy, 'MF-GP-UCB')
          F2InitPts = F1InitPts(1:numF2InitPts, :);
          F2InitVals = zeros(numF2InitPts, 1);
          F2InitTimes = zeros(numF2InitPts, 1);
          for i = 1:numF2InitPts
            F2InitVals(i) = mfFunc.evalAtFidel(2, F2InitPts(i,:));
            F2InitTimes(i) = cputime - params.initTime;
          end
          initPts{2} = F2InitPts;
          initVals{2} = F2InitVals;
          initTimes{2} = F2InitTimes;
          maxF1F2Diff = max(0, max( F2InitVals - F1InitVals(1:numF2InitPts) ));
        end
      % ComplexMFVersion ends here ==================================================

      case 'SimpleVersion'

        if numFidels > 1 & ( strcmp(params.acqStrategy, 'MF-GP-UCB') | ...
                             strcmp(params.acqStrategy, 'Illus-MF-GP-UCB') )
          numF2InitPts = floor(0.5*initBudget/mfFunc.costs(end));
          remBudget = initBudget - numF2InitPts*mfFunc.costs(2);
          numF1InitPts = min([200, 10*numFidels*numDims, ...
                                    floor(remBudget/mfFunc.costs(1))]);
          numF2InitPts = min(numF1InitPts, numF2InitPts);

          F1InitPts = boGetInitPts(bounds, numF1InitPts);
          F1InitVals = zeros(numF1InitPts, 1);
          F1InitTimes = zeros(numF1InitPts, 1);
          F2InitPts = F1InitPts(1:numF2InitPts, :);
          F2InitVals = zeros(numF2InitPts, 1);
          F2InitTimes = zeros(numF2InitPts, 1);
          for i = 1:numF1InitPts
            F1InitVals(i) = mfFunc.evalAtFidel(1, F1InitPts(i,:));
            F1InitTimes(i) = cputime - params.initTime;
          end
          for i = 1:numF2InitPts
            F2InitVals(i) = mfFunc.evalAtFidel(2, F2InitPts(i,:));
            F2InitTimes(i) = cputime - params.initTime;
          end
          initPts{1} = F1InitPts;
          initVals{1} = F1InitVals;
          initTimes{1} = F1InitTimes;
          initPts{2} = F2InitPts;
          initVals{2} = F2InitVals;
          initTimes{2} = F2InitTimes;

        else
          numFMInitPts = floor(initBudget/mfFunc.costs(end));
          FMInitPts = boGetInitPts(bounds, numFMInitPts);
          FMInitVals = zeros(numFMInitPts, 1);
          FMInitTimes = zeros(numFMInitPts, 1);
          for i = 1:numFMInitPts
            FMInitVals(i) = mfFunc.evalAtFidel(numFidels, FMInitPts(i,:)); 
            FMInitTimes(i) = cputime - params.initTime;
          end
          initPts{numFidels} = FMInitPts;
          initVals{numFidels} = FMInitVals;
          initTimes{numFidels} = FMInitTimes;
        end


      % SimpleVersion ends here =====================================================

      otherwise
        error('Unknown Init Strategy');
      end

    % Save to params
    params.initPts = initPts;
    params.initVals = initVals;
    params.initTimes = initTimes;
  else
    fprintf('Using the given Initialisation.\n');
  end
  initPts = params.initPts;
  initVals = params.initVals;
  initTimes = params.initTimes;

  % Compute the cost incurred for initialisation
  initCost = 0;
  numAllInitFuncEvals = 0;
  numIndInitFuncEvals = zeros(numFidels, 1);
  initEvalFidels = zeros(0, 1);
  for i = 1:numFidels
    initCost = initCost + size(initVals{i}, 1) * mfFunc.costs(i);
    numIndInitFuncEvals(i) = size(initVals{i}, 1);
    numAllInitFuncEvals = numAllInitFuncEvals + numIndInitFuncEvals(i);
    initEvalFidels = [initEvalFidels; i * ones(numIndInitFuncEvals(i), 1)];
  end
  params.initCost = initCost;
  params.numAllInitFuncEvals = numAllInitFuncEvals;
  params.numIndInitFuncEvals = numIndInitFuncEvals;

  % Set some stuff up to be returned as history
  params.initEvalFidels = initEvalFidels;
  params.initEvalVals = cell2mat(initVals);
  params.initEvalPts = cell2mat(initPts);
  params.initEvalTimes = cell2mat(initTimes);
  instRewards = params.initEvalVals .* mfFunc.costs(params.initEvalFidels);
  params.initCumRewards = cumsum(instRewards);

  fprintf('Initialisation used %s: %d points, %d cost, %0.2f time.\n', ...
    mat2str(numIndInitFuncEvals), numAllInitFuncEvals, initCost, ...
    params.initEvalTimes(end));

  % Some statistics we will use later on.
  stdY = std(params.initEvalVals);
  meanY = mean(params.initEvalVals);
  maxY = max(params.initEvalVals);
  minY = min(params.initEvalVals);
  rangeY = maxY - minY;
  params.rangeY = rangeY;

  % The not-so-important hyper parameters
  if ~isfield(params, 'optPtStdThreshold') | isempty(params.optPtStdThreshold)
    params.optPtStdThresholdCoeff = 0.01;
  end

  % Hyper Parameters for MF optimisation.
  % ============================================================================
  if numFidels > 1 & ( strcmp(params.acqStrategy, 'MF-GP-UCB')  | ...
                        strcmp(params.acqStrategy, 'Illus-MF-GP-UCB') )
    % zeta
    if ~isfield(params, 'fixZeta'), % Don't fix zeta by default
      params.fixZeta = false;
    end
    if ~isfield(params, 'zetas')
      if ~isfield(params, 'zeta')
        if ~exist('maxF1F2Diff', 'var'), maxF1F2Diff = 0;
        end
        params.zeta = ceil(1e4 * 2*maxF1F2Diff)/1e4;
      end
        params.zetas = params.zeta * ( (numFidels-1):-1:1 )';
    end
    % To fix gamma or not
    if ~isfield(params, 'fixGamma'), % Don't fix gamma by default,
      params.fixGamma = false;
    end
    % Initial value of gamma
    if ~isfield(params, 'gammas')
      if ~isfield(params, 'gamma')
        params.gamma = ceil(1e4 * 0.01 * rangeY)/1e4;
      end
      params.gammas = params.gamma * ones(numFidels-1, 1);
    end
  else
    params.zetas = [];
    params.gammas = [];
  end

  % Hyper parameters for the GP.
  % ============================================================================
  if ~isfield(params, 'bwRange')
    params.bwRange = sqrt(numDims) * [1e-3, 10] * mean(bounds(:,2) - bounds(:,1));
  end
  if ~isfield(params, 'scaleRange')
    params.scaleRange = [1 10] * rangeY;
  end
  if ~isfield(params, 'gpMeanFuncs'),
    if strcmp(params.acqStrategy, 'MF-GP-UCB') | ...
      strcmp(params.acqStrategy, 'Illus-MF-GP-UCB') | ...
      strcmp(params.acqStrategy, 'GP-UCB') 
      priorMeanVal = maxY + 2*rangeY; % set it to be a large value. works best for UCB.
    else
      priorMeanVal = meanY; % Do the usual thing. works best for EI and the rest.
    end
    params.gpMeanFuncs = repmat({@(t) priorMeanVal * ones(size(t,1), 1)}, ...
      numFidels, 1);
  end
  if ~isfield(params, 'gpNoiseVars'),
    params.gpNoiseVars = repmat( 1e-4 * stdY^2, numFidels, 1);
  end
  if ~isfield(params, 'multGPLearnStrategy'),
%     params.multGPLearnStrategy = 'jointLearn'; % seems to work best for all.
    params.multGPLearnStrategy = 'sepLearn'; % seems to work best for all.
  end
  if ~isfield(params, 'gpDiRectOpts')
    params.gpDirectOpts.maxevals = 200;
  end

  % Parameters for the initial GP
  % -----------------------------------------------------------------------------
  if ~isfield(params, 'initGPBw'), 
    params.initGPBw = params.bwRange(1) * ones(numFidels,1);
  end
  if ~isfield(params, 'initGPScale'), 
    params.initGPScale = std( cell2mat(initVals) ) * ones(numFidels, 1);
  end
  if ~isfield(params, 'initGPFuncHs')
    initFuncHs = cell(numFidels, 1);
    for f = 1:numFidels
      [~,~,~,currFuncH] = GPRegression([], initPts{f}, initVals{f}, ...
        params.initGPBw(f), params.initGPScale(f), params.gpMeanFuncs{f}, ...
        params.gpNoiseVars(f) );
      initFuncHs{f} = currFuncH;
    end
    params.initGPFuncHs = initFuncHs;
  end

  % Finally, is the budget time based ?
  if ~isfield(params, 'budgetType')
    params.budgetType = 'givenCost';
  end

end

