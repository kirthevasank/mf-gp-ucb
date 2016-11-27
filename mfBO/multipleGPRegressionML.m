function [funcHs, hyperParams] = multipleGPRegressionML(Xtrs, Ytrs, hyperParams)
% A function which learns the same kernel for multiple Gaussian processes using
% data from all GPs. 
% We will use a Gaussian kernel of the form  scale * exp(-||x-x'||^2/2*bw^2)
% Inputs
%   Xtrs: a numGPs x 1 cell array where each element is the X matrix for that GP.
%   Ytrs: a numGPs x 1 cell array where each element is the Y vector for that GP.
%   hyperParams: a struct which (optionally) has the following fields.
%     meanFuncs: a numGPs x 1 cell array of function handles for the prior GP mean.
%                If just a single function, we will use it for all GPs.
%     noisVars: a numGPs x 1 vector of positive scalars for the GP noise. If just a
%             single value, we will use it for all GPs.
%     bwRange: The range (lower and upper limits) for optimising the bandwidth via GP
%              marginal likelihood maximisation.
%     scaleRange: The range for choosing the scale for the kernel.
% Outputs
%   funcHs: a numGPs x 1 of function handles which output the posterior mean and
%     variance of the estimated GPs.
%   hyperParams: a struct which contains the optimal hyper parameters.

  % Check if we have only one dataset
  if isa(Xtrs, 'double'), Xtrs = {Xtrs}; end
  if isa(Ytrs, 'double'), Ytrs = {Ytrs}; end

  numGPs = numel(Xtrs);
  numDims = size(Xtrs{1}, 1);

  % Book-keeping Hyper-parameters
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFAULT_SAFE = true; 
  DEFAULT_SAFE = false; % I missed some bugs due 2 default values. Disabling for now.
  if DEFAULT_SAFE
%     fprintf('**** Checking for Default values !!\n');
    if ~exist('hyperParams', 'var') | isempty(hyperParams)
      hyperParams = struct; 
    end

    % Mean functions --------------------------------------------------------
    if ~isfield(hyperParams, 'meanFuncs')
      fprintf('** Using Default Mean Functions. ***');      
      meanY = mean(cell2mat(Ytrs));
      hyperParams.meanFuncs = @(t) meanY * ones( size(t, 1), 1);
    end
    if size(hyperParams.meanFuncs, 1) == 1,
      meanFuncs = cell(numGPs, 1);
      for i = 1:numGPs, meanFuncs{i} = hyperParams.meanFuncs; 
        end
      hyperParams.meanFuncs = meanFuncs;
    end

    % Noise parameters ------------------------------------------------------
    if ~isfield(hyperParams, 'noiseVars')
      fprintf('** Using Default Noise variances. ***');      
      stdY = std(cell2mat(Ytrs));
      hyperParams.noiseVars = (0.05 * stdY)^2;
    end
    if size(hyperParams.noiseVars, 1) == 1
      hyperParams.noiseVars = hyperParams.noiseVars * ones(numGPs, 1);
    end

    % Range for the bandwidth -----------------------------------------------
    if ~isfield(hyperParams, 'bwRange') | isempty(hyperParams.bwRange)
      fprintf('** Using Default BW ranges. ***');      
      stdX = norm(std( cell2mat(Xtrs) ));
      hyperParams.bwRange = [0.01 100] * stdX;
    end

    % Range for the scale ---------------------------------------------------
    if ~isfield(hyperParams, 'scaleRange') | isempty(hyperParams.scaleRange)
      fprintf('** Using Default scale ranges. ***');      
      if ~exist('stdY', 'var'), stdY = std(cell2mat(Ytrs));
      end
      hyperParams.scaleRange = [0.01 100] * stdY;
    end

    if ~isfield(hyperParams, 'diRectOpts') | isempty(hyperParams.diRectOpts)
      fprintf('** Using Default DiRect Opt params. ***');      
      diRectOpts.maxevals = 200;
      hyperParams.diRectOpts = diRectOpts;
    end

    if ~isfield(hyperParams, 'multGPLearnStrategy') | ...
      isempty(hyperParams.multGPLearnStrategy),
      fprintf('** Using Default Multiple GP Learn Strategy. ***');      
      hyperParams.multGPLearnStrategy = 'usePrevious';
    end
  end
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Bring all hyperparameters to workspace --------------------------------
  meanFuncs = hyperParams.meanFuncs;
  noiseVars = hyperParams.noiseVars;
  bwRange = hyperParams.bwRange;
  scaleRange = hyperParams.scaleRange;
  diRectOpts = hyperParams.diRectOpts;
  multGPLearnStrategy = hyperParams.multGPLearnStrategy;

  % Book-keeping Hyper-parameters ends here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Optimise for bw and scale using DiRect.
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  diRectBounds = [log(bwRange); log(scaleRange)];

  switch multGPLearnStrategy


    case {'usePrevious', 'jointLearn'}

      % if usePrevious ---------------------------------------------------------
      if strcmp(hyperParams.multGPLearnStrategy, 'usePrevious')
        % Use the previous GPs to initialise the mean of the following GP.
        prelimFuncHs = cell(numGPs, 1);
        bwOpt = zeros(numGPs, 1); 
        scaleOpt = zeros(numGPs, 1);
        for i = 1:(numGPs-1)
          nlmlF = @(t) normMargLikelihood( exp(t(1)), exp(t(2)), Xtrs{i}, Ytrs{i}, ...
                        meanFuncs{i}, noiseVars(i));
          [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOpts);
          bwOpt(i) = exp(optParams(1));
          scaleOpt(i) = exp(optParams(2));
          [~, ~, ~, currFuncH] = GPRegression([], Xtrs{i}, Ytrs{i}, ...
            bwOpt(i), scaleOpt(i), meanFuncs{i}, noiseVars(i));
          prelimFuncHs{i} = currFuncH;
          if i < numGPs
            meanFuncs{i+1} = prelimFuncHs{i};
          end
        end
      end
      % end if usePrevious ---------------------------------------------------

      anlmlF = @(t) allNormMargLikelihood( exp(t(1)), exp(t(2)), Xtrs, Ytrs, ...
        meanFuncs, noiseVars);
      [~, optParams] = diRectWrap(anlmlF, diRectBounds, diRectOpts);
      bwOpt = exp(optParams(1)) * ones(numGPs, 1);
      scaleOpt = exp(optParams(2)) * ones(numGPs, 1);
%       fprintf('Picked bw/scale = %0.5f (%0.3f, %0.3f), %0.5f (%0.3f, %0.3f)\n', ...
%         bwOpt(1), bwRange(1), bwRange(2), ...
%         scaleOpt(1), scaleRange(1), scaleRange(2) );


    case 'sepLearn'
      % Learn each GP Separately 
      bwOpt = zeros(numGPs, 1);
      scaleOpt = zeros(numGPs, 1);
      for i = 1:numGPs
        if i==1 | (size(Ytrs{i},1) >= min(10, size(Ytrs{i-1},1))) 
          nlmlF = @(t) normMargLikelihood( exp(t(1)), exp(t(2)), Xtrs{i}, Ytrs{i}, ...
                        meanFuncs{i}, noiseVars(i));
          [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOpts);
          bwOpt(i) = exp(optParams(1));
          scaleOpt(i) = exp(optParams(2));
        else
          bwOpt(i) = bwOpt(i-1); 
          scaleOpt(i) = scaleOpt(i-1);
        end
      end

    otherwise
      error('Unknown Multiple GP Learning strategy');
  end
  % Optimising for hyper-parameters ends here. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Finally fit the GPs and get the Function Handles.
  funcHs = cell(numGPs, 1);
  for i = 1:numGPs
    [~, ~, ~, currFuncH] = GPRegression([], Xtrs{i}, Ytrs{i}, ...
      bwOpt(i), scaleOpt(i), meanFuncs{i}, noiseVars(i));
    funcHs{i} = currFuncH;
  end

  fprintf('Picked BWs: %s, (%0.4f,%.4f) Scales: %s, (%0.4f,%.4f).\n', ...
    mat2str(round(1e4*bwOpt)/1e4), bwRange(1), bwRange(2), ...
    mat2str(round(1e4*scaleOpt)/1e4), scaleRange(1), scaleRange(2));
  hyperParams.bw = bwOpt;
  hyperParams.scale = scaleOpt;

end


% This iterates through normMargLikelihood to obtain the totla marginal likelihood of
% all data.
function allNlml = allNormMargLikelihood(bw, scale, Xs, Ys, meanFuncs, noiseVars)
  numGPs = numel(Xs);
  allNlml = 0;
  for i = 1:numGPs
    if size(Xs{i}, 1) > 2
      currNlml = normMargLikelihood(bw, scale, Xs{i}, Ys{i}, meanFuncs{i}, ...
                                    noiseVars(i));
      allNlml = allNlml + currNlml; 
    end
  end
end

