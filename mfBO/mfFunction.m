classdef mfFunction < handle

  properties 
    numFidels; % The number of fidelities.
    fidelFunctions; % A numFidelsx1 cell array containing each function.
    bounds; % The (rectangular) domain for the function.
    numDims; % The number of dimensions.
    costs; % A numFidels x 1 vector of the costs.
    zeta; % The maximum difference between the functions at each fidelity.
    hfMaxVal; % The maximum value of the highest fidelity function.
    hfMaxPt; % The maximiser of the highest fidelity function.

  end % end properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods

    % Constructor
    function obj = mfFunction(fidelFunctions, bounds, costs, zeta, hfMaxPt,hfMaxVal);
      obj.numFidels = numel(fidelFunctions);
      obj.fidelFunctions = fidelFunctions;
      obj.bounds = bounds;
      obj.numDims = size(bounds, 1);
      obj.costs = costs;
      % Set eta
      if exist('zeta', 'var'), obj.zeta = zeta;
      else, obj.zeta = [];
      end
      if exist('hfMaxPt', 'var'),
        obj.hfMaxPt = hfMaxPt;
        obj.hfMaxVal = hfMaxVal;
      else,
        obj.hfMaxPt = [];
        obj.hfMaxVal = [];
      end
      
    end

    % Evaluating at a certain fidelity.
    function val = evalAtFidel(obj, fidel, x)
      val = obj.fidelFunctions{fidel}(x);
    end

  end % end methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

