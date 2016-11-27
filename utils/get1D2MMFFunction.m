function [ret, retSF] = get1D2MMFFunction(params)

  if ~exist('params', 'var') | isempty(params)
    params = struct;
  end

  f1 = @(t) - 70*(t+0.01).* (t-0.31).* (t+0.51).* (t-0.71).* (t-0.98) + 0.05;
  f2 = @(t) - 70*(t-0).* (t-0.35).* (t+0.55).* (t-0.65).* (t-0.98);

  % Obtain eta
  th = linspace(0,1,10000)';
  F1 = f1(th);
  F2 = f2(th);
  eta = max( F2 - F1 );
  [hfMaxVal, hfMaxIdx] = max(F2);
  hfMaxPt = th(hfMaxIdx);
  costs = [1;10];

  % compute the maximum value.
  ret = mfFunction({f1; f2}, [0 1], costs, eta, hfMaxPt, hfMaxVal);

  % Now just return a single fidelity version too.
  retSF = mfFunction({f2}, [0 1], costs(2), [], hfMaxPt, hfMaxVal);

end 

