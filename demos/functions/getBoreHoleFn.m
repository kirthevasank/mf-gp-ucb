function [mff, sff] = getBoreHoleFn(numFidels, numDims)
% numFidels isn't really needed now. Keeping it since I might want to change it in 
% the future.

  bhBounds = [0.05 0.15; ...
              100, 50000; ...
              63070, 115600; ...
              990, 1110; ...
              63.1, 116; ...
              700, 820; ...
              1120, 1680; ...
              9855, 12045];

  numFidels = 2;
  funcHs = {@(x) bhF1(x, bhBounds); @(x)bhF2(x, bhBounds)};
  hfMaxPt = [];
  hfMaxVal = 309.523221;
  bounds = repmat([0 1], 8, 1);
  costs = (10.^(0:1))';

  mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
  sff = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
                   hfMaxPt, hfMaxVal);

end


function vals = bhF2(pt, bhBounds)
% Borehole fidelity 2 function

  x = mapToOrigBounds(pt, bhBounds);
  rw = x(:,1);
  r  = x(:,2);
  Tu = x(:,3);
  Hu = x(:,4);
  Tl = x(:,5);
  Hl = x(:,6);
  L  = x(:,7);
  Kw = x(:,8);
  
  intTerm = 2* L .* Tu ./ (log(r./rw) .* rw.^2 .* Kw);
  vals = 2*pi * Tu .* (Hu - Hl) ./ ( log(r./rw) .* (1 + intTerm + Tu./Tl) );

end


function vals = bhF1(pt, bhBounds)
% Borehole fidelity 1 function

  x = mapToOrigBounds(pt, bhBounds);
  rw = x(:,1);
  r  = x(:,2);
  Tu = x(:,3);
  Hu = x(:,4);
  Tl = x(:,5);
  Hl = x(:,6);
  L  = x(:,7);
  Kw = x(:,8);
  
  intTerm = 2* L .* Tu ./ (log(r./rw) .* rw.^2 .* Kw);
  vals = 5 * Tu .* (Hu - Hl) ./ ( log(r./rw) .* (1.5 + intTerm + Tu./Tl) );

end

