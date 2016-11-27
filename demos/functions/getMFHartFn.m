function [mff, sff] = getMFHartFn(numFidels, numDims)

  alpha0 = [1.0 1.2 3.0 3.2]';
  delta = 0.1 * [0.01 -0.01 -0.1 0.1]';
%   delta = [0.01 -0.01 -0.1 0.1]';

  funcHs = cell(numFidels, 1);
  costs = round(10.^(0:(numFidels-1)))';
  bounds = repmat([0 1], 6, 1);

  for i = 1:numFidels
    alphai = alpha0 + delta * (numFidels - i);
    funcHs{i} = @(x) hartman6Wrap(x, alphai);
  end

  % The highest fidelity maximum
  hfMaxPt = [0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573];
  hfMaxVal = funcHs{numFidels}(hfMaxPt);

  % A multifidelity function and a single fidelity function 
  mff = mfFunction(funcHs, bounds, costs, [], hfMaxPt, hfMaxVal);
  sff = mfFunction({funcHs{numFidels}}, bounds, costs(numFidels), [], ...
                   hfMaxPt, hfMaxVal);

end


function vals = hartman6Wrap(X, alphai)
  numPts = size(X, 1);
  vals = zeros(numPts, 1);
  for i = 1:numPts
    vals(i) = hart6sc( X(i,:), alphai);    
  end
end


function [y] = hart6sc(xx, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HARTMANN 6-DIMENSIONAL FUNCTION, RESCALED
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, x3, x4, x5, x6]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha = [1.0, 1.2, 3.0, 3.2]';
A = [10, 3, 17, 3.50, 1.7, 8;
     0.05, 10, 17, 0.1, 8, 14;
     3, 3.5, 1.7, 10, 17, 8;
     17, 8, 0.05, 10, 0.1, 14];
P = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
               2329, 4135, 8307, 3736, 1004, 9991;
               2348, 1451, 3522, 2883, 3047, 6650;
               4047, 8828, 8732, 5743, 1091, 381];

outer = 0;
for ii = 1:4
  inner = 0;
  for jj = 1:6
    xj = xx(jj);
    Aij = A(ii, jj);
    Pij = P(ii, jj);
    inner = inner + Aij*(xj-Pij)^2;
  end
  new = alpha(ii) * exp(-inner);
  outer = outer + new;
end

y = outer;

end

