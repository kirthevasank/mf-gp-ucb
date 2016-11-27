function [simRegrets, cumRegrets, cumAvgRegrets, cumCosts, infoDetail] = ...
  getSimCumRegrets(mfFunc, history, params)

  if ~exist('params', 'var'), params = struct();
  end

  hfTrueMaxVal = mfFunc.hfMaxVal;

  [hfMaxVals, cumRewards, cumAvgRewards, cumCosts] = ...
    getSimCumRewards(mfFunc, history, params);

  if isempty(hfTrueMaxVal)
%     fprintf('True maximum information not available. Returning rewards instead.\n');
    infoDetail = 'rewards';
    simRegrets = hfMaxVals;
    cumRegrets = cumRewards;
    cumAvgRegrets = cumAvgRewards;

  else
    % Obtain the regrets
    infoDetail = 'regrets';
    simRegrets = max(hfTrueMaxVal - hfMaxVals, 0);
    cumRegrets = max(hfTrueMaxVal * cumCosts - cumRewards, 0);
    cumAvgRegrets = max(hfTrueMaxVal - cumAvgRewards, 0);
  end

end

