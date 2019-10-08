function vecLDD = getLocallyDynamicDerivative(vecTimes,matValues,dblPeakScale)
	%getLocallyDynamicDerivative Returns locally dynamic derivative. Syntax:
	%   vecLDD = getLocallyDynamicDerivative(vecTimes,vecValues,dblPeakScale)
	%		- vecTimes is a vector with time points corresponding to vecValues
	%		- vecValues is a vector with values for which a derivative is computed
	%
	%dblPeakScale is optional (default: four sd)
	%
	%This functions is a subfunction for getSmoothDeriv()
	%
	%
	%Version history:
	%1.0 - October 3 2019
	%	Created by Jorrit Montijn
	
	%check inputs
	if ~exist('dblScale','var') || isempty(dblPeakScale)
		dblPeakScale = 4;
	end
	
	%assign variables
	vecTimes = vecTimes(:);
	matDeriv = diff(matValues,1,1);
	matLocalDeriv = nan(size(matDeriv));
	for intFold=1:size(matDeriv,2)
		%get data
		vecDeriv = matDeriv(:,intFold);
		%calculate
		dblMuD = mean(vecDeriv);
		dblSdD = std(vecDeriv);
		intNumS = numel(vecDeriv);
		dblLowThresh = dblMuD - dblPeakScale*dblSdD;
		dblHighThresh = dblMuD + dblPeakScale*dblSdD;
		vecLDD = zeros(size(vecDeriv));
		%run through all points
		parfor intS=1:intNumS
			%get points left and right of current point
			vecLeftD = vecDeriv(intS:-1:max(1,intS-intNumS));
			vecRightD = vecDeriv(intS:min(intNumS,intS+intNumS));
			%pad with zeros
			vecLeftD = cat(1,vecLeftD,zeros(intNumS+1-numel(vecLeftD),1));
			vecRightD = cat(1,vecRightD,zeros(intNumS+1-numel(vecRightD),1));
			%average bidirectional cumulative sum
			vecCumSumL = cumsum(vecLeftD);
			vecCumSumR = cumsum(vecRightD);
			vecCS = (vecCumSumL+vecCumSumR)/2;
			%get closest threshold, or otherwise max(abs) value
			intIdx=find(vecCS > dblHighThresh | vecCS < dblLowThresh,1);
			if isempty(intIdx)
				[dummy,intIdx] = max(abs(vecCS));
				dblVal = vecCS(intIdx);
			else
				dblVal = vecCS(intIdx);
			end
			%find time difference and divide value by dt
			intMinS = max(1,intS-intIdx);
			intMaxS = min(intNumS,intS+intIdx);
			dblDiffT = (vecTimes(intMaxS)-vecTimes(intMinS));
			vecLDD(intS) = dblVal/dblDiffT;
		end
		matLocalDeriv(:,intFold) = vecLDD;
	end
	%mean across folds
	if size(matLocalDeriv,2)>1
		vecLDD = mean(matLocalDeriv,2)./std(matLocalDeriv,[],2);
	end
end

