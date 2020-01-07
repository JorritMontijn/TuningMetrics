function [vecSmoothDeriv,sOptionalOutputs] = getSmoothDeriv(vecSpikeT,intPeaks,intSmoothSd,dblBinSize,dblPeakScale)
	%getSmoothDeriv Returns binned&smoothed locally dynamic derivative (LDD). Syntax:
	%   vecSmoothDeriv = getSmoothDeriv(vecSpikeT,intPeaks,dblSmoothSd,dblBinSize,dblPeakScale)
	%Required input:
	%	- vecSpikeT is a vector with spike times relative to event times
	%
	%Optional inputs:
	%	- intPeaks: number of peaks to return [default: 3]
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of bins) [default: 3]
	%	- dblBinSize: bin size in seconds [default: 1/1000]
	%	- dblPeakScale: critical value for locally dynamic derivative [default: 4]
	%
	%Outputs:
	%	- vecSmoothDeriv; Binned & smoothed LDD
	%	- sOptionalOutputs; structure with fields:
	%		- vecBinT; time vector for values in vecSmoothDeriv
	%		- vecPeakTimes; time in seconds for peaks
	%		- vecPeakValues; values at peaks
	%		- vecPeakWidths; widths of peaks
	%		- vecPeakProminences; prominences of peaks
	%		- vecDerivT; time vector for unsmoothed & unbinned LDD
	%		- vecLDD; unsmoothed & unbinned LDD
	%
	%Version history:
	%1.0 - October 3 2019
	%	Created by Jorrit Montijn
	
	%% set default values
	if ~exist('intPeaks','var') || isempty(intPeaks)
		intPeaks = 3;
	end
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 3;
	end
	if ~exist('dblBinSize','var') || isempty(dblBinSize)
		dblBinSize = 1/1000;
	end
	if ~exist('dblPeakScale','var') || isempty(dblPeakScale)
		dblPeakScale = 4;
	end
	
	%% get difference from uniform
	vecFracs = linspace(0,1,numel(vecSpikeT))';
	vecLinear = (vecSpikeT./max(vecSpikeT));
	vecDiff = vecFracs - vecLinear;
	vecDiff = vecDiff - mean(vecDiff);
	
	%% get locally dynamic derivative
	vecLDD = getLocallyDynamicDerivative(vecSpikeT,vecDiff,dblPeakScale);
	
	%% binning&smooth
	%build time vectors
	vecDerivT = vecSpikeT(1:(end-1)) + diff(vecSpikeT)/2;
	dblMaxT = roundi(max(vecDerivT),3,'floor');
	vecBinT = 0:dblBinSize:dblMaxT;
	%bin
	[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecDerivT,vecLDD,vecBinT);
	vecUseVals = find(vecCounts);
	%interpolate missing values
	vecInterpVals = interp1([-dblBinSize vecBinT(vecUseVals) dblMaxT+dblBinSize],[0 vecMeans(vecUseVals)' 0],vecBinT);
	%smooth
	vecFilt = normpdf(-2*(intSmoothSd):2*intSmoothSd,0,intSmoothSd);
	vecFilt = vecFilt./sum(vecFilt);
	vecSmoothDeriv = conv(vecInterpVals,vecFilt,'same');
	%find peaks
	[vecPeakValuesPos,vecPeakLocsPos,vecPeakWidthsPos,vecPeakProminencesPos] = findpeaks(vecSmoothDeriv);
	[vecPeakValuesNeg,vecPeakLocsNeg,vecPeakWidthsNeg,vecPeakProminencesNeg] = findpeaks(-vecSmoothDeriv);
	%concatenate
	vecPeakProminences = cat(2,vecPeakProminencesPos,-vecPeakProminencesNeg);
	vecPeakLocs = cat(2,vecPeakLocsPos,vecPeakLocsNeg);
	vecPeakWidths = cat(2,vecPeakWidthsPos,vecPeakWidthsNeg);
	vecPeakValues = cat(2,vecPeakValuesPos,-vecPeakValuesNeg);
	vecPeakTimes = vecBinT(vecPeakLocs);
	%sort by prominence
	[vecPeakProminences,vecReorder] = sort(abs(vecPeakProminences),'descend');
	vecPeakTimes = vecPeakTimes(vecReorder);
	vecPeakWidths = vecPeakWidths(vecReorder);
	vecPeakValues = vecPeakValues(vecReorder);
	
	%% build output
	if nargout > 1
		sOptionalOutputs = struct;
		sOptionalOutputs.dblOnset = vecPeakTimes(1);
		sOptionalOutputs.vecBinT = vecBinT;
		sOptionalOutputs.vecDerivT = vecDerivT;
		sOptionalOutputs.vecPeakTimes = vecPeakTimes(1:intPeaks);
		sOptionalOutputs.vecPeakValues = vecPeakValues(1:intPeaks);
		sOptionalOutputs.vecPeakWidths = vecPeakWidths(1:intPeaks)*dblBinSize;
		sOptionalOutputs.vecPeakProminences = vecPeakProminences(1:intPeaks);
		sOptionalOutputs.vecLDD = vecLDD;
	end
	
end

