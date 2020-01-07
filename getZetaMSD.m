function [dblZETA,sOptionalOutputs] = getZeta(vecSpikeTimes,vecEventStarts,intPlot,dblUseMaxDur,intResampNum,boolVerbose,intPeaks,intSmoothSd,dblMinScale)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZeta,sOptionalOutputs] = getZeta(vecSpikeTimes,vecEventStarts,intPlot,dblUseMaxDur,intResampNum,boolVerbose)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- intPlot: integer, plotting switch (0=none, 1=traces, 2=raster plot) (default: 0)
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- boolVerbose: boolean, switch to print messages
	%
	%	output:
	%	- dblZeta; FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%	- sOptionalOutputs; structure with fields:
	%		- dblZ; uncorrected peak z-score 
	%		- dblP; p-value corresponding to zeta
	%		- dblHzD; Cohen's D based on mean-rate stim/base difference
	%		- dblHzP; p-value based on mean-rate stim/base difference
	%		- vecSpikeT: timestamps of spike times (corresponding to vecZ)
	%		- vecZ; z-score for all time points corresponding to vecSpikeT
	%		- vecRealMSprime; multi-scale derivative of vecCDD
	%		- vecCDD; cumulative density difference
	%		- vecPeakProminences; height in z-score of peaks in vecZ
	%		- vecPeakLocs; index vector for peaks in vecZ
	%		- vecPeakWidths; width in spikes for peaks in vecZ
	%		- vecPeakValues; z-scores of peaks in vecZ
	%		- vecPeakTimes; time in trial for peaks in vecZ
	%		- vecPeakEnergy; energy under peaks in vecZ
	%
	%Version history:
	%0.9 - June 27 2019
	%	Created by Jorrit Montijn
	%1.0 - September 24 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - November 18 2019
	%	New procedure based on multi-scale derivative [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%% set default values
	if ~exist('intPeaks','var') || isempty(intPeaks)
		intPeaks = 3;
	end
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 5;
	end
	if ~exist('dblBase','var') || isempty(dblBase)
		dblBase = 1.5;
	end
	if ~exist('dblMinScale','var') || isempty(dblMinScale)
		dblMinScale = round(log(1/1000) / log(dblBase));
	end
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = false;
	end
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(vecEventStarts(:,1)));
	end
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%calculate stim/base difference?
	boolActDiff = false;
	dblHzD = nan;
	if size(vecEventStarts,2) == 2
		boolActDiff = true;
	end
	
	%% run normal
	%get data
	[vecRealMSprime,sOptReal] = getMultiScaleDeriv(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale);
	vecSpikeT = sOptReal.vecSpikeT;
	
	%% run bootstraps2
	cellRandDiff = cell(1,intResampNum);
	cellRandT = cell(1,intResampNum);
	for intResampling=1:intResampNum
		%% msg
		if boolVerbose
			fprintf('Now at resampling %d/%d\n',intResampling,intResampNum);
		end
		%% get random subsample
		vecStimUseOnTime = vecEventStarts(:,1) + 2*dblUseMaxDur*(rand(size(vecEventStarts(:,1)))-0.5);
		
		%get data
		[vecRandMSprime,sOptRand] = getMultiScaleDeriv(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur,intSmoothSd,dblMinScale);
	
		
		%assign data
		cellRandDiff{1,intResampling} = vecRandMSprime;
		cellRandT{1,intResampling} = sOptRand.vecSpikeT;
	end
	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	%define plots
	vecRandSd = cellfun(@std,cellRandDiff);
	vecRandMean = cellfun(@mean,cellRandDiff);
	vecZ = ((vecRealMSprime-mean(vecRandMean))./mean(vecRandSd));
	
	%% find peaks
	%find peaks
	[vecPeakValuesPos,vecPeakLocsPos,vecPeakWidthsPos,vecPeakProminencesPos] = findpeaks(vecZ);
	[vecPeakValuesNeg,vecPeakLocsNeg,vecPeakWidthsNeg,vecPeakProminencesNeg] = findpeaks(-vecZ);
	%concatenate
	vecPeakPosNeg = cat(1,ones(size(vecPeakProminencesPos)),zeros(size(vecPeakProminencesNeg)));
	vecPeakProminences = cat(1,vecPeakProminencesPos,-vecPeakProminencesNeg);
	vecPeakLocs = cat(1,vecPeakLocsPos,vecPeakLocsNeg);
	vecPeakWidths = cat(1,vecPeakWidthsPos,vecPeakWidthsNeg);
	vecPeakValues = cat(1,vecPeakValuesPos,-vecPeakValuesNeg);
	vecPeakTimes = vecSpikeT(vecPeakLocs);
	vecPeakEnergyEst = vecPeakWidths.*vecPeakValues;
	
	%top 3 prominences
	[dummy,vecUsePromPeaks] = findmax(abs(vecPeakProminences),intPeaks);
	
	%top 3 values
	[dummy,vecUseValPeaks] = findmax(abs(vecPeakValues),intPeaks);
	
	%top 3 energy
	[dummy,vecUseNrgPeaks] = findmax(abs(vecPeakEnergyEst),intPeaks);
	
	%get unique peak list
	vecUsePeaks = unique(cat(2,vecUsePromPeaks,vecUseValPeaks,vecUseNrgPeaks));
	vecUsePeaks(isnan(vecUsePeaks)) = [];
	if isempty(vecUsePeaks)
		[dummy,vecUsePeaks] = max(abs(vecZ));
	end
	vecPosPeak = vecPeakPosNeg(vecUsePeaks);
	
	%% calculate energy per peak
	%build overall ordering
	vecBaseThresh = find(abs(vecZ) < 1);
	vecPeakEnergy = zeros(numel(vecUsePeaks),1);
	for intPeak=1:numel(vecUsePeaks)
		intUsePeakLoc = vecPeakLocs(vecUsePeaks(intPeak));
		vecDists = vecBaseThresh - intUsePeakLoc;
		intLocLeft = max(vecDists(vecDists<0))+intUsePeakLoc;
		intLocRight = min(vecDists(vecDists>0))+intUsePeakLoc;
		vecPeakEnergy(intPeak) = sum(vecZ(intLocLeft:intLocRight));
	end
	indRem = (vecPeakEnergy > 0 & vecPosPeak == 0) | (vecPeakEnergy < 0 & vecPosPeak == 1);
	%remove local minima
	vecUseEnergy = vecPeakEnergy;
	vecUseEnergy(indRem) = [];
	vecUsePeaks(indRem) = [];
	%get values
	vecUsePeakVals = vecPeakValues(vecUsePeaks);
	
	%order
	[vecUseEnergy,vecReorder] = sort(abs(vecUseEnergy),'descend');
	vecUsePeaks = vecUsePeaks(vecReorder);
	intKeepPeaks = min([numel(vecUsePeaks) intPeaks]);
	vecKeepPeaks = vecUsePeaks(1:intKeepPeaks);
	vecUseEnergy = vecUseEnergy(1:intKeepPeaks);
	
	%define temporal anomaly value
	[dblZ,intZ] = max(abs(vecPeakValues(vecKeepPeaks)));
	dblZ = vecPeakValues(vecKeepPeaks(intZ));
	dblCorrectionFactor = 1;
	dblZETA = dblZ*dblCorrectionFactor;
	dblP=1-(normcdf(abs(dblZETA))-normcdf(-abs(dblZETA)));
	
	if boolActDiff
		%% calculate mean-rate difference
		%pre-allocate
		intMaxRep = size(vecEventStarts,1);
		vecStimHz = zeros(intMaxRep,1);
		vecBaseHz = zeros(intMaxRep,1);
		dblMedianBaseDur = median(vecEventStarts(2:end,1) - vecEventStarts(1:(end-1),2));
		
		%go through trials to build spike time vector
		for intEvent=1:intMaxRep
			%get times
			dblStartT = vecEventStarts(intEvent,1);
			dblStopT = dblStartT + dblUseMaxDur;
			dblPreT = dblStartT - dblMedianBaseDur;
			
			% build trial assignment
			vecStimHz(intEvent) = sum(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT)/(dblStopT - dblStartT);
			vecBaseHz(intEvent) = sum(vecSpikeTimes < dblStartT & vecSpikeTimes > dblPreT)/dblMedianBaseDur;
		end
		
		%get metrics
		dblHzD = abs(mean(vecStimHz - vecBaseHz)) / ( (std(vecStimHz) + std(vecBaseHz))/2);
		[h,dblHzP]=ttest(vecStimHz,vecBaseHz);
	end
	
	%% plot
	if intPlot
		%plot maximally 50 traces
		intPlotIters = min([size(cellRandDiff,2) 50]); 
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot == 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes,vecEventStarts(:,1));
			xlabel('Time from event (s)');
			ylabel('Trial #');
			title('Spike raster plot');
			fixfig;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		[vecPEPMean,vecPEPSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.005:dblUseMaxDur,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecPEPMean,vecPEPSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		
		
		subplot(2,3,3)
		cla;
		hold all
		matC = lines(2);
		plot(vecSpikeT,sOptReal.vecLinear,'Color',[0.5 0.5 0.5]);
		plot(vecSpikeT,sOptReal.vecFracs,'Color',lines(1));
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		title(sprintf('Cum. dens. real (blue), linear (grey)'));
		fixfig
		
		subplot(2,3,4)
		cla;
		plot(vecSpikeT,sOptReal.vecDiff,'Color',lines(1));
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		title(sprintf('Cum. dens. diff real/linear'));
		fixfig
	
		subplot(2,3,5)
		cla;
		hold all
		for intRandIter=1:intPlotIters
			plot(cellRandT{intRandIter},cellRandDiff{intRandIter},'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealMSprime,'Color',lines(1));
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		title(sprintf('MS'' real (blue), shuff (grey)'));
		fixfig
		
		subplot(2,3,6)
		plot(vecSpikeT,vecZ)
		xlabel('Time from event (s)');
		ylabel('Z-score over shuffled');
		hold on
		scatter(vecSpikeT(vecPeakLocs(vecKeepPeaks(~ismember(1:numel(vecKeepPeaks),intZ)))),vecPeakValues(vecKeepPeaks(~ismember(1:numel(vecKeepPeaks),intZ))),'kx');
		scatter(vecSpikeT(vecPeakLocs(vecKeepPeaks(intZ))),dblZ,'rx');
		text(vecSpikeT(vecPeakLocs(vecKeepPeaks)),vecPeakValues(vecKeepPeaks),cellfun(@num2str,vec2cell(roundi(vecPeakValues(vecKeepPeaks),2)),'UniformOutput',false))
		hold off
		
		if boolActDiff
			title(sprintf('ZETA=%.2f (p=%.3f),d(Hz)=%.2f (p=%.3f)',dblZETA,dblP,dblHzD,dblHzP));
		else
			title(sprintf('ZETA=%.2f (p=%.3f)',dblZETA,dblP));
		end
		fixfig
		
	end
	
	%% build optional output structure
	if nargin > 1
		sOptionalOutputs = struct;
		sOptionalOutputs.dblZ = dblZ;
		sOptionalOutputs.dblP = dblP;
		if boolActDiff
			sOptionalOutputs.dblHzD = dblHzD;
			sOptionalOutputs.dblHzP = dblHzP;
		end
		sOptionalOutputs.vecSpikeT = vecSpikeT;
		sOptionalOutputs.vecZ = vecZ;
		sOptionalOutputs.vecRealMSprime = vecRealMSprime;
		sOptionalOutputs.vecCDD = sOptReal.vecDiff;
		%assign peaks
		sOptionalOutputs.vecPeakProminences = vecPeakProminences(vecKeepPeaks);
		sOptionalOutputs.vecPeakLocs = vecPeakLocs(vecKeepPeaks);
		sOptionalOutputs.vecPeakWidths = vecPeakWidths(vecKeepPeaks);
		sOptionalOutputs.vecPeakValues = vecPeakValues(vecKeepPeaks);
		sOptionalOutputs.vecPeakTimes = vecPeakTimes(vecKeepPeaks);
		sOptionalOutputs.vecPeakEnergy = vecUseEnergy;
	end
end

