function [dblZeta,sOptionalOutputs] = getZeta(vecSpikeTimes,vecEventStarts,intPlot,dblUseMaxDur,intResampNum,boolVerbose)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZ,sOptionalOutputs] = getZeta(vecSpikeTimes,vecEventStarts,intPlot,dblUseMaxDur,intResampNum,boolVerbose)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- intPlot: integer, plotting switch (0=none, 1=traces, 2=raster plot) (default: 0)
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 25)
	%	- boolVerbose: boolean, switch to print messages
	%
	%	output:
	%	- dblZeta; FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%	- sOptionalOutputs; structure with fields:
	%		- dblZ; uncorrected peak z-score 
	%		- dblP; p-value corresponding to zeta
	%		- dblHzD; Cohen's D based on mean-rate stim/base difference
	%		- dblHzP; p-value based on mean-rate stim/base difference
	%		- vecInterpT: timestamps of interpolated z-scores
	%		- vecZ; z-score for all time points corresponding to vecInterpT
	%		- vecPeaksZ; z-scores of peaks in vecZ
	%		- vecPeaksIdxT; index vector for peaks in vecZ
	%		- vecPeaksTime; time in trial for peaks in vecZ
	%		- vecPeaksWidths; width in spikes for peaks in vecZ
	%		- vecPeaksProminences; height in z-score of peaks in vecZ
	%		- vecRealDiff: real offset of spikes relative to uniform rate
	%		- matRandDiff; matrix of shuffled runs with offset to uniform
	%
	%Version history:
	%0.9 - June 27 2019
	%	Created by Jorrit Montijn
	%1.0 - September 24 2019
	%	New procedure to determine statistical significance [by JM]
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%calculate stim/base difference?
	boolActDiff = false;
	dblHzD = nan;
	if size(vecEventStarts,2) == 2
		boolActDiff = true;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(vecEventStarts(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%% prepare interpolation points
	intMaxRep = size(vecEventStarts,1);
	[vecEventPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,vecEventStarts(:,1));
	indUseSpikes = vecEventPerSpike > 0 & vecEventPerSpike <= intMaxRep & vecTimePerSpike>0 & vecTimePerSpike < dblUseMaxDur;
	vecUseSpikeTimes = vecTimePerSpike(indUseSpikes);
	vecSpikeT = unique(sort(vecUseSpikeTimes,'ascend'));
	intSpikes = numel(vecSpikeT);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);

	%% run bootstraps2
	matRandDiff = nan(intSpikes,intResampNum);
	for intResampling=1:intResampNum
		%% msg
		if boolVerbose
			fprintf('Now at resampling %d/%d\n',intResampling,intResampNum);
		end
		%% get random subsample
		vecStimUseOnTime = vecEventStarts(:,1) + 2*median(diff(vecEventStarts(:,1)))*rand(size(vecEventStarts(:,1)));
		
		%get temp offset
		[vecRandDiff,vecRandFrac,vecRandFracLinear] = ...
			getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
	
		%assign data
		matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
	end

	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	%define plots
	vecRandMean = nanmean(matRandDiff,2);
	vecRandSd = nanstd(matRandDiff,[],2);
	vecZ = ((vecRealDiff-vecRandMean)./vecRandSd);
	
	%get max & min values
	[vecPosVals,vecPosPeakLocs,vecPosPeakWidth,vecPosPeakHeight]= findpeaks(vecZ,'MinPeakDistance',numel(vecZ)/10);
	[vecNegVals,vecNegPeakLocs,vecNegPeakWidth,vecNegPeakHeight]= findpeaks(-vecZ,'MinPeakDistance',numel(vecZ)/10);
	vecAllVals = cat(1,vecPosVals,-vecNegVals);
	vecAllPeakLocs = cat(1,vecPosPeakLocs,vecNegPeakLocs);
	vecAllPeakWidths = cat(1,vecPosPeakWidth,vecNegPeakWidth);
	vecAllPeakProminences = cat(1,vecPosPeakHeight,vecNegPeakHeight);
	
	%remove minor peaks
	indRem = (vecAllPeakProminences./sum(vecAllPeakProminences) < 1/20);
	vecAllVals(indRem) = [];
	vecAllPeakLocs(indRem) = [];
	vecAllPeakWidths(indRem) = [];
	vecAllPeakProminences(indRem) = [];
	%reorder
	[vecAllPeakLocs,vecReorder] = sort(vecAllPeakLocs,'ascend');
	vecAllVals = vecAllVals(vecReorder);
	vecAllPeakWidths = vecAllPeakWidths(vecReorder);
	vecAllPeakProminences = vecAllPeakProminences(vecReorder);
	vecAllPeakTimes = vecSpikeT(vecAllPeakLocs);
	
	%find highest peak and retrieve value
	if isempty(vecAllVals)
		[dummy,intLoc]= max(abs(cat(1,vecPosVals,vecNegVals)));
		veTempPeakLocs = cat(1,vecPosPeakLocs,vecNegPeakLocs);
		intInterpLoc = veTempPeakLocs(intLoc);
	else
		[dummy,intLoc]= max(abs(vecAllVals));
		intInterpLoc = vecAllPeakLocs(intLoc);
	end
	dblMaxZTime = vecSpikeT(intInterpLoc);
	dblZ = vecZ(intInterpLoc);
	dblCorrectionFactor = 2/3.5;
	dblZeta = dblZ*dblCorrectionFactor;
	dblP=1-(normcdf(abs(dblZeta))-normcdf(-abs(dblZeta)));
	
	if boolActDiff
		%% calculate mean-rate difference
		%pre-allocate
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
		intPlotIters = min([size(matRandDiff,2) 50]); 
		
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
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.1:dblUseMaxDur,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		
		subplot(2,3,3)
		plot(vecSpikeT,vecRealFrac)
		hold on
		plot(vecSpikeT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
		title(sprintf('Real data'));
		xlabel('Time from event (s)');
		ylabel('Fractional position of spike in trial');
		fixfig
		
		subplot(2,3,4)
		hold on
		plot(vecSpikeT,vecRealDiff);
		xlabel('Time  from event (s)');
		ylabel('Offset of data from linear (frac pos)');
		title(sprintf('Real diff data/baseline'));
		fixfig
		
		subplot(2,3,5)
		cla;
		hold all
		for intOffset=1:intPlotIters
			plot(vecSpikeT,matRandDiff(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealDiff,'Color',lines(1));
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		title(sprintf('1-%d, diff data/baseline',intPlotIters));
		fixfig
		
		subplot(2,3,6)
		plot(vecSpikeT,vecZ);
		hold on
		scatter(vecAllPeakTimes,vecAllVals,'kx');
		scatter(dblMaxZTime,dblZ,'rx');
		hold off
		xlabel('Time from event (s)');
		ylabel('Z-score');
		title(sprintf('Zeta=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZeta,dblP,dblHzD,dblHzP));
		fixfig
	end
	
	%% build optional output structure
	if nargin > 1
		sOptionalOutputs = struct;
		sOptionalOutputs.dblZ = dblZ;
		sOptionalOutputs.dblP = dblP;
		sOptionalOutputs.dblHzD = dblHzD;
		sOptionalOutputs.dblHzP = dblHzP;
		sOptionalOutputs.vecInterpT = vecSpikeT;
		sOptionalOutputs.vecZ = vecZ;
		sOptionalOutputs.vecPeaksZ = vecAllVals;
		sOptionalOutputs.vecPeaksIdxT = vecAllPeakLocs;
		sOptionalOutputs.vecPeaksTime = vecAllPeakTimes;
		sOptionalOutputs.vecPeaksWidths = vecAllPeakWidths;
		sOptionalOutputs.vecPeaksProminences = vecAllPeakProminences;
		sOptionalOutputs.vecRealDiff = vecRealDiff;
		sOptionalOutputs.matRandDiff = matRandDiff;
	end
end

