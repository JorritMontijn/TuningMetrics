function [dblZeta,sZETA,sMSD] = getZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,intPlot,boolVerbose)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZeta,sZETA,sMSD] = getZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,intPlot,boolVerbose)
	%	input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- intPlot: integer, plotting switch (0=none, 1=traces, 2=raster plot) (default: 0)
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
	%		- vecPeaksValues; z-scores of peaks in vecZ
	%		- vecPeaksIdxT; index vector for peaks in vecZ
	%		- vecPeaksTimes; time in trial for peaks in vecZ
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
	if size(vecEventStarts,2) > 2
		vecEventStarts = vecEventStarts';
	end
	if size(vecEventStarts,2) == 2
		boolActDiff = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(vecEventStarts(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = true;
	end
	
	%% prepare interpolation points
	%pre-allocate
	intMaxRep = size(vecEventStarts,1);
	cellSpikeTimesPerTrial = cell(intMaxRep,1);
	
	%go through trials to build spike time vector
	for intEvent=1:intMaxRep
		%get times
		dblStartT = vecEventStarts(intEvent,1);
		dblStopT = dblStartT + dblUseMaxDur;
		
		% build trial assignment
		cellSpikeTimesPerTrial{intEvent} = vecSpikeTimes(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT) - dblStartT;
	end
	
	%get spikes in fold
	vecSpikeT = sort(cell2vec(cellSpikeTimesPerTrial),'ascend');
	intSpikes = numel(vecSpikeT);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);

	%% run bootstraps2
	hTic = tic;
	matRandDiff = nan(intSpikes,intResampNum);
	for intResampling=1:intResampNum
		%% msg
		if boolVerbose && toc(hTic) > 5
			fprintf('Now at resampling %d/%d\n',intResampling,intResampNum);
			hTic = tic;
		end
		%% get random subsample
		vecStimUseOnTime = vecEventStarts(:,1) + 2*dblUseMaxDur*(rand(size(vecEventStarts(:,1)))-0.5);
		
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
	vecZ = ((vecRealDiff-mean(vecRandMean))./mean(vecRandSd));
	if numel(vecZ) < 3
		dblZeta = 0;
		sZETA = struct;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		return
	end
	
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
	[dummy,intPeakLoc]= max(abs(vecZ));
	dblMaxZTime = vecSpikeT(intPeakLoc);
	dblZ = vecZ(intPeakLoc);
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
			plotRaster(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
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
		%{
		subplot(2,3,4)
		hold on
		plot(vecSpikeT,vecRealDiff);
		xlabel('Time  from event (s)');
		ylabel('Offset of data from linear (frac pos)');
		title(sprintf('Real diff data/baseline'));
		fixfig
		%}
		subplot(2,3,4)
		cla;
		hold all
		for intOffset=1:intPlotIters
			plot(vecSpikeT,matRandDiff(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(vecSpikeT,vecRealDiff,'Color',lines(1));
		scatter(vecAllPeakTimes,vecRealDiff(vecAllPeakLocs),'kx');
		scatter(dblMaxZTime,vecRealDiff(intPeakLoc),'rx');
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolActDiff
			title(sprintf('Zeta=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZeta,dblP,dblHzD,dblHzP));
		else
			title(sprintf('Zeta=%.3f (p=%.3f)',dblZeta,dblP));
		end
		fixfig
		
	end
	if dblP < 0.05 || dblHzP < 0.05
		%calculate MSD if significant
		[vecMSD,sMSD] = getMultiScaleDeriv(vecSpikeTimes,vecEventStarts,dblUseMaxDur,[],[],[],intPlot,boolVerbose);
	else
		sMSD = [];
	end
	
	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZ = dblZ;
		sZETA.dblP = dblP;
		if boolActDiff
			sZETA.dblHzD = dblHzD;
			sZETA.dblHzP = dblHzP;
		end
		sZETA.vecSpikeT = vecSpikeT;
		sZETA.vecZ = vecZ;
		sZETA.vecPeaksValues = vecAllVals;
		sZETA.vecPeaksIdxT = vecAllPeakLocs;
		sZETA.vecPeaksTimes = vecAllPeakTimes;
		sZETA.vecPeaksWidths = vecAllPeakWidths;
		sZETA.vecPeaksProminences = vecAllPeakProminences;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
	end
end

