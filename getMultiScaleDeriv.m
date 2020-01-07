function [vecMSprime,sOptionalOutputs] = getMultiScaleDeriv(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase)
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
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 5;
	end
	if ~exist('dblBase','var') || isempty(dblBase)
		dblBase = 1.5;
	end
	if ~exist('dblMinScale','var') || isempty(dblMinScale)
		dblMinScale = round(log(1/1000) / log(dblBase));
	end
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(vecEventStarts(:,1)));
	end
	
	%% prepare normalized spike times
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
	
	%define largest scale
	dblMaxScale = log(max(vecSpikeT)) / log(dblBase);
	
	%% get difference from uniform
	vecFracs = linspace(0,1,numel(vecSpikeT))';
	vecLinear = (vecSpikeT./max(vecSpikeT));
	vecDiff = vecFracs - vecLinear;
	vecDiff = vecDiff - mean(vecDiff);
	
	%% get multi-scale derivative
	vecExp = dblMinScale:dblMaxScale;
	vecScale=dblBase.^vecExp;
	intScaleNum = numel(vecScale);
	matDeriv = zeros(intSpikes,intScaleNum);
	try %try parallel
		parfor intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intSpikes
				%select points within window
				dblT = vecSpikeT(intS);
				dblMinEdge = dblT - dblScale/2;
				dblMaxEdge = dblT + dblScale/2;
				intMinSpike = find(vecSpikeT > dblMinEdge,1);
				if isempty(intMinSpike),intMinSpike=1;end
				intMaxSpike = find(vecSpikeT > dblMaxEdge,1);
				if isempty(intMaxSpike),intMaxSpike=intSpikes;end
				matDeriv(intS,intScaleIdx) = (vecDiff(intMaxSpike) - vecDiff(intMinSpike))/dblScale;
			end
		end
	catch %otherwise try normal loop
		for intScaleIdx=1:intScaleNum
			dblScale = vecScale(intScaleIdx);
			
			%run through all points
			for intS=1:intSpikes
				%select points within window
				dblT = vecSpikeT(intS);
				dblMinEdge = dblT - dblScale/2;
				dblMaxEdge = dblT + dblScale/2;
				intMinSpike = find(vecSpikeT > dblMinEdge,1);
				if isempty(intMinSpike),intMinSpike=1;end
				intMaxSpike = find(vecSpikeT > dblMaxEdge,1);
				if isempty(intMaxSpike),intMaxSpike=intSpikes;end
				matDeriv(intS,intScaleIdx) = (vecDiff(intMaxSpike) - vecDiff(intMinSpike))/dblScale;
			end
		end
	end
	
	%% smoothing
	vecFilt = normpdf(-2*(intSmoothSd):2*intSmoothSd,0,intSmoothSd)';
	vecFilt = vecFilt./sum(vecFilt);
	matSmoothDeriv = conv2(matDeriv,vecFilt,'same');
	%mean
	vecMSprime = mean(matSmoothDeriv,2);
	
	if 0
		%% plot
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		subplot(2,3,1)
		sOpt = struct;
		sOpt.handleFig =-1;
		[vecPEPMean,vecPEPSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.005:dblUseMaxDur,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecPEPMean,vecPEPSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		
		subplot(2,3,2);
		plot(vecSpikeT,vecDiff)
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		title(sprintf('Cum. dens. diff real/linear'));
		fixfig
	
		subplot(2,3,3);
		imagesc(matDeriv');
		vecTicksY = get(gca,'ytick');
		cellTicksY = cellfun(@sprintf,cellfill('%.1e',size(vecTicksY)),vec2cell(vecScale(vecTicksY)),'UniformOutput',false)
		set(gca,'yticklabel',cellTicksY);
		ytickangle(gca,75)
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Spike number');
		title('Multi-scale derivatives');
		fixfig
		grid off
		
		subplot(2,3,4);
		imagesc(matSmoothDeriv');
		vecTicksY = get(gca,'ytick');
		cellTicksY = cellfun(@sprintf,cellfill('%.1e',size(vecTicksY)),vec2cell(vecScale(vecTicksY)),'UniformOutput',false);
		set(gca,'yticklabel',cellTicksY);
		ytickangle(gca,75)
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Spike number');
		title('Smoothed MS''');
		fixfig
		grid off
		
		subplot(2,3,5);
		plot(vecMSprime)
		xlabel('Spike number');
		ylabel('Multi-scale derivative');
		title(sprintf('MS'' by spike'));
		fixfig
	
		
		subplot(2,3,6);
		plot(vecSpikeT,vecMSprime)
		xlabel('Time from event (s)');
		ylabel('Multi-scale derivative');
		title(sprintf('MS'' by time'));
		fixfig
		
	end
		
	%% build output
	if nargout > 1
		sOptionalOutputs = struct;
		sOptionalOutputs.vecSpikeT = vecSpikeT;
		sOptionalOutputs.vecFracs = vecFracs;
		sOptionalOutputs.vecLinear = vecLinear;
		sOptionalOutputs.vecDiff = vecDiff;
		
		sOptionalOutputs.vecMSprime = vecMSprime;
		sOptionalOutputs.vecScale = vecScale;
		sOptionalOutputs.matSmoothMSprime = matSmoothDeriv;
		sOptionalOutputs.matMSprime = matDeriv;
		
	end
end

