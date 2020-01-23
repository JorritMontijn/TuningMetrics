function [vecMSD,sMSD] = getMultiScaleDeriv(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot,boolVerbose)
	%getMultiScaleDeriv Returns multi-scale derivative. Syntax:
	%   [vecMSD,sMSD] = getMultiScaleDeriv(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot,boolVerbose)
	%Required input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%
	%Optional inputs:
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of bins) [default: 3]
	%	- dblMinScale: minimum derivative scale in seconds [default: 1/1000]
	%	- dblBase: critical value for locally dynamic derivative [default: 4]
	%	- intPlot: integer, plotting switch (0=none, 1=plot)
	%	- boolVerbose: boolean, switch to print messages
	%
	%Outputs:
	%	- vecMSprime; Multi-scale derivative
	%	- sMSD; structure with fields:
	%		- vecMSprime;
	%		- vecSpikeT;
	%		- vecFracs;
	%		- vecLinear; 
	%		- vecDiff; 
	%		- vecScale; 
	%		- matSmoothMSprime; 
	%		- matMSprime;
	%
	%Version history:
	%1.0 - October 3 2019
	%	Created by Jorrit Montijn
	%1.1 - January 23 2019
	%	Syntax update [by JM]
	
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
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = true;
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
	vecMSD = mean(matSmoothDeriv,2);
	
	if intPlot > 0
		
		subplot(2,3,5);
		imagesc(matSmoothDeriv');
		vecTicksY = get(gca,'ytick');
		cellTicksY = cellfun(@sprintf,cellfill('%.1e',size(vecTicksY)),vec2cell(vecScale(vecTicksY)),'UniformOutput',false);
		set(gca,'yticklabel',cellTicksY);
		ytickangle(gca,75)
		ylabel(sprintf('Scale (s) (%.1es - %.1es)',vecScale(1),vecScale(end)));
		xlabel('Spike number');
		title('Smoothed multi-scale derivatives');
		fixfig
		grid off
		
		subplot(2,3,6);
		plot(vecSpikeT,vecMSD)
		xlabel('Time from event (s)');
		ylabel('Multi-scale derivative');
		title(sprintf('MS'' by time'));
		fixfig
		
	end
		
	%% build output
	if nargout > 1
		sMSD = struct;
		sMSD.vecMSprime = vecMSD;
		sMSD.vecSpikeT = vecSpikeT;
		sMSD.vecFracs = vecFracs;
		sMSD.vecLinear = vecLinear;
		sMSD.vecDiff = vecDiff;
		
		sMSD.vecScale = vecScale;
		sMSD.matSmoothMSprime = matSmoothDeriv;
		sMSD.matMSprime = matDeriv;
		
	end
end

