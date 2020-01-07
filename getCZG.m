function [dblZeta,sDeriv] = getCZG(vecAutoTimes,vecHeteroTimes,dblWindowT,intPlot,intResampNum,boolVerbose)
	%getZeta Calculates neuronal responsiveness index zeta
	%syntax: [dblZ,sOptionalOutputs] = getZeta(vecAutoTimes,vecHeteroTimes,dblUseMaxDur,intPlot,intResampNum,boolVerbose)
	%
	%Version history:
	%1.0 - September 24 2019
	%	New procedure to determine statistical significance [by JM]
	
	%% prep data
	%ensure orientation
	vecAutoTimes = vecAutoTimes(:);
	vecHeteroTimes = vecHeteroTimes(:);
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot) || intPlot == 2
		intPlot = 1;
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = true;
	end
	
	%% build test vector
	dblAddFrac=0.05;
	vecAddSpikes = vecAutoTimes(randperm(numel(vecAutoTimes),round(dblAddFrac*numel(vecHeteroTimes))));
	vecAddJitter = (2/1000)*randn(size(vecAddSpikes));
	vecHeteroTimes2 = sort(cat(1,vecHeteroTimes,vecAddSpikes+vecAddJitter));

	%% get ZETA
	dblWindowT=0.25;
	dblUseMaxDur = dblWindowT*2;
	[dblZeta,sOpt] = getZeta(vecAutoTimes,vecHeteroTimes2-dblWindowT,intPlot,dblUseMaxDur,intResampNum,boolVerbose);
	vecSpikeT = sOpt.vecSpikeT;
	
	%% get smoothed LDD
	intPeaks = 3;
	intSmoothSd = 3;
	dblBinSize = 1/1000;
	dblPeakScale = 4;
	[vecSmoothDeriv,sDeriv] = getSmoothDeriv(vecSpikeT,intPeaks,intSmoothSd,dblBinSize,dblPeakScale);
	vecPeakTimes = sDeriv.vecPeakTimes-dblWindowT;
	vecPeakValues = sDeriv.vecPeakValues;
	vecBinT = sDeriv.vecBinT-dblWindowT;
	%% plot
	if intPlot > 0
		vecAx=get(gcf,'Children');
		%build ticks
		vecTickX = [0 dblWindowT dblUseMaxDur];
		vecTickXlabel = [-dblWindowT 0 dblWindowT];
		for intAx=1:numel(vecAx)
			axes(vecAx(intAx));
			xlim([0 dblUseMaxDur]);
			set(gca,'xtick',vecTickX,'xticklabel',vecTickXlabel);
		end
		
		subplot(2,3,6)
		cla;
		plot(vecBinT,vecSmoothDeriv);
		hold on
		scatter(vecPeakTimes,vecPeakValues,'kx');
		for intPeak=1:intPeaks
			text(vecPeakTimes(intPeak),vecPeakValues(intPeak),sprintf('T=%.0fms',vecPeakTimes(intPeak)/dblBinSize),'FontSize',12);
		end
		hold off
		xlabel('Time from event (s)');
		ylabel('Derivative');
		title(sprintf('%s %s-%s,U%d/C%d',strArea,strDate,strBlock,intSU,intClust));
		fixfig
		drawnow;
		
		
		strFileName = sprintf('%s%s-%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
		%export_fig([strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.tif']);
		%print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
		
	end

