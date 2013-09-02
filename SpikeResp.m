classdef SpikeResp < handle
   
   properties
      cats, trials, binSize
      spikeTimes, trains, PSTH
      stim
      isi, isiHist,  isiBins
      spec, specFreq
      rate, rateMean, rateStd
      BARSrate, BARSconfint
      PSTHrate, tauOpt, ISE
      xc, lags
      Crr, Csr, Tsr, F, Pss, Psr, Prr, Pr1r2;
      
      yLab=''
      catNames=[]
   end
   
   methods (Access='public')
      function self = SpikeResp(spikeTimes, binSize)
         % SpikeResp(spikeTimes, binSize)
         if nargin==1
            binSize = 0.2;%ms
         end
         self.binSize = binSize;
         self.spikeTimes = spikeTimes;
         self.cats = size(self.spikeTimes,1);
         self.trials = size(self.spikeTimes,2);
         self.catNames = num2cellstr(1:self.cats);
         self.getRawTrains();
         
      end
      
      function getRawTrains(self)
         [self.PSTH, trainStdDummy, self.trains] = times2trains(self.spikeTimes, self.binSize);
      end
      
      function setStim(self, stim, stimBinSize)
         % add stim for each cat, performs resampling and error checking
         % setStim(stim, stimBinSize)
         if size(stim,1)~=self.cats
            fprintf('\n STIM has the wrong dimensions');
         else
            if stimBinSize~=self.binSize
               if length(stim)<10000
                  stim = resample(stim',100*self.binSize,100*stimBinSize)';
               else
                  stim = resample(stim',self.binSize/self.binSize,stimBinSize/self.binSize)';
               end
            end
            self.stim = stim;
         end
      end
      
      
      function srSub = getSubSetResp(self, subPart)
         %returns a SpikeResp with subPart spikeTimes deleted (set to 0)
         spikeTimesSub = self.spikeTimes;
         
         for cat = 1:self.cats
            for trial = 1:self.trials
               spikeTimesSub(cat, trial, rand(size(spikeTimesSub,3),1)>subPart)=0;
            end
         end
         srSub = SpikeResp(spikeTimesSub, self.binSize);
      end
      
      function srShift = getShiftResp(self, shift)
         %returns a SpikeResp with spikeTimes shifted by shift
         shiftSpikeTimes = self.spikeTimes;
         shiftSpikeTimes(shiftSpikeTimes>0) = shiftSpikeTimes(shiftSpikeTimes>0) + shift;
         srShift = SpikeResp(shiftSpikeTimes,self.binSize);
      end
      
      function srRand = getRandResp(self)
         self.getRawTrains();
         resp = self.trains;
         for cat = 1:self.cats
            for trial = 1:self.trials
               tmp = find(resp(cat,trial,randperm(length(resp)))>0);
               randSpikeTimes(cat,trial,1:length(tmp)) = tmp;
            end
         end
         srRand = SpikeResp(randSpikeTimes,self.binSize);
      end
      
      function srTrim = getTrimResp(self, start, stop, period )
         %returns a SpikeResp with spikeTimes trimmed
         trimSpikeTimes = self.spikeTimes;
         trimSpikeTimes = trimTimes(trimSpikeTimes,start,stop, period);
         srTrim = SpikeResp(trimSpikeTimes,self.binSize);
      end
      
      function getSpec(self)
         % PREPARE
         Fs = 1000/self.binSize;
         fftLen = ceil(size(self.trains,2)/1000)*1000;%ensure that spec is evaluated at integer frequency values
         
         % SPECTRAL ESTIMATION
         %self.spec = zeros(self.trials,self.cats,fftLen);
         for cat = 1:self.cats
            % estimate response amp by welch's method
            %     [psd freq] = getSpecWin(squeeze(rate(cat,:,:)),400,.5,Fs);
            %     freqIdxAmp = find(ismember(freq, harmFreq))
            %     harmAmp((cat-1)*trials + (1:trials),:) = sqrt(psd(:,freqIdxAmp));
            
            % estimate full response (amp+phase) using the periodogramm method
            [spec freq] = getSpec(squeeze(self.trains(cat,:,:)),Fs,fftLen);
            self.spec(cat,:,:) = spec;
            
         end
         self.specFreq = freq;
      end
      
      function getCoherence(self)
         
         self.getRawTrains();
         
         idx = 1:min(size(self.stim,2), size(self.trains,3));
         y = self.trains(:,:,idx);
         x = self.stim(:,idx);
         
         self.Pss = zeros(self.cats,1, 513);
         self.Prr = zeros(self.cats, self.trials, 513);
         self.Pr1r2 = zeros(self.cats, .5*(self.trials^2 - self.trials), 513);
         
         cnt = 0;
         for cat = 1:self.cats
            [self.Pss(cat,1,:), self.F] = pwelch(x(cat,:),[],512,1024,1000./self.binSize);
            for trial = 1:self.trials
               self.Psr(cat,trial,:) = cpsd(x(cat,:),y(cat,trial,:),[],512,1024,1000./self.binSize);
               self.Prr(cat,trial,:) = pwelch(y(cat,trial,:),[],512,1024,1000./self.binSize);
               for trial2 = trial+1:self.trials
                  cnt = cnt + 1;
                  [self.Pr1r2(cat,cnt,:), self.F] = cpsd(y(cat,trial2,:),y(cat,trial,:),[],512,1024,1000./self.binSize);
               end
            end
         end
         self.Csr = mean(self.Psr,2).*conj(mean(self.Psr,2))./(self.Pss.*mean(self.Prr,2));
         self.Crr = mean(self.Pr1r2,2).*conj(mean(self.Pr1r2,2))./(mean(self.Prr,2).*mean(self.Prr,2));
         self.Tsr = conj(mean(self.Psr,2))./self.Pss;
      end
      
      
      function getISI(self)
         for cat = 1:self.cats
            isiTmp = [];
            for trial = 1:self.trials
               st = self.spikeTimes(cat,trial,self.spikeTimes(cat,trial,:)>0);
               isiTmp = [isiTmp(:); diff(st(:))];
               
            end
            self.isi(cat,1:length(isiTmp)) = isiTmp;
            [isiHistTmp, bins] = hist(isiTmp, 0:1:150);
            self.isiHist(cat,:) = isiHistTmp;
            self.isiBins = bins;
         end
      end
      
      function getISIcorr(self)
         % determine serial correlations between ISI's following Farkhooi
         % et al. 2009
         self.getISI();
         thisISI = [];
         for cat = 1:self.cats%concatenate all ISI
            thisISI = [thisISI self.isi(cat,self.isi(cat,:)>0)];
         end
         matISI = [thisISI(3:end-1);thisISI(2:end-2);thisISI(1:end-3)]';%copy for correlation
         [rho,p] = corr(matISI,'type','Spearman');
         ISIrho = rho([2,3])       %take only rho for 1vs2 and 1vs3
         ISIrhop = p([2,3])        %and their p-values
         
      end
      
      
      function [xc,lags] = getCORR(self,maxLag)
         %[xc,lags] = getCORR(self,maxLag)
         xc = zeros(self.cats,self.trials,maxLag);
         for cat = 1:self.cats
            for trial = 1:self.trials
               a=xcov(squeeze(self.trains(cat,trial,:)),maxLag,'coeff');
               xc(cat,trial,:) = a(maxLag+2:end);
            end
         end
         lags = (1:maxLag)*self.binSize;
         self.lags = lags;
         self.xc = xc;
      end
      
      
      function getRate(self,varargin)
         stimLen = max(self.spikeTimes(:));
         stimStart = 0;
         if nargin>1
            stimLen = varargin{1};
         end
         if nargin>2
            stimStart = varargin{2};
         end
         self.rate = zeros(self.cats,self.trials);
         
         for cat = 1:self.cats
            for trial = 1:self.trials
               self.rate(cat,trial) = sum(self.spikeTimes(cat,trial,:)>stimStart & self.spikeTimes(cat,trial,:)<stimStart+stimLen);
            end
         end
         
         self.rate = self.rate./stimLen*1000;
         
         self.rateMean = mean(self.rate,2);
         self.rateStd = std(self.rate,[],2);
      end
      
      
      function getPSTHKernel(self, varargin)
         %PARAMS
         %  tauAll - [OPTIONAL] list of std of the Gaussian filter in ms
         if isempty(varargin)
            tauAll = round(2.^(1:.5:5)*self.binSize);
         else
            tauAll = varargin{1};
         end
         for tau = 1:length(tauAll)
            win = gausswin(tauAll(tau)*6,tauAll(tau));
            win = win./sum(win);
            
            for cat = 1:self.cats
               for trial = 1:self.trials
                  convTrains(cat,trial,:) = conv(squeeze(self.trains(cat,trial,:)),win,'same');
               end
               PSTH(cat,:) = squeeze(mean(squeeze(convTrains(cat,:,:))));
               ISE(tau,cat,:) = sum(bsxfun(@minus, squeeze(convTrains(cat,:,:)),PSTH(cat,:)).^2,2);
            end
         end
         MISE = mean(ISE,3);
         [val,tau] = min(MISE);
         self.tauOpt = tauAll(tau);
         self.ISE = ISE;
         win = gausswin(tauAll(tau)*6,tauAll(tau));
         win = win./sum(win);
         for cat = 1:self.cats
            for trial = 1:self.trials
               convTrains(cat,trial,:) = conv(squeeze(self.trains(cat,trial,:)),win,'same');
            end
            self.PSTHrate(cat,:) = squeeze(mean(squeeze(convTrains(cat,:,:))));
         end
      end
      
      
      function filterTrains(self,win,sizeTag)
         if nargin==2, sizeTag='same';end
         self.getRawTrains();
         if length(win)==1
            tau = win;
            x =-3*tau:self.binSize:3*tau;
            win = exp(-(x./tau).^2);%gauss;
            win = win./sum(win);
         end
         self.PSTH = [];
         for cat = 1:self.cats
            for trial = 1:self.trials
               newTrains(cat,trial,:) = conv(squeeze(self.trains(cat,trial,:)),win, sizeTag);
            end
            self.PSTH(cat,:) = squeeze(mean(newTrains(cat,:,:),2));
         end
         self.trains = newTrains;
      end
      
      
      function getPSTHBARS(self)
         for cat = 1:self.cats
            data(:,2) = sum(squeeze(self.trains(cat,:,:)),1);
            data(:,1) = (1:length(self.trains))*self.binSize;
            for i = 1:20
               barsP_mat(data,self.trials);
            end
            load -ascii samp_mu;
            %compute the confidence intervals
            for ndx = 1:size(samp_mu,2)
               self.BARSconfint(cat,1,ndx) = prctile(samp_mu(:,ndx),1);
               self.BARSconfint(cat,2,ndx) = prctile(samp_mu(:,ndx),99);
            end
            %get the mean estimate
            self.BARSrate(cat,:) = mean(samp_mu,1);
         end
      end
      
      
      function plotISI(self)
         for cat = 1:self.cats
            subplot(1,self.cats,cat)
            area(self.isiBins,self.isiHist(cat,:));
            ylabel(self.catNames);
         end
         drawnow
      end
      
      
      function h = plotDot(self, varargin)
         %plotDot([startTime endTime])
         %startTime and endTime are in ms and OPTIONAL
         if ~isempty(varargin)
            idx = varargin{1}/self.binSize;
         else
            idx = [min(self.spikeTimes(:)) length(self.PSTH)]/self.binSize;
         end
         h = [];
         hold on;
         for cat = 1:self.cats
            st = reshape(self.spikeTimes(cat,:,:),self.trials,size(self.spikeTimes,3));
            spikeIdx = find(st~=0.0 & st>idx(1) & st<idx(2));
            spikes = st(spikeIdx);
            [trialIdx j]= ind2sub(size(st),spikeIdx);
            y = trialIdx/self.trials+(1.2*cat-1);
            htmp = line([spikes';spikes'],[y';y'+1/(self.trials+2)],'Color','k');
            h = [h; htmp];
         end
         
         hold off;
         axis('tight');
         set(gca,'YLim',[-.5 1.2*self.cats+.5])
         ylabel(self.yLab)
         xlabel('time [ms]');
         set(gca,'YTick',(1:self.cats)*1.2-.5, 'YTickLabel',self.catNames);
         drawnow;
         
      end
      
      function plotPSTH(self,varargin)
         %plotPSTH([startTime endTime])
         %startTime and endTime are in ms and OPTIONAL
         if ~isempty(varargin)
            idx = (varargin{1}(1):varargin{1}(2))/self.binSize;
         else
            idx = (1:length(self.PSTH))/self.binSize;
         end
         col = [0 0 0];
         if length(varargin)>1,col = varargin{2};end
         PSTH = self.PSTH(:,idx);
         if ~iscell(self.catNames), self.catNames = num2cellstr(self.catNames);end
         h = myPlotMult(idx*self.binSize,flipud(PSTH),fliplr(flipud(self.catNames)));
         set(h,'Color',col)
      end
      
      function plotDotPSTH(self)
         self.plotDot()
         hold on
         plotPSTH = self.PSTH./max(self.PSTH(:))/1.2;
         
         for cat = 1:self.cats
            hp(cat) = plot(plotPSTH(cat,:) + cat - 1/(self.trials+2));
         end
         set(hp,'Color','r','LineWidth',1.5);
      end
      
      function [hLine,hError] = plotRate(self,varargin)
         X = 1:length(self.rateMean);
         if nargin>1, X = varargin{1};end
         [hLine,hError] = myErrorBar(X,self.rateMean,self.rateStd);
      end
      
      
      
      function srFlat = flatten(self)
         % concatenates all cats and trials of spikeTimes and stim (if
         % nonempty)
         len = length(self.stim)*self.binSize;
         count = 0;
         newTimes = [];
         for cat = 1:self.cats
            thisTimes = squeeze(self.spikeTimes(cat,:,:));
            thisTimes = thisTimes(:);
            newTimes = [newTimes(:)' thisTimes(thisTimes>0)' + len*(cat-1)];
            
         end
         srFlat = SpikeResp(reshape(newTimes,1,1,[]),self.binSize);
         
         if length(size(self.stim))>2
            for dim = 1:size(self.stim,3)
               newStim(:,dim) = reshape(self.stim(:,:,dim)',1,[]);
            end
            srFlat.setStim(newStim,self.binSize);
         else
            srFlat.setStim(reshape(self.stim',1,[]),self.binSize);
         end
         
      end
      
      function [stReshape,catReshape, trialReshape] = reshapeSpikeTimes(self)
         %makes spikeTimes 2D, returning cat and trial indices
         stReshape = SpikeResp(...
            reshape(permute(self.spikeTimes,[2,1,3]), self.cats*self.trials,[]),...
            self.binSize);
         catReshape = ceil((1:cats*trials)/trials);
         trialReshape = repmat(1:trials,1,cats);
      end
      
      function export2HDF5(self,fileName)
         %make stim 2D for saving
         tmpTimes = reshape(permute(self.spikeTimes,[2,1,3]),self.cats*self.trials,[]);
         hdf5write(fileName,'/spikeTimes',tmpTimes)
         hdf5write(fileName,'/trials',self.trials,'WriteMode', 'append')
         hdf5write(fileName,'/cats',self.cats,'WriteMode', 'append')
         hdf5write(fileName,'/stim',self.stim,'WriteMode', 'append')
      end
      
      function exportTxt(self, fileName)
         if strcmp(self.catNames{1},'')
            self.catNames = num2cellstr(1:self.cats);
         end
         if self.cats==1 || self.trials==1
            catSpikes = squeeze(self.spikeTimes);
            catStim   = squeeze(self.stim);
            save([fileName '_spikes.txt'], 'catSpikes','-ASCII')
            save([fileName '_stim.txt'], 'catStim'  ,'-ASCII')
         else
            for cat = 1:self.cats
               catSpikes = squeeze(self.spikeTimes(cat,:,:));
               catStim   = self.stim(cat,:);
               save([fileName '_spikes_' self.catNames{cat} '.txt'], 'catSpikes','-ASCII')
               save([fileName '_stim_'   self.catNames{cat} '.txt'], 'catStim'  ,'-ASCII')
            end
         end
         
      end
   end
   
end
