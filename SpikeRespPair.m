classdef SpikeRespPair < handle
   
   properties
      sr1, sr2
      cellLabels, catLabels
      binSize, cats, trials, minLen
      spikeCount_CC, spikeCount_p
      maxLag, spikeRate_lags, spikeRate_C, spikeRate_C_sem
      spikeTime_lags, spikeTime_C, spikeTime_C_sem
      spikeTime_Raw_lags, spikeTime_C_Raw, spikeTime_C_Raw_sem
      spikeTime_Shuffle_lags, spikeTime_C_Shuffle, spikeTime_C_Shuffle_sem
      compoundSTA, compoundSTAraw, compoundNL,
   end
   
   methods (Access='public')
      function self = SpikeRespPair(sr1,sr2,cellLabels)
         %usage:
         %  SpikeRespPair(sr1,sr2,cellLabels)
         %args:
         %  sr1, sr2   - SpikeResp pobjects
         %  cellLabels - 2-cell array with each cell's label
         self.sr1 = sr1;
         self.sr2 = sr2;
         self.binSize = self.sr1.binSize;
         self.cats = self.sr1.cats;
         self.trials = self.sr1.trials;
         self.cellLabels = cellLabels;
         
         self.sr1.getRawTrains();
         self.sr2.getRawTrains();
         
         self.minLen =  min(size(sr1.trains,3),size(sr2.trains,3));
         
      end
      
      function spikeCountCorr(self)
         % estimates the spike-count correlation for each cat
         %
         % PLOT WITH:  plotSpikeCountCorr()
         %
         % usage:
         %     SpikeResponsePair.spikeCountCorr();
         % results are written to
         %     spikeCount_CC (for the coefficient of correlation, 1 x cats)
         %     spikeCount_p (for the p-value for the coefficients, 1 x cats)
         self.sr1.getRate();
         self.sr2.getRate();
         for cat = 1:self.cats
            [cc,p] = corrcoef(self.sr1.rate(cat,:), self.sr2.rate(cat,:));
            self.spikeCount_CC(cat) = cc(2);
            self.spikeCount_p(cat) =  p(2);
         end
      end
      
      
      function spikeRateCorr(self, maxLag)
         % estimates the correlation between the time-dependent rate
         % functions (trial-averaged PSTH) for each cat. This is very close
         % to the shuffle-corrector
         %
         % PLOT WITH:  plotSpikeRateCorr()
         %
         % usage:
         %     SpikeResponsePair.spikeRateCorr([maxLag=50]);
         %     with maxLag given in ms
         % results are written to
         %     spikeRate_C    - for the coefficient of correlation between each cell's PSTH,
         %                      cats x 2*maxLag+1 x 4 (where the last dimensions means all
         %                      combinations, i.e. 1vs1, 1vs2, 2vs1, 2vs2) , 1 x cats)
         %     spikeRate_lags - a vector of time lags - for plotting purposes)
         if nargin==1
            maxLag = 50;
         end
         maxLag = round(maxLag/self.binSize);
         self.spikeRate_C = [];% to avoid dimension mistmatch errors when changing maxLag
         self.sr1.getRawTrains();
         self.sr2.getRawTrains();
         for cat = 1:self.cats
            [self.spikeRate_C(cat,:,:),self.spikeRate_lags] = ...
               xcov([self.sr1.PSTH(cat,1:self.minLen); self.sr2.PSTH(cat,1:self.minLen)]',maxLag,'coeff');
         end
         self.spikeRate_C = self.spikeRate_C * self.binSize/1000;
      end
      
      
      function spikeTimeCorr(self, maxLag)
         % estimates the correlation between the spike times for each cat
         % via trial-wise correlation of the binary spike trains
         %
         % PLOT WITH:  plotSpikeTimeCorr()
         %
         % usage:
         %     SpikeResponsePair.spikeTimeCorr([maxLag=50ms]);
         %     with maxLag given in ms
         % results are written to
         %     spikeRate_C    - for the coefficient of correlation between
         %                      each cell's PSTH, cats x 2*maxLag+1 x 4
         %                      (where the last dimensions means all
         %                      combinations, i.e. 1vs1, 1vs2, 2vs1, 2vs2)
         %     spikeRate_lags - a vector of time lags - for plotting
         %     purposes)
         if nargin==1
            maxLag = 50;
         end
         maxLag = round(maxLag/self.binSize);
         spikeTime_C_single = [];% to avoid dimension mistmatch errors when changing maxLag
         % first make sure that we work on a vanilla (that is not smoothed) rate
         self.sr1.getRawTrains();
         self.sr2.getRawTrains();
         % than correlate all trials
         for cat = 1:self.cats
            for trial = 1:self.trials
               train1 = self.sr1.trains(cat,trial,1:self.minLen);
               train2 = self.sr2.trains(cat,trial,1:self.minLen)   ;
               [spikeTime_C_single(cat,trial,:,1:4),self.spikeTime_lags] = ...
                  xcov([train1(:), train2(:)],maxLag,'coeff');
            end
         end
         spikeTime_C_single = spikeTime_C_single * self.binSize/1000;
         % now take the mean over trials
         self.spikeTime_C = reshape(mean(spikeTime_C_single,2),self.cats,2*maxLag+1,4);
         self.spikeTime_C_sem = reshape(std(spikeTime_C_single,[],2),self.cats,2*maxLag+1,4)/sqrt(self.trials);
      end
      
      function spikeTimeCorrClassic(self, maxLag)
         % estimates the correlation between the spike times for each cat
         % via trial-wise correlation of the binary spike trains
         %
         % PLOT WITH:  plotSpikeTimeCorr()
         %
         % usage:
         %     SpikeResponsePair.spikeTimeCorr([maxLag=50ms]);
         %     with maxLag given in ms
         % results are written to
         %     spikeRate_C    - for the coefficient of correlation between
         %                      each cell's PSTH, cats x 2*maxLag+1 x 4
         %                      (where the last dimensions means all
         %                      combinations, i.e. 1vs1, 1vs2, 2vs1, 2vs2)
         %     spikeRate_lags - a vector of time lags - for plotting
         %     purposes)
         if nargin==1
            maxLag = 50;
         end
         maxLag = round(maxLag/self.binSize);
         spikeTime_C_singleShuffle = [];% to avoid dimension mistmatch errors when changing maxLag
         spikeTime_C_singleRaw = [];% to avoid dimension mistmatch errors when changing maxLag
         % first make sure that we work on a vanilla (that is, non-smoothed) rate
         self.sr1.getRawTrains();
         self.sr2.getRawTrains();
         % then, correlate all trials
         for cat = 1:self.cats
            cnt = 0;
            for trial1 = 1:self.trials
               for trial2 = trial1:self.trials
                  train1 = self.sr1.trains(cat,trial1,1:self.minLen);
                  train2 = self.sr2.trains(cat,trial2,1:self.minLen);
                  if trial1 == trial2 % raw cross-correlogram
                     [spikeTime_C_singleRaw(cat,trial1,:,1:4),self.spikeTime_Raw_lags] = ...
                        xcov([train1(:),train2(:)],maxLag,'coeff');
                  else % shuffle corrector
                     cnt = cnt + 1;
                     [spikeTime_C_singleShuffle(cat,cnt,:,1:4),self.spikeTime_Shuffle_lags] = ...
                        xcov([train1(:),train2(:)],maxLag,'coeff');
                  end
               end
            end
         end
         spikeTime_C_singleRaw = spikeTime_C_singleRaw * self.binSize/1000;
         spikeTime_C_singleShuffle = spikeTime_C_singleShuffle * self.binSize/1000;
         % now take the mean over trials
         self.spikeTime_C_Raw = reshape(mean(spikeTime_C_singleRaw,2),self.cats,2*maxLag+1,4);
         self.spikeTime_C_Shuffle = reshape(mean(spikeTime_C_singleShuffle,2),self.cats,2*maxLag+1,4);
         self.spikeTime_C_Raw_sem = reshape(std(spikeTime_C_singleRaw,[],2),self.cats,2*maxLag+1,4)/sqrt(self.trials);
         self.spikeTime_C_Shuffle_sem = reshape(std(spikeTime_C_singleShuffle,[],2),self.cats,2*maxLag+1,4)/sqrt(self.trials);
         
      end
      
      %%     S . T . A .   F U N C T I O N S
      
      
      function getCompoundSTA(self, maxDel)
         % estimates STA's and nonlinearities based on compound events for the cell pair.
         % Compounds events are coinicident (both cells firing - 11),
         % anti-coincident (only one cell is firing - 10 or 01), silence (neither
         % cell is firing - 00) or independent events (1x or x1). The window size for two spikes being (anti-)
         % coincidence is defined via maxDel. STA's and NL's are computed
         % for each cat
         %
         % PLOT WITH:  plotCompound()
         %
         % usage:
         %     SpikeResponsePair.getCompoundSTA(maxDel);
         %     with maxDel given in bins
         % results are written to
         %     compoundSTA - STA's, time x cats x cmps (cmps in the order: 11, 00, 10, 01, 1x, x1)
         %     compoundNL  - NL's associated with each STA,
         %                   prj's x 2 x cats x cmps (2nd dim: 1 is prj
         %                   values (x) and 2 is NL (y)
         
         % construct delayed compound trains
         compoundTrains = zeros(self.cats,self.minLen-2*maxDel,6);
         delPSTH1 =  self.sr1.PSTH(:,maxDel  + 1 :self.minLen- (maxDel));% crop max del at both ends
         
         for del = -maxDel:maxDel
            delPSTH2 =  self.sr2.PSTH(:,maxDel - del + 1 :self.minLen - (maxDel + del));% shift
            compoundTrains(:,:,1) = compoundTrains(:,:,1) + (delPSTH1 & delPSTH2);% both firing
            compoundTrains(:,:,2) = compoundTrains(:,:,2) + (~delPSTH1 & ~delPSTH2);% none firing
            compoundTrains(:,:,3) = compoundTrains(:,:,3) + (delPSTH1 & ~delPSTH2);% #1 firing, #2 not
            compoundTrains(:,:,4) = compoundTrains(:,:,4) + (~delPSTH1 & delPSTH2);% #1 not firing, #2 does
            compoundTrains(:,:,5) = compoundTrains(:,:,5) + (delPSTH1);% #1 isolated but with delay
            compoundTrains(:,:,6) = compoundTrains(:,:,6) + (delPSTH2);% #2 isolated but with delay
         end
         % for each compound, compute category-wise sta
         for cmp = 1:size(compoundTrains,3)
            for cat = 1:self.cats
               spikeTimes = [];
               catTimes = find(compoundTrains(cat,:,cmp))*self.binSize;
               spikeTimes(:,1,1:length(catTimes)) = catTimes;
               
               clear srTmp sfTmp
               srTmp = SpikeResp(spikeTimes,1);
               srTmp.setStim(self.sr1.stim(cat,:,:),1);
               sfTmp = FeaturesSTC(srTmp,64);
               sfTmp.getFeat();
               sfTmp.getNonlinearity1D(1,1:4,32);
               
               self.compoundSTAraw(:,cat,cmp)  = sfTmp.STA;%the STA
               self.compoundSTA(:,cat,cmp)  = sfTmp.STA./abs(norm(sfTmp.STA));%noralized STA
               self.compoundNL(:,1,cat,cmp) = squeeze(sfTmp.psr1D(2,2,:));% the x axis (prj values)
               self.compoundNL(:,2,cat,cmp) = squeeze(sfTmp.psr1D(2,1,:));% the y axis (NL)
            end
         end
      end
      
      
      %%     P L O T   F U N C T I O N S
      
      
      function plotSpikeCountCorr(self)
         % plots the results of spikecountcorr
         colors = colormap('jet');%
         for cat = 1:self.cats
            h(cat) = subplot(1,self.cats,cat);%keep track of each subplot's handle for later adjustement of axes limits
            hold on
            for t = 1:size(self.sr1.rate,2)-1% do this to let each trial have its color
               plot(self.sr1.rate(cat, t),self.sr2.rate(cat, t),'.','Color',colors(round(1.5*t),:),'MarkerSize',8);
               plot([self.sr1.rate(cat, t) self.sr1.rate(cat, t+1)],[self.sr2.rate(cat, t) self.sr2.rate(cat, t+1)],'-','Color',colors(round(1.5*t),:),'Linewidth',.5);
            end
            plot(self.sr1.rate(cat, t+1),self.sr2.rate(cat, t+1),'.','Color',colors(round(1.5*t),:),'MarkerSize',8);
            hold off
            % annotate the plot...
            title({self.catLabels{cat}, ['r^2=' num2str(self.spikeCount_CC(cat)^2,2) ' | p=' num2str(self.spikeCount_p(cat),2)]})
            axis('square')
            if cat==1,ylabel([self.cellLabels{2} ' rate [Hz]']),end
            if cat==round(self.cats/2)
               xlabel([self.cellLabels{1} ' rate [Hz]'])
               title({[self.cellLabels{1} '-' self.cellLabels{2} ' covariations'], self.catLabels{cat},...
                  ['r^2=' num2str(self.spikeCount_CC(cat)^2,2) ' | p=' num2str(self.spikeCount_p(cat),2)]})
            end
         end
         set(h,'color','none',...
            'XLim',[min(self.sr1.rate(:)) max(self.sr1.rate(:))],...
            'YLim',[min(self.sr2.rate(:)) max(self.sr2.rate(:))])
      end
      
      
      function plotCorr(self,lags,C,yLim)
         % the master plotter
         for cat = 1:self.cats
            h(cat) = subplot(1,self.cats,cat);%keep track of each subplot's handle for later adjustement of axes limits
            plot(lags*self.binSize, squeeze(C(cat,:,:)).*1000/self.binSize,'LineWidth',1.5)
            axis('tight')
            % annotate the plot...
            title(self.catLabels{cat})
            if cat==1,ylabel('corr [Hz]'),end
            if cat==round(self.cats/2)
               xlabel('lag [ms]')
            end
         end
         legend({[self.cellLabels{1} self.cellLabels{1}], [self.cellLabels{1} self.cellLabels{2}], [self.cellLabels{2} self.cellLabels{1}], [self.cellLabels{2} self.cellLabels{2}]})
         legend('boxoff')
         if nargin==3, yLim = [min(C(:)*1000/self.binSize) max(C(:)*1000/self.binSize)];end
         set(h,'color','none','YLim',yLim)
      end
      
      function plotSpikeRateCorr(self)
         % plot results of the spikeRateCorr, with ylimit chosen as to
         % avoid the peak at zero lag
         self.plotCorr(self.spikeRate_lags, self.spikeRate_C, ...
            [min(self.spikeRate_C(:))*1000/self.binSize...
            max(max(max(self.spikeRate_C(:,abs(self.spikeRate_lags)>1/self.binSize,:),[],2)))*1000/self.binSize]);
      end
      
      function plotSpikeTimeCorr(self)
         % plot results of the spikeTimeCorr, with ylimit chosen as to
         % avoid the peak at zero lag
         self.plotCorr(self.spikeTime_lags,self.spikeTime_C,...
            [min(self.spikeTime_C(:))*1000/self.binSize ...
            max(max(max(self.spikeTime_C(:,abs(self.spikeTime_lags)>1/self.binSize,:),[],2)))*1000/self.binSize]);
      end
      
      function plotSpikeRateAndTimeCorr(self)
         % the master plotter
         lags = self.spikeTime_lags;
         corrLabels = {[self.cellLabels{1} self.cellLabels{1}], [self.cellLabels{1} self.cellLabels{2}], [self.cellLabels{2} self.cellLabels{1}], [self.cellLabels{2} self.cellLabels{2}]};
         for cat = 1:self.cats
            for crr = 1:4
               h = mySubPlot(self.cats,4,cat,crr);%keep track of each subplot's handle for later adjustement of axes limits
               plot(lags*self.binSize, squeeze(self.spikeTime_C(cat,:,crr)).*1000/self.binSize,'Color',[.75 0 0],'LineWidth',1.5)
               hold on
               plot(lags*self.binSize, squeeze(self.spikeRate_C(cat,:,crr)).*1000/self.binSize,'k','LineWidth',1.5)
               hold off
               axis('tight')
               if cat==1,
                  title(corrLabels{crr});
               end
               yLim = [min(min(min(self.spikeTime_C(cat,abs(self.spikeTime_lags)>1/self.binSize,crr),[],2)))*1000/self.binSize ...
                  max(max(max(self.spikeTime_C(cat,abs(self.spikeTime_lags)>1/self.binSize,crr),[],2)))*1000/self.binSize];
               if ~isnan(yLim)
                  set(h,'color','none','YLim',yLim)
               end
               if cat==self.cats
                  xlabel('lag [ms]')
               end
               if crr==1,ylabel(self.catLabels{cat}),end
            end
            % annotate the plot...
            
         end
         legend({'spike','PSTH'})
         legend('boxoff')
      end
      
      
      function plotCompound(self)
         cmps = size(self.compoundSTA,3);
         cmpLabels = {'11', '00', '10', '01', '1x' ,'x1'};
         T = -size(self.compoundSTA,1):self.binSize:-self.binSize;
         cols = lines(cmps);
         for cmp = 1:cmps
            % the STA's
            h1(cmp) = mySubPlot(2,cmps,1,cmp);
            plot(T, squeeze(self.compoundSTA(:,:,cmp)),'LineWidth',1.5)
            axis('tight')
            % annotate the plot...
            title(cmpLabels{cmp})
            %if cmp==1,ylabel('dB?'),end
            if cmp==round(cmps/2)
               xlabel('time [ms]')
            end
            
            % the nonlinearities
            h2(cmp) = mySubPlot(2,cmps,2,cmp);
            plot(squeeze(self.compoundNL(:,1,:,cmp)),squeeze(self.compoundNL(:,2,:,cmp)),'LineWidth',1.5)
            axis('tight')
            % annotate the plot...
            if cmp==1,ylabel('rate [kHz]'),end
            if cmp==round(cmps/2)
               xlabel('prj on STA')
            end
         end
         legend(self.catLabels)
         legend('boxoff')
         set(h1,'color','none','YLim',[min(self.compoundSTA(:)) max(self.compoundSTA(:))])
         set(h2,'color','none','YLim',[min(min(min(squeeze(self.compoundNL(:,2,:,:))))) max(max(max(squeeze(self.compoundNL(:,2,:,:)))))])
      end
      
      
      
   end
end


