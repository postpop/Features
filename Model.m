classdef Model < handle
   
   properties
      fit
      fitTau, tauAll
      whichFeat
      srTest, sfTest, sf
      pred, resp, trains
      freq, Csr, Crr, I_lb, I_ub
      TAU, MODE, nBins
   end
   
   methods (Abstract)
      pred = runModel(self)
   end
   
   
   
   methods (Access='public')
      
      function self = Model(srTest)
         
      end
      
      %% quantify model fit
      function fit = CC(self)
         % computes simple crosscorrelation coefficient between measured
         % and predicted response
         fit = corrcoef(self.pred, self.resp);
         fit = fit(2);
         self.fit = fit;
      end
      
      function fit = CCdebias(self)
         % computes debiased crosscorrelation coefficient between measured
         % and predicted response after Hsu et al. 2004
         
         resp1 = reshape(mean(self.trains(:,1:2:end,:),2),self.srTest.cats,size(self.trains,3));%times2trains(self.srTest.spikeTimes(:,1:2:end,:), self.srTest.binSize);
         resp2 = reshape(mean(self.trains(:,2:2:end,:),2),self.srTest.cats,size(self.trains,3));%times2trains(self.srTest.spikeTimes(:,2:2:end,:), self.srTest.binSize);
         len = min(length(resp1),length(resp2));
         resp1 = resp1(1:len);
         resp2 = resp2(1:len);
         [cx,lag] = xcorr((resp1+resp2)/2,self.pred,20);
         [val,idx] = max(cx);
         resp1 = resp1(21+lag(idx):end-20);
         resp2 = resp2(21+lag(idx):end-20);
         pred = self.pred(21:end-20-lag(idx));
         len = min([length(resp1),length(resp2),length(pred)]);
         resp1 = resp1(1:len);
         resp2 = resp2(1:len);
         pred = pred(1:len);
         
         ntrials_proper =  self.srTest.trials;
         
         cc_two_halves_best1 = diag(corrcoef(resp1, resp2),1).^2;
         temp=(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_best1))/2;
         cc_two_halves_best=1./(temp+1);
         cc_two_halves_best=sqrt(cc_two_halves_best);
         
         cc_spike_pre_best1 = diag(corrcoef((resp1+resp2)/2,pred),1).^2;
         cc_spike_pre_best=cc_spike_pre_best1*(1+sqrt(1/cc_two_halves_best1))/(-ntrials_proper+ntrials_proper*sqrt(1/cc_two_halves_best1)+2);
         cc_spike_pre_best=sqrt(cc_spike_pre_best);
         
         fit = cc_spike_pre_best./cc_two_halves_best;
         self.fit = fit;
         
      end
      
      function fit = CCdebias1(self)
         % computes debiased crosscorrelation coefficient between measured
         % and predicted response after Petersen et al. 2008
         resp1 = reshape(mean(self.trains(:,1:2:end,:),2),self.srTest.cats,size(self.trains,3));%times2trains(self.srTest.spikeTimes(:,1:2:end,:), self.srTest.binSize);
         resp2 = reshape(mean(self.trains(:,2:2:end,:),2),self.srTest.cats,size(self.trains,3));%times2trains(self.srTest.spikeTimes(:,2:2:end,:), self.srTest.binSize);
         len = min(length(resp1),length(resp2));
         resp1 = resp1(1:len);
         resp2 = resp2(1:len);
         [cx,lag] = xcorr((resp1+resp2)/2,self.pred,20);
         [val,idx] = max(cx);
         resp1 = resp1(21+lag(idx):end-20);
         resp2 = resp2(21+lag(idx):end-20);
         pred = self.pred(21:end-20-lag(idx));
         len = min([length(resp1),length(resp2),length(pred)]);
         resp1 = resp1(1:len);
         resp2 = resp2(1:len);
         pred = pred(1:len);
         
         cR1R2 = cov(resp1,resp2);
         cRP = cov((resp1+resp2)/2,pred);
         fit = cRP(2)/sqrt(cR1R2(2)*cRP(4));
         self.fit = fit;
         
      end
      
      function mse = MSE(self)
         % computes mse error between measured and predicted response
         mse = mean(mean((self.pred - self.resp).^2));
      end
      
      function fit = VAR(self)
         % computes fraction of variance explained by the model,
         % see Pillow, 2005
         fit = 1-mean((self.pred-self.resp).^2/(self.resp-mean(self.resp)).^2);
         self.fit = fit;
      end
      
      function fit = MSEdebias(self)
         % from Machens et al. 2004 - NOT WORKING
         trials = self.srTest.trials;
         for trial = 1:trials
            train = self.trains(1,trial,:);
            sE(trial,:) = mean((train(:)'- self.pred).^2);% MSE
            sR(trial,:) = mean(train(:).^2);
         end
         sE = mean(sE);
         sR = mean(sR);
         sN = trials/(trials-1) * (sR - mean(self.resp.^2));
         
         self.fit = (sR - sE)/(sR - sN);
         fit = self.fit;
      end
      
      function fit = coherence(self)
         nfft = 512;
         trials = self.srTest.trials;
         trains = squeeze(self.srTest.trains(1,:,1:length(self.pred)));
         [PSmtx, self.freq]=mtcsd([self.pred; trains]',nfft,1000/self.srTest.binSize,[],[],[],[],6);
         Crr = zeros(trials*(trials-1)/2,length(self.freq));
         Psr = zeros(trials,length(self.freq));
         Prr = Psr;
         Prs = Psr;
         Pss = PSmtx(:,1,1);Pss = Pss(:);
         
         count = 0;
         for r1 = 2:trials+1
            P11 = PSmtx(:,r1,r1);P11 = P11(:);
            Prr(r1-1,:) = P11(:);
            
            Ps1 = PSmtx(:,r1,1);
            P1s = PSmtx(:,1,r1);
            Psr(r1-1,:) = P1s(:);
            Prs(r1-1,:) = Ps1(:);
            
            for r2 = r1+1:trials+1
               count = count+1;
               P12 = PSmtx(:,r1,r2);P12 = P12(:);
               P21 = PSmtx(:,r2,r1);P21 = P21(:);
               P22 = PSmtx(:,r2,r2);P22 = P22(:);
               Crr(count,:) = ((P12.*P21)./(P11.*P22));
            end
         end
         
         Psr = mean(Psr)';
         Prs = mean(Prs)';
         Prr = mean(Prr)';
         self.Csr = (Psr.*(Prs))./(Prr.*Pss);
         self.Crr = mean(sqrt(Crr),1);
         self.I_lb = -log2(1-self.Csr);
         self.I_ub = -log2(1-self.Crr);
         self.fit = sum(self.I_lb(1:63));
         fit = self.fit;
      end
      
      function fit = likelihood(self)
         % from Pillow et al.
         logPred = log2(self.pred);
         logPred(isinf(logPred)) = 0;% avoid inf errors!!
         self.fit = nansum(logPred(self.resp>0)) - self.srTest.binSize/1000*nansum(self.pred);
         fit = self.fit;
      end
      
      function [fit] = info(self)
         
         % estimate raw info
         rate = self.pred;
         rMean = nanmean(rate);
         rNorm = rate./rMean;
         rawInfo = nanmean(rNorm.*log2(rNorm));
         
         % bias
         lens = [0.5:.05:0.9];
         biasInfo = zeros(size(lens));
         for l = 1:length(lens)
            len = lens(l);
            rate = self.pred(:,1:floor(len*size(self.pred,2)));
            rMean = nanmean(rate);
            rNorm = rate./rMean;
            biasInfo(l) = nanmean(rNorm.*log2(rNorm));
         end
         lenFit = fitPoly(1./lens,biasInfo,1,[0 1./lens]);
         biasLen = lenFit(1);%diff(lenFit(1:2));
         if biasLen<0, biasLen=0;end
         self.fit = rawInfo - biasLen;
         fit = self.fit;
         
         % from Fairhall et al. 2006
%          try
%             fullInfo = infoSpike(self.srTest);
%          catch ME
%             disp(ME)
%             fullInfo = rawInfo;
%          end
                    
         
%          if length(self.whichFeat)==1
%             %get raw info
%             self.sf.getNonlinearity1D(self.MODE, self.whichFeat, self.nBins, 0);
%             prs = squeeze(self.sf.psr1D(self.whichFeat,3,:));
%             pr  = squeeze(self.sf.psr1D(self.whichFeat,4,:));
%             rawInfo = nansum(prs.*log2(prs./pr));
%             %estimate bias
%             [val,idx] = min(abs(self.sf.V(2:end)));
%             biasFeat = idx+[1];
%             self.sf.getNonlinearity1D(self.MODE, biasFeat, self.nBins, 0);
%             prs = squeeze(self.sf.psr1D(biasFeat,3,:));
%             biasInfo = nansum(prs.*log2(prs./pr));
%          elseif length(self.whichFeat)==2
%             %get raw info
%             self.sf.getNonlinearity2D(self.MODE, self.whichFeat, self.nBins, 0);
%             prs = squeeze(self.sf.psr2D(3,:,:))+eps;
%             pr  = squeeze(self.sf.psr2D(4,:,:))+eps;
%             rawInfo = sum(sum(prs.*log2(prs./pr)));
%             %estimate bias
%             [val,idx] = min(abs(self.sf.V(2:end)));
%             biasFeat = idx+[1 2 ];
%             tmpFeat = self.sf.feat;
%             self.sf.feat(:,biasFeat(1)) = zeros(self.sf.n,1);self.sf.feat(10,biasFeat(1)) = 1;
%             self.sf.feat(:,biasFeat(2)) = zeros(self.sf.n,1);self.sf.feat(end-10,biasFeat(2)) = 1;
%             self.sf.getNonlinearity2D(self.MODE, biasFeat, self.nBins, 0);
%             self.sf.feat = tmpFeat;
%             prs = squeeze(self.sf.psr2D(3,:,:))+eps;
%             self.sf.getNonlinearity2D(self.MODE, biasFeat, self.nBins, 0);
%             biasInfo = sum(sum(prs.*log2(prs./pr)));
%          else
%             disp('only allowed for 1D and 2D STE based models!!')
%             rawInfo = [];
%             biasInfo = [];
%          end
%          modelInfo = (rawInfo - biasInfo);
%          
%          self.fit = modelInfo/fullInfo;
%          fit = self.fit;
      end
      
      function fit = SNR(self)
         %from Geffen et al. 2009
         residual = mean((self.pred - self.resp).^2);%MSE
         signal = mean((self.resp - mean(self.resp)).^2);%signal variance
         %noise = mean(mean(bsxfun(@minus, squeeze(self.srTest.trains(1,:,1:length(self.resp))),self.resp).^2));%trial2trial var
         fit = 1-residual./signal;
         self.fit = fit;
      end
      
      function fit = metric(self)
         len = min(length(self.srTest.trains),length(self.pred));
         self.tauAll = 2.^(-2:8);
         for tau = 1:length(self.tauAll)
            self.srTest.filterTrains(self.tauAll(tau));
            for cat = 1:self.srTest.cats
               testTrains = squeeze(self.srTest.trains(cat,:,1:len));
               dist = pdist([self.pred(1:len);testTrains]);
               noise(tau,cat) = mean(dist(self.srTest.trials + 1:end));%mean inter-Response variability
               error(tau,cat) = mean(dist(1:self.srTest.trials));
            end
         end
         self.fitTau = error./noise;
         self.fit = min(mean(self.fitTau,2));
         fit = self.fit;
      end
      
      %% optimize nonlinearity parameters
      function [optModel, modelPerf, optMODE, optnBins] = optimizeParams(self, bins)
         % optimize model fit through the nonlinearity by changing number
         % of bins, estimation method
         %USAGE
         %  [modelPerf] = optimizeParamsTau(self, bins)
         %PARAMS
         % bins     - array of number of bins for the nonlin
         %RETURNS
         % modelPerf   - model performance as a function of bins and method
         
         %created 20090707 JC
         mEdObj = metaclass(self);
         modelPerf = zeros(3,length(bins));
         for md = [1 2 4]
            for nBin = 1:length(bins)
               tmpModel = eval([mEdObj.Name '(self.srTest,self.sf,self.whichFeat,md,bins(nBin),self.TAU)']);
               tmpModel.runModel;
               modelPerf(md,nBin) = tmpModel.CCdebias;
            end
         end
         [v1,i1] = max(modelPerf(:));
         [x1,y1] = ind2sub(size(modelPerf),i1);
         optMODE = x1;
         optnBins = bins(y1);
         optModel = eval([mEdObj.Name '(self.srTest,self.sf,self.whichFeat,optMODE, optnBins,self.TAU)']);
         optModel.runModel();
      end
      
      
      function [optModel, modelPerf, selfPerf] = optimizeParamsTau(self, bins,tauAll)
         % optimize model fit through the nonlinearity by changing number
         % of bins, estimation method, and smoothing of psth.
         %USAGE
         %  [modelPerf, selfPerf] = optimizeParamsTau(self, bins,tauAll)         %
         %PARAMS
         % bins     - number of bins for the nonlin
         % tauAll   - width of the gaussian filter with which to smooth
         %            the psth
         %RETURNS
         % modelPerf   - model performance as a function of tau
         % selfPerf    - full performance for each tau (over methods to
         %               estimate nonlin and #bins
         
         %created 20090707 JC
         %          self.TAU = 0;
         %          self.srTest.getRawTrains();
         %          [tmpModel, selfPerf method nBin] = self.optimizeParams(bins);
         %          for tauIdx = 1:length(tauAll)
         %             self.TAU = tauAll(tauIdx);
         %             self.srTest.getRawTrains();
         %             if self.TAU>0
         %                self.srTest.filterTrains(self.TAU);
         %             end
         %             self.trains = self.srTest.trains;
         %             self.resp = self.srTest.PSTH;
         %             self.srTest.getRawTrains();
         %             self.runModel();
         %             modelPerf(tauIdx) = self.CCdebias;
         %          end
         %          [val,idx] = max(modelPerf);%fid optimal model
         %          tau = tauAll(idx);
         %          optMODE = method;
         %          optnBins = nBin;
         %          optTAU = tau;
         %          %%
         %          mEdObj = metaclass(self);
         %          optModel = eval([mEdObj.Name '(self.srTest,self.sf,self.whichFeat,optMODE,optnBins,optTAU)']);
         %          optModel.runModel();
         %
         
         % OLD AND WORKING
         mEdObj = metaclass(self);
         for tauIdx = 1:length(tauAll)
            self.TAU = tauAll(tauIdx)
            self.srTest.getRawTrains();
            [tmpModel, selfPerf(tauIdx,:,:) method(tauIdx) nBin(tauIdx)] = self.optimizeParams(bins);
            modelPerf(tauIdx) = tmpModel.CC;
         end
         [val,idx] = max(modelPerf);%fing optimal model
         tau = tauAll(idx);
         optMODE = method(idx);
         optnBins = nBin(idx);
         optTAU = tau;
         selfPerf = squeeze(selfPerf(idx,:,:));
         %%
         optModel = eval([mEdObj.Name '(self.srTest,self.sf,self.whichFeat,optMODE,optnBins,optTAU)']);
         optModel.runModel();
      end
      
      
      function spikeTimes = getTrainsPoisson(self,trials)
         %spikeTimes = getTrainsPoisson(self,trials)
         %  generates spike trains from the estimated rate according to a
         %  Poisson process
         %PARAMS:
         %  trials
         %RETURNS:
         %  spiketimes - cats x trials x times
         if self.srTest.binSize~=1
            error('sampling rate MUST be 1000Hz! binsize MUST be 1ms!')
         end
         %estimate absolute refractoriness from ISI hist
         self.srTest.getISI
         refractAbsolute = 1;%self.srTest.isiBins(find(self.srTest.isiHist>0,1,'first'));
         
         spikeTimes = zeros(1,trials,0);
         for trial = 1:trials
            thisTimes = generatePoisson(self.pred,refractAbsolute);
            spikeTimes(1,trial,1:length(thisTimes)) = thisTimes;
         end
      end
      
      function spikeTimes = getTrainsLIF(self,trials)
         %spikeTimes = getTrainsPoisson(self,trials)
         % NOT IMPLEMENTED YET!!!
         %  generates spike trains by feeding the estimated rate into a LIF
         %PARAMS:
         %  trials
         %RETURNS:
         %  spiketimes - cats x trials x times
         if self.srTest.binSize~=1
            error('not implemented yet')
         end
         spikeTimes = zeros(1,trials,0);
         %         for trial = 1:trials
         %             thisTimes =
         %             spikeTimes(1,trial,1:length(thisTimes)) = thisTimes;
         %          end
      end
      
      function plotPred(self)
         hold on
         plot(-self.srTest.stim,'k')
         plot(self.resp,'r')
         plot(self.pred,'g')
         hold off
      end
      
   end
end
