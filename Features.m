classdef Features < handle
   
   properties
      sr,      n
      SSraw,      Resp
      
      feat, V
      randFeat, randVal, randSTA
      subFeat,  subVal, subSTA, subNonlin1D
      asyFeat, subs, asyError, asyNonlin1D, asyVal
      
      psr1D,psr2D
   end
   
   methods (Abstract)
      getFeat(self);
   end
   
   methods (Access='public')
      
      %1. get spike triggered ensemble
      function getSTE(self)
         % NEW: concatenates stimuli if they are multidimensional (as
         % indiciated by self.sr.stim having a third dimension
         endIdx = min(length(self.sr.stim), length(self.sr.PSTH));
         self.sr.getRawTrains();
         if length(size(self.sr.stim))>2% multidimensional stimuli (e.g. for STRF)
            self.SSraw = zeros(self.sr.cats*endIdx,self.n*size(self.sr.stim,3));
         else % simple 1D (e.g. audio) stims
            self.SSraw = zeros(self.sr.cats*endIdx,self.n);
         end
         self.Resp = zeros(size(self.SSraw,1),1);
         lastEndIdx = 0;
         for cat = 1:self.sr.cats
            if length(size(self.sr.stim))>2
               for dim = 1:size(self.sr.stim,3)
                  SScat(:,:,dim) = makeStimRows(self.sr.stim(cat,:,dim)',self.n,1:endIdx);
               end
            else
               SScat = makeStimRows(self.sr.stim(cat,:)',self.n,1:endIdx);
            end
            self.SSraw(lastEndIdx + (1:endIdx),:,:) = reshape(SScat,[],self.n*size(self.sr.stim,3));
            clear SScat
            thisResp = self.sr.PSTH(cat,1:endIdx)';
            thisResp(1:self.n) = 0;
            self.Resp(lastEndIdx + (1:endIdx)) = thisResp;
            clear thisResp
            lastEndIdx = lastEndIdx + endIdx;
         end
         self.Resp = reshape(self.Resp,1,[])';
         %self.Resp = self.Resp./(sum(self.Resp)+eps);%normalize Resp
      end
      
      function normalizeFeatures(self)
         %normalizes Features to norm 1 and to have a positive extremum
         for fea = 1:size(self.feat,2)
            self.feat(:,fea) = self.feat(:,fea)/norm(self.feat(:,fea));%normalize feats
            [val, idx] = max(abs(self.feat(:,fea)));
            self.feat(:,fea) = self.feat(:,fea)*sign(self.feat(idx,fea));
         end
         if ~isempty(self.subFeat)
            for run = 1:size(self.subFeat,1)
               for fea = 1:size(self.subFeat,3)
                  self.subFeat(run,:,fea) = self.subFeat(run,:,fea)/norm(self.subFeat(run,:,fea));%normalize feats
                  [val, idx] = max(abs(self.subFeat(run,:,fea)));
                  self.subFeat(run,:,fea) = self.subFeat(run,:,fea)*sign(self.subFeat(run,idx,fea));
               end
            end
         end
      end
      
      function getFeatRand(self, varargin)
         if nargin==1
            runs = 10;
         else
            runs = varargin{1};
         end
         randFeat = zeros(runs,size(self.feat,1),size(self.feat,2));
         randVal = zeros(runs,size(self.feat,1)+1);
         % from subset
         fprintf('\n\t')
         oldResp = self.Resp;
         oldFeat = self.feat;
         oldV = self.V;
         parfor run = 1:runs
            sfRand = self;
            sfRand.Resp = self.Resp(randperm(length(self.Resp)));
            sfRand.getFeat();
            randFeat(run,:,:) = sfRand.feat;
            randSTA(run,:) = sfRand.STA;
            randVal(run,:) = sfRand.V;
         end
         fprintf('!\n')
         self.Resp = oldResp;
         self.feat=oldFeat;
         self.V=oldV;
         self.randFeat = randFeat;
         self.randVal = randVal;
         self.randSTA = randSTA;
      end
      
      function getFeatSub(self,subset, varargin)
         %getFeatSub(subset, runs=10)
         if nargin==2
            runs = 10;
         else
            runs = varargin{1};
         end
         
         self.subFeat = zeros(runs,size(self.feat,1),size(self.feat,2));
         self.subVal = zeros(runs,size(self.feat,1)+1);
         % from subset
         thisSR = self.sr;
         for run = 1:runs
            srSub = thisSR.getSubSetResp(subset);
            srSub.setStim(thisSR.stim,thisSR.binSize);
            sfSub = FeaturesSTC(srSub,self.n);
            sfSub.getFeat();
            self.subFeat(run,:,:) = sfSub.feat;
            self.subSTA(run,:) = sfSub.STA;
            sfSub.getNonlinearity1D(2,1:4,32,4);
            self.subNonlin1D(run,:,:,:) = squeeze(sfSub.psr1D(1:4,1:2,:));
            self.subVal(run,:) = sfSub.V;
         end
         
      end
      
      function getNonlinearity1D(self, MODE,varargin)
         %getNonlinearity1D(MODE,[dims, nBins, TAU])
         % ARGS:
         %  MODE - 1 hist, 2 kde, 3 nix, 4 RoG
         %  dims
         %  nBins
         nBins = 48;
         if nargin>=4
            nBins = varargin{2};
         end
         resp = self.Resp;
         if nargin>=5
            tau = varargin{3};
            if tau>0
               x =-3*tau:self.sr.binSize:3*tau;
               win= exp(-(x./tau).^2);%gauss;
               win = win./sum(win);
               resp = conv(resp,win,'same');
            end
         end
         if sum(resp)~=1
            resp = resp/sum(resp);
         end
         
         psr1D = zeros(size(self.feat,2),4,nBins);
         dims = 1:size(self.feat,2);
         if nargin>=3 && not(isempty(varargin{1}))
            dims = varargin{1};
         end
         
         for fea = 1:length(dims)
            featNorm(:,dims(fea)) = self.feat(:,dims(fea))/norm(self.feat(:,dims(fea)));%normalize feats
         end
         psF = self.SSraw*featNorm;%project all stim onto features p(s)
         %psF = bsxfun(@times, psF, var(self.SSraw(:,1))./var(psF));%normalize asin Nagel et al 2006
         for fea = 1:length(dims)
            ps = psF(:,dims(fea));
            maxPrj = 3.5*std(ps);%max(psF)+.0000001;
            minPrj = -3.5*std(ps);%min(psF)+.0000001;
            bins = linspace(minPrj,maxPrj,nBins);
            nps = normpdf(bins,mean(ps),std(ps));
            
            %% histogram
            if MODE == 1 % raw histogram
               %hps  = histcw(ps,ones(size(resp))./length(resp),bins);
               hpfs = histcw(ps,resp,bins);
               %hps = hps(1:end-1);
               hpfs = hpfs(1:end-1);
            elseif MODE==2 % kernel density estimation
               ps(ps<minPrj & ps>maxPrj) = [];
               % new mode
               %hps  = ksdensity(ps',bins);
               hps  = evaluate(kde(double(ps'),'rot',[] ,'Epan'),double(bins));
               %hpfs = ksdensity(ps',bins,'weights',resp');
               hpfs = evaluate(kde(double(ps'),'rot',double(resp'),'Epan'),double(bins));
               %version from STC paper
               %[bw,hpfs,x2]  = kde(ps(resp>0),4*nBins, [minPrj maxPrj]);
               %hpfs = interp1(x2,hpfs,bins,'nearest','extrap');
            elseif MODE==3 % splines
               [val, idx] = sort(ps);
               hpfs = csaps(double(ps(idx)),double(resp(idx)),.99,bins);
               hps = ones(size(hpfs));
            elseif MODE==4 % adaptive bin width 
               [bins, hpfs] = tfAdapt(double(ps),double(resp),nBins);
               hpfs = hpfs';
               nps = ones(size(hpfs));               
            else
               %% fit gaussian
               Mhpfs = resp'*ps;%STA
               Shpfs = sqrt(sum(resp'*ps.^2 - Mhpfs^2));
               hpfs = normpdf(bins,Mhpfs,Shpfs);
            end
            
            %normalize distributions
            hps = nps./sum(nps);%p(s)
            % hps = hps./sum(hps);%p(s)
            hpfs = hpfs./sum(hpfs);%p(f|s)
            
            psr1D(dims(fea),1,:) = hpfs./hps;%divide p(f|s)/p(s)
            psr1D(dims(fea),2,:) = bins;
            psr1D(dims(fea),3,:) = hpfs;
            psr1D(dims(fea),4,:) = hps;
            
         end
         self.psr1D = psr1D;
      end
      
      function varargout = getNonlinearity2D(self, MODE,dim,varargin)
         %getNonlinearity1D(MODE,dims,[nBins])
         
         varargout{1} =[];
         
         nBins =32;
         if nargin>3
            nBins = varargin{1};
         end
         resp = self.Resp;
         if nargin>4
            tau = varargin{2};
            if tau>0
               x =-3*tau:self.sr.binSize:3*tau;
               win= exp(-(x./tau).^2);%gauss;
               win = win./sum(win);
               resp = conv(resp,win,'same');
            end
         end
         for fea = 1:length(dim)
            featNorm(:,fea) = self.feat(:,dim(fea))/norm(self.feat(:,dim(fea)));%normalize feats
         end
         psF = self.SSraw*featNorm;%project all stim onto features p(s)
         
         maxPrj = 3.5*std(psF);%max(psF)+.0000001;
         minPrj = -3.5*std(psF);%min(psF)+.0000001;
         bins1 = linspace(minPrj(1),maxPrj(1),nBins);
         bins2 = linspace(minPrj(2),maxPrj(2),nBins);
         [X1,X2] = meshgrid(bins1,bins2);
         mpsF = mean(psF);mpsF(isnan(mpsF))=0;
         cpsF = cov(psF); cpsF(isnan(cpsF))=0;
         nps = mvnpdf([X1(:) X2(:)],mpsF,cpsF);
         nps = reshape(nps,length(bins1),length(bins2))';
         
         if MODE==1% raw histogram
            hpfs= flipud(hist3(psF(resp>0,:),'Edges',{bins1,bins2}));
         elseif MODE==2%kernel density estimation
            hps  = evaluate(kde(double(psF'),'rot',[] ,'Epan'),double([X1(:) X2(:)]'));
            hps  = reshape(hps,[nBins,nBins]);
            hpfs = evaluate(kde(double(psF'),'rot',double(resp'),'Epan'),double([X1(:) X2(:)]'));
            hpfs = reshape(hpfs,[nBins,nBins]);
            
            %[tmp,hpfs] = kde2d(fliplr(psF(resp>0,:)),nBins,minPrj, maxPrj);
            hpfs(hpfs<0)=0;
         elseif MODE==3
            %[~, idx] = sort(psF);
            %hpfs = csaps({double(psF(idx,1)) double(psF(idx,2))},double(resp(idx)),.99,double([X1(:) X2(:)]'));
            %hps = ones(size(hpfs));
            
         else %ratio of gaussians
            hps = nps;
            Mhpfs = resp'*psF;
            psFCenter = bsxfun(@minus,psF,Mhpfs);%psF-repmat(Mhpfs,size(psF,1),1);
            Shpfs = psFCenter'*bsxfun(@times,psFCenter,resp);%psFCenter'*(psFCenter.*repmat(resp,1,length(dim)));
            Shpfs = .5*(Shpfs + Shpfs');
            hpfs = mvnpdf([X1(:) X2(:)],Mhpfs,Shpfs);
            hpfs = reshape(hpfs,length(bins1),length(bins2))';
            if nargout>0, varargout{1} = Mhpfs;varargout{2} = Shpfs;end
         end
         a = sqrt(X1.^2 + X2.^2 );
         outIdx = find(a>max(maxPrj));
         hpfs(outIdx)=0;
         %hps = nps;
         hps(outIdx)=0;
         hps = (hps./sum(hps(:)))';%p(s)
         hpfs = (hpfs./sum(hpfs(:)))';%p(f|s)
         
         if MODE~=3
            psr2D(1,:,:) = hpfs./hps;%divide p(f|s)/p(s)
         else
            psr2D(1,:,:) = hpfs;%direct ratio estimator
         end
         psr2D(2,1,:) = bins1;
         psr2D(2,2,:) = bins2;
         psr2D(2,3,1) = 1;%mean(std(psF));
         psr2D(3,:,:) = hpfs;
         psr2D(4,:,:) = hps;
         
         self.psr2D = psr2D;
      end
      
      function plotFeat(self, N)
         
         for n = 1:N
            mySubPlot(2,N,1,n);
            plot(self.feat(:,n),'k','LineWidth',1.5)
            title(num2str(self.V(n),2))
            axis('tight')
            
            mySubPlot(2,N,2,n);
            hold on
            plot(squeeze(self.psr1D(n,2,:)), squeeze(self.psr1D(n,3,:)),'Color',[.6 .6 .6])
            plot(squeeze(self.psr1D(n,2,:)), squeeze(self.psr1D(n,4,:)),'k')
            plot(squeeze(self.psr1D(n,2,:)), squeeze(self.psr1D(n,1,:))./nansum(squeeze(self.psr1D(n,1,:))),'r','LineWidth',1.5)
            hold off
            axis('tight')
         end
      end
      
      function plotFeatVal(self, N)
         
         mySubPlot(1,N+1,1,1);
         plot(sort(self.V,'descend'),'.-k','LineWidth',1.5,'MarkerSize',18);
         axis('tight')
         title('eigenvalues')
         ylabel('\sigma')
         t = -size(self.feat,1)+1:1:0;
         for n = 1:N
            mySubPlot(1,N+1,1,n+1);
            plot(t,self.feat(:,n),'k','LineWidth',1.5)
            title(num2str(self.V(n),2))
            axis('tight')
         end
      end
      
      function plotSubRandFeat(self,N)
         %USAGE plotSubRandFeat(N)
         %plots subset and random features as well as the values
         subplot(1,N+2,[1 2])
         hold on
         h1 = plot(abs(self.randVal)','.','Color',[.6 .6 .6]);
         h2 = plot(abs(self.subVal)','.','Color',[1 .5 .2]);
         axis('tight')
         h3 = plot(abs(self.V),'.k');
         
         sigV = find(abs(self.V)>prctile(abs(self.randVal),99)');
         plot(sigV,abs(self.V(sigV)),'or');
         hold off
         legend([h1(1),h2(1),h3(1)],{'RND','SUB','REAL'},'Location','NorthEast');
         
         for f=1:N
            subplot(1,N+2,2+f)
            hold on
            errorbar(mean(self.subFeat(:,:,f)),...
               std(self.subFeat(:,:,f))/sqrt(size(self.subFeat,1)),...
               'k','LineWidth',1.0);
            axis('tight')
            hold off
         end
      end
      
      function plotNonlin2D(self)
         X1 = self.psr2D(2,1,:)/3.5;X1 = double(X1(:)/self.psr2D(2,3,1)/3.5);
         X2 = self.psr2D(2,1,:)/3.5;X2 = double(X2(:)/self.psr2D(2,3,1)/3.5);
         
         subplot(1,3,1)
         myPcolor(X1, X2, double(squeeze(self.psr2D(4,:,:))));
         title('p(s)');colorbar(),axis('square')
         set(gca,'XTick',1:8:200,'XTickLabel',round(X1(1:8:end)))
         set(gca,'YTick',1:8:200,'YTickLabel',round(X1(1:8:end)))
         
         subplot(1,3,2)
         myPcolor(X1, X2, double(squeeze(self.psr2D(3,:,:))));
         title('p(s|r)');colorbar(),axis('square')
         set(gca,'XTick',1:8:200,'XTickLabel',round(X1(1:8:end)))
         set(gca,'YTick',1:8:200,'YTickLabel',round(X1(1:8:end)))
         subplot(1,3,3)
         myPcolor(X1, X2, double(squeeze(self.psr2D(1,:,:))));
         title('p(r|s)');colorbar(),axis('square')
         set(gca,'XTick',1:8:200,'XTickLabel',round(X1(1:8:end)))
         set(gca,'YTick',1:8:200,'YTickLabel',round(X1(1:8:end)))
         
      end
      
      function [sigFeat, sigV] = checkSignificanceByVal(self)
         
         %          self.generateRandomFeatures(100);
         %          % check for significant STA
         %          srf = self.randFeat(:,:,1);
         %          ci = mean(abs(prctile(srf(:),[.5 99.5])));
         %          if sum(abs(self.feat(:,1))>ci)>10/self.n
         %             sigFeat = self.feat(:,1);
         %          end
         %
         %          % assume a single signficant eigVal
         %          ci = prctile(self.randVal(:),[.5 99.5]);
         %          sigFeat = [sigFeat,self.feat(:,1+find(self.V<ci(1) | self.V>ci(2)))];
         %          sigV    = self.V(self.V<ci(1) | self.V>ci(2) );
         
      end
      
      function getAsympError(self)
         
         self.subs = 1./(10:-1:2);
         runs = 10;
         fprintf('\n')
         for sub = 1:length(self.subs)
            fprintf('.')
            self.getFeatSub(self.subs(sub),runs);
            self.asyFeat(sub,:,:,:) = self.subFeat;
            self.asyNonlin1D(sub,:,:,:,:) = self.subNonlin1D;
            self.asyVal(sub,:,:) = self.subVal;
            for run = 1:runs
               cc = corrcoef(self.feat(:,1),self.asyFeat(sub,run,:,1));
               self.asyError(sub,run) = cc(2);
            end
         end
         fprintf('!\n')
      end
      
      function plotAsympError(self)
         plot(log(self.subs),1-mean(self.asyError,2),'.-k')
      end
      
      function st = saveStruct(self)
         fn = fieldnames(self);
         for f = 1:length(fn)
            if any(strcmp(fn{f}, {'SSraw','RRraw','RRrawDelay', 'basisPrj'})) || isempty(self.(fn{f}))
               continue
            end
            st.(fn{f}) = self.(fn{f});
         end
      end
      
   end
end
