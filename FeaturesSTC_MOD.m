classdef FeaturesSTC_MOD < Features
   
   properties
      STA, RTA
      STC, RTC
   end
   
   
   
   methods (Access='public')
      
      function self = FeaturesSTC_MOD(sr,n)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.n = n;
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      %2a. STA/STC
      function [feat, V] = getFeat(self)
         % estimate features for each cat independently
         self.STC = zeros(self.n);
         self.RTC = zeros(self.n);
         self.STA = zeros(1,self.n);
         self.RTA = zeros(1,self.n);
         for cat = 1:self.sr.cats
            catTimes = self.sr.spikeTimes(cat,:,:);
            sr = SpikeResp(catTimes,self.sr.binSize);
            sr.stim = self.sr.stim(cat,:);
            weight(cat) = sum(sr.spikeTimes(:)>0);
            sf = FeaturesSTC(sr,self.n);
            sf.getFeat();
            self.STA = self.STA + weight(cat).*sf.STA;
            self.RTA = self.RTA + sf.RTA;
            self.STC = self.STC + weight(cat).*sf.STC;
            self.RTC = self.RTC + sf.RTC;
            
         end
         % normalize
         self.STA = self.STA./sum(weight);
         self.RTA = self.RTA./self.sr.cats;
         self.STC = self.STC./sum(weight);
         self.RTC = self.RTC./self.sr.cats;
         %sum up stc matrices and perform eigenanalysis
         [feat,V] = eig(self.STC - self.RTC);% get eigs of the difference matrix
         self.feat = [self.STA',feat];
         self.V = [0;diag(V)];
         
         for f = 2:size(self.feat,2)
            [val, idx] = max(abs(self.feat(:,f)));
            self.feat(:,f) = self.feat(:,f) * sign(self.feat(idx,f));
         end
      end
      
      function getNonlinearity1D(self, MODE,varargin)
          %getNonlinearity1D(MODE,[dims],[nBins])
         nBins = 48;
         if nargin>=4
            nBins = varargin{2};
         end         
         tau=0;
         if nargin>=5
            tau = varargin{3};         
         end
         dims = 1:size(self.feat,2);
         if nargin>=3 && not(isempty(varargin{1}))
            dims = varargin{1};
         end
         
         for cat = 1:self.sr.cats            
            catTimes = self.sr.spikeTimes(cat,:,:);
            sr = SpikeResp(catTimes,self.sr.binSize);
            sr.stim = self.sr.stim(cat,:);
            weight(cat) = sum(sr.spikeTimes(:)>0);
            sf = FeaturesSTC(sr,self.n);
            sf.getFeat();
            sf.getNonlinearity1D(MODE,dims,nBins,tau);  
            
            
         end
         
         
      end
      
   end
   
end
