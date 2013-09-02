classdef FeaturesGLMsparse < Features
   
   properties
      stim, resp, spikeIdx
      basis, whitener
      fit
      pred, perf
      
   end
   
   
   methods (Access='public')
      
      function self = FeaturesGLMsparse(sr, n)
         if isa(sr,'SpikeResp')
            self.n = n;
            self.sr = sr;
            self.sr.getRawTrains();

            self.getSTE();
            self.Resp = squeeze(sum(self.sr.trains,2));
            self.basis = get1DLaplacianPyramidBasis(n,4,.5,2.5);
            self.SSraw = [self.SSraw*self.basis];
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      %2a. STA/STC
      function getFeat(self)
         %%
%          self.RRraw = makeStimRows(self.Resp, dt*self.n+delay);
%          self.RRraw = self.RRraw(:,1:end-delay);
%          self.RRraw = resample(self.RRraw',1,dt,100)';
         %%
         X = self.SSraw;
         self.whitener = diag(1./std(X,[],1));
         X = X*self.whitener; %Whiten to standard deviation = 1 (X*B*D)
         U = ones(length(self.Resp),1);
         y = self.Resp;
         
         %self.fit = cvglmfitsparseprior(y,X,U,getcvfolds(length(y),5),'modeltype','logisticr','modelextra',self.sr.trials);
         self.fit = cvglmfitsparseprior(y,X,U,getcvfolds(length(y),5),'modeltype','ls');
         feat = self.basis*self.whitener*self.fit.w;
         self.feat = [self.fit.u; feat];
         self.pred = X*self.fit.w + U*self.fit.u;
         self.perf = rsq(self.pred, self.Resp);
      end
   end
end
