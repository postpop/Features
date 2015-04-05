classdef FeaturesGLMsparse < Features
   
   properties
      stim, resp, spikeIdx
      basis, whitener
      fit
      pred, perf
      RRraw, RRrawDelay, RRrawDeltaT
      h0, h1
      
   end
   
   
   methods (Access='public')
      
      function self = FeaturesGLMsparse(sr, n)
         if isa(sr,'SpikeResp')
            self.n = n;
            self.sr = sr;
            self.sr.getRawTrains();
            
            self.RRrawDelay = 2;
            self.RRrawDeltaT = 3;
            self.getSTE();
            self.Resp = squeeze(sum(self.sr.trains,2));
            %self.basis = get1DLaplacianPyramidBasis(n,4,.5,2.5);
            self.basis = get1DLaplacianPyramidBasis(n,4,1,2.5);
            self.SSraw = [self.SSraw*self.basis];
            self.RRraw = makeStimRows(self.Resp, self.RRrawDeltaT*self.n+self.RRrawDelay);
            self.RRraw = self.RRraw(:,1:self.RRrawDeltaT:end-self.RRrawDelay);
            self.RRraw = self.RRraw*self.basis; 

         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      %2a. STA/STC
      function getFeat(self)
         %%
         X = [self.SSraw self.RRraw];
         self.whitener = diag(1./std(X,[],1));
         X = X*self.whitener; %Whiten to standard deviation = 1 (X*B*D)
         U = ones(length(self.Resp),1);
         y = self.Resp;
         
         %self.fit = cvglmfitsparseprior(y,X,U,getcvfolds(length(y),5),'modeltype','logisticr','modelextra',self.sr.trials);
         self.fit = cvglmfitsparseprior(y,X,U,getcvfolds(length(y),5),'modeltype','ls');
         feat = [self.basis self.basis]*self.whitener*self.fit.w;% unwhiten and project onto basis
         self.feat = [self.fit.u; feat];% add bias term
         self.h0 = self.fit.u;% bias term
         self.h1 = feat;% add bias term
         self.pred = X*self.fit.w + U*self.fit.u;% predict
         self.perf = rsq(self.pred, self.Resp);% performance
      end
   end
end

