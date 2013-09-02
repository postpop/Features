classdef FeaturesMLquadBasisHist < FeaturesML
   
   
   methods (Access='public')
      
      function self = FeaturesMLquadBasisHist(sr,n)
         self = self@FeaturesML(sr,n);
         self.RRrawDelay = 5;
         self.getSTEML()
      end
      
      function getSTEML(self)
         self.getSTE();             % get standard STE
         self.RRraw = makeStimRows(self.Resp, self.n+self.RRrawDelay);
         self.RRraw = self.RRraw(:,1:end-self.RRrawDelay);
         self.getBasis2D();         % get basis
         %self.getBasis1D();
         %self.respBasis = self.basis1D;
         self.getBasis1D();
         self.prj2Basis();          % prj all stc onto that basis
         
         self.SSraw = [ones(size(self.SSraw(:,1)))'; (self.SSraw*self.basis1D')'; self.basisPrj; (self.RRraw*self.basis1D')']';
      end
      
      
      
      function coef2kernel(self, varargin)
         % either works on the internal regression coefficients or on those
         % provided by the argument
         if nargin>1
            h = varargin{1};
         else
            h = self.h;
         end
         self.h0 = h(1);
         self.h1 = h(2:size(self.basis1D,1)+1)'*self.basis1D;
         %%
         self.h2 = zeros(size(self.basis2D,2));
         hCoef = h(size(self.basis1D,1)+2:size(self.basis1D,1)+size(self.basis2D,1)+1);
         for b = 1:length(hCoef)
            self.h2 = self.h2 + squeeze(self.basis2D(b,:,:))*hCoef(b);
         end
         self.hHist = h(size(self.basis1D,1)+size(self.basis2D,1)+2:end)'*self.basis1D;

      end
      
      function kernel2coef(self,varargin)
         % either works on the internal regression coefficients or on those
         % provided by the argument
         if nargin>1
            h2 = varargin{1};
         else
            h2 = self.h2;
         end
         
         hCoef = zeros(size(self.basis2D,1),1);
         for b = 1:size(self.basis2D,1)
            thisBas = squeeze(self.basis2D(b,:,:));
            hCoef(b) = sum(sum(thisBas.*h2));
         end
         self.h(size(self.basis1D,1)+2:end) = hCoef;
      end
      
      
   end
   
end
