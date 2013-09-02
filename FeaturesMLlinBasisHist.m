classdef FeaturesMLlinBasisHist < FeaturesML
   
   
   methods (Access='public')
      
      function self = FeaturesMLlinBasisHist(sr,n)
         self = self@FeaturesML(sr,n);
         self.getSTEML()
      end
      
      function getSTEML(self)
         delay = 0;
         dt = 1;
         self.getSTE();             % get standard STE
         self.RRraw = makeStimRows(self.Resp, dt*self.n+delay);
         %self.RRraw = self.RRraw(:,1:end-delay);
         self.RRraw = resample(self.RRraw',1,dt,100)';
         self.getBasis1D();
         tmp = (self.RRraw*self.basis1D')';
         tmp(end,:) = 0;
         self.SSraw = [ones(size(self.SSraw(:,1)))'; (self.SSraw*self.basis1D')' ; tmp]';
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
         self.hHist = h(size(self.basis1D,1)+2:end)'*self.basis1D;
      end
      
      function kernel2coef(self,varargin)
         % either works on the internal regression coefficients or on those
         % provided by the argument
         if nargin>1
            h1 = varargin{1};
         else
            h1 = self.h1;
         end
         
         self.h(2:end) = h1;
      end
      
      
      
   end
   
end
