classdef FeaturesMLlinBasis < FeaturesML
   
   
   methods (Access='public')
      
      function self = FeaturesMLlinBasis(sr,n)
         self = self@FeaturesML(sr,n);
         self.getSTEML()
      end
      
      function getSTEML(self)
         self.getSTE();             % get standard STE
         self.getBasis1D();
         
         self.SSraw = [ones(size(self.SSraw(:,1)))'; (self.SSraw*self.basis1D')' ]';
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
         self.h1 = h(1 + (1:size(self.basis1D,1)))'*self.basis1D;
      end
      
      function kernel2coef(self,varargin)
         % either works on the internal regression coefficients or on those
         % provided by the argument
         if nargin>1
            h1 = varargin{1};
         else
            h1 = self.h1;
         end
         
         self.h(1 + (1:size(self.basis1D,1)),1) = h1*self.basis1D';
      end
      
      
   end
   
end
