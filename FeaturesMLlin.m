classdef FeaturesMLlin < FeaturesML
   
   
   methods (Access='public')
      
      function self = FeaturesMLlin(sr, n)
         self = self@FeaturesML(sr,n);
         self.getSTEML()
                  
      end
      
      function getSTEML(self)
         self.getSTE();             % get standard STE
         
         self.SSraw = [ones(size(self.SSraw(:,1)))'; self.SSraw']';
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
         self.h1 = h(2:self.n+1)';
         %%
         self.h2 = [];
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
