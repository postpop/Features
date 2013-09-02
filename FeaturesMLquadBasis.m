classdef FeaturesMLquadBasis < FeaturesML
   
   
   methods (Access='public')
      
      function self = FeaturesMLquadBasis(sr, n)
         self = self@FeaturesML(sr,n);
         self.getSTEML()
      end
      
      function getSTEML(self, varargin)
         self.getSTE();             % get standard STE
         self.getBasis2D();         % get basis
         self.getBasis1D();
         self.prj2Basis();          % prj all stc onto that basis
         
         self.SSraw = [ones(size(self.SSraw(:,1)))'; (self.SSraw*self.basis1D')'; self.basisPrj]';
%          RTA = mean(self.SSraw,1);% prior/raw STA
%          self.SSraw = bsxfun(@minus,self.SSraw,RTA);%center dist
%          
%          if nargin>1
%          %proj. out filter and find second filter orth. to first one
%             flt1 = normalize(varargin{1});
%             %self.SSraw = self.SSraw-self.SSraw*flt1'*flt1;
%             self.Resp = self.Resp - varargin{2};
%          end
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
         hCoef = h(size(self.basis1D,1)+2:end);
         for b = 1:length(hCoef)
            self.h2 = self.h2 + squeeze(self.basis2D(b,:,:))*hCoef(b);
         end
         self.h2 = self.h2 + tril(self.h2,-1)';
         
      end
      
      function kernel2coef(self,varargin)
         % either works on the internal kernel or on that
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
