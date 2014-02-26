classdef FeaturesML < Features
   
   properties
      h, h0, h1, h2, hHist
      regH, regFitInfo
      RRraw, RRrawDelay
      pred, perf
      predTest, perfTest
      basis1D, basis2D, basisPrj
      khat
   end
   
   methods (Abstract)
      getSTEML(self)
      coef2kernel(self, varargin)
      kernel2coef(self, varargin)
   end
   
   methods (Access='public')
      
      function self = FeaturesML(sr, n)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.n = n;
            self.getSTE();
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      function getFeat(self, varargin)
         % getFeat([mode=0, trainIdx=all idx])
         regMode = 0;
         if nargin>1
            regMode = varargin{1};
         end
         trainIdx = true(size(self.Resp));
         if nargin>2
            trainIdx = varargin{2};
         end
         %% regression
         switch regMode
            case {1, 2} % ridge or lasso
%                opts = statset('UseParallel','always');
%                % parameter ALPHA interpolates between ridge (0) and lasso (1)
%                if regMode == 1, alpha = 0.001; else alpha = 1; end
%                [self.regH, self.regFitInfo] = lasso(self.SSraw(trainIdx,:), self.Resp(trainIdx), 'Options', opts, 'alpha', alpha, 'NumLambda', 64, 'CV', 4);
%                lamOpt = self.regFitInfo.Index1SE;
%                self.h = self.regH(:,lamOpt);
               self.khat = runRidgeOnly(self.SSraw(trainIdx,:), self.Resp(trainIdx), size(self.SSraw,2), 1);
               self.h = self.khat;
            case 3 % ALDsf
               self.khat = runALD(self.SSraw(trainIdx,:), self.Resp(trainIdx), size(self.SSraw,2), 1);
               self.h = self.khat.khatSF;
            otherwise % standard least squares
               self.h = pinv(self.SSraw(trainIdx,:))*self.Resp(trainIdx);
         end
         
         % rearrange kernels to get terms corr. to diff. orders
         self.coef2kernel();
         % prediction
         self.pred = self.SSraw*self.h;
         self.perf = rsq(self.pred(trainIdx), self.Resp(trainIdx));
         self.perfTest = rsq(self.pred(~trainIdx), self.Resp(~trainIdx));
         self.feat = self.h;
      end
      
      function getBasis1D(self)
         xs = size(self.SSraw,2);
         sigma = 2;
         X = 1:xs;
         nBasis = self.n/2;%self.n/2;
         
         self.basis1D = zeros(nBasis, size(self.SSraw,2));
         cnt = 0;
         xpts = linspace(0,xs,nBasis);
         for xidx = 1:length(xpts)
            x = xpts(xidx);
            mpsF = x;
            cpsF = sigma;
            nps = normpdf(X,mpsF,cpsF);
            cnt = cnt+1;
            self.basis1D(cnt,:) = nps;
         end
      end
      
      function getBasis2D(self)
         n = size(self.SSraw,2);
         xs = n;
         ys = n;
         sigma = 3;
         [X1, X2] = meshgrid(1:xs,1:ys);
         nBasis = self.n/2;
         
         self.basis2D = zeros((nBasis^2 + nBasis)/2, n, n);
         cnt = 0;
         xpts = linspace(0,xs,nBasis);
         ypts = linspace(0,ys,nBasis);
         for xidx = 1:length(xpts)
            for yidx = xidx:length(ypts)
               x = xpts(xidx);
               y = ypts(yidx);
               mpsF = [x,y];
               cpsF = diag([sigma sigma]);
               nps = mvnpdf([X1(:) X2(:)],mpsF,cpsF);
               nps = reshape(nps,xs,ys)';
               %nps = nps + nps';
               cnt = cnt+1;
               self.basis2D(cnt,:,:) = nps;
            end
         end
      end
      
      function prj2Basis(self)
         basis = self.basis2D;
         n = size(self.SSraw,2);
         SSraw = self.SSraw;
         basisPrj = zeros(size(basis,1),size(SSraw,1));
         %%
         for b = 1:size(basis,1)
            bas = reshape(basis(b,:,:), n, n);
            if isempty(strfind(lower(class(self)),'basis'))
               basisPrj(b,:) = sum(SSraw'.*(bas*SSraw'));
            else
               % do this only for indices where the basis is nonzero
               % FIX: turned of for delta basis
               xidx = sum(abs(bas.^2))>eps;
               yidx = sum(abs(bas.^2),2)>eps;
               basisPrj(b,:) = sum(SSraw(:,yidx)'.*(bas(yidx,xidx)*SSraw(:,xidx)'));
            end
         end
         %%
         self.basisPrj = basisPrj;
      end
      
      function getDeltaBasis2D(self)
         n = self.n;
         nBasis = (n.^2+n)/2;
         self.basis2D = zeros(nBasis, self.n, self.n);
         cnt = 0;
         template = zeros(n);
         for xidx = 1:n
            for yidx = xidx:n
               nps = template;
               nps(xidx,yidx) = 1;
               %nps = nps + nps':
               cnt = cnt+1;
               self.basis2D(cnt,:,:) = nps;
            end
         end
         
      end
      
   end
   
end
