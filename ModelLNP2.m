classdef ModelLNP2 < Model
   
   properties
      nonLinearity, NLbins
      linFilter1,linFilter2
      drive1,drive2
   end
   
   methods (Access='public')
      function self = ModelLNP2(srTest,sf,whichFeat,varargin)
         %self = ModelLNP2(srTest,sf,whichFeat,[MODE, nBins,TAU])
         self.MODE = 2;
         self.nBins = 16;
         self.TAU = 0;
         if ~isempty(varargin)
            self.MODE = varargin{1};
            self.nBins = varargin{2};
            if nargin>5
               self.TAU = varargin{3};
            end
         end
         self.whichFeat = whichFeat;
         self.sf = sf;
         self.sf.getNonlinearity2D(self.MODE,whichFeat(1:2),self.nBins,self.TAU);
         self.nonLinearity = double(mean(self.sf.sr.trains(:))*squeeze(self.sf.psr2D(1,:,:)));
         self.nonLinearity(isnan(self.nonLinearity))=0;
         self.nonLinearity(self.nonLinearity(:)<0)=0;
         
         self.NLbins = squeeze(self.sf.psr2D(2,2,:));
         
         self.linFilter1 = self.sf.feat(:,whichFeat(1))./norm(self.sf.feat(:,whichFeat(1)));
         self.linFilter2 = self.sf.feat(:,whichFeat(2))./norm(self.sf.feat(:,whichFeat(2)));
         
         self.srTest = srTest;
         self.srTest.getRawTrains();
         if self.TAU>0
            self.srTest.filterTrains(self.TAU);
         end
         self.trains = reshape(self.srTest.trains,1,self.srTest.trials,[]);
         self.resp = reshape(self.srTest.PSTH,1,[]);
         self.srTest.getRawTrains();
         self.sfTest = FeaturesSTC(srTest,self.sf.n);
      end
      
      function runModel(self)
         
         self.drive1 = self.sfTest.SSraw*(self.linFilter1);
         %self.drive1(self.drive1>nanmax(self.NLbins)) = nanmax(self.NLbins);
         %self.drive1(self.drive1>nanmin(self.NLbins)) = nanmin(self.NLbins);
         self.drive2 = self.sfTest.SSraw*(self.linFilter2);
         %self.drive2(self.drive2>nanmax(self.NLbins)) = nanmax(self.NLbins);
         %self.drive2(self.drive2>nanmin(self.NLbins)) = nanmin(self.NLbins);
         
         self.pred = interp2(self.NLbins,self.NLbins,...
            self.nonLinearity,self.drive1,self.drive2,'spline')';
         self.pred(self.pred>max(self.nonLinearity(:))) = max(self.nonLinearity(:));
         [cx,lag] = xcorr(self.resp,self.pred,10);
         [val,idx] = max(cx);
         self.resp = self.resp(11+lag(idx):end-10);
         self.pred = self.pred(11:end-10-lag(idx));
         self.pred(self.pred<0)=0;
         self.trains = self.trains(:,:,11+lag(idx):end-10);
         
         maxLen = min(size(self.resp,2),size(self.pred,2));
         self.pred = self.pred(1:maxLen);
         self.resp = self.resp(1:maxLen);
                  
      end
      
      
   end
end
