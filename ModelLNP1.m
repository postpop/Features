classdef ModelLNP1 < Model
   
   properties
      
      nonLinearity, NLbins
      linFilter
      drive
   end
   
   methods (Access='public')
      function self = ModelLNP1(srTest,sf,whichFeat,varargin)         
         %self = ModelLNP1(srTest,sf,whichFeat,[MODE, nBins,TAU])         
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
         self.sf.getNonlinearity1D(self.MODE,whichFeat,self.nBins);
         self.nonLinearity = double(mean(self.sf.sr.trains(:))*squeeze(self.sf.psr1D(whichFeat,1,:)));%double(mean(srTest.trains(:))*squeeze(self.sf.psr1D(whichFeat,1,:)));
         self.nonLinearity(isnan(self.nonLinearity))=0;
         self.NLbins = squeeze(self.sf.psr1D(whichFeat,2,:));
         self.linFilter = self.sf.feat(:,whichFeat)./norm(self.sf.feat(:,whichFeat));
         
         self.srTest = srTest;    
         self.srTest.getRawTrains();       
         if self.TAU>0
            self.srTest.filterTrains(self.TAU);
         end
         self.trains = reshape(self.srTest.trains,1,self.srTest.trials,[]);
         self.resp = reshape(self.srTest.PSTH,1,[]);
         self.srTest.getRawTrains();       
         %self.sfTest = FeaturesSTC(srTest,self.sf.n);        
         self.sfTest = eval([class(self.sf) '( srTest, self.sf.n )']);
         try
            self.sfTest.getSTEML();
         end
      end
      
      function runModel(self)
         self.drive = self.sfTest.SSraw*(self.linFilter);      
         [~, uniBinIdx] = unique(self.NLbins);
         self.pred = interp1(self.NLbins(uniBinIdx),self.nonLinearity(uniBinIdx),self.drive,'pchip')';                  
         [cx,lag] = xcorr(self.resp,self.pred,10);
         [~,idx] = max(cx);
         self.resp = self.resp(11+lag(idx):end-10);
         self.pred = self.pred(11:end-10-lag(idx));
         self.drive= self.drive(11:end-10-lag(idx));
         self.pred(self.pred<0)=0;
         self.trains = self.trains(:,:,11+lag(idx):end-10);
         
         maxLen = min(size(self.resp,2),size(self.pred,2));
         self.pred = self.pred(1:maxLen); 
         self.resp = self.resp(1:maxLen); 
      end
      
      
      
   end
end
