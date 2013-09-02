classdef FeaturesiSTAC < Features

   properties

   end



   methods (Access='public')

      function self = FeaturesiSTAC(sr,n)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.n = n;
            self.getSTE();
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end

      %2a. STA/STC
      function getFeat(self,varargin)
         % getFeat([nFeat=3])
         %PARAMS:
         %  nFeat  - number of features to search for
         nFeat = 3;
         if nargin>1
            nFeat = varargin{1};
         end
         sfSTC = FeaturesSTC(self.sr,self.n);
         sfSTC.getFeat();

         [self.feat, V] = compiSTAC(double(sfSTC.STA'), double(sfSTC.STC'), double(sfSTC.RTA'), double(sfSTC.RTC'), nFeat);
         self.V = [V(1); diff(V)];
         clear sfSTC;
      end

   end

end
