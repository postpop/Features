classdef SpikeRespCont < SpikeResp
   
   methods (Access='public')
      function self = SpikeRespCont(spikeTimes, binSize)
         % SpikeResp(spikeTimes, binSize)
         if nargin==1
            binSize = 0.2;%ms
         end
         self = self@SpikeResp(spikeTimes, binSize);
         
         %          self.binSize = binSize;
         %          self.spikeTimes = spikeTimes;
         %          self.cats = size(self.spikeTimes,1);
         %          self.trials = size(self.spikeTimes,2);
         %          self.catNames = num2cellstr(1:self.cats);
         %          self.getRawTrains();
         
      end
      
      function getRawTrains(self)
         self.trains = self.spikeTimes;
         tmp = mean(self.spikeTimes,2);
         self.PSTH = tmp(:,:);
         %[self.PSTH, trainStdDummy, self.trains] = times2trains(self.spikeTimes, self.binSize);
      end
      
   end
   
end
