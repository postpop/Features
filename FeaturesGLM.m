classdef FeaturesGLM < Features
   
   properties
      stim, resp, spikeIdx
      bias, stimFilter, respFilter
      nsp, delta, RResp
      feat0
      kbas, ihbas
      kbasLen, kbasNum, ihbasLen, ihbasNum
      SSrawBas, RespBas
      
   end
   
   
   methods (Access='public')
      
      function self = FeaturesGLM(sr, n)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.sr.filterTrains(2);
            
            % prepare basis for stimulus filter
            kbasprs.nh = 8;  % Number of basis vectors to use
            kbasprs.hdt = 1; % Spacing for time sampling  (set to integer)
            kbasprs.endpoints = [2 38];  % location of first and last bump peaks.
            kbasprs.b = 8;  % nonhlinearity of bump spacing: small values-> strongly nonlinear; large values -> linear spacing
            [kt,kbas] = makeRaisedCosBasis(kbasprs); % generate basis
            self.kbas = flipud(kbas);
            nkt = length(kt);
            
            % project onto basis
            [self.kbasLen, self.kbasNum] = size(self.kbas);
            self.n = self.kbasLen;
            self.getSTE();
            self.SSrawBas = double(self.SSraw*self.kbas);%project stim onto basis
            
            % prepare basis for response filter
            ihbasprs.ncols = 6;  % Number of basis vectors for post-spike kernel
            ihbasprs.hpeaks = [2 16];  % Peak location for first and last vectors
            ihbasprs.b = .5;  % How nonlinear to make spacings
            ihbasprs.absref = 1; % absolute refractory period
            [iht,self.ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,1);
            
            % project onto basis
            [self.ihbasLen, self.ihbasNum] = size(self.ihbas);
            self.RespBas = double(makeStimRows(self.Resp, self.ihbasLen)) * self.ihbas ;
            
            % prepare everything else
            self.nsp = sum(self.Resp>0);
            self.spikeIdx = find(self.Resp>0);
            self.spikeIdx = self.spikeIdx(10:end-10);
            self.delta = self.sr.binSize/1000;
            
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      %2a. STA/STC
      function [feat, V] = getFeat(self)
         
         
         % initial conditions
         bias = 1;
         stimFilter =  randn(self.kbasNum,1);%double(mean(self.SSrawBas(self.spikeIdx,:)))';%randn(size(self.SSrawBas,2),1);
         stimFilter = stimFilter./norm(stimFilter);
         respFilter = randn(self.ihbasNum,1)/100;respFilter = 0*respFilter./norm(respFilter);
         param0 = [bias, stimFilter', respFilter'];
         % remember initial conditions
         self.feat0(1,1) = param0(1);
         self.feat0(1:self.kbasLen,2) = param0(2:size(self.SSrawBas,2)+1)*self.kbas';
         self.feat0(1:self.ihbasLen,3) = param0(size(self.SSrawBas,2)+1 + (1:self.ihbasNum))*self.ihbas';
                  
         % optimize!!!
         options = optimset('GradObj','off','Display','iter','TolFun',1e-12,'TolX',1e-12);
         [params,fval,exitflag] = fminunc(@loglike,param0,options);
         %[params,fval,exitflag,output] = minFunc(@loglike,param0',options)
         
         self.feat(1,1) = params(1);
         self.feat(1:self.kbasLen,2) = params(2:size(self.SSrawBas,2)+1)*self.kbas';
         self.feat(1:self.ihbasLen,3) = params(size(self.SSrawBas,2)+1 + (1:self.ihbasNum))*self.ihbas';
         self.V = fval;
         
         function [ll, grad] = loglike(filters)
            etol = 1e-100;
            bias = filters(1);
            stimFilter = filters(2:size(self.SSrawBas,2)+1)';
            %stimFilter = stimFilter./norm(stimFilter);
            respFilter = filters(size(self.SSrawBas,2)+1 + (1:self.ihbasNum))';
            % project stim on stimfilter
            stimFilt = self.SSrawBas*stimFilter;%  conv(self.stim(:),stimFilter,'same');%
            % project resp on respfilter
            respFilt = self.RespBas*respFilter; %conv(self.resp(:),respFilter,'same');%
            % predict rate function
            lambda = exp(stimFilt + respFilt + bias);
            lambda(lambda<etol) = etol;
            lambda(lambda>1000) = 1000;
            
            % model log likelihood
            ll = -(sum(log(lambda(self.spikeIdx))) - self.delta * sum(lambda));
            if isinf(ll)
               disp('inf')
            end
            if isnan(ll)
               disp('nan')
            end
            if nargout==2
               lamb = makeStimRows(lambda, self.n);
               %compute gradients
               gradStim = sum(self.SSrawBas(self.spikeIdx,:))' - ...
                  self.delta * self.SSrawBas'*lambda;
               gradResp = sum(self.RespBas(self.spikeIdx,:))' - ...
                  self.delta * self.RespBas'*lambda;
               gradBias = self.nsp - ...
                  self.delta * sum(lambda);
               grad = -[gradBias; gradStim; gradResp];%[gradBias; gradStim];%
            end
         end
      end
      
      
   end
   
end
