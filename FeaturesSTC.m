classdef FeaturesSTC < Features
   
   properties
      STA, RTA
      STC, RTC
      orthoSTA
   end
   
   
   methods (Access='public')
      
      function self = FeaturesSTC(sr,n)
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.n = n;
            self.orthoSTA = 1;
            self.getSTE();
         else
            disp('ERROR: arg #1 is not of class SpikeResponse');
         end
      end
      
      %2a. STA/STC
      function [feat, V] = getFeat(self, varargin)
         N = size(self.SSraw,1);
         self.RTA = sum(self.SSraw,1)/N;% prior/raw STA
         self.RTC = (self.SSraw'*self.SSraw)/N - self.RTA'*self.RTA;
         
         SSrawCenter = bsxfun(@minus,self.SSraw,self.RTA);%center dist
         self.STA = self.Resp'*self.SSraw;%STA
         
         if nargin<2 
            normSTA = self.STA/norm(self.STA);%normalize STA
            %SSrawCenter = bsxfun(@minus,SSrawCenter,self.STA);%sub sta
            chunkSize = 10000;
            if size(SSrawCenter,1)>chunkSize% chunk
               chunks = [1:chunkSize:size(SSrawCenter,1)-1 size(SSrawCenter,1)];
               for chk = 1:length(chunks)-1
                  chunkIdx = chunks(chk):chunks(chk+1);
                  if self.orthoSTA
                     SSrawCenter(chunkIdx,:) = SSrawCenter(chunkIdx,:)- SSrawCenter(chunkIdx,:)*normSTA'*normSTA;%orthsta
                  end
                  chunkSTC(chk,:,:) = SSrawCenter(chunkIdx,:)'*bsxfun(@times,SSrawCenter(chunkIdx,:),self.Resp(chunkIdx,:));%STC
               end
               self.STC = squeeze(sum(chunkSTC,1));
            else
               if self.orthoSTA
                  SSrawCenter = SSrawCenter-SSrawCenter*normSTA'*normSTA;%orthsta
               end
               self.STC = SSrawCenter'*bsxfun(@times,SSrawCenter,self.Resp);%STC
            end
            [feat,V] = eig(self.STC - self.RTC);% get eigs of the difference matrix
         else
            feat = ones(length(self.STA));
            feat(:,1) = self.STA./norm(self.STA);
            V = diag(ones(1,length(self.STA)));
         end
         self.feat = [self.STA'./norm(self.STA'),feat];
         self.V = [0;diag(V)];
         
%          for f = 2:size(self.feat,2)
%             [val, idx] = max(abs(self.feat(:,f)));
%             self.feat(:,f) = self.feat(:,f) * sign(self.feat(idx,f));
%          end
      end
      
      function [feat, V] = dejitterSTE(self)
         %% estimate initial model
         self.getFeat();
         %% initial jitter distribution
         jitterSig = 4;
         jitterMax = 4*sig;
         jitterDist = normpdf(-jitterMax:.1:jitterMax,0,jitterSig);
         for run = 1:10
            % run model
            m2 = ModelLNP1(self.sr,self,2,1,16);
            m2.runModel;
            % determine jitter dist
            
            
            % estimate new model incorporating jitter
            
            % evaluate model
            
         end
         
         
      end
      
      
   end
   
end
