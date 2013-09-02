classdef FeaturesMID < Features
   
   properties
      pathToExe
      prefix
      dims
   end
   
   
   methods (Access='public')
      
      function self = FeaturesMID(sr,n)
         % estimate filters using the "Maximally Informative Dimensions"
         % (MID) by Tatyana Sharpee. Calls a LINUX binary downloaded from http://cnl-t.salk.edu/Download/
         %USAGE
         %  self = FeaturesMID(sr,n,dirName)
         %PARAMS
         %  sr      - SpikeResp object
         %  n       - length of the filter
         if isa(sr,'SpikeResp')
            self.sr = sr;
            self.n = n;
            self.dims = 1;
            self.getSTE();
            self.prefix = 'mid';
            self.pathToExe = '/Users/janclemens/Dropbox/code/MID/';
         else
            disp('ERROR: arg #1 is not of class SpikeResponse or other error!');
         end
         
      end
      
      %MID
      
      function getFeat(self)
         %% prepare
         respFile = [self.prefix 'resp.dat'];
         stimFile = [self.prefix 'stim.dat'];
         
         % write STIM and SPIKE files
         respModel = uint8(reshape(self.Resp,1,[])');%flatten response
         fid = fopen(respFile,'w+t');
         fprintf(fid,'%u\n',respModel);
         fclose(fid);
         stimModel = double(reshape(self.SSraw',1,[])');
         fid = fopen(stimFile,'wb');
         fwrite(fid,stimModel,'double');
         fclose(fid);
         %% adjust params
         param = textread([self.pathToExe 'params.xml'],'%s', 'delimiter','+');
         param{4} = [param{4}(1:end-11) respFile '"/>'];
         param{19} = [param{19}(1:end-11) stimFile '"/>'];
         param{22} = [param{22}(1:end-5) int2str(size(self.SSraw,2)) '"/>'];
         param{25} = [param{25}(1:end-6) int2str(size(self.SSraw,2)) '"/>'];
         param{end-2} = [param{end-2}(1:end-7) self.prefix '"/>'];
         fid = fopen([self.prefix 'param.xml'], 'wt');
         fprintf(fid, '%s\n', param{:});
         fclose(fid);
         
         %% call MID routine
         if self.dims==1 % 1D
            eval(['!' self.pathToExe 'mid1d 1 ' self.prefix 'param.xml 1'])
         else % ND
            % get initial cond's
            param{end-5} = [param{end-5}(1:end-7) '100"/>'];
            fid = fopen([self.prefix 'paramQuick.xml'], 'wt');
            fprintf(fid, '%s\n', param{:});
            fclose(fid);
            for dim = 1:self.dims
               eval(['!' self.pathToExe 'mid1d ' int2str(dim) ' ' self.prefix 'paramQuick.xml 1'])
            end
            % do full ND optimization
            eval(['!' self.pathToExe 'midnd ' int2str(self.dims) ' ' self.prefix 'param.xml 1'])
         end
         self.readFeat()
      end
      
      function readFeat(self)
         % read results
         if self.dims==1 %1D
            fid=fopen([self.prefix '-1D-n1-v1-p1.bst'], 'rb');
            self.feat(:,1) = fread(fid, 'double');
            fclose(fid);
         else  %ND
            for dim = 1:self.dims
               fid=fopen([self.prefix '-ND-n' int2str(dim) '-v1-p1.bst'], 'rb');
               self.feat(:,dim) = fread(fid, 'double');
               fclose(fid);
            end
         end
      end
      
   end
end
