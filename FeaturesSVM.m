clearvars
load('C:\Eigene Dateien\data.MAT\R_20081022_sgl\TN1\RepLng.mat')
%load('C:\Eigene Dateien\data.MAT\R_20081215_sgl\BSN1\RepLng.mat')
%load('C:\Eigene Dateien\data.MAT\R_20081209_sgl\AN1\RepLng.mat')
%% Simple Cell Model
scFilt = diff(gausswin(20,6));
plot(scFilt)
stim = generateCarrier(5, 1000, 1, 100);
drive = conv(stim,scFilt,'same');
respQud = drive.^2;
respLin = drive.^2;
respLin(drive<0)=0;

respLin = respLin/max(respLin);
respLin(respLin>=.1) = 1;
respLin(respLin< .1) = 0;
respQud = respQud/max(respQud);
respQud(respQud>=.1) = 1;
respQud(respQud< .1) = 0;

clf
plot(respQud(1:100),'b')
hold on
plot(respLin(1:100),'--r')
hold off
stimRows = double(makeStimRows(stim,64));
cAll = 2.^(-10:4:10);
gAll = 2.^(-10:4:10);
rAll = 0;

[accuracyL, predlabelL,HL,fitScapeL,decisionValuesL,modelsL] = svmGridLinear(stimRows, respLin, cAll);
%[accuracyP, predlabelP,HP,fitScapeP,decisionValuesP,modelsP] = svmGridPoly(stimRows, respLin, cAll,gAll);



%% basis
myPyr = get1DLaplacianPyramidBasis(64,4,.5,2.5);

%%
sfLng = FeaturesSTC(srLng,64);
stim = sfLng.SSraw;
resp = sfLng.Resp;
resp(resp>0)=1;
resp(resp<=0)=0;

save('C:\Eigene Dateien\code.ext\volterra\volt_stim.dat','stim','-ASCII')
save('C:\Eigene Dateien\code.ext\volterra\volt_resp.dat','resp','-ASCII')
%%
sfRep = FeaturesSTC(srRep,64);
stim = sfRep.SSraw;
resp = sfRep.Resp;
resp(resp>0)=1;
resp(resp<=0)=0;

save('C:\Eigene Dateien\code.ext\volterra\volt_stimRep.dat','stim','-ASCII')
save('C:\Eigene Dateien\code.ext\volterra\volt_respRep.dat','resp','-ASCII')


%% prepare data
sfLng = FeaturesSTC(srLng,64);
sfLng.getFeat();
sfRep = FeaturesSTC(srRep,64);
stimFit = sfLng.SSraw;
respFit = double(makeStimRows(sfLng.Resp,64));
stimVal = sfRep.SSraw;
respVal = double(makeStimRows(sfRep.Resp,64));
%%
N = 10000;
stimFit = stimFit(1:N,:);
respFit = respFit(1:N,:);
stimVal = stimVal;
respVal = respVal;

resp = sfLng.Resp(1:N);
resp(resp>0)=1;
resp(resp<=0)=-1;

X = [stimFit*myPyr];
%U = [respFit*myPyr];
%%
cAll = 2.^(-10:4:10);
gAll = 2.^(-10:4:10);
rAll = 0;
[accuracyL, predlabelL,HL,fitScapeL,decisionValuesL] = svmGridLinear(X, resp, cAll);
save('SVMlinear','accuracyL', 'predlabelL','HL','fitScapeL','decisionValuesL')
[accuracyP, predlabelP,HP,fitScapeP,decisionValuesP] = svmGridPoly(X, resp, cAll,gAll);
save('SVMpoly','accuracyP', 'predlabelP','HP','fitScapeP','decisionValuesP')
%[accuracyS, predlabelS,HS,fitScapeS,decisionValuesS] = svmGrid(X, resp, cAll,gAll);
