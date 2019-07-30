%%  1: Load MAT-File with two voxel timeseries

load sample.mat

NV = length(timeMR);
rLinear = [1:NV]';

RETR_RespRegr=func_RETR_Resp_regressors(resp,3,Fs); NR_RETRresp=size(RETR_RespRegr,2); RETR_RespRegr = RETR_RespRegr(ind_BOLD,:);
RETR_CardRegr=func_RETR_Card_regressors(time,PPGlocs,3);   NR_card=size(RETR_CardRegr,2);  RETR_CardRegr = RETR_CardRegr(ind_BOLD,:);

regr=[RETR_RespRegr,RETR_CardRegr,ones(NV,1),rLinear];

ind_resp = 1:NR_RETRresp;
ind_card=max(ind_resp+1):max(ind_resp+NR_card);

%% 2: Fit RETROICOR regressors on the first voxel 
% y1 is a timeseries w/ a strong cardiac artifact  -------------------------------------

voxel = y1;
B = regr\voxel;
yPred_full = regr*B;
yPred_card = regr(:,ind_resp)*B(ind_resp);
yPred_resp = regr(:,ind_card)*B(ind_card);

r_RETRresp = corr(voxel,yPred_card)
r_RETRcard = corr(voxel,yPred_resp)

figure('Position', [19         890        2351         420])
plot(timeMR, voxel), hold on
plot(timeMR, yPred_full)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Voxel time-series affected by cardiac pulsatility')
xlim([500 800])

%% 3: Fit RETROICOR regressors on the second voxel 
% y2 is a timeseries w/ a strong respiratory artifact  -------------------------------------

voxel = y2;
B = regr\voxel;
yPred_full = regr*B;
yPred_card = regr(:,ind_resp)*B(ind_resp);
yPred_resp = regr(:,ind_card)*B(ind_card);

r_RETRresp = corr(voxel,yPred_card)
r_RETRcard = corr(voxel,yPred_resp)

figure('Position', [19         890        2351         420])
plot(timeMR, voxel), hold on
plot(timeMR, yPred_full)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Voxel time-series affected by breathing-related motion')
xlim([500 800])

