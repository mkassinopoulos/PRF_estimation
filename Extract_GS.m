%% Extract the GS from fMRI data
% This script extracts 4 variants of GS. Specifically, it extracts the 
% mean timeseries from the whole brain (WB), gray matter (GM),
% white matter (WM) and the cerebrospinal fluid (CSF).
% In addition, this script removes the motion realignment parameters
% as well as the cardiac and respiratory  RETROICOR regressors from 
% each of the 4 variants of GS. In total, 8 variants of GS are extracted
% which are stored in the MAT-File 'GS.mat'.

%% 1: Define directories for necessary data

clear, clc, close all
task_list = {'Rest1_LR','Rest1_RL', 'Rest2_LR', 'Rest2_RL'};

%  Set the following parameters !!

subject_id = '118730';
task_id = task_list{4};   

% -----------------------------------------

path_scan = ['../Data/S',subject_id,'/',task_id,'/'];
path_func = [path_scan,'/fMRI/func.nii.gz'];
path_mask = [path_scan,'/fMRI/brainmask_fs.2.nii.gz'];
path_CSF = [path_scan,'/fMRI/Mask_CSF.nii.gz'];
path_GM = [path_scan,'/fMRI/Mask_GM.nii.gz'];
path_WM = [path_scan,'/fMRI/Mask_WM.nii.gz'];
path_movement = [path_scan,'/fMRI/Movement_Regressors_dt.txt'];
path_output = [path_scan,'/fMRI/GS.mat'];

%%  2: Load niftii files

nii_func = load_nii(path_func); img = nii_func.img;
nii_mask = load_nii(path_mask);   img_mask = nii_mask.img;
nii_CSF = load_nii(path_CSF);   img_CSF = nii_CSF.img;
nii_GM = load_nii(path_GM);   img_GM = nii_GM.img;
nii_WM = load_nii(path_WM);  img_WM = nii_WM.img;

[NX, NY, NZ, NV] = size(img);

%% 3: Copy all voxel timeseries from the whole-brain in a 2-d matrix

NVX = 0;
CordAll=zeros(NX*NY*NZ,3);
for x=1:NX
    for y=1:NY
        for z=1:NZ
            if img_mask(x,y,z)==1
                voxel = squeeze(img(x,y,z,:));
                voxel_mean = mean(voxel); 
                if   voxel_mean > 100
                    NVX=NVX+1;
                    CordAll(NVX,1)=x;
                    CordAll(NVX,2)=y;
                    CordAll(NVX,3)=z;
                    img(x,y,z,:) = detrend(voxel,'linear') + voxel_mean; 
                end                
            end
        end
    end
end
CordAll=CordAll(1:NVX,:);

img_col=zeros(NV,NVX);
for vox=1:NVX
    x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3); 
    img_col(:,vox) = squeeze(img(x,y,z,:));
end

%% 4: Extract the mean timeseries from WB, GM, WM and CSF


idx_GM=zeros(NVX,1); N_GM=0;
idx_WM=zeros(NVX,1); N_WM=0;
idx_CSF=zeros(NVX,1); N_CSF=0;
for vox=1:NVX
    x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
    if img_GM(x,y,z)>0.5
        N_GM=N_GM+1;
        idx_GM(N_GM)=vox;
    elseif img_CSF(x,y,z)>0.8
        N_CSF=N_CSF+1;
        idx_CSF(N_CSF)=vox;
    elseif img_WM(x,y,z)>0.9
        N_WM=N_WM+1;
        idx_WM(N_WM)=vox;
    end
end

idx_GM=idx_GM(1:N_GM,1);
idx_WM=idx_WM(1:N_WM,1);
idx_CSF=idx_CSF(1:N_CSF,1);
fprintf('   GM:  %3.1f%%  \n',100*N_GM/NVX);
fprintf('   WM:  %3.1f%%  \n',100*N_WM/NVX);
fprintf('   CSF: %3.1f%%  \n',100*N_CSF/NVX);

CSF = mean(img_col(:,idx_CSF)')';
GM = mean(img_col(:,idx_GM)')';
WM = mean(img_col(:,idx_WM)')';
WB = mean(img_col')';

fprintf(' ----- Correlation ----- \n')
fprintf('CSF-GM: %3.1f%%  \n',100*corr(CSF,GM))
fprintf('WM-GM: %3.1f%%  \n',100*corr(WM,GM))
fprintf('CSF-WM: %3.1f%%  \n',100*corr(CSF,WM))
fprintf('WB-GM: %3.1f%%  \n',100*corr(WB,GM))

TR = 0.72;   timeMR = 0 :TR: (NV-1)*TR;

plot(timeMR,[(GM),(WM),(CSF)])
xlabel('Time (s)')
legend('Gray matter','White matter','Cerebrospinal fluid')

%% 5:  Create regressor matrix with motion parameters and RETROICOR regressors

movRegr = load(path_movement);
load(path_physio);

ind_BOLD = find(trig==1);
addpath('RETROICOR')
RETR_RespRegr = func_RETR_Resp_regressors(resp,3,Fs); NR_RETRresp=size(RETR_RespRegr,2); RETR_RespRegr = RETR_RespRegr(ind_BOLD,:);
RETR_CardRegr = func_RETR_Card_regressors(time,PPGlocs,3);   NR_card=size(RETR_CardRegr,2);  RETR_CardRegr = RETR_CardRegr(ind_BOLD,:);

regr = [RETR_RespRegr, RETR_CardRegr];    regr =  detrend(regr,'linear');
rLinear = [1:NV]';
regr = [regr, rLinear, ones(NV,1)];

%% 6: Remove confounds from the 4 GS variants extracted in Section 4 and save variables

B = regr\GM;  GM_clean = GM - regr*B;
B = regr\WM;  WM_clean = WM - regr*B;
B = regr\CSF;  CSF_clean = CSF - regr*B;
B = regr\WB;  WB_clean = WB - regr*B;

save(path_output,'WB','GM','WM','CSF','WB_clean','GM_clean','WM_clean','CSF_clean')






























