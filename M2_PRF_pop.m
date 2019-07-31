
clear, close all, clc
%%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from MAT-File

load('../Data/HCP_41_subjects_phys_GS.mat')

[NP, nScans] = size(HR_all);
resp_all = zscore(resp_all);

Ts_10 = 0.1 ;                                                       % Sampling period in seconds
time_10 = 0:Ts_10:(NP-1)*Ts_10;
timeMR = time_10(ind_BOLD_10);

RF_all = zeros(NP, nScans);
for sc = 1:nScans
 
    subject = scans_41_subjects{sc,1};
    task = scans_41_subjects{sc,2};
    fprintf('Scan: %d,    Subject: %s,    Task: %s    \n', sc, subject, task)
    
    resp = resp_all(:,sc);
    resp_s = smooth(resp,10*1.5) ;
    RF=diff(resp_s); RF=[0;RF(:)]; RF = RF.^2;
    RF_all(:,sc) = RF;       
end


%% 2: Estimate PRF parameters

ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',100,'Display','iter','UseParallel',1);   % Display: iter
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'UseParallel',true,'MaxIterations',100,'MaxFunctionEvaluations',3000 ,'OptimalityTolerance',1e-8 );

PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5,  -1.1,   -2.6  ];  
lb(1:8)=0;  ub(1:8)=20; ub(2:2:8)=3;   lb(9:10) = -inf; ub(9:10)=inf;

h = @(P) func_M2_PRF_pop(P,Ts_10,HR_all,RF_all,ind_BOLD_10,GS_all);

tic
% Uncomment the following line if you want to use  Genetic Algorithm
% (GA). GA may yield better fit with the cost of longer computational time.
% PRF_par = ga(h,length(lb),[],[],[],[],lb,ub,[],[],ga_opts);
PRF_par = fmincon(h,PRF_par,[],[],[],[],lb,ub,[],options);
fprintf('Minutes elapsed: %3.1f  \n',  toc/60)

[obj_function,CRF,RRF, r_PRF_pop ] = h(PRF_par);

fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF & RRF (HR & RF): r=%3.2f  \n',mean(r_PRF_pop))

%%  3: Plot PRF curves 

t_IR = 0:Ts_10:(length(CRF)-1)*Ts_10;

figure('Position',[ 1000         754         613         584])
subplot(2,1,1)
plot(t_IR,CRF,'LineWidth',4), grid on
title('Cardiac Response Function (CRF_{pop}) ')
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
xlim([0 60])

subplot(2,1,2)
plot(t_IR,RRF,'LineWidth',4), grid on
title('Respiration response function (RRF_{pop}) ')
xlim([0 60])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')


%%  4: Apply PRF model (timeseries and PRF curves)  ----------------

%  Set the following parameters !!

sc = 140;     % choose a scan (sc) from 1-164

% -----------------------------------------

subject = scans_41_subjects{sc,1};
task = scans_41_subjects{sc,2};
fprintf('Scan: %d,    Subject: %s,    Task: %s    \n', sc, subject, task)

GS = zscore(GS_all(:,sc)); HR = zscore(HR_all(:,sc)); RF = zscore(RF_all(:,sc));
NV = length(GS);

HR_conv = conv(HR,CRF); HR_conv_MR = HR_conv(ind_BOLD_10);
RF_conv = conv(RF,RRF); RF_conv_MR = RF_conv(ind_BOLD_10);

xPhys = [HR_conv_MR, RF_conv_MR];     xPhys = detrend(xPhys,'linear');
regr = [ones(NV,1), xPhys];

B = regr\GS;     yPred = regr*B;
r_PRF(1) = corr(GS, yPred);

yPred_card = regr(:,2)*B(2);  r_PRF(2) = corr(yPred_card,GS);
yPred_resp = regr(:,3)*B(3);  r_PRF(3) = corr(yPred_resp,GS);

HR_conv = B(2)*HR_conv;   HR_conv = HR_conv(1:length(HR));   HR_conv_MR = HR_conv(ind_BOLD_10);
RF_conv = B(3)*RF_conv;     RF_conv = RF_conv(1:length(RF));     RF_conv_MR = RF_conv(ind_BOLD_10);

fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF (HR): %3.2f  \n',r_PRF(2))
fprintf('RRF (RF): %3.2f  \n',r_PRF(3))
fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF(1))


%%  5: Plot output of PRF model (timeseries and PRF curves)  

%  Set the following parameters !!

smoothPar = 5;
fontTitle = 20;
fontLabels = 8;
fontTxt = 16;
lineW = 3;
yl1 = -5.3; yl2 = 5.5;

% -----------------------------------------

screenSize = get(0,'ScreenSize'); xL = screenSize(3); yL = screenSize(4);
figure
set(gcf, 'Position', [0.2*xL 0.2*yL  0.6*xL 0.6*yL ]);
set(gcf, 'Position', [0.1*xL 0.1*yL  0.8*xL 0.8*yL ]);

ax1 = subplot(5,3,1:2);
plot(time_10,HR)
ylabel('HR (bpm)')
title(sprintf('Heart rate (HR; %2.0f±%1.0f bpm )',mean(HR),std(HR)))

ax6 = subplot(5,3,[3,6]);
plot(t_IR,CRF,'LineWidth',4), grid on
title('Cardiac Response Function (CRF_{pop}) ')
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
xlim([0 60])

ax2 = subplot(5,3,4:5);
h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2=plot(timeMR,yPred_card,'LineWidth', lineW);
legend([h1,h2],'Global signal', 'X_{HR}')
title('BOLD fluctuations due to changes in HR')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(2)) ,'FontSize',fontTxt,'FontWeight','bold')
ylabel('Amplitude (a.u.)')
ylim([yl1, yl2])
legend('boxoff')


ax3 = subplot(5,3,7:8);
h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2=plot(timeMR,yPred,'LineWidth',lineW);
title('Full model')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(1)) ,'FontSize',fontTxt,'FontWeight','bold')
ylabel('Amplitude (a.u.)')
legend([h1,h2],'Global signal','X_{FM}')
ylim([yl1, yl2])
legend('boxoff')

ax4 = subplot(5,3,10:11);
h1 = plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2 = plot(timeMR,yPred_resp,'LineWidth',lineW);
title('BOLD fluctuations due to changes in respiration')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(3)) ,'FontSize',fontTxt,'FontWeight','bold')
legend([h1,h2],'Global signal','X_{RF}'), legend('boxoff')
ylabel('Amplitude (a.u.)')
ylim([yl1, yl2])


ax7 = subplot(5,3,[12,15]);
plot(t_IR,RRF,'LineWidth',4), grid on
title('Respiration response function (RRF_{pop}) ')
xlim([0 60])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

ax5 = subplot(5,3,13:14);
plot(time_10,RF,'LineWidth',2), hold on
title('Respiratory Flow (RF)')
ylabel('RF (a.u.)')
xlabel('Time (s)')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
xlim([timeMR(1) timeMR(end)])

ax_list = [ax1,ax2,ax3,ax4,ax5,ax6,ax7];
for ax=ax_list
    subplot(ax)
    ax.XGrid = 'on';
    ax.GridAlpha=0.7;
    ax.GridLineStyle='--';
    ax.FontSize = 16;
    ax.FontWeight = 'bold';    
end



%%   6: Create matrix of Physiological Regressors for the General linear Model 

figure('Position', [ 316         673        1849         483])

plot(timeMR, xPhys)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')

subject = scans_41_subjects{sc,1};
task = scans_41_subjects{sc,2};
title(sprintf('Physiological regressors to be included in the General Linear Model for scan %s (%s) ', subject, task),'Interpreter','none')











