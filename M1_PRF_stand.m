
clear, close all, clc
%%  1:  Load and display physiological variables (heart rate and respiration) and global signal (GS) from MAT-File

%  Set the following parameters !!

sc = 140;     % choose a scan (sc) from 1-164

% -----------------------------------------

load('../Data/HCP_41_subjects_phys_GS.mat')

GS = GS_all(:,sc);  HR=HR_all(:,sc); resp=zscore(resp_all(:,sc));
Ts_10 = 0.1 ;    Fs_10 = 10;                                                   % Sampling period in seconds
time_10 = 0:Ts_10:(length(HR)-1)*Ts_10;
timeMR = time_10(ind_BOLD_10);

figure('Position',[543         425        1588         792])

ax1 = subplot(3,1,1);
plot(time_10,HR)
title('Heart rate (HR)')
ylabel('HR (bpm)')

ax2 = subplot(3,1,2);
plot(time_10,resp)
title('Respiration (HR)')
ylabel('Amplitude (a.u.)')

ax3 = subplot(3,1,3);
plot(timeMR,GS);
title('Global signal (GS)')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')

linkaxes([ax1,ax2,ax3],'x')
xlim([0,max(time_10)])

%% 2: Extract respiration volume per time (RVT) from respiration

figure('Position',[543         425        1588         792])

ax1 = subplot(3,1,1);
time = time_10;
[pks,loc] = findpeaks(resp,time,'MinPeakDistance',2,'MinPeakHeight',0.2);
plot(time,resp), hold on
plot(loc,pks,'o')
respUpp = interp1([0,loc,time(end)],[pks(1),pks',pks(end)],time_10);
plot(time_10,respUpp)
title('Respiration')

[pks,loc] = findpeaks(-resp,time,'MinPeakDistance',2,'MinPeakHeight',0.2); pks=-pks;
plot(loc,pks,'x')
respLow = interp1([0,loc,time(end)],[pks(1),pks',pks(end)],time_10);
plot(time_10,respLow)
ylabel('Amplitude (a.u.)')

BR = 60./diff(loc);
time_BR = [time(1),(loc(2:end)-loc(1:end-1))/2+loc(1:end-1),time(end)];
BR = interp1(time_BR,[BR(1),BR,BR(end)],time_10);

ax2 = subplot(3,1,2);
plot(time_10,BR)
title(sprintf('Breathing rate (BR): %3.1f ±%3.1f ',mean(BR),std(BR)))
ylabel('BR (rpm)')

ax3 = subplot(3,1,3);
RVT = ((respUpp-respLow).*BR)';
plot(time_10,zscore(RVT))
title('Respiration volume per time (RVT)')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')

linkaxes([ax1,ax2,ax3],'x')
xlim([0,max(time_10)])

%% 3: Apply standard PRF model

t_IR = 0:Ts_10:60;
RRF = 0.6*t_IR.^(2.1).*exp(-t_IR/1.6)-0.0023*t_IR.^3.54.*exp(-t_IR/4.25);
RRF  = RRF/max(RRF);
CRF = 0.6*t_IR.^2.7.*exp(-t_IR/1.6)-(16/(sqrt(2*pi*9)))*exp(-(t_IR-12).^2/18);
CRF  = CRF/max(CRF);

NV = length(GS);

HR = smooth(HR,6*Fs_10);
HR_conv = conv(HR,CRF);
RVT_conv = conv(RVT,RRF);

xPhys = [HR_conv(ind_BOLD_10),RVT_conv(ind_BOLD_10)];   xPhys = detrend(xPhys,'linear');
regr = [ones(NV,1),xPhys];

B = regr\GS;     yPred = regr*B;

r_PRF(1) = corr(yPred,GS);
yPred_card = regr(:,2)*B(2);  r_PRF(2) = corr(yPred_card,GS);
yPred_resp = regr(:,3)*B(3);  r_PRF(3) = corr(yPred_resp,GS);

fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF (HR): %3.2f  \n',r_PRF(2))
fprintf('RRF (RVT): %3.2f  \n',r_PRF(3))
fprintf('CRF & RRF (HR & RVT): %3.2f  \n',r_PRF(1))



%%  4: Plot output of PRF model (timeseries and PRF curves)  

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
title('Cardiac Response Function (CRF_{stand}) ')
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
legend([h1,h2],'Global signal','X_{RVT}'), legend('boxoff')
ylabel('Amplitude (a.u.)')
ylim([yl1, yl2])


ax7 = subplot(5,3,[12,15]);
plot(t_IR,RRF,'LineWidth',4), grid on
title('Respiration response function (RRF_{stand}) ')
xlim([0 60])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

ax5 = subplot(5,3,13:14);
plot(time_10,RVT,'LineWidth',2), hold on
title('Respiration Volume per Time (RVT)')
ylabel('RVT (a.u.)')
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



%%   5: Create matrix of Physiological Regressors for the General linear Model

figure('Position', [ 316         673        1849         483])

plot(timeMR, xPhys)
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')

subject = scans_41_subjects{sc,1};
task = scans_41_subjects{sc,2};
title(sprintf('Physiological regressors to be included in the General Linear Model for scan %s (%s) ', subject, task),'Interpreter','none')










