%%     Preprocess_Phys.m

%   Use this script to:
%  1. Extract volume triggers for fMRI acquisition
%  2. Preprocess cardiac and respiratory signal and extract physiological variables (e.g. heart rate) 

%    =========================================
%% 1: Load physiological data for a particular scan

clear, clc, close all
task_list = {'Rest1_LR','Rest1_RF', 'Rest2_LR', 'Rest2_RL'};

%  Set the following parameters !!

subject_id = '118730';
task_id = task_list{4};    % list

%  ---------------------------------------
% path_phys: the path where the 'Physio_log.txt' is located. The figures
% and physiological variables extracted with this script are saved in the
% same path.

path_phys = ['../Data/S',subject_id,'/',task_id,'/'];
x = load([path_phys,'Physio_log.txt']);

resp = zscore(x(:,2));
cardiac = zscore(x(:,3));

% Create timeline at the frequency sampling Fs = 400 Hz that the recordings
% were acquried

N = length(resp);
Fs = 400; Ts=1/Fs;  time=0:Ts:(N-1)*Ts;

%    =========================================
%% 2: Find the triggers (or timepoints) that the fMRI volumes were acquired 

% The triggers in each scan of the resting-state fMRI data in the HCP
% should sum up to 1200

trig_orig = x(:,1);
tmp = diff(trig_orig);
loc = find(tmp>0);
trig = zeros(size(trig_orig));
trig(1) = 1;
trig(loc) = 1;

fprintf('Number of triggers found: %d  \n\n', sum(trig))

figure('Position', [240         800        2200         500])

plot(time',[trig*.5,resp,cardiac])
xlabel('Time (s)')
ylabel('Amplitude')
legend('Trigger','Respiration','Cardiac')
xlim([60 90])

%    =========================================
%% 3: Preprocess the cardiac signal and extract the heart rate (HR)

% The cardiac signal is first band-pass filtered and, subsequently, the
% peaks in the signal are identified. Based on the times of the peaks, the
% HR is estimated. As the cardiac signal is often noisy, there are some
% steps that we can take in order to improve the quality of the extracted
% HR.

% First, we have to specify the 'minPeak' variable which denotes the
% minimum time interval in seconds between peaks needed for the identification of the
% peaks. Normal values for 'minPeak' vary from 0.50 to 0.90 seconds. This value should be
% adjusted by the user for each scan separately based on visual inspection

% Then, we have to specify some parameters related to the correction for
% outliers in the extracted HR based on visual inspection. Typically, the
% 'filloutliers_ThresholdFactor' has to be somewhere between 3-20 in order
% to correct for outliers in noisy epochs while also retain real abrupt changes
% in HR. By zooming in a time interval with an abrupt change in HR, based
% on the quality of the cardiac signal and the time of peaks we can be more
% confident whether this abrupt change is real or artefactual.

%  Set the following parameters !!

minPeak = 0.55;
filloutliers_window = 30*Fs;     % given in number of samples
filloutliers_ThresholdFactor = 20;     %  normal range between 3-20
%  ---------------------------------------

Fs_10 = 10; Ts_10 = 1/Fs_10; time_10 = time(1):Ts_10:time(end);

f1 = 0.3; f2 = 10;     [filt_b,filt_a] = butter(2,[f1,f2]/(Fs/2));
cardiac_filt = filtfilt(filt_b,filt_a,cardiac);
[pks,PPGlocs] = findpeaks(cardiac_filt,time,'MinPeakDistance',minPeak);

figure('Position', [240         800        2200         500])

ax1 = subplot(2,1,1);
plot(time,cardiac_filt), hold on
b = zeros(1,length(PPGlocs));
for i = 1:length(PPGlocs)
    [~,b(i)] = min(abs(PPGlocs(i)-time));
end
u = zeros(N,1);
u(b) = 1;
plot(time,u)
legend({'PPG', 'Peaks in PPG'})
title('Photoplethysmogram (PPG)')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')

ax2 = subplot(2,1,2);
HR = 60./diff(PPGlocs);
time_HR = [time(1),(PPGlocs(2:end)-PPGlocs(1:end-1))/2+PPGlocs(1:end-1),time(end)];
HR_raw_400 = interp1(time_HR,[HR(1),HR,HR(end)],time);
HR_filloutl_400 = filloutliers(HR_raw_400,'linear','movmedian',filloutliers_window,'ThresholdFactor',filloutliers_ThresholdFactor);
HR = interp1(time,HR_filloutl_400,time_10); HR=HR(:);

plot(time,HR_raw_400), hold on
plot(time_10, HR)
title(sprintf('Heart rate (%d±%d bpm)',round(mean(HR)),round(std(HR))))
legend('Raw','Corrected for outliers')

linkaxes([ax1,ax2],'x')
xlim([time(1),time(end)])

savefig([path_phys,'\Heart_rate.fig'])

%    =========================================
%% 4: Preprocess the respiratory signal

% The respiratory signal is first low-pass filtered and, subsequently,
% linearly detrended and corrected for outliers. The 
% parameters for the outlier correction have to be adjusted for each scan
% separately based on visual inspection.

%  Set the following parameters !!

filloutliers_window = 0.3*Fs;     % given in number of samples
filloutliers_ThresholdFactor = 0.2;     
%  ---------------------------------------

close all

f1=0.01; f2=5; [filt_b,filt_a] = butter(2,[f1,f2]/(Fs/2));

resp = zscore(x(:,2));
resp = detrend(resp,'linear');
resp1 = filloutliers(resp,'linear','movmedian',filloutliers_window,'ThresholdFactor',filloutliers_ThresholdFactor);
resp2 = zscore(filter(filt_b,filt_a,resp1));


figure('Position', [240         800        2200         500])
plot(time,resp), hold on
plot(time,resp1)
plot(time,resp2)


%    =========================================
%% 5: Extract breathing-related variables

% After the respiratory signal is preprocessed, it gets downsampled at 10
% Hz, and the following variables are extracted: breathing rate (BR),
% respiration volume per time (RVT) and respiratory flow (RF).

resp = resp2;
resp_10 = interp1(time,resp,time_10);

figure('Position', [90, 130, 2300, 1200])

ax1 = subplot(5,2,1);
hist(resp,100), title('Histogram of respiration')

ax2 = subplot(5,2,3:4);
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

BR = 60./diff(loc);
time_BR = [time(1),(loc(2:end)-loc(1:end-1))/2+loc(1:end-1),time(end)];
BR = interp1(time_BR,[BR(1),BR,BR(end)],time_10);

ax3 = subplot(5,2,5:6);
plot(time_10,BR)
title(sprintf('Breathing rate (BR): %3.1f ±%3.1f ',mean(BR),std(BR)))
ylabel('BR (rpm)')

ax4 = subplot(5,2,7:8);
RVT = ((respUpp-respLow).*BR)';

plot(time,zscore(resp)), hold on, plot(time_10,zscore(RVT))
legend('Respiration','Respiration volume per time (RVT)')

resp_s = smooth(resp_10,10*1.5) ;
RF = diff(resp_s); RF=[0;RF(:)]; RF = RF.^2;

ax5 = subplot(5,2,9:10);
plot(time_10,RF,'g'),   
ylim([0 10*std(RF)])
title('Respiratory flow (RF)')
ylabel('RF (a.u.)')

linkaxes([ax2,ax3,ax4,ax5],'x')
xlim([time(1),time(end)])

savefig([path_phys,'\Respiration.fig'])


%    =========================================
%% 6: Save extracted variables to the file 'Physio_and_triggers.mat'

TR = 0.72;   % Repetition time (TR) of fMRI acquisition in seconds
save([path_phys,'Physio_and_triggers.mat'],'trig','time','Fs','TR','PPGlocs','HR','Fs_10','resp','resp_10','BR','RVT','RF')













