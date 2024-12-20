%% Task a
% Plotting Against Time
clc;
clear;
close all;

load('ecg_data.mat');
N=length(ecg_signal); % Number of Signals
Ts=(0:N - 1)/fs; % Calculating Time Period

% Variables to determine starting and ending time
T_start=0;
T_end=10;

ecg_selected_range=ecg_signal((T_start*fs+1):(T_end*fs));

figure(1);
plot(Ts(T_start*fs+1:T_end*fs), ecg_selected_range, 'b')
xlabel('Time (Seconds)');
ylabel('Amplitude (mV)');
title('First 10 Seconds (Time Domain)');

%% Task b
% Calculating BPM using time domain in previously selected range
Max_Amplitude=0.8;
Min_Distance=0.4;

% Finding Peaks
[peaks, loc]=findpeaks(ecg_selected_range,fs,'MinPeakHeight',Max_Amplitude,'MinPeakDistance',Min_Distance); 

% Calculates the different periods between each Heart Beat in the signal
difference=diff(loc);
average_BPM=mean(60./difference)

%% Task C
% Plotting ECG Signal in Frequency Domain
clc;
close all;

ECG_SIGNAL_SELECTED_RANGE=fft(ecg_selected_range); % Converting to Frequency Domain using FFT
N_Freq=length(ECG_SIGNAL_SELECTED_RANGE);

Freq=((-N_Freq/2):((N_Freq/2)-1))*fs/N_Freq;  % Shifting X-Axis to be centered around Zero

% Plotting ECG in Frequency domain using fftshift
plot(Freq,abs(fftshift(ECG_SIGNAL_SELECTED_RANGE)), 'b')
xlabel('Frequency (Hz)');
ylabel('Amplitude (mV)');
title('First 10 Seconds (Frequency Domain)');

%% Task D
% Designing a filter to remove Out-Of-Band Noise

clc;
close all;

load("Filter_object.mat")   % Loading the Filter
ecg_filtered=Hd.filter(ecg_selected_range); % Applying the filter

% Setting up a tiled layout to view multiple graphs
tiledlayout(2,1);
nexttile;

% Plotting the filtered time domain signal
plot(Ts(T_start*fs+1:T_end*fs), ecg_filtered)
xlabel('Time (Seconds)');
ylabel('Amplitude (mV)');
title('First 10 Seconds (Time Domain) Filtered');

% Calculating the filtered signal in Frequency Domain
ECG_SIGNAL_SELECTED_RANGE=fft(ecg_filtered);
N_Freq=length(ECG_SIGNAL_SELECTED_RANGE);
Freq=(-N_Freq/2:N_Freq/2-1)*fs/N_Freq;

% Plotting the filtered frequency domain signal
nexttile;
plot(Freq,abs(fftshift(ECG_SIGNAL_SELECTED_RANGE))/N_Freq)
xlabel('Frequency (Hz)');
ylabel('Amplitude (mV)');
title('First 10 Seconds (Frequency Domain) Filtered');

%% Task E
% Finding BPM using Frequency Domain

clc;

mag_spectrum=abs(fftshift(ECG_SIGNAL_SELECTED_RANGE));  % Calculating the Magnitude Spectrum
freq_range=find(Freq>=0.83 & Freq<=2.5);     % Selecting the location of range of frequencies of Human Heart Rate
[~,Max_Loc]=max(mag_spectrum(freq_range));   % Finding the locations of where the max amplitude is in the selected frequency range
BPM_Freq=60*abs(Freq(freq_range(Max_Loc)))   % calculating heart beat in the segment

%% Task F
% Plot the Heart Rate across time using Time Domain Signal

clc;
close all;
load("Filter_object.mat")
ecg_filtered_all=Hd.filter(ecg_signal); % Applying Filter

% Calcualting Heart Rate Across Time
[all_peaks,all_loc]=findpeaks(ecg_filtered_all,fs,'MinPeakHeight',Max_Amplitude,'MinPeakDistance',Min_Distance);
dECG=diff(all_loc);
BPM_Time=60./dECG;

% Plotting Heart Rate Across Time
plot(all_loc(1:length(BPM_Time)),smooth(BPM_Time))
title("BPM Across Time (Time Domain)");
xlabel('Time (Sec)');
ylabel('Beats Per Minut (BPM)');

%% Task G
% Plot the Heart Rate across time using Frequency Domain Signal
clc;
close all;

% Window Setup
window_length=5*fs;             % Specifying Window Size
FTT_Length=window_length*4;     % For FFT Length, Higher for better frequency resolution
window=triang(window_length);   % Type of Window, Hanning for best solution

% Converting to Spectrogram using Short Time Fourier Transform
[ECG_stft, Freq_stft, Time_stft]=stft(ecg_filtered_all, fs,"Window",window,"OverlapLength",window_length/2, "FFTLength",FTT_Length);
ecg_spectrogram=abs(ECG_stft); % Calculating Magnitude Spectrogram

% Selecting the range of frequencies of normal heart rate
Freq_range_stft=find(Freq_stft>=0.83 & Freq_stft<=2.5);
Freq_stft_selected=Freq_stft(Freq_range_stft);
ecg_spectrogram_selected_range = ecg_spectrogram(Freq_range_stft,:);

[~,freq_locs]=max(ecg_spectrogram_selected_range); % Finding the Frequencies with the highest amplitude across time
BPM_Freq_full=60*Freq_stft_selected(freq_locs); % Calculating the BPM across Time

plot(Time_stft, BPM_Freq_full, 'B');
title("BPM Across Time (STFT)");
xlabel('Time (Sec)');
ylabel('Beats Per Minut (BPM)');

%% Task H
% Finding Abnormal Heart Rate intervals
clc;

% Selecting the Interval of Lower than normal Heart Rate
Low_range=find(BPM_Time<60);
Low_interval=all_loc(Low_range);

% Selecting the Interval of Higher than normal Heart Rate
High_range=find(BPM_Time>100);
High_interval=all_loc(High_range);

% Plotting BPM Across time with highlighting the abnormal parts
plot(all_loc(1:length(BPM_Time)),smooth(BPM_Time),'b')
hold on
plot(High_interval,smooth(BPM_Time(High_range)),'ro')
plot(Low_interval,smooth(BPM_Time(Low_range)),'gx')
yline(100)
yline(60)
legend("BPM Value","Higher than Normal (>100)", "Lower than Normal (<60)")
title("BPM Across Time (Time Domain)");
xlabel('Time (Sec)');
ylabel('Beats Per Minut (BPM)');

%state driven loop where I iterate inside loops depending on the region the
%loop is on
%initializing the counter to an appropriate number for array indexing
i=1; 
while (i<=length(BPM_Time))
if(BPM_Time(i)>100);        %an abnormal heartbeat state has been detected
    first_num=all_loc(i); %keeping my first value at which the abnormality has been detected
    while((BPM_Time(i)>100) & (i<length(BPM_Time)) ) %looping until heartbeat rate stabilizes or until I get to the end of the recorded BPM
        i=i+1;
    end
    disp("critical region, BPM above 100: ["+first_num+","+all_loc(i)+"]")

elseif(BPM_Time(i)<60);  %an abnormal heart state has been detected
    first_num=all_loc(i); %keeping my first value at which the abnormality has been detected
    while((BPM_Time(i)<60) & (i<length(BPM_Time))) %looping until heartbeat rate stabilizes or until I get to the end of the recorded BPM
        i=i+1;
    end
      disp("critical region, BPM below 60: ["+first_num+","+all_loc(i)+"]")
end
i=i+1; %updating the counter at the end of the while loop
end












