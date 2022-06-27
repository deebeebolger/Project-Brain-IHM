
audionom = 'Curare_Human_Congruent_RS.wav';
auddir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',filesep);
filename = fullfile(auddir,audionom);
wavIn = miraudio(filename,'Normal', 'Sampling',44100);
sF = get(wavIn,'Sampling'); %Extract the sampling frequency. 
alldat = mirgetdata(wavIn); % Extract the actual data. 
time = 0:1/sF{1,1}:(1/sF{1,1})*length(alldat);  %Calculate the time vector

[Y, FS] = audioread(filename);
info = audioinfo(filename);

%% Load in the textfile associated with the audiofile with feedback onsets.

txtfiles = dir(strcat(auddir, '*.txt')); % Find the wavfiles.
txtfiles_nom = {txtfiles.name};
tfm = strfind(txtfiles_nom, audionom(1:end-4)); 
itxt = [cell2mat(cellfun(@isempty, tfm, 'UniformOutput',false))==0];
txtcurr = txtfiles_nom{1, itxt};

txtfulldir = fullfile(auddir, txtcurr); 
fbinfo = readtable(txtfulldir,'ReadVariableNames',false);
fbindx = zeros(length(fbinfo{:,1}),1);

% Resample the audio to match the EEG sampling rate.
Fs_new = 512;
[Numer, Denom] = rat(Fs_new/sF{1,1});
alldat_rs = resample(alldat, Numer, Denom);
tnew = 0:1/Fs_new:(1/Fs_new)*length(alldat_rs);   %New time vector
tnew = tnew(1:end-1);

fbtime = nan(length(tnew),1);


for fcnt = 1:length(fbindx)  %Mark the time points with triggers.
    
    fbindx(fcnt) = dsearchn(tnew',fbinfo{fcnt,1});
    fbtime(fbindx(fcnt)) = 10;
    
end

%% 

freqbm = equal_xbm_bands(100,10000,10); % Subroutinte to divide the frequency interval into N bands of equal width along the BM. (Chimera toolbox)
b = quad_filt_bank(freqbm, sF{1,1});   % Create a bank of fir complex filters. 
wavfilt = zeros(size(b,2),length(alldat));
wavenv = zeros(size(b,2),length(alldat));

for bcnt = 1:size(b,2)
    plot(1:size(b,1),abs(b(:,bcnt)))
    hold on
end

for k = 1:size(b,2)
	wavfilt(k,:) = fftfilt(b(:,k), alldat);        % Apply the filter bank (fir complex filters)
    wavenv(k,:) = smooth(abs(wavfilt(k,:)),1000,'lowess');   % The envelope is the abs of complex signal as real & imaginary parts in quadrature. 
    
end

figure('Color',[1 1 1], 'Position', [0 1 1600 1000]);
for pcnt = 1:size(wavfilt,1)
    subplot(10,1,pcnt)
    plot(time(1:end-1),wavfilt(pcnt,:))
    title(['Cut-off frequency as BM width: ',num2str(ceil(freqbm(pcnt))),'-',num2str(ceil(freqbm(pcnt+1))), 'Hz']);
end

figure('Color',[1 1 1], 'Position', [0 1 1600 1000]);
for ecnt = 1:size(wavenv,1)
    subplot(10,1,ecnt)
    plot(time(1:end-1), wavenv(ecnt,:),'k')   %apply smoothing for visualisation
    title(['Narrow-band Envelope: Cut-off frequency as BM width: ',num2str(ceil(freqbm(ecnt))),'-',num2str(ceil(freqbm(ecnt+1))), 'Hz']);
end

%% To get wide band envelope, sum the narrow-band envelopes.

envwband = sum(wavenv,1);

figure('Color',[1 1 1],'Position',[0 1 1800 300])
plot(time(1:end-1), envwband,'r')
title('Wide band envelope (100hz - 10kHz)')
set(gca,'YGrid','on', 'XGrid','on')


%% Carry out the same envelope extraction but on the resampled waveform.
% The audio needs to be aligned with the EEG. This implies re-sampling the envelope to
% the EEG sampling rate (512Hz)

envwband_rs = resample(envwband, Numer, Denom);

figure('Color',[1 1 1],'Position',[0 1 1800 300])
plot(tnew, envwband_rs,'r')
title('Resampled Wide band envelope (100hz - 10kHz)')
set(gca,'YGrid','on', 'XGrid','on')

figure
plot(tnew, envwband_rs)
hold on
stem(tnew,fbtime,'r')
scrollplot(gca,'WindowSizeX',5, 'MinX',0)
x = gca;
x.Title.String = audionom;
