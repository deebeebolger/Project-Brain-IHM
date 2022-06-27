% February 2020.         Programmed by: D. Bolger
% Use of function audioread to extract audio from the *.mp4 video files.
% Audio files are written to file in *.wav format with a sampling rate of
% 48kHz (the original video sampling rate). 
% *************************************************************************

vidfile = 'Reg_Human_Congru_AA075.mp4';

videoIn = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','IHMBrain-Videos',vidfile);
InfoVid = audioinfo(videoIn);

audfile = 'Reg_Human_Congru_AA075.wav';
audioIn = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','TriggerTests',audfile);

[vidAud,Fs]  = audioread(videoIn,'native');   %Extracts audio from the video
audiowrite(audioIn,vidAud,Fs);

InfoAud = audioinfo(audioIn);

sprintf('Duration of Video and extracted Audio is %fsecs and %fsecs, respectively.', InfoVid.Duration, InfoAud.Duration)

%% Resample to from 48kHz to 44.1kHz - the optimum sampling rate for carrying out speech analysis using dedicated toolbox.

audres = resample(double(vidAud),44100,Fs);

audfilers = 'List1a-audio-rs.wav';
audioIn_rs = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',audfilers);
audiowrite(audioIn_rs,audres,44100);
InfoAud_rs = audioinfo(audioIn_rs);

sprintf('Duration of Video and extracted, resampled Audio is %fsecs and %fsecs, respectively.', InfoVid.Duration, InfoAud.Duration)
