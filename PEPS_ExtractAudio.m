% February 2020.         Programmed by: D. Bolger
% Use of function audioread to extract audio from the *.mp4 video files.
% Audio files are written to file in *.wav format with a sampling rate of
% 48kHz (the original video ampling rate). 
% *************************************************************************

vidfile = 'List3d-movie.mp4';
videoIn = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Movie-Files',vidfile);
InfoVid = audioinfo(videoIn);

audfile = 'List3d-audio.wav';
audioIn = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',audfile);

[vidAud,Fs]  = audioread(videoIn,'native');   %Extracts audio from the video
audiowrite(audioIn,vidAud,Fs);

InfoAud = audioinfo(audioIn);

sprintf('Duration of Video and extracted Audio is %fsecs and %fsecs, respectively.', InfoVid.Duration, InfoAud.Duration)