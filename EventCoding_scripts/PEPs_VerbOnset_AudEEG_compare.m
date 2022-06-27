% Script to verify the extent to which the verb onsets in the auditory
% stimulus and the EEG coincide for a given EEG dataset.
close all;
clear all;

audionom = 'Curare_Human_Incongruent_RS.wav';
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

%% Now load in the EEG dataset
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset(); % Choose *.set file manually. 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

if ~contains(audionom, EEG.video_name)
    display('*****The current dataset and audio file do not correspond!!*****')
else
    display('*****Current dataset and audio file correspond*****')
    
    X = contains(string({EEG.event.type}),'Fb verbal');
    T = [EEG.event.etimes];
    VerbT = T(X);
    diffVerbT = diff(VerbT)';
    diffAudT = diff(fbinfo{:,1});
    diffcomp = [diffVerbT, diffAudT];
end

%% 