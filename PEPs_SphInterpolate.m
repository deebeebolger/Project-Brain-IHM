close all;
clear all;

sujnom = 's02';
films = {'Film1' 'Film2' 'Film4'};
datadir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Data_Preproc',sujnom);
dirsave = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Data_Preproc','Processed_Data_cont');

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session


for count = 1:size(films,2)
    currdir = fullfile(datadir,strcat(sujnom,'-',films{1,count}),filesep);
    currset = dir(strcat(currdir, '*-icarej.set'));
    
    EEG = pop_loadset('filename',currset.name,'filepath',currdir);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',EEG.setname,'gui','off');
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    %% CARRY OUT SPHERICAL SPLINE INTERPOLATION ON THE CURRENTLY LOADED DATASET
    
    title_interp=strcat(currset.name(1:end-4),'-ssinterp');
    EEG = pop_interp(EEG,EEG.chanlocs_prerej, 'spherical');  % EEGLAB spherical spline interpolation function (eeg_interp())
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(title_interp),'gui','off');
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
    cd(dirsave)
    [status, msg, msgID] = mkdir(sujnom);
    dirsave_curr = fullfile(dirsave,sujnom);
    
    EEG = pop_saveset( EEG, 'filename',char(title_interp),'filepath',dirsave_curr);  % Saves a copy of the current resampled dataset to the current directory
    EEG = eeg_checkset( EEG );
    eeglab redraw
    
end




