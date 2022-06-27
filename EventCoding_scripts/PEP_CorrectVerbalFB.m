close all
clear all

filenom = 'Curare_Human_Incongruent_RS_fb_markerchannel.txt';
fileIn = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',filenom);
XIn = readtable(fileIn);
display(XIn)

SOA = diff(XIn{:,1});
fborig = XIn{:,3};

%% Start up EEGLAB AND LOAD IN DATASET

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset(); % Choose *.set file manually. 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw


if ~contains(filenom, EEG.video_name)
    display('*****The current dataset and audio file do not correspond!!*****')
else
    display('*****Current dataset and audio file correspond*****')
    
    % Find onset times of verbal feedbacks (Fb verbal) only.
    ifbv = contains({EEG.event.type},'Fb verbal'); % Indices of verbal feedbacks.
    alllats = [EEG.event.latency];
    fbvlats = alllats(ifbv)./EEG.srate;  % verbal feedbacks in seconds
    SOA_fbv = diff(fbvlats)';
    
    soa_diff = SOA-SOA_fbv;
    S =[SOA, SOA_fbv, soa_diff];
    fbvlats = fbvlats';
    indx = find(ifbv);
    
    for i = 1:length(fbvlats)-1
        EEG.event(indx(i+1)).latency = (fbvlats(i)+SOA(i))*EEG.srate;
        EEG.event(indx(i)).feedback = fborig{i,1};
        
        if contains(EEG.event(indx(i)+1).type,'gestuel') && EEG.event(indx(i)).trialnum==EEG.event(indx(i)+1).trialnum
            EEG.event(indx(i)+1).latency = EEG.event(indx(i)).latency;
        end
        
    end
    
end

fnom = [EEG.setname,'v2'];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom),'gui','off'); % Create a new dataset for the current raw datafile
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw   

