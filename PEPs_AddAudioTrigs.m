
currfile = 'Reg_Human_Incongruent-fb-markerchannel-doctor-v2.csv';
dataIn_file = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files','Marker-Files-doctor',currfile);
DataIn = readtable(dataIn_file);
vidoi = 'Reg_Human_Incongruent';

%% load in the EEG data and add the triggers to the event's structure. 


[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

onsets_all = DataIn{:,1};

findverb_int = find(contains({EEG.event.type},'Verbal-cong'));
newonsets_time = EEG.event(findverb_int(1)).etimes + onsets_all(2:end);
firstdiff = onsets_all(2)-onsets_all(1);
firstdoc = newonsets_time(1)-firstdiff;
newonsets_time = cat(1,firstdoc,newonsets_time);

time = EEG.times;
newonset_ms = newonsets_time.*1000;
lat = dsearchn(time',newonset_ms);
count = 0;

EEG.event = [];
EEG.urevent = [];

for counter = 1:length(newonsets_time)

    EEG.event(counter).type = [];
    EEG.event(counter).etimes = [];
    EEG.event(counter).latency = [];
    EEG.event(counter).urevent = [];
    EEG.event(counter).trialnum = [];
    EEG.event(counter).feedbacks = [];

    EEG.event(counter).type = [char(DataIn{counter,3}),'-',char(DataIn{counter,4})];
    EEG.event(counter).latency = lat(counter);
    EEG.event(counter).urevent = counter;
    EEG.event(counter).etimes = newonsets_time(counter);
    EEG.event(counter).feedbacks = char(DataIn{counter,5});

    if strcmp(char(DataIn{counter,3}),'Doctor') ==1

        count = count+1;
        EEG.event(counter).trialnum = count;

    elseif strcmp(char(DataIn{counter,3}),'Patient')

        EEG.event(counter).trialnum = count;
    end

    EEG.urevent(counter).type = [char(DataIn{counter,3}),'-',char(DataIn{counter,4})];
    EEG.urevent(counter).latency = lat(counter);
    
end


savedir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021','s39',...
    'AudTrigs','s39-Film4',filesep);

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off');
EEG = eeg_checkset( EEG );

EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',savedir);
EEG = eeg_checkset( EEG );
eeglab redraw




