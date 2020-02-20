% PEPS_CorrectFB_EEGstruct()
close all
clear all
dirIn = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-PEPs','PEPs_Preprocess_Subjects','s01',filesep);
vidnum = 'Film4';
video_name = 'Curare_Human_Incongruent';
allfiles= dir(dirIn);
fileIndex = find(~[allfiles.isdir]);
fileeeg = dir(strcat(dirIn,['*',vidnum,'.set']));
filetxt = dir(strcat(dirIn,['*',vidnum,'.txt']));

tableIn = readtable(fullfile(dirIn,filetxt.name));   %Read in textfile with feedback onset times and names.

% Open EEGLAB session and load in the corresponding *.set file.
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);


EEG = pop_loadset(fullfile(dirIn,fileeeg.name));
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw


for i = 1:length(EEG.event)
    
    if i>length(tableIn{:,1})
        EEG.event(i).latency=[]
        EEG.event(i).type=[];
        EEG.event(i).urevent=[];
        EEG.event(i).trialnum=[];
    else
        
        if i==1
            
            EEG.event(i).latency = tableIn{i,1}*EEG.srate;
        else
            addon = tableIn{i-1,4}*EEG.srate;
            EEG.event(i).latency = EEG.event(i-1).latency+addon;
        end
        
        EEG.event(i).type = char(tableIn{i,2});
        EEG.event(i).urevent = i;
        EEG.event(i).trialnum = tableIn{i,3};
    end
    
end


EEG.video_name = video_name;

newtitre = [fileeeg.name(1:end-4),'-trigcorr'];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',newtitre,'gui','off');
EEG = pop_saveset( EEG, 'filename',newtitre,'filepath',dirIn);
eeglab redraw


