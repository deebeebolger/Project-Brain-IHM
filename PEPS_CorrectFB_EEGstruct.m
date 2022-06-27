% PEPS_CorrectFB_EEGstruct()
% Programmed by: Deirdre Bolger
% This corrects the current subject's event's structure based on the
% pre-defined feedback onsets in the excel file.
% The resulting *.set file has the post-fix "trigcorr" added.
%************************************************

close all
clear all
sujcurr = 's01';
vidnum = 'Film3';
video_name = 'Reg_Agent_Incongruent';

dirgen = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2',...
    'PEPs_DataPreproc_2021',sujcurr);
dirIn = fullfile(dirgen,strcat(sujcurr,'-',vidnum),filesep);
load(fullfile(dirgen,'feedback_summary_correct.mat'));
tableIn = feedbacks{1,str2double(vidnum(end))};


allfiles= dir(dirIn);
fileIndex = find(~[allfiles.isdir]);
fileeeg = dir(strcat(dirIn,['*.set']));


% Open EEGLAB session and load in the corresponding *.set file.
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);


EEG = pop_loadset(fullfile(dirIn,fileeeg.name));
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw


for i = 1:size(tableIn,1)     %length(EEG.event)
    
    if i > length(tableIn{:,1})

        EEG.event(i).latency = [];
        EEG.event(i).type=[];
        EEG.event(i).urevent=[];
        EEG.event(i).trialnum=[];
          
    else
        
        if i==1
            
            EEG.event(i).latency = dsearchn(EEG.times',tableIn{i,1}*1000);  
            EEG.event(i).etimes = ([EEG.event(i).latency]-1)/EEG.srate;
        else
            addon = (tableIn{i-1,4}+tableIn{i-1,5})*EEG.srate;
            EEG.event(i).latency = EEG.event(i-1).latency+addon;
            EEG.event(i).etimes = ([EEG.event(i).latency]-1)/EEG.srate;
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


