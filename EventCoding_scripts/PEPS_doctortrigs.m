%% Load in textfile with the doctor trigger data.
close all
clear all
choice = 'doctor';
sujcurr = 's38';
vidnum = 'Film4';
currsujfile = [sujcurr,'-',vidnum];
mainpath = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021',filesep);

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


dircurr = fullfile(mainpath,sujcurr,currsujfile,filesep);
X = dir(dircurr);
allfiles = {X.name};
Xi = contains(allfiles,'pats.set');
currfile_open = allfiles{1,Xi};

EEG = pop_loadset('filename',currfile_open,'filepath',dircurr);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

video_name = EEG.video_name;

filep = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files','Marker-Files-doctor',filesep);
filep2 = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',filesep);


filenom = [video_name,'-fb-markerchannel-doctor-RS.txt'];
fileIn = fullfile(filep,filenom);
dataIn = readtable(fileIn);

filenom2 = [video_name,'_RS-fb-markerchannel.txt'];
fileIn2 = fullfile(filep2,filenom2);
dataIn2 = readtable(fileIn2);

%% Extract the doctor triggers
idoc = contains(dataIn{:,3},'Doctor-');
ipat = ~contains(dataIn{:,3},'Doctor-');

if strcmp(choice,'doctor') == 0


    isverb = find(contains({EEG.event.type},'Verbal'));
    ispat = find(ipat);
    
    for ecnt = 1:length(isverb)
        EEG.event(isverb(ecnt)).feedbacks = dataIn{ispat(ecnt),3};
    end

%%

newtitle = [currfile_open(1:end-4),'-pats.set'];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',newtitle,'gui','off');
EEG = pop_saveset( EEG, 'filename',newtitle,'filepath',dircurr);
eeglab redraw

elseif strcmp(choice, 'doctor')

    disp('*********Adding doctor triggers to the current dataset.********')
    
    isdoc = find(idoc);
    isverb = find(contains({EEG.event.type},'Verbal')); 
    doctimes = zeros(length(isdoc)-1,1);
    doclatency = zeros(length(isdoc),1);
    docfdback = cell(length(isdoc),1);
    doctrials = zeros(length(isdoc),1);

%     doctimes = zeros(length(isverb)-1,1);
%     doclatency = zeros(length(isverb),1);
%     docfdback = cell(length(isverb),1);
%     doctrials = zeros(length(isverb),1);
    

    x = cell2mat(cellfun(@isempty, {EEG.event.feedbacks},'UniformOutput',false));

    if sum(x)>0
       
        emptyfb = find(x);
    
        for cnt = 1:length(emptyfb)
            EEG.event(emptyfb(cnt)).feedbacks = 'gesture';
        end
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',EEG.setname,'gui','off');
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',dircurr);
        eeglab redraw
    end

    
    for evcnt = 1:length(isdoc)

        diffcur = dataIn{isdoc(evcnt)+1,1} - dataIn{isdoc(evcnt),1};
        doctimes(evcnt) =  EEG.event(isverb(evcnt)).etimes - diffcur;
        doclatency(evcnt) = dsearchn(EEG.times', doctimes(evcnt)*1000);
        docfdback{evcnt,1} = dataIn{isdoc(evcnt),3};
        doctrials(evcnt) = evcnt;
    end
    

    alllats = cat(1,[EEG.event.latency]',doclatency);
    alltimes = cat(1,[EEG.event.etimes]',doctimes);
    alltrials = cat(1,[EEG.event.trialnum]',doctrials);
    allfeedbacks = cat(1,[EEG.event.feedbacks]',cellstr(string(docfdback)));
    %allfeedbacks = cat(1,[EEG.event.feedbacks]',docfdback);
    alltypes = cat(1,{EEG.event.type}',cellstr(string(docfdback)));
    %alltypes = cat(1,{EEG.event.type}',docfdback);

    [latsort,sindx] = sort(alllats);
    etimes_sort = alltimes(sindx);
    trial_sort = alltrials(sindx);
    feedback_sort = cell(length(allfeedbacks),1);
    alltypes_sort = cell(length(allfeedbacks),1);

     for fcntr = 1:length(allfeedbacks)
        feedback_sort{fcntr,1} = allfeedbacks{sindx(fcntr),1};
        alltypes_sort{fcntr,1} = alltypes{sindx(fcntr),1};
     end

     %% Add the doctor triggers to the current EEG data structure event field.

     for counter = 1:length(feedback_sort)

         EEG.event(counter).type = alltypes_sort{counter,1};
         EEG.event(counter).latency = latsort(counter);
         EEG.event(counter).etimes = etimes_sort(counter);
         EEG.event(counter).feedbacks = feedback_sort{counter,1};
         EEG.event(counter).trialnum = trial_sort(counter);
         EEG.event(counter).urevent = counter;
     end
    

    newtitle = [currfile_open(1:end-4),'-doctrigs.set'];
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',newtitle,'gui','off');
    EEG = pop_saveset( EEG, 'filename',newtitle,'filepath',dircurr);
    eeglab redraw

    

end

% etimes_test = ([EEG.event.latency]-1)/EEG.srate;
% etimes_test1 = etimes_test.*1000;

