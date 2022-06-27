
basefile_name = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021','Segmented',filesep);

allcurr_files = dir(basefile_name);
filecurr_names = {allcurr_files.name};
fcurr_names = filecurr_names([allcurr_files.isdir]);

sub_oi = '';  % Number of the current subject of interest.
torej = {'.' '..'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  %Open EEGLAB

for fcounter = 1:length(fcurr_names)
    
    if sum(strcmp(torej,fcurr_names{1,fcounter}))==0
        
        newbase = fullfile(basefile_name, fcurr_names{1,fcounter},'BLCorrected',filesep);
        currfiles = dir(newbase);
        currfiles_titles = {currfiles(~[currfiles.isdir]).name};
        
        if sum(contains(currfiles_titles, sub_oi))>0
            X = currfiles_titles(contains(currfiles_titles, sub_oi));
            EEG = pop_loadset('filename',{X{contains(X,'.set')}},'filepath',char(newbase));
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            eeglab redraw;
        end
        
    end
    
end
    
icong =  contains({ALLEEG.setname},'cong');  %Isolate the datasets with congruent trials. 
total_trials = sum([ALLEEG(icong).trials]);
cindx = find(icong);


%% This is required to make the feedback field of the events structure
% consistent before merging. 
for cint = 1:length(cindx)
    
    y = cellfun(@char, {ALLEEG(cindx(cint)).event.feedbacks},'UniformOutput',false);
    for ycnt = 1:length(y)
        ALLEEG(cindx(cint)).event(ycnt).feedbacks = y{ycnt};
    end
end

eegmerged = pop_mergeset(ALLEEG,cindx);
eegmerged.setname = ['s',sub_oi,'_Conds_all_BLC'];
eegmerged.video_name = string({ALLEEG(cindx).video_name});
eegmerged.rejchans = {ALLEEG(cindx).rejchans};
eegmerged.rejchand_indx = {ALLEEG(cindx).rejchan_indx};
eegmerged.trialnums = [ALLEEG(cindx).trials];


v = matlab.lang.makeValidName(eegmerged.setname);
eval([v,'=eegmerged'])

savedir_path = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021','Segmented','matfiles',filesep);

save(fullfile(savedir_path,[eegmerged.setname,'.mat']), '-struct', 's03_Conds_all_BLC');



    
   