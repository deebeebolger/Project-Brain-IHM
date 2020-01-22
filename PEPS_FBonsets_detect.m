% Date: 29-11-2019             Programmed by: D. Bolger
% Script to detect the onset times of the feedbacks of interest and to
% determine the modality and type of each feedback.
% Script applied in the Brain-IHM project.
% ********************************************************
%% LOAD IN THE *.bdf FILE OF THE CURRENT SUBJECT
% Loads in the *.bdf file, ensuring that 74 channels are included so that the 
% ERGO1 and ERGO2 data are loaded.

Dirbase = fullfile(filesep,'Volumes','KINGSTON','PEPS_Sujets','s02',filesep);

allfiles= dir(Dirbase);
fileIndex = find(~[allfiles.isdir]);
filenum = dir(strcat(Dirbase,'*.bdf'));                      %find all the *.bdf files in the current folder
filenom = {filenum.name};
ergsig1 = cell(1,length(filenom));
ergsig2 = cell(1,length(filenom));

% OPEN EEGLAB SESSION
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

for scnt = 1:length(filenom)
    
    fullDir = strcat(Dirbase,filenom{1,scnt});
    fnom = filenom{1,scnt}(1:end-4);

    % The following three lines is added to resolve a bug occurring when
    % opening the *.bdf file.
    x = fileparts( which('sopen') );
    rmpath(x);
    addpath(x,'-begin');

    % Opening up *.bdf file and saving as a *.set file.
    EEG = pop_biosig(fullDir, 'channels',[1:76], 'ref', [] ,'refoptions',{'keepref' 'off'} );
    
    
    ergsig1{1,scnt} = EEG.data(75,:);
    ergsig2{1,scnt} = EEG.data(76,:);
    
    EEG = pop_select( EEG,'nochannel',{'Erg1' 'Erg2'});
    EEG = eeg_checkset( EEG );
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom),'gui','off'); % Create a new dataset for the current raw datafile
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom),'filepath',Dirbase);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
  
end


%% **************************EXTRACT PHOTODIODE ONSETS****************************
fb_sessions = filenom;  
thresh_val = 0;
fs = EEG.srate;
[b,a] = butter(2,8./(fs/2));
onoffsets_ergo1 = cell(1,size(ergsig1,2));
onoffsets_ergo2 = cell(1,size(ergsig1,2));
onsets_ergo1 = cell(1,size(ergsig1,2));
onsets_ergo2 = cell(1,size(ergsig1,2));
offsets_ergo1 = cell(1,size(ergsig1,2));
offsets_ergo2 = cell(1,size(ergsig1,2));
Times_bloc = cell(length(fb_sessions),1);

for ergcnt = 1:size(ergsig1,2)
    ergcnt
    Times_bloc{ergcnt,1} = ALLEEG(ergcnt).times;
    time = Times_bloc{ergcnt,1};
    
    ergsig1_curr = ergsig1{1,ergcnt};
    ergsig2_curr = ergsig2{1,ergcnt};
    ergsig1D = detrend(ergsig1_curr,0);
    ergsig2D = detrend(ergsig2_curr,0);
    
    ergsig1Dfilt = filtfilt(b,a,double(ergsig1D));
    ergsig2Dfilt = filtfilt(b,a,double(ergsig2D));
    
    % Fit a linear curve to filtered data to correct possible trend
    p1 = polyfit(1:length(ergsig1Dfilt),ergsig1Dfilt,1);
    f1 = polyval(p1,1:length(ergsig1Dfilt));
    ergsig1Dfilt_corr = ergsig1Dfilt - f1;
    
    p2 = polyfit(1:length(ergsig2Dfilt),ergsig2Dfilt,1);
    f2 = polyval(p2,1:length(ergsig2Dfilt));
    ergsig2Dfilt_corr = ergsig2Dfilt - f2;

    %Half wave rectify the filtered signal and invert.
    ergsig1_hwr = zeros(size(ergsig1Dfilt_corr));
    ipos = find(ergsig1Dfilt_corr>0);
    ergsig1_hwr(ipos)= ergsig1Dfilt_corr(ipos).*1;
    
    ergsig2_hwr = zeros(size(ergsig2Dfilt_corr));
    ipos1 = find(ergsig2Dfilt_corr<0);
    ergsig2_hwr(ipos1) = ergsig2Dfilt_corr(ipos1).*-1;

    % Set all activity <= mean activity to zero.
    ergsig1_hwr(ergsig1_hwr<=max(ergsig1_hwr)/2) = 0; ergsig1_hwrz = ergsig1_hwr;
    ergsig2_hwr(ergsig2_hwr<=max(ergsig2_hwr)/2) = 0; ergsig2_hwrz = ergsig2_hwr;

    % Turn all trigger signals into step functions. 
     ergsig1_hwrst = ergsig1_hwrz;
     ergsig2_hwrst = ergsig2_hwrz;
     ergsig1_hwrst( ergsig1_hwrz>0) = 1;
     ergsig2_hwrst( ergsig2_hwrz>0) = 1;

    % Find onsets (diff(D_hwr == 1)) and offsets (diff(D_hwr == -1))
    erg1_diff = diff(ergsig1_hwrst);  % 1OD of D_hwr
    erg1_diff = cat(2,erg1_diff,0);
    erg2_diff = diff(ergsig2_hwrst);  
    erg2_diff = cat(2,erg2_diff,0);

    [pks_sig1,minima_sig1,locs_pks_sig1,locs_min_sig1]= CREx_peakfinder(erg1_diff); 
    [pks_sig2,minima_sig2,locs_pks_sig2,locs_min_sig2]= CREx_peakfinder(erg2_diff); 
    
    
    if locs_min_sig1(1)<locs_pks_sig1(1)
        onoffsets_sig1 = [time(locs_pks_sig1);time(locs_min_sig1(2:end))]';
        durs_sig1 = onoffsets_sig1(:,2) - onoffsets_sig1(:,1);
        onoffsets_ergo1{1,ergcnt} = cat(2,onoffsets_sig1,durs_sig1);
        onsets_all1 = nan(size(time));
        offsets_all1 = nan(size(time));
        onsets_all1(erg1_diff == pks_sig1(1)) = 0;
        offsets_all1(erg1_diff== minima_sig1(1)) = 0;
        onsets_ergo1{1,ergcnt} = onsets_all1;
        offsets_ergo1{1,ergcnt} = offsets_all1;
        
      
    else
        onoffsets_sig1 = [time(locs_pks_sig1);time(locs_min_sig1)]';
        durs_sig1 = onoffsets_sig1(:,2) - onoffsets_sig1(:,1);
        onoffsets_ergo1{1,ergcnt} = cat(2,onoffsets_sig1,durs_sig1);
        onsets_all1 = nan(size(time));
        offsets_all1 = nan(size(time));
        onsets_all1(erg1_diff == pks_sig1(1)) = 0;
        offsets_all1(erg1_diff== minima_sig1(1)) = 0;
        onsets_ergo1{1,ergcnt} = onsets_all1;
        offsets_ergo1{1,ergcnt} = offsets_all1;
    end
    
       
    if locs_min_sig2(1)<locs_pks_sig2(1)
       
        onoffsets_sig2 = [time(locs_pks_sig2);time(locs_min_sig2(2:end))]';
        durs_sig2 = onoffsets_sig2(:,2) - onoffsets_sig2(:,1);
        onoffsets_ergo2{1,ergcnt} = cat(2,onoffsets_sig2,durs_sig2);
        onsets_all2 = nan(size(time));
        offsets_all2 = nan(size(time));
        onsets_all2(erg2_diff == pks_sig2(1)) = 0;
        offsets_all2(erg2_diff== minima_sig2(1)) = 0;
        onsets_ergo2{1,ergcnt} = onsets_all2;
        offsets_ergo2{1,ergcnt} = offsets_all2;
    else
        
        onoffsets_sig2 = [time(locs_pks_sig2);time(locs_min_sig2)]';
        durs_sig2 = onoffsets_sig2(:,2) - onoffsets_sig2(:,1);
        onoffsets_ergo2{1,ergcnt} = cat(2,onoffsets_sig2,durs_sig2);
        onsets_all2 = nan(size(time));
        offsets_all2 = nan(size(time));
        onsets_all2(erg2_diff == pks_sig2(1)) = 0;
        offsets_all2(erg2_diff== minima_sig2(1)) = 0;
        onsets_ergo2{1,ergcnt} = onsets_all2;
        offsets_ergo2{1,ergcnt} = offsets_all2;
    end
    
    figure('Name',fb_sessions{1,ergcnt},'NumberTitle','off');
    subplot(2,2,1)
    plot(time,ergsig1_hwrst);
    hold on
    plot(time,onsets_ergo1{1,ergcnt},'or','MarkerFaceColor','r')
    hold on
    plot(time,offsets_ergo1{1,ergcnt},'og','MarkerFaceColor','g'); 
    set(gca,'YLim',[0 1.5])
    title('Photodiode Signal as Step-function: onsets and offsets'); 
    subplot(2,2,2)
    plot(time,ergsig1Dfilt_corr)
    hold on
    plot(time,onsets_ergo1{1,ergcnt},'or','MarkerFaceColor','r');
    hold on
    plot(time,offsets_ergo1{1,ergcnt},'og','MarkerFaceColor','g');
    title('Original photodiode signal (detrended+filtered) with onsets (red) and offsets (green) (ERG1)');
    
    subplot(2,2,3)
    plot(time,ergsig2_hwrst);
    hold on
    plot(time,onsets_ergo2{1,ergcnt},'or','MarkerFaceColor','r')
    hold on
    plot(time,offsets_ergo2{1,ergcnt},'og','MarkerFaceColor','g'); 
    set(gca,'YLim',[0 1.5])
    title('Photodiode Signal as Step-function: onsets and offsets'); 
    subplot(2,2,4)
    plot(time,ergsig2Dfilt_corr)
    hold on
    plot(time,onsets_ergo2{1,ergcnt},'or','MarkerFaceColor','r');
    hold on
    plot(time,offsets_ergo2{1,ergcnt},'og','MarkerFaceColor','g');
    title('Original photodiode signal (detrended+filtered) with onsets (red) and offsets (green) (ERG2)');

end

%% CATEGORIZE THE FEEDBACK ONSET TIMES BASED ON DURATIONS 
triginfo_file = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-PEPs','Brain_IHM_Lists.xlsx');  %Path to excel file with trigger info. 

onsets_all = cell(1,size(ergsig1,2));
fbtypes_all = cell(1,size(ergsig1,2));
lats_all = cell(1,size(ergsig1,2));
onset_data = cell(1,size(ergsig1,2));
Allergo_onsets = cell(1,size(ergsig1,2));
bloc_num = length(filenom);

trignames = {'erg1-cong' 'erg2-cong' 'erg1-incong' 'erg2-incong'};

for bcnt = 1:bloc_num
    
    % Read in information on Lists, corresponding videos and trigger-types
    % from xlsx file. 
    ix = strfind(filenom{1,bcnt},'_');
    sheetnom = filenom{1,bcnt}(ix(1)+1:ix(2)-1);
    XIn = readtable(triginfo_file,'FileType','spreadsheet', 'sheet',sheetnom, 'ReadVariableNames',1, 'ReadRowNames',1,...
        'TreatAsEmpty','NA', 'Range','A1:F5');
    display(XIn);   %print to screen a table showing the properties of current list. 
    fprintf('The current film is:\t%s\n', char(XIn{bcnt,1}));   %print to screen the current film 
    
    %% Detect the current triggers according to the video-type presented. 
    alltrigs = XIn{bcnt,2:5};
    itrig = ~isnan(XIn{bcnt,2:5});
    trigs_curr = alltrigs(itrig);
    trignom_curr = trignames(itrig);
      
    durs_erg1onset = onoffsets_ergo1{1,bcnt}(:,3);  %durations of photodiode signal for ergo1
    durs_erg2onset = onoffsets_ergo2{1,bcnt}(:,3);  %durations of photodiode signal for ergo2

    %ERGO1 (visual) to begin with...
    if length(trigs_curr)<=2
        trig_erg1incong = 400;
        trig_erg2incong = 320;
    else
        trig_erg1incong = trigs_curr(3);
        trig_erg2incong = trigs_curr(4);
    end
    
    erg1incong_diff = arrayfun(@(eg1) abs(eg1-trig_erg1incong), durs_erg1onset);
    erg1cong_diff = arrayfun(@(eg1) abs(eg1-trigs_curr(1)), durs_erg1onset);
    [erg1min, ierg1] = min([erg1cong_diff, erg1incong_diff],[],2);    % 1=congruent; 2=incongruent
    

    %ERGO2 (auditory)...
    erg2incong_diff = arrayfun(@(eg2) abs(eg2-trig_erg2incong), durs_erg2onset);
    erg2cong_diff = arrayfun(@(eg2) abs(eg2-trigs_curr(2)), durs_erg2onset);
    [erg2min, ierg2] = min([erg2cong_diff, erg2incong_diff], [],2);   % 1=congruent; 2=incongruent 
    
    % Determine the latencies for the ERGO1 and ERGO2 onsets
    time = Times_bloc{bcnt,1};
    erg1onsets = onoffsets_ergo1{:,bcnt}(:,1);
    erg1onsets_cong = erg1onsets(ierg1==1);   %ergo1 congruent onset times
    erg1onsets_incong = erg1onsets(ierg1==2); %ergo1 incongruent onset times
    
    erg2onsets = onoffsets_ergo2{:,bcnt}(:,1);
    erg2onsets_cong = erg2onsets(ierg2==1);   %ergo2 congruent onset times
    erg2onsets_incong = erg2onsets(ierg2==2); %ergo2 incongruent onset times
    
    if isempty(erg1onsets_incong) || isempty(erg2onsets_incong)
        
        dumerg1_cong = ones(length(erg1onsets_cong),1);
        dumerg2_cong = ones(length(erg2onsets_cong),1).*2;
        Allergs = [erg1onsets_cong; erg2onsets_cong];
        Alldums = [dumerg1_cong; dumerg2_cong];
        dum_erg = [Allergs, Alldums];
       [onsets_all{1,bcnt},idx] = sort(dum_erg(:,1),'ascend');
       idum_all = Alldums(idx,:);  %1=vis_cong, 2=aud_cong
    else
        
        dumerg1_cong = ones(length(erg1onsets_cong),1);
        dumerg2_cong = ones(length(erg2onsets_cong),1).*2;
        dumerg1_incong = ones(length(erg1onsets_incong),1).*3;
        dumerg2_incong = ones(length(erg2onsets_incong),1).*4;
        Allergs = [erg1onsets_cong; erg2onsets_cong; erg1onsets_incong; erg2onsets_incong];
        Alldums = [dumerg1_cong; dumerg2_cong; dumerg1_incong; dumerg2_incong];
        dum_erg = [Allergs, Alldums];
       [onsets_all{1,bcnt},idx] = sort(dum_erg(:,1),'ascend');
       idum_all = Alldums(idx,:);  %1=vis_cong, 2=aud_cong, 3=vis_incong, 4=aud_incong
        
        
    end 
    
    
    tdiff = arrayfun(@(tin) abs(tin-onsets_all{1,bcnt}),time, 'UniformOutput', false);
    idx_t = cellfun(@(l) sum(l==0), tdiff);
    Allergo_onsets = time(idx_t==1);
    fbtype = cell(length(Allergo_onsets),1);
    lats = find(idx_t);
    
    %Assign trigger names to the onset times 
    for i = 1:length(Allergo_onsets)
        
            fbtype{i,1}= trignames{1,idum_all(i)};
        
    end
    
    fbtypes_all{1,bcnt} = fbtype;
    lats_all{1,bcnt} = Allergo_onsets;
    
    onset_data{1,bcnt}.video = string(fb_sessions{1,bcnt});
    onset_data{1,bcnt}.onset_times = onsets_all{1,bcnt}';
    onset_data{1,bcnt}.latencies = lats_all{1,bcnt}';
    onset_data{1,bcnt}.types = fbtypes_all{1,bcnt}';
    
    for cntr =1:length(lats_all{1,bcnt})
        ALLEEG(bcnt).event(cntr).type = string(fbtypes_all{1,bcnt}(cntr));
        ALLEEG(bcnt).event(cntr).latency = lats(cntr);
        ALLEEG(bcnt).event(cntr).urevent = cntr;
        ALLEEG(bcnt).urevent(cntr).type = string(fbtypes_all{1,bcnt}(cntr));
        ALLEEG(bcnt).urevent(cntr).latency = lats(cntr);
    end
    
    [ALLEEG, ~, CURRENTSET] = pop_newset(ALLEEG, ALLEEG(bcnt), CURRENTSET,'setname',char(fb_sessions{1,bcnt}),'gui','off');
    [ALLEEG, EEG] = eeg_store(ALLEEG, ALLEEG(bcnt), CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fb_sessions{1,bcnt}),'filepath',Dirbase);
    eeglab redraw  
end


eeglab redraw


