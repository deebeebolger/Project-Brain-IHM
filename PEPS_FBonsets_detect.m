% Date: 29-11-2019             Programmed by: D. Bolger
% Script to detect the onset times of the feedbacks of interest and to
% determine the modality and type of each feedback.
% Script applied in the Brain-IHM project.
% ********************************************************
%% LOAD IN THE *.bdf FILE OF THE CURRENT SUBJECT
% Loads in the *.bdf file, ensuring that 74 channels are included so that the 
% ERGO1 and ERGO2 data are loaded.

Dirbase = fullfile(filesep,'Volumes','KINGSTON','PEPS_Sujets','s01',filesep);

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
fb_sessions = {'List1a-Film1' 'List1a-Film3' 'List1a-Film4'};  
thresh_val = 0;
fs = EEG.srate;
[b,a] = butter(2,8./(fs/2));
onoffsets_ergo1 = cell(1,size(ergsig1,2));
onoffsets_ergo2 = cell(1,size(ergsig1,2));
onsets_ergo1 = cell(1,size(ergsig1,2));
onsets_ergo2 = cell(1,size(ergsig1,2));
offsets_ergo1 = cell(1,size(ergsig1,2));
offsets_ergo2 = cell(1,size(ergsig1,2));

for ergcnt = 1:size(ergsig1,2)
    ergcnt
    time = ALLEEG(ergcnt).times;
    
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
    title('Original photodiode signal (detrended+filtered) with onsets (red) and offsets (green)');
    
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
    title('Original photodiode signal (detrended+filtered) with onsets (red) and offsets (green)');

end

%% CATEGORIZE THE FEEDBACK ONSET TIMES BASED ON DURATIONS 

onsets_all = cell(1,size(ergsig1,2));
fbtypes_all = cell(1,size(ergsig1,2));
lats_all = cell(1,size(ergsig1,2));
onset_data = cell(1,size(ergsig1,2));
sujnum = length(ALLEEG);


fbonset_vistypes = {'vis-cong' 'vis-incong'};
fbonset_audtypes = {'aud-cong' 'aud-incong'};
fbonset_visdurs = [160 320]; 
fbonset_auddurs = [200 400]; % duration of corresponding photodiode signal in ms
 % corresponds to the 3 columns of onoffsets_aud/vis cell arrays.


for sujcnt = 1:sujnum

    dur_aud = onoffsets_ergo1{1,sujcnt}(:,3);
    Iaud = arrayfun(@(x) nearest(fbonset_auddurs,x),dur_aud);
    fba = fbonset_audtypes(1,Iaud);

    dur_vis = onoffsets_ergo2{1,sujcnt}(:,3);
    Ivis = arrayfun(@(x) nearest(fbonset_visdurs,x),dur_vis);
    fbv = fbonset_vistypes(1,Ivis);
    
    fb_all = cat(2,fba,fbv);
    [onsets_all{1,sujcnt},idx] = sort([onoffsets_ergo1{1,sujcnt}(:,1)',onoffsets_ergo2{1,sujcnt}(:,1)'],'ascend');
    lats_all{1,sujcnt} = arrayfun(@(tin) find(time == tin),onsets_all{1,sujcnt});
    fbtypes_all{1,sujcnt} = fb_all(1,idx);
    
    onset_data{1,sujcnt}.video = string(fb_sessions{1,sujcnt});
    onset_data{1,sujcnt}.onset_times = onsets_all{1,sujcnt}';
    onset_data{1,sujcnt}.latencies = lats_all{1,sujcnt}';
    onset_data{1,sujcnt}.types = fbtypes_all{1,sujcnt}';
    
    for cntr =1:length(lats_all{1,sujcnt})
        ALLEEG(sujcnt).event(cntr).type = string(fbtypes_all{1,sujcnt}(cntr));
        ALLEEG(sujcnt).event(cntr).latency = lats_all{1,sujcnt}(cntr);
        ALLEEG(sujcnt).event(cntr).urevent = cntr;
        ALLEEG(sujcnt).urevent(cntr).type = string(fbtypes_all{1,sujcnt}(cntr));
        ALLEEG(sujcnt).urevent(cntr).latency = lats_all{1,sujcnt}(cntr);
    end
    
    [ALLEEG, ~, CURRENTSET] = pop_newset(ALLEEG, ALLEEG(sujcnt), CURRENTSET,'setname',char(fb_sessions{1,sujcnt}),'gui','off');
    [ALLEEG, EEG] = eeg_store(ALLEEG, ALLEEG(sujcnt), CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fb_sessions{1,sujcnt}),'filepath',Dirbase);
    eeglab redraw  
end





