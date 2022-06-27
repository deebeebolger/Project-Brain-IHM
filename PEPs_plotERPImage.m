%% Script PEPs_plotERPImage
% Date: 3/3/2021         Programmed by: Deirdre Bolger
%
%*************************************************************************************

[paramfile_nom, paramfile_path] = uigetfile('*.txt','Select a Parameters text-file');  %Load in the parameters textfile with ERP plotting parameters.

fid = fopen(fullfile(paramfile_path,paramfile_nom));
C = textscan(fid, '%s','Delimiter',',','CommentStyle','//');
fclose(fid);

params = C{1,1};                         %Extract the contents of the parameters file.
line1 = string(split(params{1,1},' '));  %video_names
line2 = string(split(params{2,1},' '));  %video_types
line3 = string(split(params{3,1},' '));  %video_subject
line4 = string(split(params{4,1},' '));  %feedback_type
line5 = split(params{5,1},' ');
line6 = split(params{6,1},' ');
line7 = split(params{7,1},' ');          %Defines the path to the chaninfo mat file.
chaninfo = load(line7{2,1},'chaninfo');   %loads in a 1 X 72 structure
chaninfo = chaninfo.chaninfo;

[indx1,tf1] = listdlg('PromptString',{'Select datasets for condition 1.',...
    'You can select several items at once.',''},...
    'SelectionMode','multiple','ListString',line1(2:end,1));

[indx2,tf2] = listdlg('PromptString',{'Select datasets for condition 2.',...
    'You can select several items at once.',''},...
    'SelectionMode','multiple','ListString',line1(2:end,1));

vidall = line1(2:end-1,1);
vidcond1 = vidall(indx1,1);  % Videos for condition1

if ~strcmp(vidall(indx2,1),'None')
    vidcond2 = vidall(indx2,1);  % Videos for condition2
else
    vidcond2 = [];
end

%% Access the data for the defined conditions

patients = {'Human','Human'};
condcong = {'-cong','-incong'};  % The congruous condition to select for each dataset being compared.
condnames = {[patients{1,1},'-',condcong{1,1}],[patients{1,2},'-',condcong{1,2}]};
fbk_oi = 'Enonce-lexicalise';   % Enonce-lexicalise
fbk_action = 'rej';  % or keep
chanoi = 'Cz';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  %Open EEGLAB
GrandAvg1 = cell(length(vidcond1),1);
indx2exl_v1 = cell(length(vidcond1),1);
indx2exl_v2 = cell(length(vidcond2),1);
AllData_vid1 = cell(length(vidcond1),1);
AllData_vid2 = cell(length(vidcond2),1);
vid1_trialindx = cell(length(vidcond1),1);
vid2_trialindx = cell(length(vidcond2),1);


for vidcntr = 1:length(vidcond1)
    ALLEEG = [];
    pathcurr = fullfile(line6{2,1},vidcond1(vidcntr,1),'BLCorrected',filesep);
    f = dir(fullfile(pathcurr));
    fnoms = {f.name};
    f1 = contains(fnoms,condcong{1,1});
    filenoms_all = {fnoms{f1}};
    filenoms = {filenoms_all{contains(filenoms_all,'.set')}};
    EEG = pop_loadset('filename',filenoms,'filepath',char(pathcurr));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    eeglab redraw;
    
    if vidcntr == 1
        time = ALLEEG(vidcntr).times;
        chanscurr = ALLEEG(vidcntr).chanlocs;
        chanoi_ind = find(ismember({chanscurr.labels},chanoi));
    end
    
    Dall = {ALLEEG(:).data};
    indx2exl_v1 = cell(length(filenoms),1);
    trialnum_vid1 = 0;
    tnum_ds1 = cell(length(filenoms),1);
    
    if ~strcmp(condcong{1,1},'-incong')
        for dscount = 1:length(filenoms)
            trialnum_vid1 = trialnum_vid1 + size(Dall{1,dscount}, 3);
            if isfield(ALLEEG(dscount).event,'feedback')
                if strcmp(fbk_action,'rej')
                    indx2exl_v1{dscount,1} = find(~contains([ALLEEG(dscount).event.feedback],fbk_oi));
                elseif strcmp(fbk_action,'keep')
                    indx2exl_v1{dscount,1} = find(contains([ALLEEG(dscount).event.feedback],fbk_oi));
                end
                if ~isempty(indx2exl_v1{dscount,1})
                    Dall{1,dscount} = squeeze(Dall{1,dscount}(chanoi_ind,:,indx2exl_v1{dscount,1}));
                    
                else
                    Dall{1,dscount} = squeeze(Dall{1,dscount}(chanoi_ind,:,:));
                end
            else
                disp('0Oops! Feedback field does not exist in current dataset.');
                Dall{1,dscount} = squeeze(Dall{1,dscount}(chanoi_ind,:,:));
                indx2exl_v1{dscount,1} = [];
            end
            
        end
    else
        for dscount1 = 1:length(filenoms)
            
            trialnum_vid1 = trialnum_vid1 + size(Dall{1,dscount1},3);
            Dall{1,dscount} = squeeze(Dall{1,dscount}(chanoi_ind,:,:));
               if size(Dall{1,dscount},2)>size(Dall{1,dscount},1)
                    Dall{1,dscount} = squeeze(Dall{1,dscount}(chanoi_ind,:,:))';
                end
            
        end
    end
    
    tnum_ds1 = cell2mat(cellfun(@(x) [1:size(x, 2)], Dall, 'UniformOutput',false));
    AllData_vid1{vidcntr,1} = cell2mat(Dall);
    vid1_trialindx{vidcntr,1} = tnum_ds1;
    
end   % end of vidcntr


if ~isempty(vidcond2)   % if we have defined a second condition dataset
    GrandAvg2 = cell(length(vidcond2),1);
    for vidcntr1 = 1:length(vidcond2)
        ALLEEG = [];
        pathcurr = fullfile(line6{2,1},vidcond2(vidcntr1,1),'BLCorrected',filesep);
        f = dir(fullfile(pathcurr));
        fnoms = {f.name};
        f1 = contains(fnoms,condcong{1,2});
        filenoms_all = {fnoms{f1}};
        filenoms = {filenoms_all{contains(filenoms_all,'.set')}};
        EEG = pop_loadset('filename',filenoms,'filepath',char(pathcurr));
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        eeglab redraw;
        % Need to find the grand average for each dataset here...
        
        Dall2 = {ALLEEG(:).data};
        indx2exl_v2 = cell(length(filenoms),1);
        trialnum_vid2 = 0;
        tnum_ds2 = cell(length(filenoms),1);
        
        
        if ~strcmp(condcong{1,2},'-incong')
            for dscount = 1:length(filenoms)
                if isfield(ALLEEG(dscount).event,'feedback')
                    if strcmp(fbk_action,'rej')
                        indx2exl_v2{dscount,1} = find(~contains([ALLEEG(dscount).event.feedback],fbk_oi));
                    elseif strcmp(fbk_action,'keep')
                        indx2exl_v2{dscount,1} = find(contains([ALLEEG(dscount).event.feedback],fbk_oi));
                    end
                    if ~isempty(indx2exl_v2{dscount,1})
                        Dall2{1,dscount} = squeeze(Dall2{1,dscount}(chanoi_ind,:,indx2exl_v2{dscount,1}));
                    else
                        Dall2{1,dscount} = squeeze(Dall2{1,dscount}(chanoi_ind,:,:));
                    end
                else
                    disp('0Oops! Feedback field does not exist in current dataset.');
                    Dall2{1,dscount} = squeeze(Dall2{1,dscount}(chanoi_ind,:,:));
                    indx2exl_v2{dscount,1} = [];
                end
                trialnum_vid2 = trialnum_vid2 + size(Dall2{1,dscount}, 3);
            end
        else
            
            for dscount1 = 1:length(filenoms)
                
                trialnum_vid2 = trialnum_vid2 + size(Dall2{1,dscount1},3);
                Dall2{1,dscount1} = squeeze(Dall2{1,dscount1}(chanoi_ind,:,:));
                if size(Dall2{1,dscount1},2)>size(Dall2{1,dscount1},1)
                    Dall2{1,dscount1} = squeeze(Dall2{1,dscount1}(chanoi_ind,:,:))';
                end
            end
            
        end
        tnum_ds2 = cell2mat(cellfun(@(x) [1:size(x, 2)], Dall2, 'UniformOutput',false));
        AllData_vid2{vidcntr1,1} = cell2mat(Dall2);
        vid2_trialindx{vidcntr1,1} = tnum_ds2;
        
    end % end of vidcntr1
end %end of if vidcond2

Cond1 = cell2mat(AllData_vid1');
Cond2 = cell2mat(AllData_vid2');
Tnum_cond1 = cell2mat(vid1_trialindx');
Tnum_cond2 = cell2mat(vid2_trialindx');

%% Create the two dataset to be plotted and the "sortvar" variables.

Cond1 = cell2mat(AllData_vid1');
Cond2 = cell2mat(AllData_vid2');
Tnum_cond1 = cell2mat(vid1_trialindx');
Tnum_cond2 = cell2mat(vid2_trialindx');
times = ALLEEG(1).times;
k = dsearchn(times', [150 250]');   %Find indices of time samples between defined time limits.

% Find peak amplitude in the defined time interval
pkamp_cond1 = max(abs(Cond1(k(1):k(2),:)));
pkamp_cond2 = max(abs(Cond2(k(1):k(2),:)));


figure;
erpimage(Cond1,Tnum_cond1,times, [condnames{1,1},':',chanoi], 20, 1, 'erp','cbar','vert',0,...
    'coher',[10 12 .01],'srate',ALLEEG(1).srate,'cycles',3);

figure;
erpimage(Cond2,Tnum_cond2,times, [condnames{1,2},':',chanoi], 10, 1, 'erp','cbar','vert',0,...
    'coher',[4 7 .01],'srate',ALLEEG(1).srate,'cycles',3);
