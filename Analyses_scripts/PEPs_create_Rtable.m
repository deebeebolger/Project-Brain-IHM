% Script to prepare table of EEG data for R long-form.
% Table structure:
% Channel  Subject  Group  Condition T1 T2 T3 T4...Tn

load_dir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021','Segmented','matfiles',filesep);
load_data = dir(load_dir);
dir_fnames = {load_data.name};
fois = dir_fnames(~[load_data.isdir]);
matfiles2load = fois(contains(fois,'BLC.mat'));

Groups = {'Human' 'Agent'};
Conditions = {'Congruent' 'Congruent-Incong' 'Incongruent'};
chanoi = {'Fz' 'FCz' 'Cz' 'CPz' 'Pz' 'F3' 'F4' 'FC3' 'FC4' 'C3' 'C4' 'CP3' 'CP4' 'P3' 'P4'};


% Define the table.
fid = fopen(fullfile(load_dir,'myTextFile_allsuj_v2.txt'),'a+');

%% 
for mcount = 3:length(matfiles2load)
    
    mcount
    s = matfiles2load{1,mcount}(1:3);
    curr_sub = load(fullfile(load_dir,matfiles2load{1,mcount}));
    flds = fieldnames(curr_sub);
    dataIn = curr_sub.(flds{1,1});

    % Define the indices of groups and conditions
    group1_indx = contains(curr_sub.video_name,Groups{1,1});
    group2_indx = contains(curr_sub.video_name,Groups{1,2});
    cond1_indx = contains(curr_sub.video_name,Conditions{1,1});     % Find data from congruent videos only
    cond2_indx = contains(curr_sub.video_name,Conditions{1,3});     % Find data from incongruent videos.
    cong_incong_indx = and(cond2_indx,[curr_sub.trialnums>10]);     % Congruent conditions from incongruent videos only
    incong_incong_indx = and(cond2_indx,[curr_sub.trialnums<=10]);  %Incongruent conditions only.
    allcong_indx = any([cond1_indx; cong_incong_indx]);              % All congruent conditions (incong and cong videos);
    
    group1_cond1 = all([group1_indx; cond1_indx]);                  %Human congruent from congruent videos only
    group2_cond1 = all([group2_indx; cond1_indx]);                  %Agent congruent from congruent videos only
    group1_cond2 = all([group1_indx; incong_incong_indx]);          %Human incongruent
    group2_cond2 = all([group2_indx; incong_incong_indx]);          %Agent incongruent
    group1_cond1_incong = all([group1_indx; cong_incong_indx]);     %Human congruent in incongruent videos
    group2_cond1_incong = all([group2_indx; cong_incong_indx]);     %Agent congruent in incongruent videos
    
    dataIndx = cell(1,2);
    dataIndx{1,1} = [group1_cond1;group1_cond1_incong;group1_cond2];
    dataIndx{1,2} = [group2_cond1;group2_cond1_incong;group2_cond2];
    
    % Find the grand average for each of the datasets.
    currsub_data = zeros(64,length(curr_sub.times),length(curr_sub.trialnums));
    tstart = 1;
    tend = 0;
    
    for scount = 1:size(currsub_data,3)
        
        tend = tend + curr_sub.trialnums(scount);
        currsub_data(:,:,scount) = mean(curr_sub.data(1:64,:,tstart:tend),3);
        tstart = tend+1;
        
    end
    
    % Create matrix of data for the current chanoi.
    
   
    
    for chancnt = 1:length(chanoi)
        
        chanidx = find(strcmp({curr_sub.chanlocs.labels},chanoi{1,chancnt}));
        
        chancurr = [currsub_data(chanidx,:,group1_cond1);currsub_data(chanidx,:,group1_cond1_incong);currsub_data(chanidx,:,find(group1_cond2));...
            currsub_data(chanidx,:,group2_cond1);currsub_data(chanidx,:,group2_cond1_incong); currsub_data(chanidx,:,group2_cond2)];
        
        
        for gcnt = length(Groups):-1:1
            for ccnt = 1:length(Conditions)
                    fprintf(fid,'%s\t %s\t %s\t %s\t %.4f\t',s,chanoi{1,chancnt}, Groups{1,gcnt}, Conditions{1,ccnt},chancurr(dataIndx{1,gcnt}(ccnt,:),:));
                    fprintf(fid,'\n');
            end % end of Condition loop.
        end %end of group loop
    end %end of channel loop
    
    if mcount==1
        fid1 = fopen(fullfile(load_dir,'myTimeFile.txt'),'w');
        fprintf(fid1,'%.4f\n',curr_sub.times);
        fclose(fid1);
        
    end
    
    
    
end




