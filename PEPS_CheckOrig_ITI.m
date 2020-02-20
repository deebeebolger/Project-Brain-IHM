% Date: 31-01-2020      Programmed by: D. Bolger
% function to take in feedback information from excel file and order the
% information into trial number based on calculated ITI.
% This information can be compared with the trial numbers assigned for each
% video for each subject to asses if the feedbacks are in the correct
% order.
% The calculated ITI can also serve to check the feedbacks calculated
% automatically from the photodiode signal. This is important as the
% photodiode signal is not always very clean.
%**************************************************************************

xlfile_in = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-PEPs','Brain-IHM-feedback-trigger-check.xlsx');
sheetnom = 'Curare-Agent-Incongruent';
XIn = readtable(xlfile_in,'FileType','spreadsheet', 'sheet',sheetnom, 'ReadVariableNames',1, 'ReadRowNames',0,...
        'TreatAsEmpty','NA', 'Range','C1:E65');
display(XIn);
 
if size(XIn,2)>2                % For the incongruent videos
    fbtypes = cellfun(@(str1,str2) strcat(str1,'-',str2), XIn{1:end,1}, XIn{1:end,2},'UniformOutput',false);
    T= XIn{1:end,3};     %extract the times in seconds. 
elseif size(XIn,2)==2           % For the congruent videos
    fbtypes = XIn{1:end,1};
    T= XIn{1:end,2};     %extract the times in seconds. 
end
[tsort, tint] = sort(T,'ascend');  %sort the times in ascending order
fbtypes_sort = fbtypes(tint,1);    %sort the feedback types.

%% Assign a trial number to each feedback instance on the basis of ITI. 

min_dist = 0.1;   %define the minimum distance between two feedbacks. Below this value, the feedbacks are considered as being in the same trial.
diff_time = diff(tsort);
trl_dummy = zeros(length(diff_time)+1,1);
trlcounter = 1;

for trlcnt = 1:length(tsort)

    if trlcnt==1
        t_end=trlcounter+1;
    end

    if t_end<=length(T)

        diff_curr = tsort(t_end)-tsort(trlcounter);

        if diff_curr<=min_dist
            trl_dummy(trlcounter:t_end)=trlcnt;
            trlcounter = t_end+1;
            t_end=t_end+2;

        elseif diff_curr> min_dist
            trl_dummy(trlcounter)=trlcnt;
            trlcounter = trlcounter+1;
            t_end=t_end+1;
        end

    elseif t_end>length(T) && diff_curr<=min_dist
        
            trl_dummy(end) = trl_dummy(end-1);
            
    elseif t_end>length(T) && diff_curr> min_dist

            trl_dummy(end)=trl_dummy(end-1)+1; 
    end

end

Trial_numbers = trl_dummy;
Feedbacks = fbtypes_sort;
Times = tsort;
ITI = cat(1,diff_time,nan);

%% CREATE A TABLE OF THE TIMES, FEEDBACKS AND TRIAL NUMBERS

fbtable_orig = table(Times, Feedbacks, Trial_numbers, ITI);

savepath = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-PEPs',strcat(sheetnom,'.mat'));

save(savepath,'fbtable_orig');
