function EEGIn = PEPS_ITIcalc(EEGIn)

min_dist = 0.1;   %define the minimum distance between two feedbacks. Below this value, the feedbacks are considered as being in the same trial.
diff_lat = diff([EEGIn.event.latency]);
difflat_sec = diff_lat./EEGIn.srate;
trl_dummy = zeros(length(difflat_sec)+1,1);
trlcounter = 1;

for trlcnt = 1:length(EEGIn.event)
    
    if trlcnt==1
        t_end=trlcounter+1;
    end
    
    if t_end<=length(EEGIn.event)
        
        diff_curr = (EEGIn.event(t_end).latency-EEGIn.event(trlcounter).latency)/EEGIn.srate;
        
        if diff_curr<=min_dist
            trl_dummy(trlcounter:t_end)=trlcnt;
            trlcounter = t_end+1;
            t_end=t_end+2;
            
        elseif diff_curr> min_dist
            trl_dummy(trlcounter)=trlcnt;
            trlcounter = trlcounter+1;
            t_end=t_end+1;
        end
        
    elseif t_end>length(EEGIn.event) && diff_curr<=min_dist
        
        trl_dummy(end) = trl_dummy(end-1);
        
    elseif t_end>length(EEGIn.event) && diff_curr> min_dist
        
        trl_dummy(end)=trl_dummy(end-1)+1;
    end
    
    trl_dummy(end)=trl_dummy(end-1)+1;
    
end

%% Add the trial number information to the event structure of the EEG structure.

for icnt = 1:length([EEGIn.event])
    
    EEGIn.event(icnt).trialnum = trl_dummy(icnt);
    
end


%% SAVE ALL THE INFORMATION INTO A MAT FILE AS A TABLE FOR EACH FILM

Trial_numbers = trl_dummy;
Feedbacks = {EEGIn.event.type}';
Times = [EEGIn.event.latency]'./EEGIn.srate;
ITI = cat(1,difflat_sec',nan);

fb_table = table(Times, Feedbacks, Trial_numbers, ITI);

save(fullfile(EEGIn.filepath,[EEGIn.setname,'-feedback_summary.mat']),'fb_table');

eeglab redraw

end