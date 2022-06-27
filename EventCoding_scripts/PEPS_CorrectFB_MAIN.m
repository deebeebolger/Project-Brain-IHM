
[fbtable_curr, icorr] = PEPS_CorrectFB(fbtable_orig, fb_table{1,4}, EEG.filepath, EEG.setname);

%% 
writetable(fbtable_curr,fullfile(EEG.filepath,[EEG.setname,'.txt']),'Delimiter',' ');   % Write the corrected feedback information to table and save.
writetable(fbtable_orig(fbtable_curr.incorrectindx(1:icorr),:),fullfile(EEG.filepath,[EEG.setname,'-missing.txt']),'Delimiter',' '); 