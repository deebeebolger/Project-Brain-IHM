EEGOut = PEPS_ITIcalc(EEG);

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEGOut, CURRENTSET,'setname',EEG.setname,'gui','off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',EEG.filepath);
eeglab redraw