[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset(); % Choose *.set file manually. 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% Load in the audio file corresponding to the current EEG dataset.

audnom = EEG.video_name;               % Title of corresponding video.
auddir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',filesep);
audfiles = dir(strcat(auddir,'*.wav'));
findwav=find(~[audfiles.isdir]);        
Allwavs= {audfiles(findwav).name}; 
X = contains(string(Allwavs),audnom);  % Search auddir for the file with name matching that of audnom.

txtfiles = dir(strcat(auddir,'*.txt'));
findtxt=find(~[txtfiles.isdir]);        
Alltxt= {txtfiles(findtxt).name}; 
Xtxt = contains(string(Alltxt),audnom);

diraud_full = fullfile(auddir,Allwavs{X});
dirtxt_full = fullfile(auddir,Alltxt{Xtxt});

% Read in the *.wav file (audio file)
[Y, sF] = audioread(diraud_full);
info = audioinfo(diraud_full);
time = 0:1/sF:(1/sF)*length(Y);    % Create the time vector.

% Need to resample the audio vector to match the EEG signal
Fs_new = EEG.srate;
[Numer, Denom] = rat(Fs_new/sF);
alldat_rs = resample(Y, Numer, Denom);
tnew = 0:1/Fs_new:(1/Fs_new)*length(alldat_rs);   %New time vector
time_rs = tnew(1:end-1);
alldat_rs = mean(alldat_rs,2);

% Read in the *.txt file (with onset info for audio file);
% First column gives onset times of triggers in audio file.
audtxtIn = readtable(dirtxt_full);

% Onset times of the triggers in EEG.
eegtrig_onsets = [EEG.event.latency]./EEG.srate;

%% Need to equalize the length of the time vectors of the audio file and EEG dataset.
% Load in the *.xls file containing the trigger onset information for both
% EEG and auditory data.

eegtime = EEG.times./1000;  %Express in seconds
audtime = time_rs;          %Time vector of the audio signal after resampling.

fprintf('The EEG dataset has length %4.2fsecs.\n The auditory file has length %4.2fsecs\n',eegtime(end)/1000,audtime(end));

filmz = {'Film1', 'Film2', 'Film3', 'Film4'};
filmcols = {'A5:D60';'E5:H60';'I5:L60';'M5:P60'};
etitlez = {EEG.setname, EEG.setname, EEG.setname, EEG.setname};
X1 = cell2mat(cellfun(@contains, etitlez, filmz, 'UniformOutput',false));
filmcurr = filmz{X1};  % The current film

onsetinfo = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Onset-comparisons','VerbOnset-compare.xlsx');
Onsetdata = readtable(onsetinfo,'Range',filmcols{X1},'Sheet',strcat(lower(EEG.setname(1)),EEG.setname(2:3)));

% Calculate the padding that should be added to the beginning of either the audio or the EEG
% signal. 
% Re-calculate the new trigger onset times.
% If the EEG signal is longer than the auditory signal.

toadd = zeros(ceil((eegtrig_onsets(1) - audtxtIn{:,1}(1))*EEG.srate),1);
audtrigs = audtxtIn{:,1}+(eegtrig_onsets(1) - audtxtIn{:,1}(1));
Ynew1 = cat(1,toadd,alldat_rs);
audtime_new1 = 0:1/Fs_new:(1/Fs_new)*length(Ynew1);   

% Calculate the padding that should be added to the end of the shorter
% signal.
toadd_end = zeros(ceil((eegtime(end)-audtime_new1(end))*EEG.srate)+1,1);
Ynew2 = cat(1,Ynew1,toadd_end);
audtime_new2 = 0:1/Fs_new:(1/Fs_new)*(length(Ynew2)-1);  % The length-corrected audio signal.

%% Mark the trigger time points on the audio time vector 

trigindx = zeros(length(audtrigs),1);
trigtime = nan(length(audtime_new2),1);

for fcnt = 1:length(trigindx)  
    trigindx(fcnt) = dsearchn(audtime_new2',audtrigs(fcnt));
    trigtime(trigindx(fcnt)) = .01;
    
end

Etrigindx = zeros(length(eegtrig_onsets),1);
Etrigtime = nan(length(audtime_new2),1);

for fcnt1 = 1:length(Etrigindx)  
    Etrigindx(fcnt1) = dsearchn(audtime_new2',eegtrig_onsets(fcnt1));
    Etrigtime(Etrigindx(fcnt1)) = 20;
    
end

%% Call of function to calculate the envelope of audio signal.

bbandenv = PEPs_EnvelopeCalc(Ynew2, Fs_new, audtime_new2, audnom, trigtime);

plotdata = [EEG.data(48,:); bbandenv]';

figure
% Create a placeholder for axes objects
ax = gobjects( 2, 1 );

ax(1) = subplot(2,1,1)
plot(audtime_new2,EEG.data(48,:))
hold on
stem(audtime_new2,Etrigtime,'r')
ax(2) = subplot(2,1,2)
plot(audtime_new2, bbandenv)
hold on
stem(audtime_new2,trigtime,'r')

dx=10;
CREx_scrollplot(dx, audtime_new2, ax)
