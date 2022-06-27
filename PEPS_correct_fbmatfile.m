%% Script to correct the Feedback *.mat files after verification of the trigger onsets.
%  For the PEPs project.
% November 6 2020
%**************************************
close all;

sujcurr = 's01';
vidnum = {'Film1','Film2','Film3'};
range1 = {'A3:D100','A3:D88','A3:D66'};
range2 = {'J3:J100','J3:K88','K3:J66'};

feedbacks = cell(1,length(vidnum));

for vidcnt = 1:length(vidnum)
     
    FBin = readtable(['/Users/bolger/Brain-IHM/Data_Preproc/',sujcurr,'_Trigger_check_Nov2020.xlsx'],...
        'Range',range1{1,vidcnt}, 'Sheet',vidnum{1,vidcnt});
    offset_diff = readtable(['/Users/bolger/Brain-IHM/Data_Preproc/',sujcurr,'_Trigger_check_Nov2020.xlsx'],...
        'Range',range2{1,vidcnt}, 'Sheet',vidnum{1,vidcnt});    % Difference in offset between predefined onset trigger times and those extracted automatically.
    
    
    dirgen = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','PEPs_DataPreproc_2021',sujcurr);
    dirIn = fullfile(dirgen,strcat(sujcurr,'-',vidnum{1,vidcnt}),filesep);
    fbinfo = load(fullfile(dirgen,'feedback_summary.mat'));

    
    sz = [size(FBin,1) 5];
    varTypes = {'double','string','double', 'double','double'};
    FBtable = table('Size',sz,'VariableTypes',varTypes); 
    
    FBtable{:,1} = FBin{:,1};
    typecell = cell(sz(1),1);
    
    for c = 1:sz(1)
    
        typecell{c,1} = FBin{c,2};
    end
    FBtable{:,2} = typecell;
    if strcmp(vidnum{1,vidcnt},'Film12')
        disp('************Weird one here*********************')
        FBtable{:,3} = FBin{:,4};
        FBtable{:,4} = FBin{:,5};
    else
        FBtable{:,3} = FBin{:,3};
        FBtable{:,4} = FBin{:,4};
    end
    
    FBtable{:,5} = offset_diff{:,1};
   
    feedbacks{1,vidcnt} = FBtable;

end

save(fullfile(dirgen,'feedback_summary_correct.mat'),'feedbacks')
