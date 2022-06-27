function [Mean1Out, Mean2Out,  StdDurs1Out,  StdDurs2Out] = PEPs_fbstat_calc(fnoms1, fnoms2, file_path, condtype)

D = dir(fullfile(file_path,'*.txt'));
fileall = {D.name};

if sum(contains(fnoms1,'Congruent'))>0
    
    FeedBacks1 = cell(length(fnoms1),1);
    Durs1 = cell(length(fnoms1),1);
    
    for fcount1 = 1:length(fnoms1)
        
        indx = contains(fileall,fnoms1{fcount1,1});
        DataIn = readtable(fullfile(file_path,fileall{1,indx}));
        indx_excl = ~contains(DataIn{:,3},'-lexicalise');
        FeedBacks1{fcount1,1} = DataIn{indx_excl,3};
        Durs1{fcount1,1} = DataIn{indx_excl,2}-DataIn{indx_excl,1};
        
    end
    
    AllFBks1 = cat(1,FeedBacks1{1,:},FeedBacks1{2,:},FeedBacks1{3,:});
    AllDurs1 = cat(1,Durs1{1,:}, Durs1{2,:}, Durs1{3,:});
    MeanDurs1_cong = mean(AllDurs1);
    StdDurs1_cong = std(AllDurs1);
    Mean1Out = MeanDurs1_cong;
    StdDurs1Out = StdDurs1_cong;
    
elseif sum(contains(fnoms1,'Incongruent'))>0
    
    FeedBacks1_cong = cell(length(fnoms1),1);
    Durs1_cong = cell(length(fnoms1),1);
    FeedBacks1_incong = cell(length(fnoms1),1);
    Durs1_incong = cell(length(fnoms1),1);
    
    for fcount1 = 1:length(fnoms1)
        
        indx = contains(fileall,fnoms1{fcount1,1});
        DataIn = readtable(fullfile(file_path,fileall{1,indx}));
        indx_excl = ~contains(DataIn{:,3},'-lexicalise');
        FB = DataIn{indx_excl,3};
        FB_onset = DataIn{indx_excl,1};
        FB_offset = DataIn{indx_excl,2};
        
        % Isolate the congruent trials only.
        indx_cong = contains(FB,'Cong');
        Feedbacks1_cong{fcount1,1} = FB(indx_cong,1);
        Durs1_cong{fcount1,1} = FB_offset(indx_cong,1) - FB_onset(indx_cong,1);
        
        % Isolate the incongruent trials only.
        indx_incong = contains(FB,'Incong');
        Feedbacks1_incong{fcount1,1} = FB(indx_incong,1);
        Durs1_incong{fcount1,1} = FB_offset(indx_incong,1) - FB_onset(indx_incong,1);
        
    end
    
    AllFBks1_cong = cat(1,FeedBacks1_cong{1,:},FeedBacks1_cong{2,:},FeedBacks1_cong{3,:});
    AllFBks1_incong = cat(1,FeedBacks1_incong{1,:},FeedBacks1_incong{2,:},FeedBacks1_incong{3,:});
    AllDurs1_cong = cat(1,Durs1_cong{1,:}, Durs1_cong{2,:}, Durs1_cong{3,:});
    AllDurs1_incong = cat(1,Durs1_incong{1,:}, Durs1_incong{2,:}, Durs1_incong{3,:});
    MeanDurs1_cong = mean(AllDurs1_cong);
    MeanDurs1_incong = mean(AllDurs1_incong);
    StdDurs1_cong = std(AllDurs1_cong);
    StdDurs1_incong = std(AllDurs1_incong);
    
    if strcmpi(condtype{1,1},'-incong')
        Mean1Out = MeanDurs1_incong;
        StdDurs1Out = StdDurs1_incong;
    elseif strcmpi(condtype{1,1},'-cong')
        Mean1Out = MeanDurs1_cong;
        StdDurs1Out = StdDurs1_cong;
    end
    
end

% *************** Second condition*******************

if sum(contains(fnoms2,'Congruent'))>0
    
    FeedBacks2 = cell(length(fnoms2),1);
    Durs2 = cell(length(fnoms2),1);
    
    for fcount2 = 1:length(fnoms2)
        
        indx = contains(fileall,fnoms2{fcount2,1});
        DataIn = readtable(fullfile(file_path,fileall{1,indx}));
        indx_excl = ~contains(DataIn{:,3},'-lexicalise');
        FeedBacks2{fcount2,1} = DataIn{indx_excl,3};
        Durs2{fcount2,1} = DataIn{indx_excl,2}-DataIn{indx_excl,1};
        
    end
    
    AllFBks2 = cat(1,FeedBacks2{1,:},FeedBacks2{2,:},FeedBacks2{3,:});
    AllDurs2 = cat(1,Durs2{1,:}, Durs2{2,:}, Durs2{3,:});
    MeanDurs2_cong = mean(AllDurs2);
    StdDurs2_cong = std(AllDurs2);
    Mean2Out = MeanDurs2_cong;
    StdDurs2Out = StdDurs2_cong;
    
elseif sum(contains(fnoms2,'Incongruent'))>0  %Needs to be finished
    
    FeedBacks2_cong = cell(length(fnoms2),1);
    Durs2_cong = cell(length(fnoms2),1);
    FeedBacks2_incong = cell(length(fnoms2),1);
    Durs2_incong = cell(length(fnoms2),1);
    
    
    for fcount2 = 1:length(fnoms2)
        
        
        indx = contains(fileall,fnoms2{fcount2,1});
        DataIn = readtable(fullfile(file_path,fileall{1,indx}));
        indx_excl = ~contains(DataIn{:,3},'-lexicalise');
        FB = DataIn{indx_excl,3};
        FB_onset = DataIn{indx_excl,1};
        FB_offset = DataIn{indx_excl,2};
        
        % Isolate the congruent trials only.
        indx_cong = contains(FB,'Cong');
        Feedbacks2_cong{fcount2,1} = FB(indx_cong,1);
        Durs2_cong{fcount2,1} = FB_offset(indx_cong,1) - FB_onset(indx_cong,1);
        
        % Isolate the incongruent trials only.
        indx_incong = contains(FB,'Incong');
        Feedbacks2_incong{fcount2,1} = FB(indx_incong,1);
        Durs2_incong{fcount2,1} = FB_offset(indx_incong,1) - FB_onset(indx_incong,1);
    end
    
    AllFBks2_cong = cat(1,FeedBacks2_cong{1,:},FeedBacks2_cong{2,:},FeedBacks2_cong{3,:});
    AllFBks2_incong = cat(1,FeedBacks2_incong{1,:},FeedBacks2_incong{2,:},FeedBacks2_incong{3,:});
    AllDurs2_cong = cat(1,Durs2_cong{1,:}, Durs2_cong{2,:}, Durs2_cong{3,:});
    AllDurs2_incong = cat(1,Durs2_incong{1,:}, Durs2_incong{2,:}, Durs2_incong{3,:});
    MeanDurs2_cong = mean(AllDurs2_cong);
    MeanDurs2_incong = mean(AllDurs2_incong);
    StdDurs2_cong = std(AllDurs2_cong);
    StdDurs2_incong = std(AllDurs2_incong);
    
    if strcmp(condtype{1,2},'-incong')
        Mean2Out = MeanDurs2_incong;
        StdDurs2Out = StdDurs2_incong;
    elseif strcmp(condtype{1,2},'-cong')
        Mean2Out = MeanDurs2_cong;
        StdDurs2Out = StdDurs2_cong;
    end
    
end




end


