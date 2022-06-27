%% Script PEPs_plotERP
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
condcong = {'-Cong','-Incong'};  % The congruous condition to select for each dataset being compared.
condcongFB = {'-cong','-incong'};
condnames = {[patients{1,1},'-',condcong{1,1}],[patients{1,2},'-',condcong{1,2}]};
fbk_oi = 'Doctor';   % Enonce-lexicalise
fbk_action = 'rej';  % or keep

% Call of function to calculate the mean and standard-deviation of feedback
% durations for the videos in question.

fpath = fullfile(filesep,'Users','bolger','Brain-IHM','Audio-Files',filesep);
[meandur1, meandur2,  std_durs1,  std_durs2] = PEPs_fbstat_calc(vidcond1, vidcond2, fpath, condcongFB);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  %Open EEGLAB
GrandAvg1 = cell(length(vidcond1),1);
indx2exl_v1 = cell(length(vidcond1),1);
indx2exl_v2 = cell(length(vidcond1),1);
vid1_trialindx = cell(length(vidcond1),1);
vid2_trialindx = cell(length(vidcond2),1);


for vidcntr = 1:length(vidcond1)
    ALLEEG = [];
    pathcurr = fullfile(line6{2,1},vidcond1(vidcntr,1),'BLCorrected',filesep);
    f = dir(fullfile(pathcurr));
    fnoms = {f.name};
    f1 = contains(fnoms,condcong{1,1});
    filenoms_all = {fnoms{f1}};
    filenoms = {filenoms_all{contains(filenoms_all,'clean.set')}};
    EEG = pop_loadset('filename',filenoms,'filepath',char(pathcurr));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    eeglab redraw;

    % Need to find the grand average for each dataset here...
    if vidcntr == 1
        time = ALLEEG(vidcntr).times;
        chanscurr = ALLEEG(vidcntr).chanlocs;
    end

    Dall = {ALLEEG(:).data};
    indx2exl_v1 = cell(length(filenoms),1);

    if ~strcmpi(condcong{1,1},condcong{1,2})
        for dscount = 1:length(filenoms)
            if isfield(ALLEEG(dscount).event,'feedback')
                if strcmp(fbk_action,'rej')
                    indx2exl_v1{dscount,1} = find(~contains([ALLEEG(dscount).event.feedback],fbk_oi));
                elseif strcmp(fbk_action,'keep')
                    indx2exl_v1{dscount,1} = find(contains([ALLEEG(dscount).event.feedback],fbk_oi));
                end
                if ~isempty(indx2exl_v1{dscount,1})
                    Dall{1,dscount} = Dall{1,dscount}(:,:,indx2exl_v1{dscount,1});
                else
                    Dall{1,dscount} = Dall{1,dscount};
                end
            else
                disp('0Oops! Feedback field does not exist in current dataset.');
                Dall{1,dscount} = Dall{1,dscount};
                indx2exl_v1{dscount,1} = [];
            end

        end
    end

    tnum_ds1 = cell2mat(cellfun(@(x) [1:size(x, 3)], Dall, 'UniformOutput',false));
    vid1_trialindx{vidcntr,1} = tnum_ds1;
    Dmean_curr = cell2mat(cellfun(@(x) mean(x,3),Dall,'UniformOutput',false));
    Dmean_curr = reshape(Dmean_curr,[size(Dall{1,1},1),size(Dall{1,1},2),length(ALLEEG)]);
    GrandAvg1{vidcntr,1} = Dmean_curr;
    if vidcntr==1
        DMCond1 = Dmean_curr;
    elseif vidcntr>1
        DMCond1 = cat(3,DMCond1,Dmean_curr);
    end
end


if ~isempty(vidcond2)   % if we have defined a second condition dataset

    GrandAvg2 = cell(length(vidcond2),1);


    for vidcntr1 = 1:length(vidcond2)
        ALLEEG = [];
        pathcurr = fullfile(line6{2,1},vidcond2(vidcntr1,1),'BLCorrected',filesep);
        f = dir(fullfile(pathcurr));
        fnoms = {f.name};
        f1 = contains(fnoms,condcong{1,2});
        filenoms_all = {fnoms{f1}};
        filenoms = {filenoms_all{contains(filenoms_all,'clean.set')}};
        EEG = pop_loadset('filename',filenoms,'filepath',char(pathcurr));
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        eeglab redraw;
        % Need to find the grand average for each dataset here...


        Dall2 = {ALLEEG(:).data};
        indx2exl_v2 = cell(length(filenoms),1);

        if ~strcmpi(condcong{1,2},condcong{1,2})
            for dscount = 1:length(filenoms)
                if isfield(ALLEEG(dscount).event,'feedback')
                    if strcmp(fbk_action,'rej')
                        indx2exl_v2{dscount,1} = find(~contains([ALLEEG(dscount).event.feedback],fbk_oi));
                    elseif strcmp(fbk_action,'keep')
                        indx2exl_v2{dscount,1} = find(contains([ALLEEG(dscount).event.feedback],fbk_oi));
                    end
                    if ~isempty(indx2exl_v2{dscount,1})
                        Dall2{1,dscount} = Dall2{1,dscount}(:,:,indx2exl_v2{dscount,1});
                    else
                        Dall2{1,dscount} = Dall2{1,dscount};
                    end
                else
                    disp('0Oops! Feedback field does not exist in current dataset.');
                    Dall2{1,dscount} = Dall2{1,dscount};
                    indx2exl_v2{dscount,1} = [];
                end
 
            end
        end

        tnum_ds2 = cell2mat(cellfun(@(x) [1:size(x, 3)], Dall2, 'UniformOutput',false));
        vid2_trialindx{vidcntr1,1} = tnum_ds2;
        Dmean_curr2 = cell2mat(cellfun(@(x) mean(x,3),Dall2,'UniformOutput',false));
        Dmean_curr2 = reshape(Dmean_curr2,[size(Dall2{1,1},1),size(Dall2{1,1},2),length(ALLEEG)]);
        GrandAvg2{vidcntr1,1} = Dmean_curr2;
        if vidcntr1==1
            DMCond2 = Dmean_curr2;
        elseif vidcntr1>1
            DMCond2 = cat(3,DMCond2,Dmean_curr2);
        end
    end
end

Tnum_cond1 = cell2mat(vid1_trialindx');
Tnum_cond2 = cell2mat(vid2_trialindx');

%% Find the grand average activity over all subjects for each video and each condition.

Data2plot = cell(1,2);
Data2plot{1,1} = mean(DMCond1,3); % The mean over all the subjects for each condition
Data2plot{1,2} = mean(DMCond2,3);

Data2plot_subs = cell(1,2);
Data2plot_subs{1,1} = DMCond1;
Data2plot_subs{1,2} = DMCond2;

%% Define the channels of interest
[Eindx,Etf] = listdlg('PromptString',{'Select electrodes to plot.',...
    'You can select several electrodes at once.',''},...
    'SelectionMode','multiple','ListString',{chanscurr.labels});  %Eindx gives the indices of the selected channels

%% Define the configuration of the ERP plot.

prompt = 'How do you want me to configure the ERP data: "electrode" or "grid" or "grid_intconf"?';
dlgtitre = 'Specify ERP Layout';
dims = [1 80];
answr = inputdlg(prompt,dlgtitre,dims);

%% PLOTTING

allchans = length(Eindx);  % The number of channels to plot;
notmtchans = Eindx;         % The indices of the channels to plot.
chanoi = {chanscurr(Eindx).labels};
t = time;
colrs=[ones(1,3).*[0.4 0.15 0.15];ones(1,3).*[0.4 0.4 0.4];ones(1,3).*[0 0.6 0.5];ones(1,3).*[0.2 0.1 0.8]];

if strcmp(answr{1,1},'electrode')

    % SET UP ELECTRODE CONFIGURATION

    hndl=figure; set(hndl,'PaperUnits','normalized','Position',[680 417 727 681],'Color',[1 1 1]);
    orient portrait; axis ('normal')

    xvals=zeros(allchans,1);
    yvals=zeros(allchans,1);
    pwidth    = 0.8;     % 0.75, width and height of plot array on figure
    pheight   = 0.8;
    axwidth=0.07;
    axheight=0.08;

    %Read in channel locations file
    [elocs,titres,theta,rads,inds]=readlocs(chaninfo(notmtchans));
    channoms= strvcat(chaninfo(1:allchans).labels);
    Th=pi/180*theta;           %convert degrees to radians

    %Convert from polar to cartesian
    [ycart,xcart]=pol2cart(Th,rads);
    xvals(notmtchans) = ycart;
    yvals(notmtchans) = xcart;

    %Find the positions of all the channels
    mtchans = setdiff(1:allchans,notmtchans);        %find the channels indices common to both
    allchans_sqrt = floor(sqrt(allchans))+1;

    for i=1:length(mtchans)

        xvals(mtchans(i))=0.5+0.2*floor((i-1)/74);  %allchans - x axes
        yvals(mtchans(i))=-0.2+mod(i-1,allchans)/74;   %allchans - y axes

    end

    channoms2=channoms(1:allchans,:);
    xvals = xvals(1:allchans);
    yvals = yvals(1:allchans);

    if length(xvals) > 1
        if length(unique(xvals)) > 1
            gcapos = get(gca,'Position'); axis off;
            xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % this recenters
            xvals = gcapos(1)+gcapos(3)/2+pwidth*xvals;  %this controls width of plot
        end
    end

    gcapos = get(gca,'Position'); axis off;
    yvals = gcapos(2)+gcapos(4)/2+pheight*yvals;  % controls height of plot


    ho=zeros(length(chanoi),1);
    sig_elecs = zeros(length(chanoi),length(t));
    sig_times = cell(length(chanoi),1);
    nosigs=0;
    pvalue_corr = cell(length(chanoi),1);
    tvals = cell(length(chanoi),1);

    %% PLOT THE DATA IN 10-20 CONFIGURATION

    axs = zeros(2,1);   %initialise the axes
    linz = zeros(2,1);  %initialise the line object
    toperm = cell(1,2);
    legnoms = cell(1,2); %initialise the legend names cell string array

    Axes = [];

    for chancnt =1:length(chanoi)

        xcenter(chancnt) = xvals(notmtchans(chancnt));
        ycenter(chancnt) = yvals(notmtchans(chancnt));
        Axes = [Axes axes('Units','normalized','Position', [ycenter(chancnt)-axheight/2 xcenter(chancnt)-axwidth/2 axheight axwidth])];
        hold on;

        for dcnt = 1:size(Data2plot,2)     % for each condition
            Dcurr = [];
            Dcurr = Data2plot{1,dcnt}(chancnt,:)';
            ploth1 = plot(time,Dcurr,'Color',colrs(dcnt,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
            axs1 = gca;
            hold on
            axs(dcnt) = axs1;
            linz(dcnt) = ploth1;
            toperm{1,dcnt} = squeeze(Data2plot_subs{1,dcnt}(chancnt,:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end

        set(axs1,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties
        title(chanoi{1,chancnt})
        [tvals_curr,pvalue_corr,permout]=plot_perm(toperm,time,3,[time(1) time(end)],'no');  %Call of function to carry out permutation t-test with fdr correction.

        set(hndl,'CurrentAxes',Axes(chancnt));
        set(Axes(chancnt),'HitTest','on','SelectionHighlight','on','UserData',{chancnt, squeeze(Data2plot_subs{1,1}(chancnt,:,:)),squeeze(Data2plot_subs{1,2}(chancnt,:,:)),...
            colrs, time, legnoms, chanoi, toperm, [meandur1 meandur2; std_durs1 std_durs2]},'NextPlot','add');
        set(Axes(chancnt),'ButtonDownFcn',@plotsinglePEPs_pe)

        if chancnt == 1
            legend(linz,legnoms,'Position',[0.1 0.75 0.2 0.1],'FontSize',14,'Box','off');
        end

    end  % end of chancnt loop
elseif strcmp(answr{1,1},'grid')

    hndl=figure; set(hndl,'PaperUnits','normalized','Position',[680 417 727 681],'Color',[1 1 1]);
    orient portrait; axis ('normal')

    axs = zeros(2,1);   %initialise the axes
    linz = zeros(2,1);  %initialise the line object
    toperm = cell(1,2);
    legnoms = cell(1,2); %initialise the legend names cell string array
    
    

    Axes = [];

    for chancnt =1:length(chanoi)

        Axes(chancnt) = subplot(4,3,chancnt);


        for dcnt = 1:size(Data2plot,2)     % for each condition
            Dcurr = [];
            Dcurr = Data2plot{1,dcnt}(chancnt,:)';
            ploth1 = plot(time,Dcurr,'Color',colrs(dcnt,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
            axs1 = gca;
            hold on
            axs(dcnt) = axs1;
            linz(dcnt) = ploth1;
            toperm{1,dcnt} = squeeze(Data2plot_subs{1,dcnt}(chancnt,:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end

        set(axs1,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties
        title(chanoi{1,chancnt})
        [tvals_curr,pvalue_corr,permout]=plot_perm(toperm,time,3,[time(1) time(end)],'no');  %Call of function to carry out permutation t-test with fdr correction.

        set(hndl,'CurrentAxes',Axes(chancnt));
        set(Axes(chancnt),'HitTest','on','SelectionHighlight','on','UserData',{chancnt, squeeze(Data2plot_subs{1,1}(chancnt,:,:)),squeeze(Data2plot_subs{1,2}(chancnt,:,:)),...
            colrs, time, legnoms, chanoi, toperm, [meandur1 meandur2; std_durs1 std_durs2]},'NextPlot','add');
        set(Axes(chancnt),'ButtonDownFcn',@plotsinglePEPs_pe)

        if chancnt == 1
            legend(linz,legnoms,'Position',[0.1 0.75 0.2 0.1],'FontSize',14,'Box','off');
        end

    end  % end of chancnt loop%
elseif strcmp(answr{1,1},'grid_intconf')



    f1 = figure;
    set(f1,'Color',[1 1 1]);

    chanoi = [];
    chanoi = {'F3' 'Fz' 'F4' 'FC3' 'FCz' 'FC4' 'C3' 'Cz' 'C4' 'P3' 'Pz' 'P4'};

    allchans = {chaninfo.labels};
    for i = 1:length(chanoi)
        chanindx(i) = find(strcmp(allchans,chanoi{1,i}));
    end

    for chancnt = 1:length(chanoi)

        toperm = cell(1,2);
        legnoms = cell(1,2);
        for dcnt = 1:size(Data2plot,2)     % for each condition
            toperm{1,dcnt} = squeeze(Data2plot_subs{1,dcnt}(chanindx(chancnt),:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end

        Axes(chancnt) = subplot(4,3,chancnt);

        axs = zeros(1,size(toperm,2));
        linz = zeros(1,size(toperm,2));

        Cond1 = squeeze(Data2plot_subs{1,1}(chanindx(chancnt),:,:));
        Cond2 = squeeze(Data2plot_subs{1,2}(chanindx(chancnt),:,:));
        fbstat = [meandur1 meandur2; std_durs1 std_durs2];

        [axs1,ploth1] = plot_ConfInt2(Cond1,time,nanmean(Cond1',1)',colrs(1,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
        hold on
        axs(1,1) = axs1;
        linz(1,1) = ploth1;

        [axs2,ploth2]= plot_ConfInt2(Cond2,time,nanmean(Cond2',1)',colrs(2,:));
        hold on
        axs(1,2) = axs2;
        linz(1,2) = ploth2;

        % Add duration annotations
        ha1 = annotation('arrow');
        ha1.Parent = f1.CurrentAxes;  % associate annotation with current axes
        % now you can use data units
        ha1.X = [0 fbstat(1,1)*1000];
        ha1.Y = [-1.5 -1.5];
        ha1.Color = 'blue';
        ha1.LineWidth = 2;

        ha2 = annotation('arrow');
        ha2.Parent = f1.CurrentAxes;  % associate annotation with current axes
        % now you can use data units
        ha2.X = [0 fbstat(1,2)*1000];
        ha2.Y = [-1.3 -1.3];
        ha2.Color = 'red';
        ha2.LineWidth = 2;

        % Add the standard deviation information of dur1
        std1_upper = fbstat(1,1)+fbstat(2,1);
        if (std1_upper*1000)>time(end)
            std1_upper = time(end)/1000;
        end
        std1_lower = fbstat(1,1)-fbstat(2,1);
        ha1b = annotation('doublearrow');
        ha1b.Parent = f1.CurrentAxes;
        ha1b.X = [std1_lower*1000 std1_upper*1000];
        ha1b.Y = [-1.5 -1.5];
        ha1b.Color = [0.3 0.3 0.3];
        ha1b.LineStyle = '--';


        % Add the standard deviation information00 of dur2
        std2_upper = fbstat(1,2)+fbstat(2,2);
        if (std2_upper*1000)>time(end)
            std2_upper = time(end)/1000;
        end
        std2_lower = fbstat(1,2)-fbstat(2,2);
        ha2b = annotation('doublearrow');
        ha2b.Parent = f1.CurrentAxes;
        ha2b.X = [std2_lower*1000 std2_upper*1000];
        ha2b.Y = [-1.3 -1.3];
        ha2b.Color = [0.3 0.3 0.3];
        ha2b.LineStyle = '--';

        set(Axes(chancnt),'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties
        Titre = title(allchans{1,chanindx(chancnt)});
        Titre.FontSize = 16;
        xlabel(gca,'Time (ms)','FontSize',12)
        ylabel(gca,'Potential (\muV)','FontSize', 12);

        if chancnt ==1
            legend(linz,legnoms,'FontSize', 14);
        end
    end
end

trialcong_num = cellfun(@(x) (size(x,2)), vid1_trialindx, 'UniformOutput', false);
trialincong_num = cellfun(@(x)(size(x,2)), vid2_trialindx, 'UniformOutput', false);
trialtotal_cong = sum(cell2mat(trialcong_num));
trialtotal_incong = sum(cell2mat(trialincong_num));

sprintf('********Trial total for condition 1 (Dall): %d ***************',trialtotal_cong);
sprintf('********Trial total for condition 2 (Dall2): %d **************', trialtotal_incong);


%%



