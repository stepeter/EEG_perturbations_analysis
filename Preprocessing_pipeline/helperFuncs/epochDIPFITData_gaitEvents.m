function [numEpchs]=epochDIPFITData_gaitEvents(EEG,pathname1,pathname2,subjnum,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix)

%Define paths
subjcode = ['WMISM_' num2str(subjnum)];
subj_basepath = [pathname1 subjcode filesep]; % Top level subject folder
% EEGsets_outpath = [subj_basepath 'EEG_sets'];
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
% epch_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group '_epch'];
mergeEpch_path = [pathname2 'STUDY_sets_' folderSuffix]; %group '_gait'];
% cd(subj_basepath)
numEpchs=[];

global ALLEEG CURRENTSET ALLCOM;

%Create epoch folder (if needed)
% if ~exist(epch_path_out)
%     mkdir(epch_path_out)
% end
if ~exist(mergeEpch_path)
    mkdir(mergeEpch_path)
end

%Load dataset before epoching
% EEG=pop_loadset([EEGsets_outpath filesep 'Merge_' group '_CAR.set']);
load([ICA_path_out filesep group '_ICA_DIPFIT_5.mat']);
EEG = update_EEG(EEG,ICA_STRUCT);

%Remove bad ICs (only keeping cortical ICs)
compInds=setdiff([1:size(EEG.icaweights,1)],ICA_STRUCT.good_comps.brain);
EEG=pop_subcomp(EEG,compInds);
EEG=eeg_checkset(EEG);

%Create trig labels
if analyzeFootDirection==0
    %Remove CW/CCW from event labels
    for i=1:length({EEG.event(1:end).type})
        bcLeft=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_Left_*']),'once');
        tcLeft=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_Left_*']),'once');
        bcRight=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_Right_*']),'once');
        tcRight=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_Right_*']),'once');
        if ~isempty(bcLeft)
            EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(9:end)];
        elseif ~isempty(tcLeft)
            EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(9:end)];
        elseif ~isempty(bcRight)
            EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(10:end)];
        elseif ~isempty(tcRight)
            EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(10:end)];
        end
    end
end

%Go through events and add numbers for timewarping
for j=1:length(condsTW)
    lastSawTC=0; lastSawBC=0; numBC=0; numTC=0;
    for i=1:length({EEG.event(1:end).type})
        bcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_' condsTW{j} '*']),'once');
        tcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_' condsTW{j} '*']),'once');
        if ~isempty(bcEv)
            if lastSawTC==1
                numBC=1; lastSawTC=0;
            else
                numBC=numBC+1;
            end
            lastSawBC=1;
            EEG.event(1,i).type=[EEG.event(1,i).type num2str(numBC)];
        elseif ~isempty(tcEv)
            if lastSawBC==1
                numTC=1; lastSawBC=0;
            else
                numTC=numTC+1;
            end
            lastSawTC=1;
            EEG.event(1,i).type=[EEG.event(1,i).type num2str(numTC)];
        end
    end
end


%Epoch parameters
CURRENTSET=1;
% Set epoch bounds (add in a little extra for ERSPs)
EpPre = epochTimes(1,1)+(-0.559);
EpPost = epochTimes(1,2)+0.559;
EEG_preEpoch=EEG; %save this for iterations in for loop
epchNums=''; timewarps=[]; warptos=[];

for p=1:length(condsTW)
    EEG=EEG_preEpoch; %use original (non-epoched) data
    EEG = pop_epoch(EEG,{[twSeq{1} condsTW{p}]},[EpPre EpPost],... %'1']},[EpPre EpPost],... %trigs(1,p),[EpPre EpPost],...
        'newname',[subjcode '_Epc'],'epochinfo','yes');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0);
    EEG=eeg_checkset(EEG,'makeur');
    EEG = pop_rmbase( EEG,[]);
    %Get rid of really noisy epochs (due to ungrounding)
    [EEG,inds]=pop_eegthresh(EEG,1,1:EEG.nbchan,-500,500,EpPre,EpPost,0,1);
    
    %Set up sequence
    trueTWseq=cell(1,length(twSeq)); nums={'1','2','1','2','3'};
    for i=1:length(twSeq)
        trueTWseq{i}=[twSeq{i} condsTW{p}]; % nums{i}];
    end
    timewarp = make_timewarp(EEG,trueTWseq,'baselineLatency',0, ...
        'maxSTDForAbsolute',3,...
        'maxSTDForRelative',3);
    
    %Check for any latencies outside of true epoch boundaries or where 1st
    %event not at 0
    tooHigh=find(timewarp.latencies(:,end)>(epochTimes(1,2)*1000));
    notZero=find(timewarp.latencies(:,1)~=0);
    badInds=unique([tooHigh; notZero]);
    timewarp.epochs(badInds)=[];
    timewarp.latencies(badInds,:)=[];
    
    timewarp.warpto = median(timewarp.latencies);
    goodepochs=sort([timewarp.epochs]);
    EEG=pop_select(EEG,'trial',goodepochs);
    EEG.timewarp.latencies = timewarp.latencies;
    EEG.warpto = timewarp.warpto;
    
    EEG=pop_editset(EEG,'setname',['Epch_' group '_' num2str(subjnum) '_' condsTW{p}]);
    EEG=eeg_checkset(EEG,'ica');
%     pop_saveset(EEG, 'filename', ['Epch_' group '_' num2str(subjnum) '_' trigs{1,p}], 'filepath', epch_path_out); %commented this out beceause only care about merged file SP 4/26/17
    EEG_epched(p)=EEG; %save set to be merged later
    epchNums=[epchNums ' ' num2str(length(timewarp.epochs))];
    numEpchs=[numEpchs length(timewarp.epochs)];
    
    %Save timewarp info for merging sets
    if isfield(EEG,'warpto')
        timewarps=[timewarps; EEG.timewarp.latencies]; %concatenate timewarp latencies
        warptos=[warptos; EEG.warpto];
    end
end

%Now merge these saved sets
if length(EEG_epched) > 1
    EEG = [];
    EEG = pop_mergeset(EEG_epched, 1:length(EEG_epched));
    if isfield(EEG,'warpto')
        EEG.warpto=median(warptos);
        EEG.timewarp.latencies=timewarps; %save median timewarp latencies
    end
    EEG.setname = ['Merge_epch_' group '_' num2str(subjnum)];
    EEG = eeg_checkset(EEG); % always checkset 
    EEG=update_EEG(EEG,ICA_STRUCT,1);
    EEG=pop_subcomp(EEG,compInds);
    EEG=eeg_checkset(EEG);
    pop_saveset(EEG, 'filename', ['Merge_epch_' group '_' num2str(subjnum)],...
        'filepath', mergeEpch_path);
else
    display('No need to merge datasets...')
end
disp(['Finished ' subjcode '!']);
disp(['Epochs remaining are: ' epchNums]);
