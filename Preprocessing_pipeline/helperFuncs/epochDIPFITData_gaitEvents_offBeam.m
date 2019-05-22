function [numEpchs]=epochDIPFITData_gaitEvents_offBeam(EEG,pathname,subjnum,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix)

%Define paths
subjcode = ['WMISM_' num2str(subjnum)];
% subj_basepath = [pathname subjcode filesep]; % Top level subject folder
% EEGsets_outpath = [subj_basepath 'EEG_sets'];
% ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
% epch_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group '_epch'];
mergeEpch_path = [pathname 'STUDY_sets_' folderSuffix]; %group '_gait'];
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
% load([ICA_path_out filesep group '_ICA_DIPFIT_5.mat']);
% EEG = update_EEG(EEG,ICA_STRUCT);
% 
% %Remove bad ICs (only keeping cortical ICs)
% compInds=setdiff([1:size(EEG.icaweights,1)],ICA_STRUCT.good_comps.brain);
% EEG=pop_subcomp(EEG,compInds);
EEG=eeg_checkset(EEG);

% %Create trig labels
% if analyzeFootDirection==0
%     %Remove CW/CCW from event labels
%     for i=1:length({EEG.event(1:end).type})
%         bcLeft=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_Left_*']),'once');
%         tcLeft=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_Left_*']),'once');
%         bcRight=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_Right_*']),'once');
%         tcRight=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_Right_*']),'once');
%         if ~isempty(bcLeft)
%             EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(9:end)];
%         elseif ~isempty(tcLeft)
%             EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(9:end)];
%         elseif ~isempty(bcRight)
%             EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(10:end)];
%         elseif ~isempty(tcRight)
%             EEG.event(1,i).type=[EEG.event(1,i).type(1:3) EEG.event(1,i).type(10:end)];
%         end
%     end
% end

% % EEGevOrig=EEG.event;
% % %First, remove all events that aren't gait events
% % badInds=[];
% % for i=1:length({EEG.event(1:end).type})
% %     bcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_*']),'once');
% %     tcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_*']),'once');
% %     if isempty(bcEv) && isempty(tcEv)
% %         badInds=[badInds i];
% %     end
% % end
% % EEG.event(badInds)=[];
% % 
% % %Now, define events
% % for j=1:length(condsTW)
% %     %First, remove all events that aren't gait events
% %     lastSawTC=0; lastSawBC=0; numBC=0; numTC=0; EEG.event=EEGevOrig; badInds=[];
% %     for i=1:length({EEG.event(1:end).type})
% %         bcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_*']),'once');
% %         tcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_*']),'once');
% %         if isempty(bcEv) && isempty(tcEv)
% %             badInds=[badInds i];
% %         end
% %     end
% %     EEG.event(badInds)=[];
% % end
% %         if ~isempty(bcEv)
% %             if lastSawTC==1
% %                 numBC=1; lastSawTC=0;
% %             else
% %                 numBC=numBC+1;
% %             end
% %             lastSawBC=1;
% % %             EEG.event(1,i).type=[EEG.event(1,i).type num2str(numBC)];
% %         elseif ~isempty(tcEv)
% %             if lastSawBC==1
% %                 numTC=1; lastSawBC=0;
% %                 EEG.event(1,i).type=[EEG.event(1,i).type num2str(numTC)];
% %             else
% %                 numTC=numTC+1;
% %             end
% %             lastSawTC=1;
% % %             EEG.event(1,i).type=[EEG.event(1,i).type num2str(numTC)];
% %         end
% %     end
% % end

%Go through events and add numbers for timewarping
for j=1:length(condsTW)
    lastSawTC=0; lastSawBC=0; numBC=0; numTC=0;
    for i=1:length({EEG.event(1:end).type})
        bcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['BC_*_' condsTW{j} '*']),'once');
        tcEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['TC_*_' condsTW{j} '*']),'once');
        if ~isempty(bcEv)
            if lastSawTC==1
                numBC=1; lastSawTC=0;
            else
                numBC=numBC+1;
            end
            lastSawBC=1;
%             EEG.event(1,i).type=[EEG.event(1,i).type num2str(numBC)];
        elseif ~isempty(tcEv)
            if lastSawBC==1
                numTC=1; lastSawBC=0;
                EEG.event(1,i).type=[EEG.event(1,i).type num2str(numTC)];
            else
                numTC=numTC+1;
            end
            lastSawTC=1;
%             EEG.event(1,i).type=[EEG.event(1,i).type num2str(numTC)];
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
    EEG = pop_epoch(EEG,{[twSeq{3} 'Right_' condsTW{p} '1'],[twSeq{3} 'Left_' condsTW{p} '1']},[EpPre EpPost],... %'1']},[EpPre EpPost],... %trigs(1,p),[EpPre EpPost],...
        'newname',[subjcode '_Epc'],'epochinfo','yes');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0);
    EEG=eeg_checkset(EEG,'makeur');
    EEG = pop_rmbase( EEG,[]);
%     %Get rid of really noisy epochs (due to ungrounding)
%     [EEG,inds]=pop_eegthresh(EEG,1,1:EEG.nbchan,-500,500,EpPre,EpPost,0,1);
    
    %Set up sequence
    trueTWseq1=cell(1,length(twSeq)); trueTWseq2=cell(1,length(twSeq));
    nums={'','','1',''}; %{'1','2','1','2','3'};
    sides={'Left_','Right_','Left_','Right_','Left_'};
    for i=1:length(twSeq)
        trueTWseq1{i}=[twSeq{i} sides{i} condsTW{p} nums{i}];
        trueTWseq2{i}=[twSeq{i} sides{i+1} condsTW{p} nums{i}];
    end
    timewarp1 = make_timewarp(EEG,trueTWseq1,'baselineLatency',epochTimes(1,1)*1000, ...
        'maxSTDForAbsolute',3,...
        'maxSTDForRelative',3);
    timewarp2 = make_timewarp(EEG,trueTWseq2,'baselineLatency',epochTimes(1,1)*1000, ...
        'maxSTDForAbsolute',3,...
        'maxSTDForRelative',3);
%     %Make sure no overlap in epochs (if so, remove from both sides)
%     inds2Rem=[];
%     for k=1:length(timewarp1.epochs)
%         if any(timewarp1.epochs(k)==timewarp2.epochs)
%             indT=find(timewarp1.epochs(k)==timewarp2.epochs);
%             timewarp2.latencies(indT,:)=[];
%             timewarp2.epochs(indT)=[];
%             inds2Rem=[inds2Rem,k];
%         end
%     end
%     timewarp1.latencies(inds2Rem,:)=[];
%     timewarp1.epochs(inds2Rem)=[];
    
    %Check for any latencies outside of true epoch boundaries or where 1st
    %event not at 0
    tooHigh=find(timewarp1.latencies(:,end)>(epochTimes(1,2)*1000));
    notZero=find(timewarp1.latencies(:,3)~=0);
    doubleEvs=find(timewarp1.latencies(:,1)==timewarp1.latencies(:,2));
    badInds=unique([tooHigh; notZero; doubleEvs]);
    timewarp1.epochs(badInds)=[];
    timewarp1.latencies(badInds,:)=[];
    
    %Check for any latencies outside of true epoch boundaries or where 1st
    %event not at 0
    tooHigh=find(timewarp2.latencies(:,end)>(epochTimes(1,2)*1000));
    notZero=find(timewarp2.latencies(:,3)~=0);
    doubleEvs=find(timewarp2.latencies(:,1)==timewarp2.latencies(:,2));
    badInds=unique([tooHigh; notZero; doubleEvs]);
    timewarp2.epochs(badInds)=[];
    timewarp2.latencies(badInds,:)=[];

    %Merge timewarps
    timewarp.latencies=cat(1,timewarp1.latencies,timewarp2.latencies);
    timewarp.epochs=sort(unique(cat(2,timewarp1.epochs,timewarp2.epochs)));
    
    if length(timewarp.epochs)~=(length(timewarp1.epochs)+length(timewarp2.epochs))
        error('Epochs still somehow overlapping!');
    end
    
    timewarp.warpto = median(timewarp.latencies);
    goodepochs=sort([timewarp.epochs]);
    EEG=pop_select(EEG,'trial',goodepochs);
    EEG.timewarp.latencies = timewarp.latencies;
    EEG.warpto = timewarp.warpto;
    
    %Now change gait events so no foot side
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
    %Do it for EEG.epoch as well (not like it matters...)
    EEG=eeg_checkset(EEG,'eventconsistency');
    
    EEG=pop_editset(EEG,'setname',['Epch_' group '_' num2str(subjnum) '_' condsTW{p}]);
%     EEG=eeg_checkset(EEG,'ica');
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
%     EEG=update_EEG(EEG,ICA_STRUCT,1);
%     EEG=pop_subcomp(EEG,compInds);
%     EEG=eeg_checkset(EEG);
    pop_saveset(EEG, 'filename', ['Merge_epch_' group '_' num2str(subjnum)],...
        'filepath', mergeEpch_path);
else
    display('No need to merge datasets...')
end
disp(['Finished ' subjcode '!']);
disp(['Epochs remaining are: ' epchNums]);
