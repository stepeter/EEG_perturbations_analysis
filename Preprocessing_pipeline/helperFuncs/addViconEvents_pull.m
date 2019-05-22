%Load EEG, Vicon events, and Vicon sync
pathname='/media/stepeter/New Volume/VR_rotations/Data/';
subjnum=2;
Vicon_sample_rate=1000; %Hz
syncFreq=0.5; %Hz
group='all';

%Load data
cd(pathname);
EEGsets_outpath=[pathname 'WMISM_' num2str(subjnum) filesep 'EEG_sets/'];
EEG=pop_loadset('filename', ['Merge_' group '_CAR.set'], 'filepath', EEGsets_outpath);
load([EEGsets_outpath 'Events/sync.mat']);
load([EEGsets_outpath 'Events/finalPeaks.mat']);

switch group
    case 'mismatches'
        error('No pull perturbations for these trials!');
    case 'pull'
        csv_sync_ch_idx=[1:3 13:15];
        break_inds=[3 6];
        trials={'Stn','Stn','Stn','Wlk','Wlk','Wlk'};
        %Vicon file should be ordered: Stand, Walk
    case 'train'
        csv_sync_ch_idx=[4:12];
        break_inds=[4 7 10]; 
        trials={'Pre','Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3','Post'};
        %Vicon file should be ordered: Pre, Train1, Train2, Train3, Post
    case 'all'
        csv_sync_ch_idx=[4:12]; EEGTrialnames={'Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3'}; 
        EEG_inds=[14:22]; %order of trials in EEG
        trials={'Stn','Stn','Stn','Wlk','Wlk','Wlk','SVZ','SVZ','SVZ','WVZ','WVZ','WVZ','Pre','Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3','Post'};
    otherwise
        error('I do not recognize this group name!');
end

% %Check consistency with fix_sync() %don't think this is necessary (off by
% 1-2 msec)
% for i=csv_sync_ch_idx
%     [new_eeg_trig_latencies,new_csv_trig_latencies,EEG,new_biomech] = ...
%         fix_sync_ViconEvents(EEG, sync, i, csv_sample_rate, syncFreq);
% end

sync_Rise_EEG=[];
%Find rising sync latencies
for i=1:length({EEG.event(1:end).type})
    syncEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['Sync Rising' '*']),'once'); %[trigPrefix 'CCW_*']),'once');
    if ~isempty(syncEv)
        sync_Rise_EEG=[sync_Rise_EEG EEG.event(1,i).latency];
    end
end

sync_Fall_EEG=[];
%Find falling sync latencies
for i=1:length({EEG.event(1:end).type})
    syncEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['Sync Falling' '*']),'once'); %[trigPrefix 'CCW_*']),'once');
    if ~isempty(syncEv)
        sync_Fall_EEG=[sync_Fall_EEG EEG.event(1,i).latency];
    end
end

% if subjnum==1
%     bound=[];
%     for i=1:length({EEG.event(1:end).type})
%         bEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['boundary' '*']),'once'); %[trigPrefix 'CCW_*']),'once');
%         if ~isempty(bEv)
%             bound=[bound EEG.event(1,i).latency];
%         end
%     end
% end

%Find breaks in EEG sync
eeg_inds_rise=find(diff(sync_Rise_EEG)>(1.2*EEG.srate/syncFreq)); %Find breaks in EEG sync
eeg_inds_fall=find(diff(sync_Fall_EEG)>(1.2*EEG.srate/syncFreq)); %Find breaks in EEG sync

%Check start, end, and around breaks for sync signals that aren't 2 secs in
%between (RISE)
if (sync_Rise_EEG(end)-sync_Rise_EEG(end-1))<(EEG.srate/syncFreq-5) || (sync_Rise_EEG(end)-sync_Rise_EEG(end-1))>(EEG.srate/syncFreq+5)
    sync_Rise_EEG(end)=[];
end
for i=length(eeg_inds_rise):-1:1
    %After break
    if (sync_Rise_EEG(eeg_inds_rise(i)+2)-sync_Rise_EEG(eeg_inds_rise(i)+1))<(EEG.srate/syncFreq-5) || (sync_Rise_EEG(eeg_inds_rise(i)+2)-sync_Rise_EEG(eeg_inds_rise(i)+1))>(EEG.srate/syncFreq+5)
        sync_Rise_EEG(eeg_inds_rise(i)+1)=[];
    end
    %Before break
    if (sync_Rise_EEG(eeg_inds_rise(i))-sync_Rise_EEG(eeg_inds_rise(i)-1))<(EEG.srate/syncFreq-5) || (sync_Rise_EEG(eeg_inds_rise(i))-sync_Rise_EEG(eeg_inds_rise(i)-1))>(EEG.srate/syncFreq+5)
        sync_Rise_EEG(eeg_inds_rise(i))=[];
        eeg_inds_rise(i)=eeg_inds_rise(i)-1;
    end
end
if (sync_Rise_EEG(2)-sync_Rise_EEG(1))<(EEG.srate/syncFreq-5) || (sync_Rise_EEG(2)-sync_Rise_EEG(1))>(EEG.srate/syncFreq+5)
    sync_Rise_EEG(1)=[];
end

%Check start, end, and around breaks for sync signals that aren't 2 secs in
%between (FALL)
if (sync_Fall_EEG(end)-sync_Fall_EEG(end-1))<(EEG.srate/syncFreq-5) || (sync_Fall_EEG(end)-sync_Fall_EEG(end-1))>(EEG.srate/syncFreq+5)
    sync_Fall_EEG(end)=[];
end
for i=length(eeg_inds_fall):-1:1
    %After break
    if (sync_Fall_EEG(eeg_inds_fall(i)+2)-sync_Fall_EEG(eeg_inds_fall(i)+1))<(EEG.srate/syncFreq-5) || (sync_Fall_EEG(eeg_inds_fall(i)+2)-sync_Fall_EEG(eeg_inds_fall(i)+1))>(EEG.srate/syncFreq+5)
        sync_Fall_EEG(eeg_inds_fall(i)+1)=[];
    end
    %Before break
    if (sync_Fall_EEG(eeg_inds_fall(i))-sync_Fall_EEG(eeg_inds_fall(i)-1))<(EEG.srate/syncFreq-5) || (sync_Fall_EEG(eeg_inds_fall(i))-sync_Fall_EEG(eeg_inds_fall(i)-1))>(EEG.srate/syncFreq+5)
        sync_Fall_EEG(eeg_inds_fall(i))=[];
        eeg_inds_fall(i)=eeg_inds_fall(i)-1;
    end
end
if (sync_Fall_EEG(2)-sync_Fall_EEG(1))<(EEG.srate/syncFreq-5) || (sync_Fall_EEG(2)-sync_Fall_EEG(1))>(EEG.srate/syncFreq+5)
    sync_Fall_EEG(1)=[];
end

%Find breaks in EEG sync (now that removed some)
eeg_inds_rise=find(diff(sync_Rise_EEG)>(1.2*EEG.srate/syncFreq)); %Find breaks in EEG sync
eeg_inds_fall=find(diff(sync_Fall_EEG)>(1.2*EEG.srate/syncFreq)); %Find breaks in EEG sync

% if subjnum==1
%     for i=2:length(bound)
%         ind_rise=find(sync_Rise_EEG<bound(i),1,'last');
%         ind_fall=find(sync_Fall_EEG<bound(i),1,'last');
%         eeg_inds_rise=[eeg_inds_rise ind_rise];
%         eeg_inds_fall=[eeg_inds_fall ind_fall];
%     end
%     eeg_inds_rise=sort(eeg_inds_rise,2,'ascend');
%     eeg_inds_fall=sort(eeg_inds_fall,2,'ascend');    
% end

if (length(eeg_inds_rise)+1)~=length(trials) %length(csv_sync_ch_idx)
    error('EEG rise and CSV breaks are not equal!');
end
if (length(eeg_inds_fall)+1)~=length(trials) %length(csv_sync_ch_idx)
    error('EEG fall and CSV breaks are not equal!');
end

inds_Prior=zeros(2,length(csv_sync_ch_idx));
for i=1:length(eeg_inds_rise)
    inds_Prior(1,i+1)=eeg_inds_rise(i);
    inds_Prior(2,i+1)=eeg_inds_fall(i);
end

%If no errors, find Vicon sync rising and falling edges
ViconSyncInds=cell(2,length(csv_sync_ch_idx));
for i=1:length(csv_sync_ch_idx)
    diffSync=diff(sync{2,csv_sync_ch_idx(i)});
    ViconSyncInds{1,i}=find(diffSync>2); %rising edge
    ViconSyncInds{2,i}=find(diffSync<-2); %falling edge
end

%Check start, end for Vicon sync signals that aren't 2 secs in
%between (FALL)
for i=1:2
    for j=1:size(ViconSyncInds,2)
        if (ViconSyncInds{i,j}(end)-ViconSyncInds{i,j}(end-1))<(Vicon_sample_rate/syncFreq-5) || (ViconSyncInds{i,j}(end)-ViconSyncInds{i,j}(end-1))>(Vicon_sample_rate/syncFreq+5)
            ViconSyncInds{i,j}(end)=[];
        end
        if (ViconSyncInds{i,j}(2)-ViconSyncInds{i,j}(1))<(Vicon_sample_rate/syncFreq-5) || (ViconSyncInds{i,j}(2)-ViconSyncInds{i,j}(1))>(Vicon_sample_rate/syncFreq+5)
            ViconSyncInds{i,j}(1)=[];
        end
    end
end

%Check that number of rising and falling edges are equal between Vicon and EEG
totalEEGInds=[inds_Prior [length(sync_Rise_EEG); length(sync_Fall_EEG)]];
totalEEGInds=diff(totalEEGInds')'; notEq=zeros(2,size(ViconSyncInds,2));
for i=1:2
    for j=1:size(ViconSyncInds,2)
        totalViconInds(i,j)=length(ViconSyncInds{i,j});
        if totalEEGInds(i,EEG_inds(j))~=totalViconInds(i,j)
            notEq(i,j)=1;
        end
    end
end



%Go through each Vicon trial and find closest sync to each and difference,
%and add to EEG events
N=length(EEG.event); q=0;
for i=1:length(csv_sync_ch_idx)
    %Left pulls are 2*i-1
    for j=1:length(finalPeaks{3,(2*csv_sync_ch_idx(i)-1)})
        ev_Left=finalPeaks{3,(2*csv_sync_ch_idx(i)-1)}(j);
        [Y_left_rise,I_left_rise]=min(abs(ViconSyncInds{1,i}-ev_Left));
        [Y_left_fall,I_left_fall]=min(abs(ViconSyncInds{2,i}-ev_Left));
%         Y_left_rise=-(1-notEq(1,i))*Y_left_rise;
%         Y_left_fall=-(1-notEq(2,i))*Y_left_fall;
        if Y_left_rise<=Y_left_fall
            ev_diff_Vicon_Left=ev_Left-ViconSyncInds{1,i}(I_left_rise); %difference btwn sync and event in Vicon
            EEG.event(1,N+q).latency=sync_Rise_EEG(I_left_rise+inds_Prior(1,EEG_inds(i)))+...
                (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
            EEG.event(1,N+q).type=['Pull_Left_' EEGTrialnames{i}];
        else
            ev_diff_Vicon_Left=ev_Left-ViconSyncInds{2,i}(I_left_fall); %difference btwn sync and event in Vicon
            EEG.event(1,N+q).latency=sync_Fall_EEG(I_left_fall+inds_Prior(2,EEG_inds(i)))+...
                (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
            EEG.event(1,N+q).type=['Pull_Left_' EEGTrialnames{i}];   
        end
        if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
            error('Event is out of bounds!');
        end
        q=q+1;
    end
    
    %Right pulls are 2*i
    for j=1:length(finalPeaks{3,(2*csv_sync_ch_idx(i))})
        ev_Right=finalPeaks{3,(2*csv_sync_ch_idx(i))}(j);
        [Y_right_rise,I_right_rise]=min(abs(ViconSyncInds{1,i}-ev_Right));
        [Y_right_fall,I_right_fall]=min(abs(ViconSyncInds{2,i}-ev_Right));
%         Y_right_rise=-(1-notEq(1,i))*Y_right_rise;
%         Y_right_fall=-(1-notEq(2,i))*Y_right_fall;
        if Y_right_rise<=Y_right_fall
            ev_diff_Vicon_Right=ev_Right-ViconSyncInds{1,i}(I_right_rise); %difference btwn sync and event in Vicon
            EEG.event(1,N+q).latency=sync_Rise_EEG(I_right_rise+inds_Prior(1,EEG_inds(i)))+...
                (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
            EEG.event(1,N+q).type=['Pull_Right_' EEGTrialnames{i}];
        else
            ev_diff_Vicon_Right=ev_Right-ViconSyncInds{2,i}(I_right_fall); %difference btwn sync and event in Vicon
            EEG.event(1,N+q).latency=sync_Fall_EEG(I_right_fall+inds_Prior(2,EEG_inds(i)))+...
                (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
            EEG.event(1,N+q).type=['Pull_Right_' EEGTrialnames{i}];
        end
        if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
            error('Event is out of bounds!');
        end
        q=q+1;
    end
end

EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

trialsCompressed=unique(EEGTrialnames); %unique(trials);
if strcmp(group,'pull') || strcmp(group,'all')
    %Check event latencies
    load('/media/stepeter/New Volume/VR_rotations/Data/template_pullEvs.mat');
%     startInds=[7.07 3.14; 4.32 11.79]; %L_stand, R_stand; L_walk, R_walk
    RMSE_vals=cell(2,length(trialsCompressed));
    for j=1:length(trialsCompressed)
        
        latL=[]; latR=[]; latL_eegEvInds=[]; latR_eegEvInds=[];
        for i=1:length({EEG.event(1:end).type})
            if strcmp(EEG.event(1,i).type,['Pull_Left_' trialsCompressed{j}])
                latL=[latL (EEG.event(1,i).latency-1)/EEG.srate];
                latL_eegEvInds=[latL_eegEvInds i];
            elseif strcmp(EEG.event(1,i).type,['Pull_Right_' trialsCompressed{j}])
                latR=[latR (EEG.event(1,i).latency-1)/EEG.srate];
                latR_eegEvInds=[latR_eegEvInds i];
            end
        end
        
%         %Checking absolute expected error (can also remove an extra event
%         %if really high error)
%         template_diffs=[];
%         for i=1:length(template_pullEvs(2*j-1,:));
%             template_diffs(i,:)=template_pullEvs(2*j-1,i)-template_pullEvs(2*j-1,:);
%         end
%         ind2remove=1; RMSE_vals{1,j}=[];
%         while (any(RMSE_vals{1,j}>20) || isempty(RMSE_vals{1,j})) && ind2remove<=length(latL)
%             RMSE_vals{1,j}=[];
%             if length(latL)~=length(template_diffs)
%                 inds=setdiff(1:length(latL),ind2remove);
%                 ind2remove=ind2remove+1;
%             else
%                 inds=1:length(template_diffs);
%             end
%             latL_cut=latL(inds);
%             for i=1:length(latL_cut)
%                 latDiffs=latL_cut(i)-latL_cut;
%                 minRMSE=100000;
%                 for k=1:length(template_diffs)
%                     RMSE_val=sum(abs(latDiffs-template_diffs(k,:)));
%                     if RMSE_val<minRMSE
%                         minRMSE=RMSE_val;
%                     end
%                 end
%                 RMSE_vals{1,j}(i,:)=minRMSE;
%             end
%             
%         end
%         %Remove the bad event
%         if ind2remove<=length(latL) && ind2remove>1
%             latL=latL(setdiff(1:length(latL),ind2remove-1));
%             EEG.event(latL_eegEvInds(ind2remove-1))=[];
%         end
%         if any(RMSE_vals{1,j}>20)
%             error('Error still too high! Manually inspect');
%         end
%         
%         %Checking absolute expected error (can also remove an extra event
%         %if really high error)
%         template_diffs=[];
%         for i=1:length(template_pullEvs(2*j,:));
%             template_diffs(i,:)=template_pullEvs(2*j,i)-template_pullEvs(2*j,:);
%         end
%         ind2remove=1; RMSE_vals{2,j}=[];
%         while (any(RMSE_vals{2,j}>20) || isempty(RMSE_vals{2,j})) && ind2remove<=length(latR)
%             RMSE_vals{2,j}=[];
%             if length(latR)~=length(template_diffs)
%                 inds=setdiff(1:length(latR),ind2remove);
%                 ind2remove=ind2remove+1;
%             else
%                 inds=1:length(template_diffs);
%             end
%             latR_cut=latR(inds);
%             for i=1:length(latR_cut)
%                 latDiffs=latR_cut(i)-latR_cut;
%                 minRMSE=100000;
%                 for k=1:length(template_diffs)
%                     RMSE_val=sum(abs(latDiffs-template_diffs(k,:)));
%                     if RMSE_val<minRMSE
%                         minRMSE=RMSE_val;
%                     end
%                 end
%                 RMSE_vals{2,j}(i,:)=minRMSE;
%             end
%             
%         end
%         %Remove bad event
%         if ind2remove<=length(latR) && ind2remove>1
%             latR=latR(setdiff(1:length(latR),ind2remove-1));
%             EEG.event(latR_eegEvInds(ind2remove-1))=[];
%         end
%         if any(RMSE_vals{2,j}>20)
%             error('Error still too high! Manually inspect');
%         end
% %         latL=latL-min(latL)+startInds(j,1);
% %         latR=latR-min(latR)+startInds(j,2);
    end
end

EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

%Quick tally of EEG events added
for j=1:length(trialsCompressed)
    countL=0; countR=0;
    for i=1:length({EEG.event(1:end).type})
        if strcmp(EEG.event(1,i).type,['Pull_Left_' trialsCompressed{j}])
            countL=countL+1;
        elseif strcmp(EEG.event(1,i).type,['Pull_Right_' trialsCompressed{j}])
            countR=countR+1;
        end
    end
    disp(['# Pull Left ' trialsCompressed{j} ': ' num2str(countL)]);
    disp(['# Pull Right ' trialsCompressed{j} ': ' num2str(countR)]);
end

%Save new EEG struct (save copy of original version with a different name)
% pop_saveset(EEG, 'filename', ['Merge_' group '_CAR'], 'filepath', EEGsets_outpath,'savemode','twofiles');