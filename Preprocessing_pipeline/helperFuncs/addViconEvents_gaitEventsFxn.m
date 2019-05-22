function [EEGout]=addViconEvents_gaitEventsFxn(pathname,subjnum,group,EEG)
    %Load EEG, Vicon events, and Vicon sync
    % pathname='/media/stepeter/New Volume/VR_rotations/Data/';
    % subjnum=3;
    Vicon_sample_rate=100; %Hz (100 for mocap; 1000 for pulls)
    sync_sample_rate=1000; %Hz
    syncFreq=0.5; %Hz
%     group='all';

    %Load data
%     cd(pathname);
    EEGsets_outpath=[pathname 'WMISM_' num2str(subjnum) filesep 'EEG_sets/'];
%     EEG=pop_loadset('filename', ['Merge_' group '_CAR.set'], 'filepath', EEGsets_outpath);
    load([EEGsets_outpath 'Events/sync.mat']);
    % load([EEGsets_outpath 'Events/finalPeaks.mat']);

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
            csv_sync_ch_idx=[9:17]; %[4:12]; 
            EEGTrialnames={'Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3'}; 
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
    
    if subjnum==29
        sync_Fall_EEG(1)=[];
    end
    
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
    
    %Subject-specific cleaning of sync breaks
    if subjnum==9
        eeg_inds_rise(21)=[]; %removes ind 2026
    elseif subjnum==11
        eeg_inds_fall(24)=[]; %removes ind 2213
        eeg_inds_fall(11)=[]; %removes ind 1086
    elseif subjnum==12
        eeg_inds_rise(18)=[]; %removes ind 1747
        eeg_inds_rise(13)=[]; %removes ind 1260
        eeg_inds_rise=sort([eeg_inds_rise 1411]); %add in nd 1411
    elseif subjnum==13
        eeg_inds_rise=[eeg_inds_rise 1392]; %add in ind 1392
        eeg_inds_rise=sort(eeg_inds_rise);
    elseif subjnum==20
        eeg_inds_fall(18)=[]; %removes ind 1739
    elseif subjnum==25
        eeg_inds_fall(19)=[]; %removes ind 1810
    elseif subjnum==28
        eeg_inds_rise(21)=[]; %removes ind 1947
        eeg_inds_fall(22)=[]; %removes ind 1952
        eeg_inds_fall(15)=[]; %removes ind 1365
    elseif subjnum==29
        eeg_inds_rise(19)=[]; %removes ind 1874
    elseif subjnum==31
        eeg_inds_rise=sort([eeg_inds_rise 359 360]); %add in dummy inds
        eeg_inds_fall=sort([eeg_inds_fall 359 360]); %add in dummy inds
        eeg_inds_rise(19)=[]; %removes ind 1577
    elseif subjnum==30
        eeg_inds_fall(13)=eeg_inds_fall(13)+1;
    elseif subjnum==1
        eeg_inds_rise=sort([eeg_inds_rise 1223 1022 1222]); %1283]); %1901]);
        eeg_inds_fall=sort([eeg_inds_fall 1223 1222]); %1283]); %1901]);
    end
    
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
    if any(csv_sync_ch_idx~=9:17)
        error('This function is specifically set up for training trials only!');
    end
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
            if subjnum==1 && any(csv_sync_ch_idx==9) && j==1
                ViconSyncInds{i,j}=[];
            else
                if (ViconSyncInds{i,j}(end)-ViconSyncInds{i,j}(end-1))<(sync_sample_rate/syncFreq-5) || (ViconSyncInds{i,j}(end)-ViconSyncInds{i,j}(end-1))>(sync_sample_rate/syncFreq+5)
                    ViconSyncInds{i,j}(end)=[];
                end
                if (ViconSyncInds{i,j}(2)-ViconSyncInds{i,j}(1))<(sync_sample_rate/syncFreq-5) || (ViconSyncInds{i,j}(2)-ViconSyncInds{i,j}(1))>(sync_sample_rate/syncFreq+5)
                    ViconSyncInds{i,j}(1)=[];
                end
            end
        end
    end

    %Check that number of rising and falling edges are equal between Vicon and EEG
%     totalEEGInds=[inds_Prior [length(sync_Rise_EEG); length(sync_Fall_EEG)]];
    q=0;
    for i=EEG_inds
        q=q+1;
        eegFinalSyncInds{1,q}=sync_Rise_EEG((eeg_inds_rise(i-1)+1):eeg_inds_rise(i))';
        eegFinalSyncInds{2,q}=sync_Fall_EEG((eeg_inds_fall(i-1)+1):eeg_inds_fall(i))';
    end
    %Removing from Vicon inds or dummy EEG inds so pass check (should only
    %remove from end of EEG inds here, otherwise need to change
    %eeg_inds_rise or eeg_inds_fall
    if subjnum==1
        eegFinalSyncInds{1,1}=[]; eegFinalSyncInds{2,1}=[];
        eegFinalSyncInds{1,8}(end-1:end)=[];
        eegFinalSyncInds{2,3}(end)=[]; eegFinalSyncInds{2,4}(end)=[]; 
        eegFinalSyncInds{2,8}(end-1:end)=[];
    elseif subjnum==4
        eegFinalSyncInds{1,8}(end)=[];
    elseif subjnum==5
        ViconSyncInds{2,6}(1)=[];
    elseif subjnum==8
        ViconSyncInds{2,6}(end)=[];
        ViconSyncInds{1,9}(end)=[];
        ViconSyncInds{2,9}(end-1:end)=[];
    elseif subjnum==9
        ViconSyncInds{1,8}(end)=[];
    elseif subjnum==10
        ViconSyncInds{2,9}(43)=[];
    elseif subjnum==12
        ViconSyncInds{1,5}(end)=[];
        ViconSyncInds{1,3}(end-1:end)=[];
        ViconSyncInds{2,3}(end)=[];
        eegFinalSyncInds{1,1}(end)=[];
    elseif subjnum==13
        eegFinalSyncInds{1,1}(end)=[];
    elseif subjnum==15
        ViconSyncInds{2,4}(19)=[];
        ViconSyncInds{2,7}(94)=[];
    elseif subjnum==16
        ViconSyncInds{2,3}(51)=[];
    elseif subjnum==20
        ViconSyncInds{2,5}(end)=[];
        eegFinalSyncInds{2,9}(end)=[];
    elseif subjnum==24
        ViconSyncInds{1,3}(end)=[];
        ViconSyncInds{2,3}(end)=[];
    elseif subjnum==25
        ViconSyncInds{2,6}(end)=[];
    elseif subjnum==28
        ViconSyncInds{1,8}(end)=[];
        ViconSyncInds{2,8}(end)=[];
        ViconSyncInds{2,2}(end)=[];
    elseif subjnum==29
        ViconSyncInds{2,3}(20)=[];
        eegFinalSyncInds{2,6}([97 99])=[];
        eegFinalSyncInds{1,6}(end)=[];
        eegFinalSyncInds{2,6}(end)=[];
        eegFinalSyncInds{1,9}(96:116)=[];
        eegFinalSyncInds{2,9}(95:112)=[];
    elseif subjnum==31
        ViconSyncInds{2,6}(8)=[];
        ViconSyncInds{1,6}(end)=[];
        ViconSyncInds{1,3}(16)=[];
    elseif subjnum==32
        ViconSyncInds{2,4}(56)=[];
    elseif subjnum==33
        ViconSyncInds{2,6}(80)=[];
    end
%     totalEEGInds=diff(totalEEGInds')'; notEq=zeros(2,size(ViconSyncInds,2));
    for i=1:2
        for j=1:size(ViconSyncInds,2)
%             totalViconInds(i,j)=length(ViconSyncInds{i,j});
            if length(ViconSyncInds{i,j})~=length(eegFinalSyncInds{i,j}) %totalEEGInds(i,EEG_inds(j))~=totalViconInds(i,j)
                error(['Sync signals do not line up for subject ' num2str(subjnum) ': ' num2str(i) ' and ' num2str(j) '!']);%notEq(i,j)=1;
            end
            %Convert Vicon sync to same sampling rate of Vicon signal (gait events
            %or pull)
            ViconSyncInds{i,j}=ViconSyncInds{i,j}*Vicon_sample_rate/sync_sample_rate;
        end
    end
    
    
    
    %Go through each Vicon trial and find closest sync to each and difference,
    %and add to EEG events
    N=length(EEG.event); q=1;
    for i=1:length(csv_sync_ch_idx)
        if subjnum==1 && any(csv_sync_ch_idx(i)==[9 2]) %isempty(eegFinalSyncInds{1,i}) || isempty(ViconSyncInds{1,i}) %
        elseif subjnum==31 && any(csv_sync_ch_idx(i)==[18:20])
        else
            %Load data from gait events
            load([pathname 'trainGait' filesep 'LFOOT_' num2str(i) '_' num2str(subjnum) '.mat']);
            load([pathname 'trainGait' filesep 'RFOOT_' num2str(i) '_' num2str(subjnum) '.mat']);

            %Left HS on beam (BC)
            for j=1:length(LFOOT.HS.ON)
                ev_Left=double(LFOOT.HS.ON(j));
                [Y_left_rise,I_left_rise]=min(abs(ViconSyncInds{1,i}-ev_Left));
                [Y_left_fall,I_left_fall]=min(abs(ViconSyncInds{2,i}-ev_Left));
        %         Y_left_rise=-(1-notEq(1,i))*Y_left_rise;
        %         Y_left_fall=-(1-notEq(2,i))*Y_left_fall;
                if Y_left_rise<=Y_left_fall
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{1,i}(I_left_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_left_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['BC_Left_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{2,i}(I_left_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_left_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['BC_Left_' EEGTrialnames{i}];   
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end

            %Left HS off beam (TC)
            for j=1:length(LFOOT.HS.OFF)
                ev_Left=double(LFOOT.HS.OFF(j));
                [Y_left_rise,I_left_rise]=min(abs(ViconSyncInds{1,i}-ev_Left));
                [Y_left_fall,I_left_fall]=min(abs(ViconSyncInds{2,i}-ev_Left));
        %         Y_left_rise=-(1-notEq(1,i))*Y_left_rise;
        %         Y_left_fall=-(1-notEq(2,i))*Y_left_fall;
                if Y_left_rise<=Y_left_fall
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{1,i}(I_left_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_left_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['TC_Left_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{2,i}(I_left_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_left_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['TC_Left_' EEGTrialnames{i}];   
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end

            %Right HS on beam (BC)
            for j=1:length(RFOOT.HS.ON)
                ev_Right=double(RFOOT.HS.ON(j));
                [Y_right_rise,I_right_rise]=min(abs(ViconSyncInds{1,i}-ev_Right));
                [Y_right_fall,I_right_fall]=min(abs(ViconSyncInds{2,i}-ev_Right));
        %         Y_right_rise=-(1-notEq(1,i))*Y_right_rise;
        %         Y_right_fall=-(1-notEq(2,i))*Y_right_fall;
                if Y_right_rise<=Y_right_fall
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{1,i}(I_right_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_right_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['BC_Right_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{2,i}(I_right_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_right_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['BC_Right_' EEGTrialnames{i}];
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end

            %Right HS off beam (TC)
            for j=1:length(RFOOT.HS.OFF)
                ev_Right=double(RFOOT.HS.OFF(j));
                [Y_right_rise,I_right_rise]=min(abs(ViconSyncInds{1,i}-ev_Right));
                [Y_right_fall,I_right_fall]=min(abs(ViconSyncInds{2,i}-ev_Right));
        %         Y_right_rise=-(1-notEq(1,i))*Y_right_rise;
        %         Y_right_fall=-(1-notEq(2,i))*Y_right_fall;
                if Y_right_rise<=Y_right_fall
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{1,i}(I_right_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_right_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['TC_Right_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{2,i}(I_right_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_right_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['TC_Right_' EEGTrialnames{i}];
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end
        end
    end

    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

    trialsCompressed=unique(EEGTrialnames); %unique(trials);
   

    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

    %Quick tally of EEG events added
    for j=1:length(trialsCompressed)
        countLon=0; countRon=0; countLoff=0; countRoff=0;
        for i=1:length({EEG.event(1:end).type})
            if strcmp(EEG.event(1,i).type,['BC_Left_' trialsCompressed{j}])
                countLon=countLon+1;
            elseif strcmp(EEG.event(1,i).type,['BC_Right_' trialsCompressed{j}])
                countRon=countRon+1;
            elseif strcmp(EEG.event(1,i).type,['TC_Left_' trialsCompressed{j}])
                countLoff=countLoff+1;
            elseif strcmp(EEG.event(1,i).type,['TC_Right_' trialsCompressed{j}])
                countRoff=countRoff+1;
            end
        end
        disp(['# BC Left ' trialsCompressed{j} ': ' num2str(countLon)]);
        disp(['# BC Right ' trialsCompressed{j} ': ' num2str(countRon)]);
        disp(['# TC Left ' trialsCompressed{j} ': ' num2str(countLoff)]);
        disp(['# TC Right ' trialsCompressed{j} ': ' num2str(countRoff)]);
    end
    
    EEGout=EEG;
    %Save new EEG struct (save copy of original version with a different name)
    % pop_saveset(EEG, 'filename', ['Merge_' group '_CAR'], 'filepath', EEGsets_outpath,'savemode','twofiles');
end