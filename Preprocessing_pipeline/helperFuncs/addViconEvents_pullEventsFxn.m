function [EEGout]=addViconEvents_pullEventsFxn(pathname,subjnum,group)
    %Load EEG, Vicon events, and Vicon sync
    % pathname='/media/stepeter/New Volume/VR_rotations/Data/';
    % subjnum=3;
    Vicon_sample_rate=1000; %Hz (100 for mocap; 1000 for pulls)
    sync_sample_rate=1000; %Hz
    syncFreq=0.5; %Hz
%     group='all';

    %Load data
    cd(pathname);
    EEGsets_outpath=[pathname 'WMISM_' num2str(subjnum) filesep 'EEG_sets/'];
    EEG=pop_loadset('filename', ['Merge_' group '_CAR_v2.set'], 'filepath', EEGsets_outpath); 
    load([EEGsets_outpath 'Events/sync.mat']);
    load([EEGsets_outpath 'Events/finalPeaksOnset.mat']); %finalPeaksOnset - inds at onset of peak; finalPeaks - inds at peak

    switch group
        case 'mismatches'
            error('No pull perturbations for these trials!');
        case 'pull'
            csv_sync_ch_idx=[1:3 13:15];
%             break_inds=[3 6];
            trials={'Stn','Stn','Stn','Wlk','Wlk','Wlk'};
            %Vicon file should be ordered: Stand, Walk
        case 'train'
            csv_sync_ch_idx=[4:12];
%             break_inds=[4 7 10]; 
            trials={'Pre','Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3','Post'};
            %Vicon file should be ordered: Pre, Train1, Train2, Train3, Post
        case 'all'
            csv_sync_ch_idx=[3:5 18:20 6:8  21:23];%[9:17]; %[4:12]; 
            EEGTrialnames={'Stn','Stn','Stn','Wlk','Wlk','Wlk','SVZ','SVZ','SVZ','WVZ','WVZ','WVZ'}; %{'Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3'}; 
            EEG_inds=1:12; % %order of trials in EEG
            trials={'Stn','Stn','Stn','Wlk','Wlk','Wlk','SVZ','SVZ','SVZ','WVZ','WVZ','WVZ'}; %,'Pre','Tr1','Tr1','Tr1','Tr2','Tr2','Tr2','Tr3','Tr3','Tr3','Post'};
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
%     if subjnum==9
%         eeg_inds_rise(21)=[]; %removes ind 2026
%     elseif subjnum==11
%         eeg_inds_fall(24)=[]; %removes ind 2213
%         eeg_inds_fall(11)=[]; %removes ind 1086
%     elseif subjnum==12
%         eeg_inds_rise(18)=[]; %removes ind 1747
%         eeg_inds_rise(13)=[]; %removes ind 1260
%         eeg_inds_rise=sort([eeg_inds_rise 1411]); %add in nd 1411
%     elseif subjnum==13
%         eeg_inds_rise=[eeg_inds_rise 1392]; %add in ind 1392
%         eeg_inds_rise=sort(eeg_inds_rise);
%     elseif subjnum==20
%         eeg_inds_fall(18)=[]; %removes ind 1739
%     elseif subjnum==25
%         eeg_inds_fall(19)=[]; %removes ind 1810
%     elseif subjnum==28
%         eeg_inds_rise(21)=[]; %removes ind 1947
%         eeg_inds_fall(22)=[]; %removes ind 1952
%         eeg_inds_fall(15)=[]; %removes ind 1365
%     elseif subjnum==29
%         eeg_inds_rise(19)=[]; %removes ind 1874
%     elseif subjnum==31
%         eeg_inds_rise=sort([eeg_inds_rise 359 360]); %add in dummy inds
%         eeg_inds_fall=sort([eeg_inds_fall 359 360]); %add in dummy inds
%         eeg_inds_rise(19)=[]; %removes ind 1577
%     elseif subjnum==30
%         eeg_inds_fall(13)=eeg_inds_fall(13)+1;
%     elseif subjnum==1
%         eeg_inds_rise=sort([eeg_inds_rise 1223 1022 1222]); %1283]); %1901]);
%         eeg_inds_fall=sort([eeg_inds_fall 1223 1222]); %1283]); %1901]);
%         trials(end-1)=[]; %EEG_inds(end)=[];
%     end

    if subjnum==1
        eeg_inds_rise=sort([eeg_inds_rise 1023]);
    elseif subjnum==11
        eeg_inds_fall(end-1)=[];
    elseif subjnum==31
        eeg_inds_rise=sort([eeg_inds_rise 359 360]); %add in dummy inds
        eeg_inds_fall=sort([eeg_inds_fall 359 360]); %add in dummy inds
    end
    
%     figure; subplot(2,1,1);
%     Dsync_Rise_EEG=diff(sync_Rise_EEG); Dsync_Fall_EEG=diff(sync_Fall_EEG);
%     hold on; plot(Dsync_Rise_EEG); scatter(eeg_inds_rise,Dsync_Rise_EEG(eeg_inds_rise),'r'); hold off;
%     subplot(2,1,2);
%     hold on; plot(Dsync_Fall_EEG); scatter(eeg_inds_fall,Dsync_Fall_EEG(eeg_inds_fall),'r'); hold off;
    
    if (length(eeg_inds_rise)+1)~=length(trials) %length(csv_sync_ch_idx)
        error('EEG rise and CSV breaks are not equal!');
    end
    if (length(eeg_inds_fall)+1)~=length(trials) %length(csv_sync_ch_idx)
        error('EEG fall and CSV breaks are not equal!');
    end

    

    %If no errors, find Vicon sync rising and falling edges
%     if any(csv_sync_ch_idx~=[3:8 18:23]) %9:17)
%         error('This function is specifically set up for stand/walk trials only!');
%     end
    
%     if subjnum==1
%         csv_sync_ch_idx(end)=[];
%     end
    
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
        if i==1
            eegFinalSyncInds{1,q}=sync_Rise_EEG(1:eeg_inds_rise(i))';
            eegFinalSyncInds{2,q}=sync_Fall_EEG(1:eeg_inds_fall(i))';
        elseif i==12
            eegFinalSyncInds{1,q}=sync_Rise_EEG((eeg_inds_rise(i-1)+1):end)';
            eegFinalSyncInds{2,q}=sync_Fall_EEG((eeg_inds_fall(i-1)+1):end)';
        else
            eegFinalSyncInds{1,q}=sync_Rise_EEG((eeg_inds_rise(i-1)+1):eeg_inds_rise(i))';
            eegFinalSyncInds{2,q}=sync_Fall_EEG((eeg_inds_fall(i-1)+1):eeg_inds_fall(i))';
        end
    end
    
    %Removing from Vicon inds or dummy EEG inds so pass check (should only
    %remove from end of EEG inds here, otherwise need to change
    %eeg_inds_rise or eeg_inds_fall)
    if subjnum==4
        ViconSyncInds{2,9}(42)=[]; ViconSyncInds{2,6}(79)=[];
    elseif subjnum==1
        eegFinalSyncInds{1,10}(end)=[]; eegFinalSyncInds{1,11}(1)=[]; %Removed first event based on rise/fall order
        eegFinalSyncInds{2,12}(end)=[];
        eeg_inds_rise(10)=eeg_inds_rise(10)+1;
    elseif subjnum==3
        ViconSyncInds{1,1}(87)=[];
    elseif subjnum==5
        ViconSyncInds{1,4}(end)=[];
        eegFinalSyncInds{1,6}(45:end)=[];
        eegFinalSyncInds{2,6}(45:end)=[];
    elseif subjnum==6
        ViconSyncInds{2,10}(61)=[];
        ViconSyncInds{1,12}(49)=[];
        eegFinalSyncInds{1,10}(1)=[]; %Removed first event based on rise/fall order
        eeg_inds_rise(9)=eeg_inds_rise(9)+1;
    elseif subjnum==7
        eegFinalSyncInds{1,4}(1)=[]; %Removed first event based on rise/fall order
        eeg_inds_rise(3)=eeg_inds_rise(3)+1;
    elseif subjnum==8
        ViconSyncInds{2,1}(1)=[];
        eegFinalSyncInds{2,4}(end)=[];
    elseif subjnum==11
        eegFinalSyncInds{2,11}=sort([eegFinalSyncInds{2,11}; 583245.5]);
    elseif subjnum==12
        eegFinalSyncInds{2,12}(end)=[];
    elseif subjnum==13
        ViconSyncInds{2,7}(1)=[];
        ViconSyncInds{2,9}(46)=[];
        ViconSyncInds{2,11}(46)=[];
    elseif subjnum==17
        ViconSyncInds{2,2}(77)=[];
    elseif subjnum==24
        eegFinalSyncInds{1,7}(1)=[];
        eeg_inds_rise(6)=eeg_inds_rise(6)+1;
        ViconSyncInds{1,12}(end)=[];
        ViconSyncInds{2,12}(end)=[];
    elseif subjnum==25
        ViconSyncInds{2,12}(53)=[];
        eegFinalSyncInds{2,10}(116)=[];
    elseif subjnum==26
        ViconSyncInds{2,10}(60)=[];
    elseif subjnum==28
        ViconSyncInds{2,7}(10)=[];
    elseif subjnum==29
        eegFinalSyncInds{2,3}(98:end)=[];
    elseif subjnum==30
        eegFinalSyncInds{1,7}(1)=[];
        eeg_inds_rise(6)=eeg_inds_rise(6)+1;
    elseif subjnum==32
        ViconSyncInds{2,7}(24)=[];
        ViconSyncInds{2,7}(1)=[];
    elseif subjnum==33
        ViconSyncInds{1,8}(18)=[];
    elseif subjnum==31
        %Trials 5 and 6 don't exist in EEG data
        ViconSyncInds{1,4}(64:end)=[];
        ViconSyncInds{2,4}(64:end)=[];
        ViconSyncInds{1,5}=[]; eegFinalSyncInds{1,5}=[];
        ViconSyncInds{2,5}=[]; eegFinalSyncInds{2,5}=[];
        ViconSyncInds{1,6}=[]; eegFinalSyncInds{1,6}=[];
        ViconSyncInds{2,6}=[]; eegFinalSyncInds{2,6}=[];
        finalPeaks{3,25}(19:end)=[]; finalPeaks{3,26}(13:end)=[]; %remove excess events not in EEG data
        finalPeaks{3,27}=[]; finalPeaks{3,28}=[]; finalPeaks{3,29}=[]; finalPeaks{3,30}=[];
    end
%     ViconSyncInds
%     eegFinalSyncInds
    
    inds_Prior=zeros(2,length(csv_sync_ch_idx));
    for i=1:length(eeg_inds_rise)
        inds_Prior(1,i+1)=eeg_inds_rise(i);
        inds_Prior(2,i+1)=eeg_inds_fall(i);
    end
    teve=1;
    
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
    pullIndsL=[1:2:6 25:2:30]; pullIndsR=[2:2:6 26:2:30];
    N=length(EEG.event); q=1;
    for i=1:6 %length(csv_sync_ch_idx) %only care about pull trials
        if subjnum==31 && any(i==5:6)
        else
            %Load data from peak events
%             load([pathname 'trainGait' filesep 'LFOOT_' num2str(i) '_' num2str(subjnum) '.mat']);
%             load([EEGsets_outpath 'Events/finalPeaks.mat']);
            
            %Left pull event
            for j=1:length(finalPeaks{3,pullIndsL(i)})
                ev_Left=double(finalPeaks{3,pullIndsL(i)}(j));
                [Y_left_rise,I_left_rise]=min(abs(ViconSyncInds{1,i}-ev_Left));
                [Y_left_fall,I_left_fall]=min(abs(ViconSyncInds{2,i}-ev_Left));
                if Y_left_rise<=Y_left_fall
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{1,i}(I_left_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_left_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['L_pull_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Left=ev_Left-ViconSyncInds{2,i}(I_left_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_left_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Left/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['L_pull_' EEGTrialnames{i}];   
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end
            
            
            %Right pull event
            for j=1:length(finalPeaks{3,pullIndsR(i)})
                ev_Right=double(finalPeaks{3,pullIndsR(i)}(j));
                [Y_right_rise,I_right_rise]=min(abs(ViconSyncInds{1,i}-ev_Right));
                [Y_right_fall,I_right_fall]=min(abs(ViconSyncInds{2,i}-ev_Right));
                if Y_right_rise<=Y_right_fall
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{1,i}(I_right_rise); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Rise_EEG(I_right_rise+inds_Prior(1,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['R_pull_' EEGTrialnames{i}];
                else
                    ev_diff_Vicon_Right=ev_Right-ViconSyncInds{2,i}(I_right_fall); %difference btwn sync and event in Vicon
                    EEG.event(1,N+q).latency=sync_Fall_EEG(I_right_fall+inds_Prior(2,EEG_inds(i)))+...
                        (ev_diff_Vicon_Right/Vicon_sample_rate)*EEG.srate;
                    EEG.event(1,N+q).type=['R_pull_' EEGTrialnames{i}];   
                end
                if EEG.event(1,N+q).latency<1 || EEG.event(1,N+q).latency>EEG.pnts
                    error('Event is out of bounds!');
                end
                q=q+1;
            end

        end
    end

    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

    trialsCompressed=unique(EEGTrialnames(1:6)); %unique(trials);
   

    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset

    %Quick tally of EEG events added
    for j=1:length(trialsCompressed)
        countpullL=0; countpullR=0;
        for i=1:length({EEG.event(1:end).type})
            if strcmp(EEG.event(1,i).type,['L_pull_' trialsCompressed{j}])
                countpullL=countpullL+1;
            elseif strcmp(EEG.event(1,i).type,['R_pull_' trialsCompressed{j}])
                countpullR=countpullR+1;
            end
        end
        disp(['# Left pulls ' trialsCompressed{j} ': ' num2str(countpullL)]);
        disp(['# Right pulls ' trialsCompressed{j} ': ' num2str(countpullR)]);
    end
    
    EEGout=EEG;
    %Save new EEG struct (save copy of original version with a different name)
    % pop_saveset(EEG, 'filename', ['Merge_' group '_CAR'], 'filepath', EEGsets_outpath,'savemode','twofiles');
end