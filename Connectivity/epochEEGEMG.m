%%Epoch datasets for group SIFT and add in EMG data
savePath='/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/'; %GroupConn_aveScalp/newData_noEMG/';
RootData='/usr/local/VR_connectivity/Data/'; 
group='all';
epochTimes=[-1 2];
trigPrefix='M_on_'; %event label prefix
trigPrefix1='M_on_CW_';
trigPrefix2='M_on_CCW_';
load('/usr/local/VR_connectivity/Data/GroupConn_aveScalp/subjectInds.mat');
%subjectInds
subjectInds=[1:7 9:13 15:17 19:20 22:33];
for i=subjectInds
    [EEGout]=addViconEvents_pullEventsFxn_aveScalp(RootData,i,group,savePath); %get pull events
    
    %Create pull sets (with EMG added)
    trigs={'L_pull_Stn','R_pull_Stn'};
    if length(trigs)<=1
        analyzePullDirection=0;
    else
        analyzePullDirection=1;
    end
    epochDIPFITDataPulls_aveScalp(EEGout,trigs,RootData,i,group,epochTimes,analyzePullDirection);
    
    trigs={'L_pull_Wlk','R_pull_Wlk'};
    if length(trigs)<=1
        analyzePullDirection=0;
    else
        analyzePullDirection=1;
    end
    epochDIPFITDataPulls_aveScalp(EEGout,trigs,RootData,i,group,epochTimes,analyzePullDirection);
    
    %Create viz sets (with EMG added)
    trigs={'M_on_CW_SVZ','M_on_CCW_SVZ'};
    if length(trigs)<=1
        analyzeMismDirection=0;
    else
        analyzeMismDirection=1;
    end
    epochDIPFITDataViz_aveScalp(EEGout,trigs,i,group,trigPrefix,analyzeMismDirection,epochTimes,trigPrefix1,trigPrefix2)
    
    trigs={'M_on_CW_WVZ','M_on_CCW_WVZ'};
    if length(trigs)<=1
        analyzeMismDirection=0;
    else
        analyzeMismDirection=1;
    end
    epochDIPFITDataViz_aveScalp(EEGout,trigs,i,group,trigPrefix,analyzeMismDirection,epochTimes,trigPrefix1,trigPrefix2)
    
    
%     %Do gait events (HS and TO)
%     [EEGout]=addViconEvents_walkPullGaitEventsFxn_aveScalp(RootData,i,group,savePath);
%     
%     %Create gait sets (with EMG added)
%     trigs={'LHS','RHS'}; %'L_pull_Stn','R_pull_Stn' %'pull_Stn'
%     if length(trigs)<=1
%         analyzePullDirection=0;
%     else
%         analyzePullDirection=1;
%     end
%     epochDIPFITDataWalkPullGait_aveScalp(EEGout,trigs,RootData,i,group,epochTimes,analyzePullDirection);
%     
%     trigs={'LTO','RTO'}; %'L_pull_Wlk','R_pull_Wlk' %'pull_Wlk'
%     if length(trigs)<=1
%         analyzePullDirection=0;
%     else
%         analyzePullDirection=1;
%     end
%     epochDIPFITDataWalkPullGait_aveScalp(EEGout,trigs,RootData,i,group,epochTimes,analyzePullDirection);
end