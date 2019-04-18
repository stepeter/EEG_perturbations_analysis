%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_connectivity/Code/EEGconn/'));
    startUpEEGLAB('close');
end

%Fieldtrip will mess DIPFIT up, so make sure its removed from path
rmpath(genpath('/home/stepeter/Documents/eeglab13_0_1b_octave_rem/external/fieldtrip-partial'));

%Specify parameters
subjsToAnalyze=[1:13 15:17 19:20 22:33];
group='all'; %'mismatches';
epochTimes=[-1 2]; %[-4 4]; %2];  %-0.3 0.8
trigPrefix='M_on_'; %event label prefix
trigPrefix1='M_on_CW_';
trigPrefix2='M_on_CCW_';
analyzeMismDirection=0; %1; %1 - keep CW and CCW labels on events; 0 - lump together
analyzeFootDirection=1;
trigs={'M_on_WVZ'}; %'M_on_CW_SVZ','M_on_CCW_SVZ','M_on_CW_WVZ','M_on_CCW_WVZ'}; %{'M_on_SVZ','M_on_WVZ'};
% Add EEGLAB path
RootData='/media/stepeter/Local_Data/VR_connectivity/Data/';

%Specify whether to run DIPFIT, check DIPFIT, and/or epoch
runDip=0; %1 - run DIPFIT, 0 - no DIPFIT
checkDip=0; %1 - check all results (will break to manually allow this), 0 - don't check results
runEpochMism=0; %1 - epoch data, 0 - no epoching
runEpochPull=1; %1 - epoch data, 0 - no epoching
% runTimewarp=1; %1 - timewarp data, 0 - no timewarping
addgaitEvsEpch=0; %1 - add in gait events, epoch, and TW, 0 - don't do this

%%Run DIPFIT
if runDip==1
    for frodo=subjsToAnalyze
        run_DIPFIT(frodo,RootData,group);
    end
    disp('Running DIPFIT completed!');
end

%%Check DIPFIT
if checkDip==1
    disp('Check DIPFIT by running cells in checkDIPFIT.m (Change subjnum!)');
    edit checkDIPFIT
    return;
%     for frodo=subjsToAnalyze
%         checkDIPFIT(frodo,group,RootData)
%     end
%     disp('Checking DIPFIT completed!');
end

%%Epoch Data for pulls
if runEpochPull==1
    for frodo=subjsToAnalyze
        trigs={'pull_Stn','pull_Wlk'}; %'pull_Wlk'}; %'L_pull_Stn','R_pull_Stn','L_pull_Wlk','R_pull_Wlk'}; %
        if length(trigs)<=2
            analyzePullDirection=0;
        else
            analyzePullDirection=1;
        end
        savePathSuffix='pullsOnset_noDir_v2'; %'Walk_4sift'; %
        [EEGout]=addViconEvents_pullEventsFxn(RootData,frodo,group);
        epochDIPFITDataPulls(EEGout,trigs,RootData,frodo,group,epochTimes,savePathSuffix,analyzePullDirection);
    end
    disp('Epoching data completed!');
end

%%Epoch Data for mismatches
if runEpochMism==1
    for frodo=subjsToAnalyze
        savePathSuffix='WVZ_4sift';
        epochDIPFITData(trigs,RootData,frodo,group,trigPrefix,analyzeMismDirection,epochTimes,trigPrefix1,trigPrefix2,savePathSuffix);
    end
    disp('Epoching data completed!');
end

%%Add gait events and epoch
condsTW={'Tr1','Tr3'}; %,'Tr2','Tr3'};
twSeq={'TC_','TC_','BC_','BC_'}; %{'BC_','BC_','TC_','TC_','TC_'}; %{'TC_','TC_','BC_','BC_'}; %{'BC_','BC_','BC_'}; %
folderSuffix=[group '_gaitON_justTrain']; allEpchs=zeros(length(subjsToAnalyze),length(condsTW)); q=0;
if addgaitEvsEpch==1 && runEpochMism==0
    for frodo=subjsToAnalyze
        EEG=[]; ALLEEG=[];
        [EEGout]=addViconEvents_gaitEventsFxn(RootData,frodo,group);
        
        %Only keep first and last 5 minutes
        maxLat=0; minLat=1e10;
        for i=1:length(EEGout.event)
            if strcmp(EEGout.event(i).type(1:3),'BC_') || strcmp(EEGout.event(i).type(1:3),'TC_')
                if EEGout.event(i).latency<minLat
                    minLat=EEGout.event(i).latency;
                end
                if EEGout.event(i).latency>maxLat
                    maxLat=EEGout.event(i).latency;
                end
            end
        end
        maxLat=maxLat+128; minLat=minLat-128; %readjust for removing frames
        EEGout=pop_select(EEGout,'point',[minLat (76799+minLat); (maxLat-76799) maxLat]);
        
%         [numEpchs]=epochDIPFITData_gaitEvents(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
%         [numEpchs]=epochDIPFITData_gaitEvents_offBeam(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
        [numEpchs]=epochDIPFITData_gaitEvents_onBeam(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
        q=q+1; allEpchs(q,:)=numEpchs;
    end
    disp('Epoching gait data completed!');
end
disp('Step 3 Done :)');
disp('Data has been run through DIPFIT, checked, epoched, and/or timewarped.');
disp('Proceed to step 4 to create a study from multiple subjects.');
