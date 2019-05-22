%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_rotationsDenoise/Code/EEG/'));
    startUpEEGLAB('close');
end

%Specify parameters
subjsToAnalyze=[1:13 15:17 19:20 22:33];
group='all'; %'mismatches';
epochTimes=[-.5 2];
trigPrefix='M_on_'; %event label prefix
trigPrefix1='M_on_CW_';
trigPrefix2='M_on_CCW_';
analyzeMismDirection=0; %1 - keep CW and CCW labels on events; 0 - lump together
analyzeFootDirection=1;
% Add EEGLAB path
RootData='/media/stepeter/Local_Data/VR_rotationsDenoise/Data/';

%Specify whether to run DIPFIT, check DIPFIT, and/or epoch
runDip=0; %1 - run DIPFIT, 0 - no DIPFIT
checkDip=0; %1 - check all results (will break to manually allow this), 0 - don't check results
runEpoch=1; %1 - epoch data, 0 - no epoching
% runTimewarp=1; %1 - timewarp data, 0 - no timewarping
addgaitEvsEpch=0; %1 - add in gait events, epoch, and TW, 0 - don't do this

%%Run DIPFIT (automatic, using DIPFIT2 EEGLAB plugin)
if runDip==1
    for frodo=subjsToAnalyze
        run_DIPFIT(frodo,RootData,group);
    end
    disp('Running DIPFIT completed!');
end

%%Check DIPFIT (manual rejection of non-neural components)
if checkDip==1
    disp('Check DIPFIT by running cells in checkDIPFIT.m (Change subjnum!)');
    edit checkDIPFIT
    return;
%     for frodo=subjsToAnalyze
%         checkDIPFIT(frodo,group,RootData)
%     end
%     disp('Checking DIPFIT completed!');
end

%%Epoch Data (automatic)
if runEpoch==1
    for frodo=subjsToAnalyze
        epochDIPFITData(RootData,frodo,group,trigPrefix,analyzeMismDirection,epochTimes,trigPrefix1,trigPrefix2);
    end
    disp('Epoching data completed!');
end

%%Add gait events and epoch (extraneous code for epoching gait events)
condsTW={'Tr1','Tr2','Tr3'}; %,'Tr3'}; %
twSeq={'BC_','BC_','TC_','TC_','TC_'}; %{'BC_','BC_','TC_','TC_','TC_'}; %{'TC_','TC_','BC_','BC_'}; %{'BC_','BC_','BC_'}; %
folderSuffix=[group '_gaitOFFall_v2Denoise']; allEpchs=zeros(length(subjsToAnalyze),length(condsTW)); q=0;
if addgaitEvsEpch==1 && runEpoch==0
    for frodo=subjsToAnalyze
        EEG=[]; ALLEEG=[];
        [EEGout]=addViconEvents_gaitEventsFxn(RootData,frodo,group);
        %Pull out training and pre/post trials
        segmentBounds='boundary'; EEGin=EEGout;
        if ischar(segmentBounds) && strcmp(segmentBounds,'boundary')
            %Get boundary events
            boundary_latencies=[];
            for i=1:length(EEGin.event)
                if strcmp(EEGin.event(i).type,'boundary')
                    boundary_latencies=[boundary_latencies EEGin.event(i).latency];
                end
            end
            %     %Special condition for subject data
            %     if strcmp(EEGin.subject,'WMISM_1')
            %         boundary_latencies(end-1)=[];
            %     end
            boundary_latencies2=[1 boundary_latencies EEGin.pnts];
        else
            boundary_latencies2=segmentBounds;
        end
        EEGout2=pop_select(EEGin,'point',boundary_latencies2([5 10]));
        evs=EEGout2.event; urevs=EEGout2.urevent; EEGout2=[];
        EEGsets_outpath=[RootData 'WMISM_' num2str(frodo) filesep 'EEG_sets/' filesep 'v2'];
        EEGout=pop_loadset('filename', ['Merge_' group '_CAR_v2.set'], 'filepath', EEGsets_outpath);
        EEGout.icachansind=[]; EEGout.event=evs; EEGout.urevent=urevs;
%         EEG=EEGout;
        
        %Only keep first and last 5 minutes
%         maxLat=0; minLat=1e10;
%         for i=1:length(EEGout.event)
%             if strcmp(EEGout.event(i).type(1:3),'BC_') || strcmp(EEGout.event(i).type(1:3),'TC_')
%                 if EEGout.event(i).latency<minLat
%                     minLat=EEGout.event(i).latency;
%                 end
%                 if EEGout.event(i).latency>maxLat
%                     maxLat=EEGout.event(i).latency;
%                 end
%             end
%         end
%         maxLat=maxLat+128; minLat=minLat-128; %readjust for removing frames
%         EEGout=pop_select(EEGout,'point',[minLat (76799+minLat); (maxLat-76799) maxLat]);
        
%         [numEpchs]=epochDIPFITData_gaitEvents(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
        [numEpchs]=epochDIPFITData_gaitEvents_offBeam(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
%         [numEpchs]=epochDIPFITData_gaitEvents_onBeam(EEGout,RootData,frodo,group,analyzeFootDirection,epochTimes,twSeq,condsTW,folderSuffix);
        q=q+1; allEpchs(q,:)=numEpchs;
    end
    disp('Epoching gait data completed!');
end
disp('Step 3 Done :)');
disp('Data has been run through DIPFIT, checked, epoched, and/or timewarped.');
disp('Proceed to step 4 to create a study from multiple subjects.');
