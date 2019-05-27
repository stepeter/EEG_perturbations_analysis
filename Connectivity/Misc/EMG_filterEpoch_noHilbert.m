%% EMG Processing for connectivity analysis
clear;
load('/usr/local/VR_connectivity/Data/GroupConn_aveScalp/subjectInds.mat');
subjectInds=setdiff([1:7 9:13 15:17 19:20 22:33],[1:4 6:7 9:10 12:13 15 17 23 25:28 30 33]);
for subjnum=subjectInds %[1:13 15:17 19:20 22:33] %
    emgloadpath='/usr/local/VR_connectivity/Data/EMG/';
    signalSrate=1000; %Hz (1000 for EMG, 100 for mocap)
    savepath=[emgloadpath 'EMG4SIFT_revisions/'];

    %Start up EEGLAB if not running already
    if ~exist('ALLCOM')
        addpath(genpath('/usr/local/VR_connectivity/Code/EEGconn/'));
        startUpEEGLAB('close');
    end
    %% Load in data (already detrended)
    load([emgloadpath 'finalEMG_S' num2str(subjnum) '.mat']);

    %% 1) 1 Hz high-pass filter & full rectification (for EEG-EMG coherence)
    finalEMG_4EEG=finalEMG;
    for i=1:size(finalEMG_4EEG,2)
        %Use dummy EEG for filtering
        EEGdummy=EEG;
        EEGdummy.srate=signalSrate;
        EEGdummy.data=finalEMG_4EEG{2,i}'; EEGdummy=eeg_checkset(EEGdummy);

        %High-pass filter at 1 Hz
        EEGdummy = pop_eegfiltnew(EEGdummy,1,0);
        EEGdummy.data=abs(EEGdummy.data);
%         X=hilbert(EEGdummy.data')';
%         for jj=1:size(EEGdummy.data,1)
%             EEGdummy.data(jj,:)=cos(angle(X(jj,:)'));
%         end
        finalEMG_4EEG{2,i}=EEGdummy.data'; %abs(EEGdummy.data');
    end

%     %% 2) 20 Hz high-pass filter & full recification (for EMG analysis)
%     %Set filter values
%     lo_cutoff=20; %(Hz)
%     [b,a]=butter(4,lo_cutoff/(signalSrate/2),'high');
%     for i=1:size(finalEMG,2)
%         for j=1:size(finalEMG{2,1},2)
%             finalEMG{2,i}(:,j)=abs(filtfilt(b,a,double(finalEMG{2,i}(:,j)))); %zero-phase filtering
%         end
%     end
% %     finalEMG_4EEG=finalEMG;
    %% Find events
    EEGpath='/usr/local/VR_connectivity/Data/';
    [finalEvLatenciesEMG]=find_pullVizPerturbationsEvents(EEGpath,subjnum,'all',0);
%     [finalEvLatenciesMocap]=find_pullVizPerturbationsEvents(EEGpath,subjnum,'all',1);
    
    %Add two more instances of Walk Pull
    finalEMG_4EEG(:,13:15)=finalEMG_4EEG(:,4:6);
    finalEMG_4EEG(:,16:18)=finalEMG_4EEG(:,4:6);

    %% Epoch datasets
    epochTimes=[-1 2];
    EpPre = epochTimes(1,1)+(-0.6); %(-0.559);
    EpPost = epochTimes(1,2)+0.6; %0.559;
    ptsBelow=floor(EpPre*signalSrate);
    ptsAbove=ceil(EpPost*signalSrate);

    %Epoch 1)
    finalEMG_4EEGepoch=cell(1,6); badEpochInds=cell(1,6); oldCondInd=0;
    for i=1:length(finalEvLatenciesEMG)
        condInd=floor((i-1)/3)+1; 
        if condInd~=oldCondInd
            epCounter=0;
        end
        for j=1:length(finalEvLatenciesEMG{i})
            if ~isnan(finalEvLatenciesEMG{i}(j))
                epCounter=epCounter+1;
                lowVal=ptsBelow+round(finalEvLatenciesEMG{i}(j));
                hiVal=ptsAbove+round(finalEvLatenciesEMG{i}(j));
                if lowVal<1 || hiVal>length(finalEMG_4EEG{2,i})
                    badEpochInds{condInd}=[badEpochInds{condInd} epCounter];
                else
                    finalEMG_4EEGepoch{condInd}=cat(3,finalEMG_4EEGepoch{condInd},finalEMG_4EEG{2,i}(lowVal:hiVal,:));
                end
            end
        end
        oldCondInd=condInd;
    end
    badEpochInds
    
%     emgSrate=1000; %Hz; sampling rate of EMG signal
% %     EEG = pop_importdata('dataformat','array','nbchan',0,'data','sig1','setname','dummyEMGset','srate',emgSrate,'pnts',0,'xmin',0);
%     
%     epochTimes=[-0.5 1.5];
%     EpPre = epochTimes(1,1); %+(-0.6); %(-0.559);
%     EpPost = epochTimes(1,2); %+0.6; %0.559;
%     ptsBelow=floor(EpPre*signalSrate);
%     ptsAbove=ceil(EpPost*signalSrate);
%     %Epoch 2)
% %     finalEMG(:,13:15)=finalEMG(:,4:6);
% %     finalEMG(:,16:18)=finalEMG(:,4:6);
%     finalEMGepoch=cell(1,6); badEpochInds_noEEG=cell(1,6); oldCondInd=0;
%     finalEvLatenciesEMG(13:end)=[];
%     
%     %Load baseline walking values
%     load(['/media/stepeter/Local_Data/VR_connectivity/Data/EMG/pre_EMG_baseline/baseEMGPeakVals_S' num2str(subjnum) '.mat']);
%     baseVals=maxBaseEMGVals/1000; %keep in Volts
%     
%     %Divide out baseline value
%     for i=1:size(finalEMG,2)
%         for j=1:size(finalEMG{2,i}(:,:),2)
%             finalEMG{2,i}(:,j)=finalEMG{2,i}(:,j)/baseVals(j);
%         end
%     end
%     
%     for i=1:length(finalEvLatenciesEMG)
%         condInd=floor((i-1)/3)+1; 
%         if condInd~=oldCondInd
%             epCounter=0;
%         end
%         for j=1:length(finalEvLatenciesEMG{i})
%             if ~isnan(finalEvLatenciesEMG{i}(j))
%                 epCounter=epCounter+1;
%                 lowVal=ptsBelow+round(finalEvLatenciesEMG{i}(j));
%                 hiVal=ptsAbove+round(finalEvLatenciesEMG{i}(j));
%                 if lowVal<1 || hiVal>length(finalEMG{2,i})
%                     badEpochInds_noEEG{condInd}=[badEpochInds_noEEG{condInd} epCounter];
%                 else
%                     finalEMGepoch{condInd}=cat(3,finalEMGepoch{condInd},finalEMG{2,i}(lowVal:hiVal,:));
%                 end
%             end
%         end
%         oldCondInd=condInd;
%     end
%     badEpochInds_noEEG
%     
%     if subjnum==1
%         for i=1:length(finalEMGepoch)
%             meanfinalEMGepoch{i}=mean(finalEMGepoch{i},3);
%         end
%     else
%         for i=1:length(finalEMGepoch)
%             meanfinalEMGepoch{i}=cat(3,meanfinalEMGepoch{i},mean(finalEMGepoch{i},3));
%         end
%     end
    %% Save results
    %save([savepath 'EMGanalysis_S' num2str(subjnum) '.mat'],'finalEMG');
    save([savepath 'EMGforEEGcohere_S' num2str(subjnum) '.mat'],'finalEMG_4EEGepoch','badEpochInds');
end
% % save([emgloadpath 'justEMGnewLatencies/' 'meanfinalEMGepoch_20HPF.mat'],'meanfinalEMGepoch');
% %% Remove baseline for EMG plots
% baseIdx=1:500; meanfinalEMGepoch_baseRem=meanfinalEMGepoch;
% for j=1:4
%     baseVal=mean(meanfinalEMGepoch{j}(baseIdx,:,:));
%     meanfinalEMGepoch_baseRem{j}=meanfinalEMGepoch{j}-repmat(baseVal,[size(meanfinalEMGepoch{j},1) 1 1]);
% end
% 
% 
% %% Plot mean results
% % condNames={'Stn Pull','Wlk Pull','Stn Viz','Wlk Viz'};
% % for cond=1:4; %1-Stn Pull, 2-Wlk Pull, 3-Stn Viz, 4-Wlk Viz
% %     figure; muscles={'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};
% %     for i=1:8
% %         subplot(4,2,i);
% %         plot(-0.5:0.001:1.5,mean(squeeze(meanfinalEMGepoch_baseRem{cond}(:,i,:)),2)*1000); %multiply by 1000 to convert to mV
% %         title(muscles{i}); ylim([0 20]);
% %     end
% %     suptitle(condNames{cond});
% % end
% 
% %% Plot mean results (all together)
% figure; muscles={'LTA','LSOL','LMG','LPL','RTA','RSOL','RMG','RPL'};
% colsMean=[1 0 0; 1 0 1; 0 0 1; 0 1 1]; %[.251 .5098 .8235; .8235 .3529 .3529; .6473 .7456 .4188; .9763 .9831 .0538];
% maxY=500; %20; 
% newOrder=[1 3 5 7 2 4 6 8];
% for i=1:8
%     subplot(4,2,newOrder(i));
%     hold on;
%     for cond=1:4
%         plot(-0.5:0.001:1.5,mean(squeeze(meanfinalEMGepoch_baseRem{cond}(:,i,:)),2)*1000,'LineWidth',2,'Color',colsMean(cond,:)); %multiply by 1000 to convert to mV
%     end
%     plot([0 0],[0 maxY],'--k','LineWidth',2);
%     plot([-0.5 1.5],[100 100],':g','LineWidth',1);
%     hold off;
%     title(muscles{i}); ylim([0 maxY]);
% end
% legend('Stn Pull','Wlk Pull','Stn Viz','Wlk Viz');
