%Create theta and alpha band networks (4-8, 8-13 Hz), averaging 1st second
%of activity; then do statistical testing using bootstrap distributions
%with same significance mask as average, this will determine which edges
%are significantly different from the other conditions

tVals=[0 1];
tInds=find(EEG.CAT.Conn.erWinCenterTimes>=tVals(1) & EEG.CAT.Conn.erWinCenterTimes<=tVals(2));
fThetaInds=find(EEG.CAT.Conn.freqs>=4 & EEG.CAT.Conn.freqs<=8);
fAlphaInds=find(EEG.CAT.Conn.freqs>=8 & EEG.CAT.Conn.freqs<=13);

condS={'allSVZ','allWVZ','pullsStn','pullsWalk'};
typeCondS={'brainBrain','EMGEMG','brainEMG','EMGbrain'};

for condInds=1:4
    cond=condS{condInds};
    
    netVals.theta=cell(16,16); netVals.alpha=cell(16,16); netVals.all=cell(16,16);
    netVals_ave.theta=zeros(16,16); netVals_ave.alpha=zeros(16,16); netVals_ave.all=zeros(16,16);
    
    if any(strcmpi(typeCond,{'brainEMG','EMGbrain'}))
        load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_baseSub_EMGbrainMask.mat']); %load connectivity data
    else
        load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_baseSub.mat']);
    end
    load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_boot.mat']);
    
    for j=1:16
        for k=1:16
            baseline=[-0.5 0];
            latencies=EEG.CAT.Conn.erWinCenterTimes;
            baseidx=find(latencies>=-0.5 & latencies<=0);
            baseVals=mean(connStruct_boot{j,k}(:,:,baseidx),3);
            curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(EEG.CAT.Conn.erWinCenterTimes)]);
            bootDat_theta=curr_ersp(:,fThetaInds,tInds);
            bootDat_alpha=curr_ersp(:,fAlphaInds,tInds);
            bootDat_all=curr_ersp(:,:,tInds);
            for frodo=1:size(bootDat_theta,1)
                A=bootDat_theta(frodo,:,:);
                A(squeeze(connStruct(j,k,fThetaInds,tInds))==0)=0; %mask using average mask
                netVals.theta{j,k}=[netVals.theta{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                tmp_dat = connStruct(j,k,fThetaInds,tInds);
                netVals_ave.theta(j,k)= mean(tmp_dat(:));
                
                A=bootDat_alpha(frodo,:,:);
                A(squeeze(connStruct(j,k,fAlphaInds,tInds))==0)=0; %mask using average mask
                netVals.alpha{j,k}=[netVals.alpha{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                tmp_dat = connStruct(j,k,fAlphaInds,tInds);
                netVals_ave.alpha(j,k)= mean(tmp_dat(:));
                
                A=bootDat_all(frodo,:,:);
                A(squeeze(connStruct(j,k,:,tInds))==0)=0; %mask using average mask
                netVals.all{j,k}=[netVals.all{j,k} mean(squeeze(mean(squeeze(A),1)))]; %A(:))];
                tmp_dat = connStruct(j,k,:,tInds);
                netVals_ave.all(j,k)= mean(tmp_dat(:));
            end
        end
    end
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/plots/median/netVals_' cond '_aveTime_sbjs.mat'],'netVals');%'_noTimeAve.mat'],'netVals');
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/plots/median/netVals_' cond '_aveTime.mat'],'netVals_ave');
end