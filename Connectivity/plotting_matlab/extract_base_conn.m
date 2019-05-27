%Get average connectivity values
condS={'allSVZ','allWVZ','pullsStn','pullsWalk'};
typeCondS={'brainBrain','EMGEMG','brainEMG','EMGbrain'};
tVals=[-0.5 1.5];
tInds=find(EEG.CAT.Conn.erWinCenterTimes>=tVals(1) & EEG.CAT.Conn.erWinCenterTimes<=tVals(2));

fThetaInds=find(EEG.CAT.Conn.freqs>=4 & EEG.CAT.Conn.freqs<=8);
fAlphaInds=find(EEG.CAT.Conn.freqs>=8 & EEG.CAT.Conn.freqs<=13);
fBetaInds=find(EEG.CAT.Conn.freqs>=13 & EEG.CAT.Conn.freqs<=30);

% baseConn.theta=zeros(16,16); baseConn.alpha=zeros(16,16); baseConn.beta=zeros(16,16); baseConn.all=zeros(16,16);
baseConn.theta=cell(16,16); baseConn.alpha=cell(16,16); baseConn.beta=cell(16,16); baseConn.all=cell(16,16);
for condInds=1:4
    for typeCondInds=1:4
        cond=condS{condInds}; %'pullsStn'; %close all;
        typeCond=typeCondS{typeCondInds}; %'brainBrain'; %'brainBrain','EMGEMG','brainEMG','EMGbrain' (from -> to)

        load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_boot.mat']);
        
        
%         %Set diagonal to zero
%         for i=1:size(connStruct,1)
%             connStruct(i,i,:,:)=0;
%         end
        for i=1:size(connStruct_boot,1)
            connStruct_boot{i,i}=connStruct_boot{i,i}*0;
        end
        
        for i=1:16
            for j=1:16
                tmpVal=mean(mean(connStruct_boot{i,j}(:,fThetaInds,tInds),3),2);
%                 baseConn.theta(i,j)=mean(tmpVal(:));
                baseConn.theta{i,j}=tmpVal(:);
                tmpVal=mean(mean(connStruct_boot{i,j}(:,fAlphaInds,tInds),3),2);
%                 baseConn.alpha(i,j)=mean(tmpVal(:));
                baseConn.alpha{i,j}=tmpVal(:);
                tmpVal=mean(mean(connStruct_boot{i,j}(:,fBetaInds,tInds),3),2);
%                 baseConn.beta(i,j)=mean(tmpVal(:));
                baseConn.beta{i,j}=tmpVal(:);
                tmpVal=mean(mean(connStruct_boot{i,j}(:,:,tInds),3),2);
%                 baseConn.all(i,j)=mean(tmpVal(:));
                baseConn.all{i,j}=tmpVal(:);
            end
        end
        save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/plots/median/all_baseConn_boot_' cond '_' typeCond '.mat'],'baseConn');
    end
end