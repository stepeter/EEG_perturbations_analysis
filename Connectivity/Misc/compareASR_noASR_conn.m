%Compare brain-brain connectivity from ASR and non-ASR

corrs=zeros(2,4);
conds={'allSVZ','allWVZ','pullsStn','pullsWalk'};
for subjs=1:2
    for j=1:length(conds)
        EEG=pop_loadset('filename',['S' num2str(subjs) '.set'],'filepath',['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' conds{j}]);
        EEG_noASR=pop_loadset('filename',['S' num2str(subjs) '_' conds{j} '.set'],'filepath','/usr/local/VR_connectivity/Data/regularGroupConn_ASRtest/conns_hiOrder');
        
        X=EEG.CAT.Conn.dDTF08(1:8,1:8,:,:);
        Y=EEG_noASR.CAT.Conn.dDTF08(1:8,1:8,:,:);
        corrs(subjs,j)=corr2(X(:),Y(:));
    end
end

disp(['Mean: ' num2str(mean(corrs(:)))]);
disp(['Std: ' num2str(std(corrs(:)))]);