%Compare brain-brain connectivity from low and high model orders

corrs=zeros(29,4);
conds={'allSVZ','allWVZ','pullsStn','pullsWalk'};
for subjs=[1:7 9:13 15:17 19:20 22:33] %1:2
    for j=1:length(conds)
        EEG=pop_loadset('filename',['S' num2str(subjs) '.set'],'filepath',['/usr/local/VR_connectivity/Data/regularGroupConn/' conds{j}]);
        EEG_hiMod=pop_loadset('filename',['S' num2str(subjs) '.set'],'filepath',['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' conds{j}]);
        
        X=EEG.CAT.Conn.dDTF08(1:8,1:8,:,:);
        Y=EEG_hiMod.CAT.Conn.dDTF08(1:8,1:8,:,:);
        corrs(subjs,j)=corr2(X(:),Y(:));
    end
end
%%
clc;
for i = 1:length(conds)
    disp(conds{i});
    disp(['Mean: ' num2str(mean(corrs(:,i)))]);
    disp(['Std: ' num2str(std(corrs(:,i)))]);
end