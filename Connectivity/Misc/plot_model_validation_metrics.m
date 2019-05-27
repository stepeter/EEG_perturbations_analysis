%Extracts model validation information for autoregressive model fits and plots the information criteria results across a range of model orders.

cond='pullsWalk'; %{'allSVZ','allWVZ','pullsStn','pullsWalk'};
modelOrders=[]; residualWhiteness=[]; stability=[]; consistency=[]; param2Dp_ratio=[];
%Load each subject's merged EEG file and ICA_STRUCT; then calculate
%appropriate icaweights variable based on icawinv and icasphere
connStruct=zeros(16,16,30,244);
numConns=zeros(16,16);
connStructS=connStruct;
% cond=conds{frodo}; %'allSVZ';
for i=subjectInds
    subjcode = ['WMISM_' num2str(i)];
    subj_basepath = [RootData subjcode filesep];
    EEGsets_outpath = ['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond]; %allSVZ'; %[subj_basepath 'EEG_sets'];
    ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
    
    EEG=pop_loadset('filename', ['S' num2str(i) '.set'], 'filepath', EEGsets_outpath);
    modelOrders=[modelOrders EEG.CAT.MODEL.morder];
    residualWhiteness=[residualWhiteness mean(EEG.CAT.VALIDATION.whitestats.acf.pval)];
    stability=[stability mean(EEG.CAT.VALIDATION.stabilitystats.lambda(:))];
    consistency=[consistency mean(EEG.CAT.VALIDATION.PCstats.PC)];
    param2Dp_ratio=[param2Dp_ratio (EEG.CAT.MODEL.morder*(EEG.CAT.nbchan^2))/(EEG.CAT.trials*(EEG.CAT.MODEL.winlen*EEG.srate))];
end
clc
disp(['Model orders: ' num2str(mean(modelOrders)) ' (' num2str(std(modelOrders)) ')']);
disp(['Residual whiteness: ' num2str(mean(residualWhiteness)) ' (' num2str(std(residualWhiteness)) ')']);
disp(['Stability: ' num2str(mean(stability)) ' (' num2str(std(stability)) ')']);
disp(['Consistency: ' num2str(mean(consistency)) ' (' num2str(std(consistency)) ')']);
disp(['param2Dp_ratio: ' num2str(mean(param2Dp_ratio)) ' (' num2str(std(param2Dp_ratio)) ')']);

%%
cond='allSVZ'; %{'allSVZ','allWVZ','pullsStn','pullsWalk'};
IC=zeros(4,length(subjectInds),40);
%Load each subject's merged EEG file and ICA_STRUCT; then calculate
%appropriate icaweights variable based on icawinv and icasphere

% cond=conds{frodo}; %'allSVZ';
q=0;
for i=subjectInds
    q=q+1;
    subjcode = ['WMISM_' num2str(i)];
    subj_basepath = [RootData subjcode filesep];
    EEGsets_outpath = ['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond]; %allSVZ'; %[subj_basepath 'EEG_sets'];
    ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
    
    EEG=pop_loadset('filename', ['S' num2str(i) '.set'], 'filepath', EEGsets_outpath);
    IC(1,q,:) = mean(EEG.CAT.IC.sbc.ic,2)';
    IC(2,q,:) = mean(EEG.CAT.IC.aic.ic,2)';
    IC(3,q,:) = mean(EEG.CAT.IC.fpe.ic,2)';
    IC(4,q,:) = mean(EEG.CAT.IC.hq.ic,2)';
end
%%
figure; hold on;
cols='brgc';
for i = 1:4
    errorbar(1:40,squeeze(mean(squeeze(IC(i,:,:)))),squeeze(std(squeeze(IC(i,:,:)))),cols(i));
end
hold off;
xlim([0 40]);
legend('sbc','aic','fpe','hq');