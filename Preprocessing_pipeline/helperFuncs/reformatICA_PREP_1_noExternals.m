%Reformat ICA_PREP_1 for 128 channels without externals
extChanNums=129:133;
for i=[1:13 15:17 19:20 22:33]
    load(['/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Data/WMISM_' num2str(i) '/ICA_Stuff/files_AMICAall/all_ICA_PREP_1.mat']);
    ICA_STRUCT.good_chans=setdiff(ICA_STRUCT.good_chans,extChanNums);
    ICA_STRUCT.good_cap = ICA_STRUCT.good_chans;
    ICA_STRUCT.good_ext = [];
    save(['/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Data/WMISM_' num2str(i) '/ICA_Stuff/files_AMICAall/all_ICA_PREP_1.mat'],'ICA_STRUCT');
end

%% Load each Merge_all.set and remove 5 external channels

for i=[1:13 15:17 19:20 22:33]
    EEG=pop_loadset(['/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Data/WMISM_' num2str(i) '/EEG_sets/Merge_all.set']);
    EEG=pop_select(EEG,'nochannel',extChanNums);
    pop_saveset(EEG,'filepath',['/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Data/WMISM_' num2str(i) '/EEG_sets/'],'filename','Merge_all.set');
end