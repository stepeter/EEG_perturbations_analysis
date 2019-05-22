%Created clustDat_new variable that maps EEG signals to group-level corticomuscular clusters

%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/usr/local/VR_connectivity/Code/EEGconn/'));
    startUpEEGLAB('close');
end

EEG=[]; ALLEEG=[]; CURRENTSET=[]; ALLCOM=[]; CURRENTSTUDY = []; STUDY = [];

%Load cluster information
load('/usr/local/VR_connectivity/Data/cluster_1222018.mat');
clustDat_new = clustDat;

%Load each subject's merged EEG file and ICA_STRUCT; then calculate
%appropriate icaweights variable based on icawinv and icasphere
for i=subjectInds %1:2 %
    subjcode = ['WMISM_' num2str(i)];
    subj_basepath = [RootData subjcode filesep];
    EEGsets_outpath = [subj_basepath 'EEG_sets'];
    ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
    cd(subj_basepath)
    
    EEG=pop_loadset('filename', ['Merge_' group '_CAR_v2.set'], 'filepath', EEGsets_outpath);
    EEG=eeg_checkset(EEG); EEG2=EEG; 
    load([ICA_path_out filesep group '_ICA_DIPFIT_5.mat']);
    
    %Copy of EEG with original IC's
    EEG2.icachansind=[]; EEG2.setname=ICA_STRUCT.associated_set;
    EEG2 = update_EEG(EEG2,ICA_STRUCT);
    %Remove bad ICs (only keeping cortical ICs)
    compInds=setdiff([1:size(EEG2.icaweights,1)],ICA_STRUCT.good_comps.brain);
    EEG2=pop_subcomp(EEG2,compInds);
    EEG2=eeg_checkset(EEG2);
    
    %Now figure out which clusters this subject is missing
    icaDat2Add=[]; icaEst4_crossVal=[];
    for j=1:length(goodClusts)
        if ~any(allSubjs(clustDat_new(goodClusts(j)).sets)==i)

        else
            %Find IC from original data to use (use one with higher
            %variance if multiple ones)
            setIndVals=find(allSubjs(clustDat_new(goodClusts(j)).sets)==i);
            dipoleInds=clustDat_new(goodClusts(j)).comps(setIndVals);
            if length(dipoleInds)>1
                dipoleInds=min(dipoleInds); %min number has higher variance
                
                badDipInds = setdiff(clustDat_new(goodClusts(j)).comps(setIndVals),dipoleInds);
                rm_inds = [];
                for ppp=badDipInds
                    tmp_ind = find(clustDat_new(goodClusts(j)).comps(setIndVals)==ppp);
                    rm_ind = setIndVals(tmp_ind);
                    rm_inds = [rm_inds, rm_ind];
                    
                    %Move to cluster 2 (outlier)
                    clustDat_new(2).comps = [clustDat_new(2).comps, badDipInds];
                    clustDat_new(2).sets = [clustDat_new(2).sets, clustDat_new(goodClusts(j)).sets(rm_ind)];
                    clustDat_new(2).allinds{1} = [clustDat_new(2).allinds{1}, badDipInds];
                    clustDat_new(2).allinds{2} = [clustDat_new(2).allinds{2}, badDipInds];
                    clustDat_new(2).setinds{1} = [clustDat_new(2).setinds{1}, clustDat_new(goodClusts(j)).sets(rm_ind)];
                    clustDat_new(2).setinds{2} = [clustDat_new(2).setinds{2}, clustDat_new(goodClusts(j)).sets(rm_ind)];
                end
                %Remove from original cluster
                clustDat_new(goodClusts(j)).comps(rm_inds) = [];
                clustDat_new(goodClusts(j)).sets(rm_inds) = [];
                clustDat_new(goodClusts(j)).allinds{1}(rm_inds) = [];
                clustDat_new(goodClusts(j)).allinds{2}(rm_inds) = [];
                clustDat_new(goodClusts(j)).setinds{1}(rm_inds) = [];
                clustDat_new(goodClusts(j)).setinds{2}(rm_inds) = [];
            end
            if isempty(icaDat2Add)
                icaDat2Add=EEG2.icaact(dipoleInds,:);
            else
                icaDat2Add=cat(1,icaDat2Add,EEG2.icaact(dipoleInds,:));
            end
        end
    end
    
    %Finally, add resulting ICA data to new EEG set and save it to use
    %later
    EEGout=EEG;
    EEGout.icaweights=eye(size(icaDat2Add,1));
    EEGout.icasphere=EEGout.icaweights;
    EEGout.icaact=[]; EEGout.icawinv=[]; EEGout.icachansind=[]; EEGout.chanlocs=[];
    EEGout.nbchan=8; EEGout.data=icaDat2Add;
    EEGout=eeg_checkset(EEGout,'ica');
    
#    pop_saveset(EEGout,'filename',['S' num2str(i)],'filepath',savePath);
end
disp('Done! :)');

%Remove dipole info so gets recomputed
for i=1:length(clustDat_new)
    clustDat_new(i).dipole=[];
end