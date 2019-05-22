function epochDIPFITData_Math540mod(pathname,subjnum,filter_lo,group,trigPrefix1,trigPrefix2,cerbellGroup,analyzeMismDirection,epochTimes)

%Define paths
subjcode = ['MISM_' num2str(subjnum)];
session = 'Day_2';
subj_basepath = [pathname subjcode filesep session filesep]; % Top level subject folder
EEGsets_outpath = [subj_basepath 'EEG_sets'];
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group cerbellGroup];
[filtertag]=filterSizeStr(filter_lo);
epch_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group cerbellGroup '_epch'];
mergeEpch_path = [pathname 'STUDY_sets_' group cerbellGroup '_short'];
cd(subj_basepath)

global ALLEEG EEG CURRENTSET ALLCOM;

%Create epoch folder (if needed)
if ~exist(epch_path_out)
    mkdir(epch_path_out)
end
if ~exist(mergeEpch_path)
    mkdir(mergeEpch_path)
end

%Load dataset before epoching
EEG=pop_loadset([EEGsets_outpath filesep 'Merge_' group cerbellGroup '_HPF_' filtertag '.set']);
load([ICA_path_out filesep group cerbellGroup '_ICA_DIPFIT_5.mat']);
[EEG,ICA_STRUCT]=interpolEEG(EEG,ICA_STRUCT,0,1);
EEG = update_EEG(EEG,ICA_STRUCT);
EEG=pop_reref(EEG,[]);

%Remove bad ICs (only keeping cortical ICs)
compInds=setdiff([1:size(EEG.icaweights,1)],ICA_STRUCT.good_comps.brain);
EEG=pop_subcomp(EEG,compInds);
EEG=eeg_checkset(EEG)

%Create trig labels
if analyzeMismDirection==0 
    %Remove CW/CCW from event labels
    for i=1:length({EEG.event(1:end).type})
        ccwEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigPrefix1 'CCW_*']),'once');
        cwEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',[trigPrefix1 'CW_*']),'once');
        if ~isempty(ccwEv)
            EEG.event(1,i).type=[trigPrefix1 EEG.event(1,i).type((length(trigPrefix1)+5):end)];
        elseif ~isempty(cwEv)
            EEG.event(1,i).type=[trigPrefix1 EEG.event(1,i).type((length(trigPrefix1)+4):end)];
        end
    end
end

%Build trig labels
types=unique({EEG.event(1:end).type});
trigsSize=sum(cell2mat(cellfun(@(x) regexp(x, regexptranslate('wildcard',[trigPrefix1 '*'])),types,'UniformOutput',false)))+sum(cell2mat(cellfun(@(x) regexp(x, regexptranslate('wildcard',[trigPrefix2 '*'])),types,'UniformOutput',false)));
trigs=cell(1,trigsSize); q=1;
for i=1:length(types)
    if regexp(types{1,i}, regexptranslate('wildcard',[trigPrefix1 '*']))
        trigs{1,q}=types{1,i};
        q=q+1; %update iterator term
    end
end
for i=1:length(types)
    if regexp(types{1,i}, regexptranslate('wildcard',[trigPrefix2 '*']))
        trigs{1,q}=types{1,i};
        q=q+1; %update iterator term
    end
end

%Epoch parameters
CURRENTSET=1;
% Set epoch bounds (add in a little extra for ERSPs)
EpPre = epochTimes(1,1)+(-0.559);
EpPost = epochTimes(1,2)+0.559;
EEG_preEpoch=EEG; %save this for iterations in for loop

for p=1:length(trigs)
    EEG=EEG_preEpoch; %use original (non-epoched) data
    EEG = pop_epoch(EEG,trigs(1,p),[EpPre EpPost],...
        'newname',[subjcode '_Epc'],'epochinfo','yes');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0);
    EEG=eeg_checkset(EEG,'makeur')
    EEG = pop_rmbase( EEG,[]);
    %Get rid of really noisy epochs (due to ungrounding)
    [EEG,inds]=pop_eegthresh(EEG,1,1:EEG.nbchan,-500,500,EpPre,EpPost,0,1);
    
    EEG=pop_select(EEG,'trial',1:25); %added in to fix number of trials
    
    EEG=pop_editset(EEG,'setname',['Epch_' group cerbellGroup '_' num2str(subjnum) '_' filtertag '_' trigs{1,p}]);
    EEG=eeg_checkset(EEG,'ica')
%     pop_saveset(EEG, 'filename', ['Epch_' group cerbellGroup '_' num2str(subjnum) '_' filtertag '_' trigs{1,p}], 'filepath', epch_path_out);
    EEG_epched(p)=EEG; %save set to be merged later
end

%Now merge these saved sets
if length(EEG_epched) > 1
    EEG = [];
    EEG = pop_mergeset(EEG_epched, 1:length(EEG_epched));
    EEG.setname = ['Merge_epch_' group cerbellGroup '_' num2str(subjnum)];
    EEG = eeg_checkset(EEG); % always checkset 
    EEG=update_EEG(EEG,ICA_STRUCT,1);
    EEG=pop_subcomp(EEG,compInds);
    EEG=eeg_checkset(EEG);
    pop_saveset(EEG, 'filename', ['Merge_epch_' group cerbellGroup '_' num2str(subjnum)],...
        'filepath', mergeEpch_path);
else
    display('No need to merge datasets...')
end
disp(['Finished ' subjcode '!']);