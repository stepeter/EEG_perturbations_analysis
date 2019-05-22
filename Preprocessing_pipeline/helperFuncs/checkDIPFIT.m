%function []=checkDIPFIT(subjnum,group,RootData)
%Check results from DIPFIT and select good components
%RUN CELL BY CELL (not all at once!)
subjnum=33; %[1:13 15:17 19:20 22:33];
%Redo: 5,8,11,22,31

%File paths
subjcode = ['WMISM_' num2str(subjnum)];
subj_basepath = [RootData subjcode filesep]; % Top level subject folder
EEGsets_outpath = [subj_basepath 'EEG_sets']; % Where you want to save your .set files
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group ];

%Load data
EEG=pop_loadset([EEGsets_outpath filesep 'Merge_' group '_CAR_v2.set']);
load([ICA_path_out filesep group '_ICA_DIPFIT_3']);
EEG.icachansind=[]; EEG.setname=ICA_STRUCT.associated_set;
EEG = update_EEG(EEG,ICA_STRUCT);

if exist([ICA_path_out filesep group '_ICA_DIPFIT_5.mat'])
    error('Already did this!');
end
% %Reject frames (just for looking at power spectra)
% EEG.data(:,ICA_STRUCT.rej_frame_idx) = [];
% EEG.icaact(:,ICA_STRUCT.rej_frame_idx) = [];
% EEG.times = EEG.times(1:size(EEG.data,2));
% EEG.pnts=size(EEG.data,2);

% Definitions for dipfit
GOOD_COMPS_CAT = {'other','brain'};
RV_THRES = 0.15; %comps with RV < RV_THRES will be pre-selected

%Pick good comps based on topoplot
%preset any comps with rv < RV_THRES as selected
for i = 1:length(EEG.dipfit.model)
    if EEG.dipfit.model(i).rv < RV_THRES
        EEG.reject.gcompreject(i) = 1;
    end
end

% mod_selectcomps(EEG,[1:size(EEG.icaact,1)],0); %J Gwin modified version of pop_selectcomps

ICA_STRUCT.good_comps.all = EEG.reject.gcompreject;

% save([ICA_path_out filesep subjcode '_ICA_DIPFIT_4'], 'ICA_STRUCT')
disp('Run cell 2 now!');
%%final selection and categorize comps
disp('Mean centering component activations...');
mean_ica_act = EEG.icaact(find(EEG.reject.gcompreject == 1),:)-...
    repmat(mean(EEG.icaact(find(EEG.reject.gcompreject == 1),:),2),1,EEG.pnts);
tmp_good_comps = find(EEG.reject.gcompreject == 1);

% avoid figure for defining good or bad components
% ICA_STRUCT.good_comps.interest = tmp_good_comps;

if ~isempty(tmp_good_comps)
    for i = 1:length(tmp_good_comps)
        tmp_comp_struct(i).labels = num2str(tmp_good_comps(i));
    end
    eegplot(mean_ica_act,'submean','off','srate',EEG.srate,...
        'eloc_file',tmp_comp_struct,'events',EEG.event);
    mod_selectcomps(EEG,find(EEG.reject.gcompreject == 1),1);
    pop_dipplot( EEG,find(EEG.reject.gcompreject == 1) , 'normlen', 'on','projlines', 'on');
    uiwait(gcf);
%     eeglab_handle=gcf;
    
    %hold on until all pop_selectcomps windows are closed
%     while gcf ~= eeglab_handle
%         pause(1);
%     end
    
    %Sort selected comps into categories
    tmp = find(EEG.reject.gcompreject == 1);
    tmp_string = [];
    for i = 1:length(tmp)
        tmp_string{i} = ['Comp: ' num2str(tmp(i))];
    end
    for i = 1:length(GOOD_COMPS_CAT)
        if isempty(tmp_string)
            eval(['ICA_STRUCT.good_comps.' GOOD_COMPS_CAT{i} ' = [];'])
        else
            selected = listdlg('ListString',tmp_string,'SelectionMode','multiple',...
                'PromptString',['Select ' GOOD_COMPS_CAT{i} ' comps'],...
                'ListSize',[300 300]);
            eval(['ICA_STRUCT.good_comps.' GOOD_COMPS_CAT{i} ' = tmp(selected);'])
            tmp_string(selected) = [];
            tmp(selected) = [];
        end
    end
    ICA_STRUCT.good_comps.all = EEG.reject.gcompreject;
    
    %create a structure of just those ICs categorized above, I think this is
    %redundant, but is used in later analysis and doesn't seem worth the effort
    %to remove/change
    if exist('ICA_STRUCT.good_comps.brainstem', 'var') ~= 0
        ICA_STRUCT.good_comps.interest=[ICA_STRUCT.good_comps.brain ICA_STRUCT.good_comps.brainstem ICA_STRUCT.good_comps.other];
    else
        ICA_STRUCT.good_comps.interest=[ICA_STRUCT.good_comps.brain ICA_STRUCT.good_comps.other];
    end
end

% save
save([ICA_path_out filesep group '_ICA_DIPFIT_5.mat'], 'ICA_STRUCT');
close all;
disp('Done! :)');