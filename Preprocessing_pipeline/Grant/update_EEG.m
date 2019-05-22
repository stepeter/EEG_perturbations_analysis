function EEG = update_EEG(EEG,ICA_STRUCT,override_popup,model)
%function EEG = update_EEG(EEG,ICA_STRUCT,override_popup)
%
%this function will update the EEG structure based on the paramaters stored
%in the ICA_STRUCT. The point of this code is to reduce the need to store
%lots of copies of the same data in multiple EEG.set files. ICA_STRUCT can
%also me a .mat filename
%
%Inputs:
%   EEG: EEG structure to be updated
%   ICA_STRUCT: ICA structure. This is not an EEGLAB variable, this is a
%               structure that can be used to manage ICA results. This 
%               structure is used in many of my EEGLAB add-on functions.
%   override_popup: binary, if true then won't present a popup warning if 
%                   EEG.setname is not the same as ICA_STRUCT.associated_set
%   model: integer, if ICA_STRUCT contains multiple ICA models then this 
%          model will be automatically selected. Otherwise user will be 
%          prompted to select a model to use 
%
%Output: Updated EEG structure
%Created by J. Gwin 05/21/2009 for U.Mich Human Neuromechanics Laboratory
%
%Modifications 
%   J. Gwin 07/30/09 added optional input override_poput
%       if true then overrides popup gui default = false
%   J. Gwin 05/28/10 alternatively you can pass a filename in as ICA_STRUCT
%   J. Gwin 09/09/10 added optional input 'model'

if nargin < 3
    override_popup = false;
end

if ~isstruct(ICA_STRUCT)
    EEG.etc.ICA_STRUCT = ICA_STRUCT;
    load(ICA_STRUCT);
end

%check that setname matches the orignial set used to create ICA_STRUCT
if ~override_popup && isfield(ICA_STRUCT,'associated_set') && ~strcmp(EEG.setname,ICA_STRUCT.associated_set)
    if ~strcmp(questdlg('EEG.setname does not match ICA_STRUCT.asociated_set',...
            'Continue?','CONTINUE ANYWAY','CANCEL','CONTINUE ANYWAY') ,...
            'CONTINUE ANYWAY')
        return;
    end
end


%update EEG structure
ica_fields = fieldnames(ICA_STRUCT);
for i = 1:length(ica_fields)
    switch ica_fields{i}
        case 'good_chans'
            all_chans = [];
            for n = 1:EEG.nbchan
                all_chans(n) = EEG.chanlocs(n).urchan;
            end
            if length(intersect(all_chans,ICA_STRUCT.good_chans)) < length(ICA_STRUCT.good_chans)
                error('EEG does not contain all the ICA_STRUCT.good_chans');
                return;
            else
                [tmp rej_chan_idx] = setxor(all_chans,ICA_STRUCT.good_chans);
                if ~isempty(rej_chan_idx)
                    if ischar(EEG.data)
                        %%%this is a set with no data just info
                        EEG.nbchan = length(ICA_STRUCT.good_chans);
                        EEG.chanlocs = EEG.chanlocs(ICA_STRUCT.good_chans);
                    else
                        %%%this set has data
                        %new versions of eeglab seem to require channel names
                        %for channel rejection
                        %EEG = pop_select(EEG,'nochannel',rej_chan_idx);
                        chan_names = {EEG.chanlocs.labels};
                        EEG = pop_select(EEG,'nochannel',chan_names(rej_chan_idx));
                    end
                    del_chan = true;
                else
                    del_chan = false;
                end
            end
        case 'rej_frame_idx'
            %do nothing...for now could insert a uioption to remove frames
        case 'weights'
            if size(ICA_STRUCT.weights,3) > 1
                if ~exist('model') || model > size(ICA_STRUCT.weights,3) || model < 0 || mod(model,1) ~= 0
                    for i = 1:size(ICA_STRUCT.weights,3)
                        str(i) = {num2str(i)};
                    end
                    [model] = listdlg('PromptString','Select a model:',...
                          'SelectionMode','single','ListString',str);
                end
                EEG.icaweights = ICA_STRUCT.weights(:,:,model);
            else
                EEG.icaweights = ICA_STRUCT.weights;
            end
        case 'sphere'
            EEG.icasphere = ICA_STRUCT.sphere;
        case 'frame_rej_params'
            %do nothing
        case 'good_comps'
            if isfield(ICA_STRUCT.good_comps,'all')
                EEG.reject.gcompreject = ICA_STRUCT.good_comps.all;
            end
        case 'dipfit'
            EEG.dipfit = ICA_STRUCT.dipfit;
    end
end
%check if need to re-ref
if  isfield(ICA_STRUCT,'ref')
    if (del_chan || ~strcmp(EEG.ref,ICA_STRUCT.ref)) ...
       && (strcmp(ICA_STRUCT.ref,'averef') || strcmp(ICA_STRUCT.ref,'avgref')) ...
       && ~ischar(EEG.data)
            
        disp('Computing average reference...');
        EEG = pop_reref( EEG, [], 'refstate',0);
    end
end
%recompute ICA acts if needed
EEG.icaact = [];
EEG.icawinv = [];
EEG = eeg_checkset(EEG); %compute ICA activations


