function [ICA_STRUCT] = frame_rej_hjh(EEG,ICA_STRUCT)

%Function to select bad frames based on several criterion. Requires user
%input also user to view results
%
%INPUTS
%   EEG: EEG set to use. 
%   ICA_STRUCT: optional input. It input then new field are appended to 
%       ICA_STRUCT. If no input then new ICA_STRUCT field is created
%
%Created by: J Gwin 5/11/2009

%Modifications:


%%%%%%%%%%definitions%%%%%%%%%%
DEFAULT_MEAN_IQR_THRES = {'0.8'};
DEFAULT_PRCTILE_THRES = {'2, 0.9'};
DEFAULT_MIN_SPACE_BETWEEN_BAD_FRAMES = 50;
DEFAULT_BAD_FRAME_EVENT_BORDER_DUR = 10;

%%%%%handle inputs
if nargin < 2
    ICA_STRUCT = [];
end
if isempty(ICA_STRUCT)
    ICA_STRUCT.rej_chan = [];
    ICA_STRUCT.ref = 'averef';
end


%%%%%%%%automatic frame rejection%%%%%%%%
%NOTE: Frame rejection is for ICA only. Rejected frames are NOT removed
%from EEG.data. Indicies of rejected frames are saved in the ICA_STRUCT
%FR (Frame Reject) is a structure with parameters related to the frame
%rejection algs and is used temporarily and then discarded before running
%ICA to save space
FR.has_run_eeg_badframes = false; %note only need to run eeg_badframes once
FR.finished_frame_rej = false;
FR.all_bad_frames = [];
while ~FR.finished_frame_rej
    %%only use the mean threshold option, the percentile option never 
    %%worked well, it was just an idea. J. Gwin 12-6-10
    %FR.options = {'eeg_badframes mean thres','eeg_badframes percentile thres'};
    %FR.selected = listdlg('PromptString','Frame Rejection Method?',...
    %        'Name','Automatic Frame Rejection','ListString',FR.options,...
    %        'SelectionMode','single','ListSize',[300 150]);
    %if isempty(FR.selected)
    %    FR.method = 'Cancel';
    %else
    %    FR.method = FR.options{FR.selected};
    %end
    FR.method = 'eeg_badframes mean thres';
    switch FR.method
        case 'eeg_badframes mean thres'
            FR.thres = inputdlg('iqr average:','',1,DEFAULT_MEAN_IQR_THRES,'on');
            if isempty(FR.thres)
                FR.thres = DEFAULT_MEAN_IQR_THRES;
            end
            if ~FR.has_run_eeg_badframes
                disp('Running eeg_badframes...');
                [cleanData FR.bad_frames FR.frame_type ...
                    FR.norm_power FR.mean_norm_power] = ...
                    eeg_badframes(EEG,str2double(char(FR.thres)));
                FR.has_run_eeg_badframes = true; %only need to run once
            else
                %don't need to re-run just recompute bad_frames    
                FR.bad_frames = find(FR.mean_norm_power>str2double(char(FR.thres)));
                FR.frame_type = ones(1,EEG.pnts);
                FR.frame_type(FR.bad_frames) = 0;
            end
        case 'eeg_badframes percentile thres'
            FR.thres = inputdlg('thres,prctile:','',1,DEFAULT_PRCTILE_THRES,'on');
            if isempty(FR.thres)
                FR.thres = DEFAULT_PRCTILE_THRES;
            end
            thres_num = str2double(char(FR.thres));
            if ~FR.has_run_eeg_badframes
                disp('Running eeg_badframes...');
                [cleanData FR.bad_frames FR.frame_type ...
                    FR.norm_power FR.mean_norm_power] = ...
                    eeg_badframes(EEG,thres_num(1));
                FR.has_run_eeg_badframes = true; %only need to run once
            end
            disp('Computing percentiles...');
            [FR.rows,FR.cols] = find(FR.norm_power > thres_num(1));
            FR.bad_frame_mat = zeros(size(FR.norm_power));
            FR.bad_frame_mat(sub2ind(size(FR.norm_power),FR.rows,FR.cols))=1;
            FR.bad_frames = find(sum(FR.bad_frame_mat,2) >= ...
                EEG.nbchan*thres_num(2));
            FR.frame_type = ones(1,EEG.pnts);
            FR.frame_type(FR.bad_frames) = 0;
    end
    if ~strcmp(FR.method,'Cancel') && ~isempty(FR.bad_frames)
        %add a border to each bad frame. i.e. if the frame is bad mark as bad
        %plus-minus XX frames on either side
        disp('Adding border and min spacing to bad_frames...');
        FR.bad_frame_border_dur = inputdlg('Bad frame border duration:','',1,...
            {num2str(DEFAULT_BAD_FRAME_EVENT_BORDER_DUR)},'on');
        if isempty(FR.bad_frame_border_dur)
            FR.bad_frame_border_dur = {num2str(DEFAULT_BAD_FRAME_EVENT_BORDER_DUR)};
        end
        FR.bad_frame_border_dur = str2num(char(FR.bad_frame_border_dur));
        FR.diff_frame_type = diff(FR.frame_type);
        tmp_idx = find(FR.diff_frame_type == -1);
        for i = 1:length(tmp_idx)
            if tmp_idx(i)-FR.bad_frame_border_dur+1 < 1
                FR.frame_type([1:tmp_idx(i)]) = 0;
            else
                FR.frame_type([tmp_idx(i)-FR.bad_frame_border_dur+1:tmp_idx(i)]) = 0;
            end
        end
        tmp_idx = find(FR.diff_frame_type == 1);
        for i = 1:length(tmp_idx)
            FR.frame_type([tmp_idx(i):tmp_idx(i)+FR.bad_frame_border_dur]) = 0;
        end
        FR.frame_type(EEG.pnts:end) = [];
        FR.bad_frames = find(FR.frame_type == 0);
        %if bad frames are closer than XXX time points then fill in the spaces
        %with more bad frames
        FR.min_bad_frame_spacing = inputdlg('Minimum bad frame spacing:','',1,...
            {num2str(DEFAULT_MIN_SPACE_BETWEEN_BAD_FRAMES)},'on');
        if isempty(FR.bad_frame_border_dur)
            FR.min_bad_frame_spacing = {num2str(DEFAULT_MIN_SPACE_BETWEEN_BAD_FRAMES)};
        end
        FR.min_bad_frame_spacing = str2num(char(FR.min_bad_frame_spacing));
        FR.bad_frame_spacing = diff(FR.bad_frames);
        tmp_idx = find(FR.bad_frame_spacing < FR.min_bad_frame_spacing ...
            & FR.bad_frame_spacing > 1);
        for i = 1:length(tmp_idx)
            FR.frame_type(FR.bad_frames(tmp_idx(i))+1:FR.bad_frames(tmp_idx(i)+1)-1) = 0;
        end
        FR.bad_frames = find(FR.frame_type == 0);
        disp([num2str(length(FR.bad_frames)/EEG.pnts) ' percent of frames are bad']);
        %create bad frame events
        disp('Creating bad frame events...');
        dur = 1; %duration counter
        curr_num_marked = 0;
        if exist('bad_frame_event')
            n = length(FR.bad_frame_event)+1;
        else 
            n = 1;
        end
        FR.bad_frame_event(n).type = 'bad frame';
        FR.bad_frame_event(n).latency = FR.bad_frames(1);
        FR.bad_frame_event(n).duration = [];
        FR.bad_frame_event(n).urevent = 1;    
        for i = 2:length(FR.bad_frames)
            if FR.bad_frames(i) == FR.bad_frames(i-1)+1
                dur = dur + 1;
            else
                %apply duration to previous event
                FR.bad_frame_event(n).duration = dur;
                curr_num_marked = curr_num_marked + FR.bad_frame_event(n).duration;
                %start new event
                n = n + 1;
                FR.bad_frame_event(n).type = 'bad frame'; 
                FR.bad_frame_event(n).latency = FR.bad_frames(i);
                FR.bad_frame_event(n).urevent = n;    
                FR.bad_frame_event(n).duration = [];
                %reset duration counter
                dur = 1;
            end
        end
        %fill in the duration for the last event
        FR.bad_frame_event(n).duration = dur;
        curr_num_marked = curr_num_marked + FR.bad_frame_event(n).duration;
        %plot results
%         h = mod_eegplot(EEG.data,'command',[],'eloc_file',...
%             EEG.chanlocs,'events',FR.bad_frame_event,'srate',EEG.srate);
%         set(h,'Name','Close the figure to continue processing');
%         waitfor(h);

%         eegplot(EEG.data,'command',[],'eloc_file',...
%             EEG.chanlocs,'events',FR.bad_frame_event,'srate',EEG.srate);
        
        clear view_badframes
view_badframes(:,1) = [FR.bad_frame_event.latency]';
view_badframes(:,2) = [FR.bad_frame_event.latency]'+[FR.bad_frame_event.duration]';
view_badframes(:,3:5) = repmat([0.5 0.5 0.5],length(FR.bad_frame_event),1);
view_badframes = [view_badframes repmat(ones(1,EEG.nbchan),length(FR.bad_frame_event),1)];
           
eegplot(EEG.data,'command',[],'eloc_file',EEG.chanlocs,'winrej', view_badframes,'srate',EEG.srate);
h = gcf;
waitfor(h);
%         set(h,'Name','Close the figure to continue processing');
%         h = gcf;
%         waitfor(h);


        reject_frame_answer = questdlg('Reject selected frames?',...
                        [num2str(curr_num_marked) ' frames marked'],'Yes','No','Yes');
        if strcmp(reject_frame_answer,'Yes')
            if isfield(ICA_STRUCT, 'frame_rej_method')
                ICA_STRUCT.frame_rej_method{length(ICA_STRUCT.frame_rej_method)+1} = ...
                    [FR.method ': ' char(FR.thres)];
            else
                ICA_STRUCT.frame_rej_method = {[FR.method ': ' char(FR.thres)]};
            end
            FR.all_bad_frames = [FR.all_bad_frames FR.bad_frames];
        else
            FR.bad_frame_event(n:end) = []; %clear the events you just added
        end
    end
    rerun_answer = questdlg('Re-run Auto Frame Rejection?',...
        'Frame Rejection','Yes','No','Yes');
    if strcmp(rerun_answer,'No') || strcmp(rerun_answer,'')
        FR.finished_frame_rej = true;    
    end
end

%%%%%%%%%%%%%%%%save ICA_STRUCT%%%%%%%%%%%%%%%%%%%%%%
FR.all_bad_frames = unique(FR.all_bad_frames);
ICA_STRUCT.rej_frame_idx = FR.all_bad_frames;
 

% I didn't get to this but i want to add an option for manual frame
% rejectino using the built in eegplot function...
% %%%%%%%%%%%%%manual frame rejection%%%%%%%%%%%%%%%%%%%%%%%
% while ~finish_frame_rej
%     eegplot(tmpEEG.data)
%      h = mod_eegplot(EEG.data,'command',[],'color',color,...
%                     'eloc_file',EEG.chanlocs);
%     answer = questdlg('Run Frame Rejection Again?',...
%         'Frame Rejection','Yes','No','Yes');
%     if strcmp(answer,'No') | strcmp(answer,'')
%         finish_frame_rej = true;
%     end
% end

   