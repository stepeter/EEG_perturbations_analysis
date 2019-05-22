function [ICA_STRUCT, chan_rej_log] = chan_rej_hjh(EEG,good_chans)

%Function to run user selected channel rejection algorithms
%the algorithm options that the user selects as well as the channel
%rejection results are stored in a structure called ICA_STRUCT and saved
%in a user selected .MAT file
%
%INPUTS:
%   EEG: EEGLAB Structure (only used locally, not returned)
%   good_chans: optional input, index of channels that will not be rejected
%               these channels will be removed from the input EEG structure
%               before performing channel rejection and will automatically
%               be added to ICA_STRUCT.good_chans. Default = [];
%
%Created by: J Gwin 5/11/2009
%
%Modifications:
%J Gwin 05-13-10 - added optional input good_chans

if nargin < 2, good_chans = []; end

EEG = pop_select(EEG,'nochannel',good_chans);

%%%%%%%%automatic channel rejection%%%%%%%%%%%%%%%%%%%%
%eeg_badChannelsByCorr (by SCCN) is not a built in EEGLAB function
%note you need to iterate because with each run you find more and more
%bad channels because it is based on correlation with surrounding channels
original_chan_list = [1:EEG.nbchan];
finish_chan_rej = false;
frames = [1:EEG.pnts]; %to start with use all frames
while ~finish_chan_rej
    options = {'By Correlation','By Joint Probability (Kurt)',...
        'By Standard Deviation','By Range','Manually','From chan_rej_log'};
    selected = listdlg('PromptString','Channel Rejection Task?',...
            'Name','Automatic Channel Rejection','ListString',options,...
            'SelectionMode','single','ListSize',[300 150]);
    if isempty(selected)
        method = 'Cancel';
    else
        method = options{selected};
    end
    switch method
        case 'By Correlation'
            DEFAULT_THRES = {'0.0001'};
            thres = inputdlg('Threshold for eeg_badChannelsByCorr:',...
                '',1,DEFAULT_THRES,'on');
            if isempty(thres)
                thres = DEFAULT_THRES;
            end
            disp('Computing correlations...');
            rej_chan_idx = eeg_badChannelsByCorr(EEG,str2num(char(thres)));
            disp('Channels are being removed from temporary data set');
            disp('You will be able to accept OR reject this change later');
            chan_names = {EEG.chanlocs.labels};
            REJ_EEG = pop_select(EEG,'nochannel',chan_names(rej_chan_idx));
            finish_view = false;
        case 'By Joint Probability (Kurt)'
            DEFAULT_THRES = {'5'};
            thres = inputdlg('Threshold for pop_rejchan:',...
                '',1,DEFAULT_THRES,'on');
            if isempty(thres)
                thres = DEFAULT_THRES;
            end
            [REJ_EEG rej_chan_idx measure] = pop_rejchan(EEG,'elec',...
                [1:EEG.nbchan],'threshold',str2num(char(thres)),'norm',...
                'on','measure','kurt');
            disp('You will be able to accept OR reject this change later');
            finish_view = false;
        case 'By Standard Deviation'
            disp('Computing standard deviation...');
            try
                stdev = std(EEG.data,0,2);
            catch
                %catch memory error and compute std one channel at a time
                stdev = [];
                for i = 1:EEG.nbchan
                    disp(['Computing standard deviation for channel ' num2str(i) '...']);
                    stdev(i) = std(EEG.data(i,:));
                end
            end
            std_answer = questdlg(['Do you want to use default threshold (mean + 2 x S.D) or select a threshold?'], ...
                         'Threshold Question', ...
                         'Default', 'Manual Select','Default');
            if strcmp(std_answer,'Default')
                thres = mean(stdev)+2*std(stdev);
                thres = {num2str(round(thres*1000)/1000)};
            else
                h = figure; set(gcf,'color','w','name','Select threshold, close figure to accept')
                plot(sort(stdev),'r.');
                while gcf == h 
                    try
                        [crp,thres] = ginput(1);
                        clf(h);
                        plot(sort(stdev),'r.');
                        hold on; plot([0 EEG.nbchan],[thres thres],'k-','linewidth',2);
                    catch
                        break;
                    end
                end
                thres = {num2str(round(thres*1000)/1000)};
            end
            rej_chan_idx = find(stdev > str2num(char(thres)));
            disp('Channels are being removed from temporary data set');
            disp('You will be able to accept OR reject this change later');
            chan_names = {EEG.chanlocs.labels};
            REJ_EEG = pop_select(EEG,'nochannel',chan_names(rej_chan_idx));
            finish_view = false;
        case 'Manually'
            clear str
            for i = 1:length(EEG.chanlocs)
                str{i} = EEG.chanlocs(i).labels;
            end

            rej_chan_idx = eegchan_listdlg(EEG,'PromptString','Select Chans to Reject',...
                'SelectionMode','multiple','ListString',str);
  
            
            disp('Channels are being removed from temporary data set');
            disp('You will be able to accept OR reject this change later');
            chan_names = {EEG.chanlocs.labels};
            REJ_EEG = pop_select(EEG,'nochannel',chan_names(rej_chan_idx));
            thres = {''};
            for i = 1:length(rej_chan_idx)
                thres{1} = [thres{1} str{rej_chan_idx(i)} ' '];
            end
            finish_view = false;
        case 'By Range'
            DEFAULT_THRES = {'[30 10000]'};
            thres = inputdlg('Threshold for range (uV) [min max]:',...
                '',1,DEFAULT_THRES,'on');
            if isempty(thres)
                thres = DEFAULT_THRES;
            end
            disp('Computing range...');
            selected_range = str2num(thres{1});
            
            hr = figure;
            plot(range(EEG.data,2), 'o');
            hold on
            plot([1,EEG.nbchan], selected_range(1)*[1 1], 'r-', 'linewidth',1);
            plot([1,EEG.nbchan], selected_range(2)*[1 1], 'r-', 'linewidth',1);
            ylimits = ylim(gca);
            set(gca, 'ylim', [0 max(selected_range(2), ylimits(2))*1.1],...
                'xlim', [0 EEG.nbchan+1]);
            ylabel('Range (uV)'); xlabel('channels');
            text(260, selected_range(1), '30', 'horizontalalignment','left','color','r');
            text(260, selected_range(2), '10000', 'horizontalalignment','left','color','r');
            
            rej_chan_idx = find(range(EEG.data,2) < selected_range(1) ...
                | range(EEG.data,2) > selected_range(2));

            disp('Channels are being removed from temporary data set');
            disp('You will be able to accept OR reject this change later');
            chan_names = {EEG.chanlocs.labels};
            REJ_EEG = pop_select(EEG,'nochannel',chan_names(rej_chan_idx));
            finish_view = false;
        case 'From chan_rej_log'
            if EEG.nbchan == 256
                
                [logname pathname] = uigetfile('*.mat','Select chan_rej_log.mat to load');
                load(fullfile(pathname, logname));
                
                disp('Channels are being removed from temporary data set');
                disp('You will be able to accept OR reject this change later');
                chan2rej = [chan_rej_log{3,:}];
                REJ_EEG = pop_select(EEG,'nochannel',chan2rej);
                chan_names = {EEG.chanlocs.labels};
                [c ia rej_chan_idx] = intersect(chan2rej, chan_names);
                thres = ''; % for logging purposes
                finish_view = false;
            else
                h = helpdlg('From chan_rej_log option can only be performed on the full 256 EEG channel set');
                waitfor(h)
                finish_view = true;
            end
%         case 'Select Frames'
%             frames = inputdlg('Enter Frames','Frame Selection',1,{['[1:'...
%                 num2str(EEG.pnts) ']']});
%             if ~isempty(frames)
%                 eval(['frames = ' frames{1} ';']);
%                 EEG.data = EEG.data(:,frames);
%             end
%             finish_view = true;
        case 'Cancel'
            finish_view = true;
    end
    while ~finish_view
        
        % view topo and eeg plots        
        [crp, idx] = get_EEG_event_array(EEG,'boundary','latency');
        
        data2 = nan(size(EEG.data));
        data2(rej_chan_idx,:) = EEG.data(rej_chan_idx,:);
        
        eegplot(EEG.data,'command',[],'data2',data2,...
            'eloc_file',EEG.chanlocs,'srate',EEG.srate,...
            'events',EEG.event(idx));
        h2 = gcf;
        
        h1 = figure;
        mod_topoplot([],EEG.chanlocs,...
            'electrodes','on','emarker',{[1:EEG.nbchan],'.','k',10,1},...
            'emarker2',{rej_chan_idx,'.','r',20,1});
        
%         set(h2,'Name','Close the figure to continue processing');              
        waitfor(h2)
%         
        if ~ishandle(h2) && ishandle(h1), close(h1), end
        
        prompt = sprintf('%s\n%s','Reject channels?',...
        [num2str(EEG.nbchan-length(rej_chan_idx)) ' channels would remain']);
        answer = questdlg(prompt,...
            'Reject channels?','REJECT','Ignore','REJECT');
        if strcmp(answer, '')
            answer = 'Ignore';
        end
    
        switch answer
            case 'REJECT'
                % must find chan_names each time before replacing EEG with REJ_EEG
                % EEG iteratively changes with each rejection
                chan_names = {EEG.chanlocs.labels}; 
                EEG = REJ_EEG;
                clear REJ_EEG;
                EEG = eeg_checkset(EEG);
                disp([num2str(length(rej_chan_idx)) ' channels removed: ']);
                
                disp(chan_names(rej_chan_idx));
                finish_view = true;
                if strcmp(method,'From chan_rej_log')
                    EEG.etc.chan_rej = chan_rej_log;
                elseif ~isfield(EEG.etc, 'chan_rej')
                    EEG.etc.chan_rej{1,1} = [method ': ' char(thres)];
                    EEG.etc.chan_rej{2,1} = length(rej_chan_idx);
                    EEG.etc.chan_rej{3,1} = chan_names(rej_chan_idx);
                    EEG.etc.chan_rej{4,1} = sum([EEG.etc.chan_rej{2,:}]);
                    EEG.etc.chan_rej{5,1} = EEG.nbchan;
                else
                    [r c] = size(EEG.etc.chan_rej);
                    EEG.etc.chan_rej{1,c+1} = ...
                        [method ': ' char(thres)];
                    EEG.etc.chan_rej{2,c+1} = length(rej_chan_idx);
                    EEG.etc.chan_rej{3,c+1} = chan_names(rej_chan_idx);
                    EEG.etc.chan_rej{4,c+1} = sum([EEG.etc.chan_rej{2,:}]);
                    EEG.etc.chan_rej{5,c+1} = EEG.nbchan;
                end
                
                % save chan_rej variable to cluster in case chan_rej
                % crashes, which happens if you close eegplot before it is
                % done runnin                
                chan_rej_log = EEG.etc.chan_rej;
                save_chan_rej_log('tmp_chan_rej_log', chan_rej_log); 
                
            case 'Ignore'
                disp('Warning: No Channels Rejected');
                finish_view = true;
%                 chan_rej_log = {};
        end
    end
    
    if isfield(EEG.etc, 'chan_rej')
        num_badchans = sum([EEG.etc.chan_rej{2,:}]);
    else
        num_badchans = 0;
    end
    prompt = sprintf('%s\n%s\n%s','Run Channel Rejection Again?',...
        [num2str(num_badchans) ' channels rejected'],...
        [num2str(EEG.nbchan) ' channels remaining']);
    answer = questdlg(prompt,...
        'Automatic Channel Rejection','Yes','No','Yes');
    if strcmp(answer,'No') || strcmp(answer,'')
        finish_chan_rej = true;
%         chan_rej_log = {};
    end
end

%%%save chan_rej results to ICA_STRUCT
for n = 1:EEG.nbchan
    ICA_STRUCT.good_chans(n) = EEG.chanlocs(n).urchan;
end
if isfield(EEG.etc,'chan_rej')
    ICA_STRUCT.chan_rej_method = EEG.etc.chan_rej;
end
ICA_STRUCT.good_chans = sort([ICA_STRUCT.good_chans good_chans]);
ICA_STRUCT.associated_set = EEG.setname;
ICA_STRUCT.chan_rej_frames_used = frames;

% save text file and final chan_rej_log variable
chan_rej_log_files = dir(['chan_rej_log*.mat']);
vn = length(chan_rej_log_files);
answer = inputdlg('Save log file as?', 'Save log file as', 1, {['chan_rej_log_v' num2str(vn+1)]}); 
if ~isempty(answer)
    save_chan_rej_log(answer{1}, chan_rej_log);
    disp('**chan_rej_log saved**')
else
    disp('**chan_rej_log was not saved**')
end


