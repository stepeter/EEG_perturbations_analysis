function [EEG_raw,subj_basepath]=EEG_load_CAR_beforePREP(subjnum,RootData,numChans)


removeExternals=1; %decide if to remove externals (1 - yes; 0 - no)
%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjcode = ['MISM_' num2str(subjnum)]; %define subjcode
session = 'Day 2';

subj_basepath = [RootData subjcode filesep session filesep ]; % Top level subject folder
rawEEG_inpath = [subj_basepath 'RawEEG']; % Where your raw .bdf files are stored

% Digitized electrode location file (not necessary for LSL xdf files)
ELocs_file = [subj_basepath session filesep 'ElecLoc/MISM_MISM_' subjcode(6:end) '.sfp'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change directory to Subject's base folder
cd(subj_basepath)
global ALLEEG EEG CURRENTSET ALLCOM;

%% DETERMINE REFERENCE TYPE
MasRef = 2; % 1 = reference to mastoid / 2 = use default references
MastCh = 1; % 1 = both EXT1 and EXT2 / 2 = only Channel EXT1 / 3 = only Channel EXT2
if MasRef == 1 & MastCh == 1
    RefMs = [numChans-6 numChans-5];
elseif MasRef == 1 & MastCh == 2
    RefMs = numChans-6;
elseif MasRef == 1 & MastCh == 3
    RefMs = numChans-5;
end
%% Get BDF and CSV filenames
% Find .bdf files in directory of specified subject basepath
xdf_files = dir_list(rawEEG_inpath,'.xdf');
if isempty(xdf_files)
    error(['No .xdf files in ' rawEEG_inpath]);
    return;
end
%% Get HEAD REF from ELOCS FILE (should be embedded in xdf now)

% % Use Head models 
% if ~exist([ELocs_file(1:end-4) '_eeglabformat_final.sfp'])
%     ELocs_rawfile=ELocs_file;
%     ELocs_file=[ELocs_rawfile(1:(length(ELocs_rawfile)-4)) '_eeglabformat.sfp'];
%     if method==0
%         standardcapelocs=['/home/stepeter/HNL_Cluster/share/data3/stepeter/Steve_Code/Biosemi/cap128.sph'];
%     elseif method==1
%         standardcapelocs=['/share/data3/stepeter/Steve_Code/Biosemi/cap128.sph'];
%     end
%     collectionsys='Biosemi';
%     [ELocs_file_eeglab]=hnl_convertZebris128(ELocs_rawfile,ELocs_file, standardcapelocs,collectionsys);
% else
%     ELocs_file_eeglab = [ELocs_file(1:end-4) '_eeglabformat_final.sfp'];
% end
%% LOAD BDF FILES, IMPORT GAIT EVENTS, FILTER, AND SAVE SETS
% Open each bdf file, store rising edge events, open new events, sync
% events with CSV data, get biomechanics from CSV, store into ALLEEG
tic
%eval(['matlabpool open ' num2str(num_procs)])
%matlabpool open
for i = 1:length(xdf_files)
    EEG=[];
    
    EEG = eeg_load_xdf([rawEEG_inpath filesep xdf_files{i}]);
    EEG=eeg_checkset(EEG);
     
    %Handles case of recording from all channels in Biosemi app
    if EEG.nbchan==281
        EEG=moveElecLocs(EEG,162:169,258:265);
        EEG=pop_select(EEG,'nochannel',130:257); %remove unused channels
        if MasRef == 1
            RefMs=RefMs-length(130:257); %change mastoid channels
        end
    end
    %Chanlocs automatically added during recording
    %EEG = pop_chanedit(EEG,'load',{ELocs_file_eeglab,'filetype','autodetect'});    
    
    [EEG]=trigChan2syncEvs(EEG); %get sync signal (and remove trigger channel and AUX channels)
    if MasRef == 1
        RefMs=RefMs-1-16; %subtract of trigger channel and AUX channels
    end
    EEG=removeExtraKeyboardEvents(EEG); %remove extra keyboard events from LSL
    
    % Store trial types for event renaming
    EEG.trialtypes=MISM_condAbbrevs(xdf_files{i}(1:end-4));
    
    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset before append to ALLEEG
    
    EEG_raw(i) = EEG;
    EEG_raw(i).setname = xdf_files{i}(1:end-4);
    EEG_raw(i).subject = subjcode;
    EEG_raw(i).filename = xdf_files{i};
    EEG_raw(i).filepath = rawEEG_inpath;
    
    fprintf('set %d done\n',i)
    
    %Common-average reference
    if MasRef == 1
        EEG_raw(i)=pop_reref(EEG_raw(i),RefMs,'keepref','on');
        display('MASTOIDS ARE USED AS REFERENCE IN THIS ANALYSIS!! ')
    elseif MasRef == 2
        EEG_raw(i)=pop_reref(EEG_raw(i),[]);
    end

    if removeExternals==1
        EEG_raw(i)=pop_select(EEG_raw(i),'nochannel',(EEG_raw(i).nbchan-7):EEG_raw(i).nbchan); %remove externals
    end
end
%matlabpool close 
toc

display('Done!! :)')
end