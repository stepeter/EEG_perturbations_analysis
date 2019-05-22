function [EEG_raw,subj_basepath]=EEG_load_beforePREP(subjnum,RootData,numChans,refType,group,removeExternals)
%Loads in EEG raw files (xdf in this case) and formats them as needed, also
%performing the desired referencing ('common','mastoid_both','mastoid_left'
%'mastoid_right')
switch refType
    case 'common'
        MasRef=2; %1 = reference to mastoid / 2 = use default references
    case 'mastoid_both'
        MasRef=1;
        MastCh = 1; % 1 = both EXT1 and EXT2 / 2 = only Channel EXT1 / 3 = only Channel EXT2
    case 'mastoid_left'
        MasRef=1;
        MastCh = 2;
    case 'mastoid_right'
        MasRef=1;
        MastCh = 3;
    otherwise
        error('Unknown reference type requested!');
end


% removeExternals=1; %decide if to remove externals (1 - yes; 0 - no)
%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjcode = ['WMISM_' num2str(subjnum)]; %define subjcode

subj_basepath = [RootData subjcode filesep]; % Top level subject folder
rawEEG_inpath = [subj_basepath 'RawEEG']; % Where your raw .bdf files are stored

% Digitized electrode location file (not necessary for LSL xdf files)
ELocs_file = [subj_basepath 'ElecLoc/WMISM_WMISM_' num2str(subjnum) '.sfp'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change directory to Subject's base folder
cd(subj_basepath)
global ALLEEG EEG CURRENTSET ALLCOM;

%% DETERMINE REFERENCE TYPE
if MasRef == 1 & MastCh == 1
    RefMs = [numChans-7 numChans-6];
elseif MasRef == 1 & MastCh == 2
    RefMs = numChans-7;
elseif MasRef == 1 & MastCh == 3
    RefMs = numChans-6;
end
%% Get BDF and CSV filenames
% Find .bdf files in directory of specified subject basepath
xdf_files = dir_list(rawEEG_inpath,'.xdf');
if isempty(xdf_files)
    error(['No .xdf files in ' rawEEG_inpath]);
    return;
end

%% Get HEAD REF from ELOCS FILE (to replace locations from xdf)
[ELocs_file_eeglab]=createSFPcorrect(ELocs_file,numChans-8,[subj_basepath 'ElecLoc' filesep]);

%Find which files to load
sets2Load=mergeEEGsetsMISMCombs(group,subjnum);

%% LOAD BDF FILES, IMPORT GAIT EVENTS, FILTER, AND SAVE SETS
% Open each bdf file, store rising edge events, open new events, sync
% events with CSV data, get biomechanics from CSV, store into ALLEEG
tic
for i = 1:length(sets2Load)
    EEG=[];
    EEG = eeg_load_xdf([rawEEG_inpath filesep sets2Load{i} '.xdf']);
    EEG=eeg_checkset(EEG);
     
    %Handles case of recording from all channels in Biosemi app
    if EEG.nbchan==281
        EEG=moveElecLocs(EEG,162:169,258:265);
        EEG=pop_select(EEG,'nochannel',130:257); %remove down to 136 chans
    end 
    
    [EEG]=trigChan2syncEvs(EEG); %get sync signal (and remove trigger channel and AUX channels)
    EEG=removeExtraKeyboardEvents(EEG); %remove extra keyboard events from LSL (when key is released)
    
    %Add in correct channel locations
    EEG = pop_chanedit(EEG,'load',{ELocs_file_eeglab,'filetype','autodetect'}); 
    
    % Store trial types for event renaming
    EEG.trialtypes=MISM_condAbbrevs([sets2Load{i}]);
    
    EEG=eeg_checkset(EEG,'eventconsistency', 'makeur') %always checkset before append to ALLEEG
    
    EEG_raw(i) = EEG;
    EEG_raw(i).setname = sets2Load{i};
    EEG_raw(i).subject = subjcode;
    EEG_raw(i).filename = [sets2Load{i} '.xdf'];
    EEG_raw(i).filepath = rawEEG_inpath;
    
    fprintf('set %d done\n',i)
    
    %Remove externals if desired
    if removeExternals==1
        EEG_raw(i)=pop_select(EEG_raw(i),'nochannel',(EEG_raw(i).nbchan-7):EEG_raw(i).nbchan); %remove externals
    else
        EEG_raw(i)=pop_select(EEG_raw(i),'nochannel',[(EEG_raw(i).nbchan-7):(EEG_raw(i).nbchan-6) EEG_raw(i).nbchan]); %remove exts. 1,2,and 8
    end
end
toc
end
