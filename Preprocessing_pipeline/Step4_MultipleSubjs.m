close all; clc;
%Multiple Subjects Analysis

%Run after DIPFIT and epoching. This puts all good dipoles into a study for
%clustering and ERSP plotting.

%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define input and output paths
% If subject code for convenience
RootData = '/media/stepeter/Local_Data/VR_connectivity/Data/';
subjcode = 'WMISM_'; %don't include number
subjsToAnalyze=[1:13 15:17 19:20 22:33];
group='all';
dateTime=clock;
studyname=[group '_' num2str(dateTime(2)) num2str(dateTime(3)) num2str(dateTime(1)) '_pull_noDir_v2.study'];
trigs2Analyze={'pull_Stn','pull_Wlk','M_on_SVZ','M_on_WVZ'}; %{'L_pull_Stn','R_pull_Stn','L_pull_Wlk','R_pull_Wlk'};%{'M_on_CW_SVZ','M_on_CCW_SVZ','M_on_CW_WVZ','M_on_CCW_WVZ'};

%Options to choose for this
createStudy=1; %1 - create a study from scratch; 0 - load study
redrawAfterStudyCreated=1; %1 - eeglab redraw and break so can look at study; 0 - continue on
preclust=1; %1 - precluster study; 0 - no preclustering
clust=1; %1 - precluster study; 0 - no preclustering
precomp_nonERSPs=0; %1 - pulls up precompute gui; 0 - no gui
precomp_ersp=0; %1 - precompute ersps; 0 - don't precompute ersps
erspComp='full'; %'light' - quicker computation; 'full' - with usual parameters (takes longer)
showClusterPlotGUI=1; %1 - show cluster plot gui at end; 0 - don't show it (and clear study)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mergeEpch_path = [RootData 'STUDYsets_pullspullsOnset_noDir_v2']; %'STUDY_sets_' group 'Pulls'];
%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_connectivity/Code/EEG/'));
    startUpEEGLAB('close');
end

%Fieldtrip will mess this part up, so make sure its removed from path
rmpath(genpath('/home/stepeter/Documents/eeglab13_0_1b_octave_rem/external/fieldtrip-partial'));

%Add necessary dependent path
addpath(genpath('/media/stepeter/New Volume/Grant/common/'))

% Edit eeg_options.m file
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);

cd(RootData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY STUDY;
%% Load Sets & Create Study
%Use final epoched sets from DIPFIT script (create separate Study for each
%subject
EEG=[]; ALLEEG=[]; CURRENTSET=[]; ALLCOM=[]; CURRENTSTUDY = []; STUDY = [];

%Create Study
if createStudy==1
    % [STUDY ALLEEG]=pop_study([],[],'gui','off');
    indx=1; 
    for frodo=subjsToAnalyze
        [STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name',studyname,'commands',...
        	{{'index' indx 'load' [mergeEpch_path filesep 'Merge_epch_' group '_' num2str(frodo) '.set'] }...
        	{'inbrain' 'on' 'dipselect' 0.15}},'updatedat','on','savedat','on' );
        indx=indx+1;
        subjCell{1,frodo}=[subjcode num2str(frodo)];
    end
    %Remove empty cells from subjCell
    subjCell=subjCell(~cellfun('isempty',subjCell));
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
    %Make STUDY design
    [STUDY] = std_makedesign(STUDY, ALLEEG, 1, 'subjselect', subjCell,'variable1','type','values1', trigs2Analyze);
    
    % Save STUDY
    [STUDY ALLEEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',RootData);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

    if redrawAfterStudyCreated==1
        eeglab redraw
        return;
    end
else
    %Load STUDY (assuming it exists)
    [STUDY ALLEEG]=pop_loadstudy('filename',studyname,'filepath',RootData);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
end

%Precompute ERSPs
if precomp_ersp==1
    if TW==0
        switch erspComp
            case 'full'
                tic
                [STUDY ALLEEG] = mod_std_precomp_v10_2_5_5a(STUDY, ALLEEG, 'components','ersp','on','itc','on','erspparams',{'cycles',[3 0.5],'alpha',0.05, 'padratio',2,'savetrials','off','baseline',nan}, 'recompute','on');
                toc
            case 'light'
                tic
                %Parameters for quicker computation (still works quite well)
                [STUDY ALLEEG] = mod_std_precomp_v10_2_5_5a(STUDY, ALLEEG, 'components','ersp','on','itc','off','erspparams',{'cycles',[3 0.9],'freqs',[3 256],'nfreqs',100,'alpha',0.05,'freqscale','log','savetrials','off','baseline',nan}, 'recompute','on');
                toc
            otherwise
                error('Incorrect case for erspComp!');
        end
    elseif TW==1
        %Run ERSPs (timewarped)
        warps=zeros(length(ALLEEG),length(ALLEEG(1,1).warpto));
        for i=1:length(ALLEEG)
            warps(i,:)=ALLEEG(1,i).warpto;
        end
%         warpingvalues=median(warps);
        roundNear=50; %round numbers to the closest multiple of this value
        warpingvalues=round(median(warps)/roundNear)*roundNear;
        
        switch erspComp
            case 'full'
                tic
                [STUDY ALLEEG] = mod_std_precomp_v10_2_5_5a(STUDY, ALLEEG, 'components','ersp','on','itc','on','erspparams',{'cycles',[3 0.5],'alpha',0.05, 'padratio',2,'savetrials','off','baseline',nan,'timewarp','subject tw matrix','timewarpms', warpingvalues}, 'recompute','on');
                toc
            case 'light'
                tic
                %Parameters suggested by Makoto
                [STUDY ALLEEG] = mod_std_precomp_v10_2_5_5a(STUDY, ALLEEG, 'components','ersp','on','itc','off','erspparams',{'cycles',[3 0.9],'freqs',[3 256],'nfreqs',100,'alpha',0.05,'freqscale','log','savetrials','off','baseline',nan,'timewarp','subject tw matrix','timewarpms', warpingvalues}, 'recompute','on');
                toc
            otherwise
                error('Incorrect case for erspComp!');
        end
    end
end

% Save STUDY again
[STUDY ALLEEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',RootData);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

%Precompute other measures
if precomp_nonERSPs==1
    [STUDY ALLEEG] = pop_precomp(STUDY, ALLEEG,'components');
end

%Precluster components
if preclust==1
    [STUDY ALLEEG]=pop_preclust(STUDY,ALLEEG);
end

%Cluster components
if clust==1
    [STUDY]=pop_clust(STUDY,ALLEEG);
end

if showClusterPlotGUI==1
    [STUDY]=pop_clustedit(STUDY,ALLEEG); %pulls up cluster visualization interface
else
    clear EEG;  clear EEG_orig; clear ALLEEG; clear CURRENTSET; clear STUDY;
end

% Save STUDY final time
[STUDY ALLEEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',RootData);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

disp('Step 4 finished!')
disp('Data from multiple subjects have been combined into a STUDY.');
disp('Proceed to Step 5 to plot STUDY results.');
