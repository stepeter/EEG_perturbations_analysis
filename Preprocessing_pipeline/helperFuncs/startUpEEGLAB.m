function startUpEEGLAB(showGUI)
%%Starts up EEGLAB with all paths. Use showGUI option to decide whether to
%%leave the GUI up or automatically close it after startup
%%Input: - showGUI:  'close' to close the GUI automatically after startup
close all; clc;
method=0; %0 - local or 1 - Cluster
exp='MISM';
% Add EEGLAB path
if method==0
    %     addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b_octave_rem'))
    %     rmpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b/functions/octavefunc/signal')
    cd('/home/stepeter/Documents/eeglab13_5_4b');
    % Add necessary dependent paths
    addpath(genpath('/home/stepeter/Documents/Grant/'))
    %addpath(genpath(['/home/stepeter/HNL_Cluster/share/data3/stepeter/Projects/' exp '/Code/EEG/']))
    %addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_5_4b/plugins/PrepPipeline0.5/'));
elseif method==1
    %     addpath(genpath('/share/data3/stepeter/eeglab13_0_1b_octave_rem'))
    %     rmpath('/share/data3/stepeter/eeglab13_0_1b/functions/octavefunc/signal')
    cd('/share/data3/stepeter/eeglab13_5_4b');
    % Add necessary dependent paths
    addpath(genpath('/share/data3/stepeter/Grant/'))
    %addpath(genpath(['/share/data3/stepeter/Projects/' exp '/Code/EEG/']))
    addpath(genpath('/share/data3/stepeter/eeglab13_5_4b/plugins/PrepPipeline0.5/'));
end
eeglab
% Edit eeg_options.m file
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 0, 'option_memmapdata', 0, ...
    'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);
if nargin==1 && strcmp(showGUI,'close')
    close all;
end
end
