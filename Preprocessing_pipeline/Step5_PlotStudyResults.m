%Plots ERSPs, dipoles & centroids together, or pulls up gui
% Assumes the study is already loaded and clustered
close all; clc; %clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define input and output paths
% If subject code for convenience
RootData = '/media/stepeter/Local_Data/VR_connectivity/Data';%'/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Data';
savePlots=0; %1 - save ERSP plots; 0 - don't save
studyname='all_5222018_pull_noDir_v2';
savePath = [RootData filesep 'STUDY_results_' STUDY.filename(1:end-6) filesep];
plotERSPs=1; %1 - plot ERSPs; 0 - no ERSPs
plotDipCentroids=0; %1 - plot dipoles and centroids; 0 - don't plot
pullUpClustGUI=0; %1 - pulls up study cluster gui at end; 0 - don't pull it up
saveSTUDY=0; %save study at end of processing
loadSTUDY=0; %1 - loads in study; 0 - study already loaded
% designs=2:4;
%Specific plot parameters
clusters_to_plot=3:10; %[3:6 8]; %[3:7 9:13]; %[3:5 10 13 15:17]; %[2:8]; 
condNums2plot=1:2; %1:4; %1:2; %[1 2];
cond_labels={'pull_Stn','pull_Wlk'}; %{'L_pull_Stn','R_pull_Stn','L_pull_Wlk','R_pull_Wlk'};%{'Stand Viz CW','Stand Viz CCW','Walk Viz CW','Walk Viz CCW'}; %{'Stand_Viz','Walk_Viz'}; %%{'Stand','Walk','Pre','Train1','Train2','Train3','Post'}; %{'Stand','Walk','Pre','Train1','Train2','Train3','Post'}; %,'Balance'};
tlim = [-500 1500];
flim = [4 100];
clim = [-1 1];
basetime = [-1000 0]; %-500 0];
evPlotLines = [0 1000];
colors2use={[1 0 0], [0 1 0], [0 0 1], [1 0 1], [0 1 1], [.8549 .6471 .1255], [1 .4902 0], ... %[1 1 0], [1 .4902 0], ...
    [1 .0784 .5765], [.8 0 1], [.6 0 0], [0 0 0]};
centrDipsTogether=0; %1 - plot centroids and dipoles together; 0 - only dipoles
%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/home/stepeter/HNL_Cluster/share/data4/stepeter/WMISM/Code/EEG/'));
    startUpEEGLAB('close');
end

% Add necessary dependent paths
addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/Grant/common/'))

% Edit eeg_options.m file
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);

cd(RootData)
if ~isdir(savePath)
    mkdir(savePath)
end
if loadSTUDY==1
    %Load in study
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    [STUDY ALLEEG] = pop_loadstudy('filename', [studyname '.study'], 'filepath', RootData);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
end
%% Get ERSPs, bootstrap, and then plot them
if plotERSPs==1      
    %Plots time-frequency spectral power
ERSP_bootPlotter(STUDY,ALLEEG,tlim,flim,clim,clusters_to_plot,cond_labels,savePlots,savePath,basetime,evPlotLines,condNums2plot);
%     ERSP_bootPlotter_designs(STUDY,ALLEEG,tlim,flim,clim,clusters_to_plot,cond_labels,savePlots,savePath,basetime,evPlotLines,designs);
%     Spec_bootPlotter_designs(STUDY,ALLEEG,tlim,flim,clim,clusters_to_plot,cond_labels,savePlots,savePath,basetime,evPlotLines,designs);
end

%Plot dipoles and centroids
if plotDipCentroids==1
    diplotfig(STUDY,ALLEEG,clusters_to_plot, colors2use(1:length(clusters_to_plot)),savePath,centrDipsTogether);
end

%Pull up cluster visualization interface
if pullUpClustGUI==1
    [STUDY]=pop_clustedit(STUDY,ALLEEG);
end

% Save STUDY
if saveSTUDY==1
    [STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename',studyname,'filepath',RootData);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
end

close all;
disp('Step 5 Finished!');
disp('Study results have been plotted.');
disp('This is the end of the standard processing pipeline! :)');