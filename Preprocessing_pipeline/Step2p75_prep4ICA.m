%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/media/stepeter/Local_Data/VR_connectivity/Code/EEG/'));
%     addpath(genpath('/media/stepeter/Local_Data/VR_rotationsDenoise/Code/EEG/'));
    startUpEEGLAB('close');
%     addpath(genpath('/media/stepeter/Local_Data/VR_connectivity/Code/EEG/'));
end
close all;
%Specify parameters
%[1,6,7,8,10,12,15,16,17,31]
subjnum=33; %[3:5,9,11,13,19:20,22:30 32:33]%subject number
RootData='/media/stepeter/Local_Data/VR_connectivity/Data/';

subjcode = ['WMISM_' num2str(subjnum)];
group='all'; %'mismatches'; %'pull';
subj_basepath = [RootData subjcode filesep]; % Top level subject folder
EEGsets_outpath = [subj_basepath 'EEG_sets']; % Where you want to save your .set files
ICA_path = [subj_basepath 'ICA_Stuff'];
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
global EEG ALLEEG
%% (No) Frame rejection (AMICA already does some frame rejection and have run ASR on data)
%Load merged set
EEG=pop_loadset('filename', ['Merge_' group '_CAR_v2.set'], 'filepath', EEGsets_outpath);
load(['/media/stepeter/Local_Data/VR_rotationsDenoise/Data/' 'WMISM_' num2str(subjnum) filesep 'ICA_Stuff' filesep 'files_AMICA' group filesep group '_ICA_PREP_1']);

%Calculate number of frames needed to keep for k values of 60 and 30
target_nframes = length(ICA_STRUCT.good_chans)^2 * 60 / EEG.pnts; 
fprintf('You need to keep %.4f frames to achieve k=60\n',target_nframes)

target_nframes2 = length(ICA_STRUCT.good_chans)^2 * 30 / EEG.pnts; 
fprintf('You need to keep %.4f frames to achieve k=30\n',target_nframes2)

% Don't reject bad frames
% [ICA_STRUCT] = frame_rej_hjh(EEG,ICA_STRUCT);

%% update ICA Struct
%Rereference
ICA_STRUCT.good_cap=1:EEG.nbchan; ICA_STRUCT.good_chans=1:EEG.nbchan;
ICA_STRUCT.ref = 'average'; EEG.setname=ICA_STRUCT.associated_set;
EEG = update_EEG(EEG, ICA_STRUCT);

%Don't remove bad frames
% EEG_temp=EEG;
% EEG_temp.data(:,ICA_STRUCT.rej_frame_idx) = [];
% EEG_temp.times = EEG_temp.times(1:size(EEG_temp.data,2));

%Plot power spectrum for frame rejected data
% Vfig=figure;
% pop_spectopo(EEG_temp, 1, [0 EEG_temp.times(end)], 'EEG' , 'freq', [10 25 42], 'freqrange',[0 80],'electrodes','off');
% saveas(Vfig,[EEGsets_outpath filesep 'PowerSpectrum' group '_frame_rej'],'fig');

%Plot final channel distribution after rejections
FigCh = figure;
topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
saveas(FigCh,[ICA_path_out filesep 'Channels_' subjcode '_' group],'fig') %Tsk(TaskType).name

%Save ICA_STRUCT with bad channels and bad frames
ICA_STRUCT.filename = [group '_ICA_PREP_2.mat'];
ICA_STRUCT.filepath = ICA_path_out;
save([ICA_path_out filesep group '_ICA_PREP_2'], 'ICA_STRUCT')


%% SAVE FLOAT FILE FOR ICA

% Save float file for running ICA
tmpdata = EEG.data;
% tmpdata(:,ICA_STRUCT.rej_frame_idx) = []; %remove bad frames
disp('Converting data to double...');
tmpdata = double(tmpdata);
filename = [subjcode group '.fdt'];
data2write = tmpdata;
tmpdatafile = fullfile(ICA_path,filename);
disp(['Writing data file: ' tmpdatafile]);
floatwrite(data2write, tmpdatafile);
% %Write file to data3 as well
% tmpdatafile3 = fullfile(ICA_path3,filename);
% disp(['Writing data file: ' tmpdatafile3]);
% floatwrite(data2write, tmpdatafile3);
num_chans = size(tmpdata,1); % Number of channels used
disp(['Number of channels: ' num2str(num_chans)]);
num_frames = size(tmpdata,2); % Number of frames used
disp(['Number of frames: ' num2str(num_frames)]);
k = num_frames / num_chans^2; % K value achieved
disp(['Final k value: ' num2str(k)]);

% disp(['Data rank: ' num2str(rank(data2write))]);
%% Define parameter file for ICA

% RootCluster = '/sshfs/kin-hnl-frontend-stepeter/share/data4/stepeter/Projects/WMISM/Data';
RootFlux = '/ufrc/dferris/s.peterson/Projects/VR_connectivity/Data';
% pGran = [RootCluster '/' subjcode '/ICA_Stuff/'];
pFlux = [RootFlux '/' subjcode '/ICA_Stuff/'];

% % Save CUDAICA parameter file to subject folder
% saveCUDAICAparamfile(ICA_path,pGran,filename,num_chans,num_frames,group);
% %Create CUDAICA folder
% cudaicaPath=[subj_basepath 'ICA_Stuff' filesep 'files_CUDAICA' group];
% if ~exist(cudaicaPath)
%     mkdir(cudaicaPath)
% end

% Save AMICA .param file to subject folder
maxthreads=4; %24
doPCA=1; numPCs=90;
if doPCA==1
    k = num_frames / numPCs^2; % K value achieved
    disp(['Final k value (with PCA): ' num2str(k)]);
end
saveAMICAparamfile(ICA_path,pFlux,filename,num_chans,num_frames,doPCA,numPCs,maxthreads,group);
% saveAMICAparamfile(ICA_path3,pFlux,filename,num_chans,num_frames,doPCA,numPCs,maxthreads,group);

% Save qsub parameter file to subject folder (for running on Flux)
email='stepeter@umich.edu';
account='dferris';
procs=4; mem='8000mb'; walltime='24:00:00';
% saveQSUBfile(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group);
% saveQSUBfile(ICA_path3,pFlux,filename,email,account,procs,mem,walltime,group);
saveSLURMfile(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group)

% disp(['DATA RANK: ' num2str(rankEEG)]);

disp('Step 2 Done :)');
disp('Data has been channel/frame rejected and all parameter files have been generated.');
disp('Run ICA (AMICA or CUDAICA) by submitting .qsub or .sc file and then proceed to step 3 to fit dipoles from ICs.');
