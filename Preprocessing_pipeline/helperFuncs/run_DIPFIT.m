function run_DIPFIT(subjnum,RootData,group)

% Add EEGLAB path
% addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b_octave_rem'))
% rmpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b/functions/octavefunc/signal')

% Add necessary dependent paths
addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/Grant/common/'))

% Edit eeg_options.m file
pop_editoptions( 'option_storedisk', 1, 'option_savetwofiles', 1, ...
    'option_single', 1, 'option_memmapdata', 0, ...
    'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%% SET THESE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define input and output paths
% If subject code for convenience
% RootData = '/home/stepeter/HNL_Cluster/share/data3/stepeter/EEG_data_processing/CSR_Data/';
subjcode = ['WMISM_' num2str(subjnum)];

subj_basepath = [RootData subjcode filesep]; % Top level subject folder
EEGsets_outpath = [subj_basepath 'EEG_sets']; % Where you want to save your .set files
ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];

% Digitized electrode location file
ELocs_file = [subj_basepath 'ElecLoc/WMISM_WMISM_' subjcode(7:end) '.sfp']; % Where you want to save your ICA files

cd(subj_basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ALLEEG EEG CURRENTSET ALLCOM;

% Use Head models 
ELocs_file_eeglab = [ELocs_file(1:end-4) '_eeglabformat_final.sfp'];

% Load EEG
EEG=pop_loadset([EEGsets_outpath filesep 'Merge_' group '_CAR_v2.set']);
load([ICA_path_out filesep group '_ICA_PREP_2.mat']);

% Load  weights and spheres
icatype = 'amica';
EEG = eeg_checkset(EEG); EEG.icachansind=[];

[ICA_STRUCT.weights, ICA_STRUCT.sphere] = load_wgts_spheres(EEG, icatype, ICA_path_out);

%Reformat sphere matrix if used doPCA
diffSphereWeights=size(ICA_STRUCT.sphere,1)-size(ICA_STRUCT.weights,2);
if diffSphereWeights>0
    ICA_STRUCT.sphere=ICA_STRUCT.sphere(1:(end-diffSphereWeights),:);
end

ICA_STRUCT.ref='mast'; EEG.setname=ICA_STRUCT.associated_set;
EEG=update_EEG(EEG,ICA_STRUCT);
% EEG.icaweights=ICA_STRUCT.weights;
% EEG.icasphere=ICA_STRUCT.sphere;

dipfitelocdata='digitizer_zebris';
fiducial_labels = {'LPA', 'Nz','RPA'};

%Update EEG with ICA_STRUCT
EEG = eeg_checkset(EEG,'ica'); %create EEG.icaact
CURRENTSET = 1;
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
% eeglab redraw;

%% Dipole Fit
subj_elocsfile = [ELocs_file_eeglab(1:end-4) '.fid'];
[eeglab_pathstr, name, ext] = fileparts(which('eeglab'));
% elocfile_standard='/share/data3/hjhuang/code/elocs_standard/cap256+8neck.sph';
eloc_BEM_HDMfile=[eeglab_pathstr '/plugins/dipfit2.3/standard_BEM/standard_vol.mat'];
eloc_BEM_MRIfile=[eeglab_pathstr '/plugins/dipfit2.3/standard_BEM/standard_mri.mat'];
eloc_BEM_CHANfile=[eeglab_pathstr '/plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc'];

if strcmp(dipfitelocdata,'digitizer_zebris') == 1 
    ansdipfit = 0;
    if isfield(ICA_STRUCT,'dipfit') == 1
        [ansdipfit,ok] = listdlg('PromptString','Dipfit already exists, would you like to run dipfit anyway?',...
                    'SelectionMode','single','ListString',{'Yes'; 'No'},...
                    'InitialValue', 1,'ListSize',[400 150]);
    end
    if ansdipfit == 1 ||  isfield(ICA_STRUCT,'dipfit') == 0 
        if EEG.nbchan <= 256
            EEG =  pop_dipfit_settings(EEG,'hdmfile',eloc_BEM_HDMfile,...
                    'mrifile',eloc_BEM_MRIfile,'chanfile',eloc_BEM_CHANfile,...
                    'coordformat','MNI','chansel', [1:EEG.nbchan]);
        end
    end
end

% %Ensure calling functions from correct path
% addpath(genpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b_octave_rem'))
% rmpath('/home/stepeter/HNL_Cluster/share/data3/stepeter/eeglab13_0_1b/functions/octavefunc/signal')

tic
if ansdipfit == 1 ||  isfield(ICA_STRUCT,'dipfit') == 0
    if strcmp(dipfitelocdata,'digitizer_zebris') == 1  

        fiducials = readlocs(subj_elocsfile,'filetype','sfp');

        fiducials = rmfield(fiducials,setdiff(fieldnames(fiducials),fieldnames(EEG.chanlocs)));

        [chanlocs_out, EEG.dipfit.coord_transform] = coregister(fiducials,EEG.dipfit.chanfile,...
            'mesh',EEG.dipfit.hdmfile,'warp',fiducial_labels, 'manual', 'off') 
    end

    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    %fit dipoles
    EEG = pop_multifit(EEG, [1:EEG.nbchan] , 'threshold',100, 'plotopt',{ 'normlen', 'on'});
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    %save results to ICA_STRUCT
    ICA_STRUCT.dipfit = EEG.dipfit;
    save([ICA_path_out filesep group '_ICA_DIPFIT_3'], 'ICA_STRUCT')
end
toc
end
