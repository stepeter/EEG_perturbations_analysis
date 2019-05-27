%Matlab function used on supercomputer to compute phase randomized surrogate statistics. Performed this in batches of 50 for 200 total per subject. Each batch of 50 took many hours to run, so this can take weeks to complete without extra computing resources.

function calcSurrStats1_hiperGator(orderNum,subFolderName,subjName,jobID)
copyfile('~/Projects/VR_connectivity/Code/compileTest/eeg_optionsbackup.txt',['/scratch/local/' jobID '/mcr_cache/.mcrCache8.3/'])
%% Run surrogate statistics (very slow :/)
folderpath=['~/Projects/VR_connectivity/Data/GroupConn_aveScalp/regularGroupConn/' subFolderName];
loadpath=['~/Projects/VR_connectivity/Data/GroupConn_aveScalp/regularGroupConn/' subFolderName];
% number of null distribution samples to draw
NumSamples = 50;
%Load in data set
    EEG=pop_loadset('filepath',loadpath,'filename',['S_' subjName '.set']);
    EEG.CAT.configs.stat_surrogateGen=struct([]);

    % obtain the bootstrap distributions for each condition
    EEG = pop_stat_surrogateGen(EEG,'nogui', ...
         'modelingApproach', EEG.CAT.configs.est_fitMVAR, ...
         'connectivityModeling',EEG.CAT.configs.est_mvarConnectivity, ...
         'mode',{'PhaseRand' 'nperms' NumSamples}, ...
         'autosave',[], ...
         'verb',1);

    %Result is EEG.CAT.PConn.RPDC is (num comps x num comps x freq x time x num samples)
    %rpdcPConn=EEG.CAT.PConn.RPDC;
    %save([folderpath filesep 'phaseRand' filesep 'rpdcPConn_' condName '_' num2str(orderNum) '.mat'],'rpdcPConn');
    ddtf08PConn=EEG.CAT.PConn.dDTF08;
    save([folderpath filesep 'phaseRand' filesep 'ddtf08PConn_' condName '_' num2str(orderNum) '.mat'],'ddtf08PConn');
    %cohPConn=EEG.CAT.PConn.Coh;
    %save([folderpath filesep 'phaseRand' filesep 'cohPConn_' condName '_' num2str(orderNum) '.mat'],'cohPConn');
    %icohPConn=EEG.CAT.PConn.iCoh;
    %save([folderpath filesep 'phaseRand' filesep 'icohPConn_' condName '_' num2str(orderNum) '.mat'],'icohPConn');
    %gpdcPConn=EEG.CAT.PConn.GPDC;
    %save([folderpath filesep 'phaseRand' filesep 'gpdcPConn_' condName '_' num2str(orderNum) '.mat'],'gpdcPConn');
    %ffddtfPConn=EEG.CAT.PConn.ffDTF;
    %save([folderpath filesep 'phaseRand' filesep 'ffddtfPConn_' condName '_' num2str(orderNum) '.mat'],'ffddtfPConn');
    %ggcPConn=EEG.CAT.PConn.GGC;
    %save([folderpath filesep 'phaseRand' filesep 'ggcPConn_' condName '_' num2str(orderNum) '.mat'],'ggcPConn');
    SPConn=EEG.CAT.PConn.S;
    save([folderpath filesep 'phaseRand' filesep 'SPConn_' condName '_' num2str(orderNum) '.mat'],'SPConn');
end
