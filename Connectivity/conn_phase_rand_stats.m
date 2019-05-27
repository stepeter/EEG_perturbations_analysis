%Function used to compute phase randomization surrogate statistics for connectivity data. Performed on multiple subjects in parallel on supercomputer to drastically cut down compute time.


% orderNum,subFolderName,subjName
loadpath='/usr/local/VR_connectivity/Data/regularGroupConn/allSVZ';
for subjName=[1:7 9:13 15:17 19:20 22:33]
%     subjName=1;

    NumSamples=200; %50;
    EEG=pop_loadset('filepath',loadpath,'filename',['S' num2str(subjName) '.set']);
    EEG.CAT.configs.stat_surrogateGen=struct([]);

    % obtain the bootstrap distributions for each condition
    EEG = pop_stat_surrogateGen(EEG,'nogui', ...
        'modelingApproach', EEG.CAT.configs.est_fitMVAR, ...
        'connectivityModeling',EEG.CAT.configs.est_mvarConnectivity, ...
        'mode',{'PhaseRand' 'nperms' NumSamples}, ...
        'autosave',[], ...
        'verb',1);


    ddtf08PConn=EEG.CAT.PConn.dDTF08;
    save([loadpath filesep 'phaseRand' filesep 'ddtf08PConn' num2str(subjName) '.mat'],'ddtf08PConn');
    SPConn=EEG.CAT.PConn.S;
    save([loadpath filesep 'phaseRand' filesep 'SPConn' num2str(subjName) '.mat'],'SPConn');
end