%%Group connectivity plots
%Load all_1222018.study and then run script

%Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(genpath('/usr/local/VR_connectivity/Code/EEGconn/'));
    startUpEEGLAB('close');
end

%Plot results by placing into large connectivity structure
RootData='/usr/local/VR_connectivity/Data/';
savePath='/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/';
group='all'; usefdr = 0;

goodClusts=3:10; %clusterNums;
numChans=128;
% subjectInds=trueSubjNumbers; %subjects in this STUDY
load('/usr/local/VR_connectivity/Data/STUDYtopo.mat');
load('/usr/local/VR_connectivity/Data/cluster_1222018.mat'); %_pruned.mat');
clustDat=cluster; %cluster_1222018_pruned;
subjectInds=[1:7 9:13 15:17 19:20 22:33];

%Convert from STUDY indices to my EEG set subject numbers
trueSubjNumbers=subjectInds; q=0; counter2=0;
allSubjs=[1:13 15:17 19:20 22:33];
for i=allSubjs
    q=q+1;
    if any(subjectInds==q)
        counter2=counter2+1;
        trueSubjNumbers(counter2)=i;
    end
end

%Delete study
EEG=[]; ALLEEG=[]; CURRENTSET=[]; ALLCOM=[]; CURRENTSTUDY = []; STUDY = [];

%%
%Create variables for baseline-subtracted connectivity
conds={'allSVZ','allWVZ','pullsStn','pullsWalk'};
for frodo=1:length(conds)
    alpha = 0.01; %0.05;%Set a significance value
    %Put in data for bootstrap distribution later
    %Load each subject's merged EEG file and ICA_STRUCT; then calculate
    %appropriate icaweights variable based on icawinv and icasphere
    connStruct_boot=cell(16,16); %zeros(16,16,30,244);
    connStructS_boot=connStruct_boot;
    numConns_boot=zeros(16,16);
    cond=conds{frodo}; %'allSVZ';
    for i=subjectInds
        subjcode = ['WMISM_' num2str(i)];
        subj_basepath = [RootData subjcode filesep];
        EEGsets_outpath = ['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond]; %allSVZ'; %[subj_basepath 'EEG_sets'];
        ICA_path_out = [subj_basepath 'ICA_Stuff' filesep 'files_AMICA' group];
        cd(subj_basepath)

        EEG=pop_loadset('filename', ['S' num2str(i) '.set'], 'filepath', EEGsets_outpath);

        %Now figure out which clusters this subject is missing
        icaDat2Add=[]; icaEst4_crossVal=[]; icNums=[];
        for j=1:length(goodClusts)
            if ~any(allSubjs(clustDat(goodClusts(j)).sets)==i)
            else
                icNums=[icNums j];
            end
        end

        icNums=[icNums 9:16]; %add in muscles
        numConns_boot(icNums,icNums)=numConns_boot(icNums,icNums)+1;
        for j=1:length(icNums)
            for k=1:length(icNums)
                connStruct_boot{icNums(j),icNums(k)}(numConns_boot(icNums(j),icNums(k)),:,:)=squeeze(EEG.CAT.Conn.dDTF08(j,k,:,:));
                connStructS_boot{icNums(j),icNums(k)}(numConns_boot(icNums(j),icNums(k)),:,:)=squeeze(EEG.CAT.Conn.S(j,k,:,:));
%                 disp('Hi!');
            end
        end 
    end
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStruct_boot.mat'],'connStruct_boot');
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStructS_boot.mat'],'connStructS_boot');
    
    %Load connStruct data
    load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct.mat']);
    load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStructS.mat']);
    
    % Set basetime to NaN if you don't want to significant mask
    % Otherwise set basetime to the limits of your cycle, ex. stride
    baseline=[-0.5 0];
    latencies=EEG.CAT.Conn.erWinCenterTimes;
    baseidx=find(latencies>=baseline(1) & latencies<=baseline(2));baselines=[];
    for j=1:16
        for k=1:16
            %Calculate and subtract baseline
            baseVals=median(connStructS_boot{j,k}(:,:,baseidx),3);
            curr_ersp = connStructS_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
            curr_ersp=permute(curr_ersp,[2 3 1]);
            
            if usefdr==1
                curr_power = median(curr_ersp,3);
                p_fdr = run_bootstat_fdr(curr_ersp, curr_power, 'median(arg1,3);', 'shuffle', 'ERSP', 'both', 200, baseidx, alpha, 2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(p_fdr>=alpha) = 0;
            else
                %Bootstrap and significance mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
            end
            curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
            connStructS(j,k,:,:)=squeeze(curr_maskedersp);
        end
    end
    
    for j=1:16
        for k=1:16
            %Calculate and subtract baseline
            baseVals=median(connStruct_boot{j,k}(:,:,baseidx),3);
            curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
            curr_ersp=permute(curr_ersp,[2 3 1]);
            
            if usefdr==1
                curr_power = median(curr_ersp,3);
                p_fdr = run_bootstat_fdr(curr_ersp, curr_power, 'median(arg1,3);', 'shuffle', 'ERSP', 'both', 200, baseidx, alpha, 2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(p_fdr>=alpha) = 0;
            else
                %Bootstrap and significance mask
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
            end
            curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
            connStruct(j,k,:,:)=squeeze(curr_maskedersp);
        end
    end
    
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStruct_baseSub.mat'],'connStruct');
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStructS_baseSub.mat'],'connStructS');
end
%%
%Create baseline-subtracted connectivity masked for EMGbrain vs. brainEMG
%differences
conds={'allSVZ','allWVZ','pullsStn','pullsWalk'};
for frodo=1:length(conds)
    cond=conds{frodo};
    %Load connStruct data
    load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct.mat']);
    load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStruct_boot.mat']);
    
    % Set basetime to NaN if you don't want to significant mask
    % Otherwise set basetime to the limits of your cycle, ex. stride
    baseline=[-0.5 0];
    latencies=EEG.CAT.Conn.erWinCenterTimes;
    baseidx=find(latencies>=baseline(1) & latencies<=baseline(2));baselines=[];
    
    for j=1:16
        for k=1:16
%             %Calculate and subtract baseline
            baseVals=median(connStruct_boot{j,k}(:,:,baseidx),3);
            curr_ersp = connStruct_boot{j,k}-repmat(baseVals, [1, 1, length(latencies)]);
            curr_ersp=permute(curr_ersp,[2 3 1]);
%             
            %Bootstrap and significance mask
            if usefdr==1
                curr_power = median(curr_ersp,3);
                p_fdr = run_bootstat_fdr(curr_ersp, curr_power, 'median(arg1,3);', 'shuffle', 'ERSP', 'both', 200, baseidx, alpha, 2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(p_fdr>=alpha) = 0;
            else
                pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                    'label','ERSP','bootside','both','naccu',200,...
                    'basevect',baseidx,'alpha',alpha,'dimaccu',2);
                curr_ersp = median(curr_ersp,3);
                curr_maskedersp = curr_ersp;
                curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
            end
            
            %Mask again if corticomuscular
            if (j<9 && k>8) || (j>8 && k<9)
                curr_ersp = connStruct_boot{j,k}; %-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp=permute(curr_ersp,[2 3 1]);
                curr_ersp2 = connStruct_boot{k,j}; %-repmat(baseVals, [1, 1, length(latencies)]);
                curr_ersp2=permute(curr_ersp2,[2 3 1]);
                
                %Replace curr_ersp2 as baseline in 0.5 sec increments
                timeVals=[0 0.5 1 1.5];
                curr_ersp_orig=curr_ersp;
                for pp=1:(length(timeVals)-1)
                    replaceInds=find(latencies>=timeVals(pp) & latencies<=timeVals(pp+1));
                    curr_ersp(:,baseidx,:)=curr_ersp2(:,replaceInds,:);
                    
                    if usefdr==1
                        curr_power = median(curr_ersp,3);
                        p_fdr = run_bootstat_fdr(curr_ersp, curr_power, 'median(arg1,3);', 'shuffle', 'ERSP', 'both', 200, baseidx, alpha, 2);
                        curr_maskedersp_piece=curr_maskedersp(:,replaceInds,:);
                        curr_maskedersp_piece(p_fdr(:,replaceInds,:)>=alpha) = 0;
                    else
                        %Bootstrap and significance mask
                        pboot = bootstat(curr_ersp,'median(arg1,3);','boottype','shuffle',...
                            'label','ERSP','bootside','both','naccu',200,...
                            'basevect',baseidx,'alpha',alpha,'dimaccu',2);

                        curr_maskedersp_piece=curr_maskedersp(:,replaceInds,:);
                        curr_maskedersp_piece_copy=curr_maskedersp_piece;
                        curr_maskedersp_piece(curr_maskedersp_piece_copy > repmat(pboot(:,1),[1 size(curr_maskedersp_piece_copy,2)]) & curr_maskedersp_piece_copy < repmat(pboot(:,2),[1 size(curr_maskedersp_piece_copy,2)])) = 0;
                    end
                    curr_maskedersp(:,replaceInds,:)=curr_maskedersp_piece; %put masked piece back in
                end
            end
            curr_maskedersp=permute(curr_maskedersp,[3 1 2]);
            connStruct(j,k,:,:)=squeeze(curr_maskedersp);
        end
    end
    
    save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStruct_baseSub_EMGbrainMask.mat'],'connStruct');
%     save(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/connStruct_EMGbrainMask.mat'],'connStruct');
end
%%
condS={'allSVZ','allWVZ','pullsStn','pullsWalk'};
typeCondS={'brainBrain','EMGEMG','brainEMG','EMGbrain'};
for condInds=1:4
    for typeCondInds=1:4
        cond=condS{condInds}; %'pullsStn'; %close all;
        typeCond=typeCondS{typeCondInds}; %'brainBrain'; %'brainBrain','EMGEMG','brainEMG','EMGbrain' (from -> to)

        maxVal4Plot=[]; %allSVZ: brainBrain: .0013, EMGEMG: .00066, brainEMG: .0049, EMGbrain: .0002
                        %allWVZ: brainBrain: .00088, EMGEMG: .001, brainEMG: .00041, EMGbrain: .00034
                        %pullsStn: brainBrain: .0019, EMGEMG: .0026, brainEMG: .0005, EMGbrain: .00036
                        %pullsWalk: brainBrain: .0018, EMGEMG: .0026, brainEMG: .00053, EMGbrain: .00072
        if any(strcmpi(typeCond,{'brainEMG','EMGbrain'}))
            load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_baseSub_EMGbrainMask.mat']);
        else
            load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStruct_baseSub.mat']);
        end
        load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStructS_baseSub.mat']);

        %Set diagonal to zero
        for i=1:size(connStruct,1)
            connStruct(i,i,:,:)=0;
        end


        baseline=[-0.5 0];
        latencies=EEG.CAT.Conn.erWinCenterTimes;
        baseidx=find(latencies>=baseline(1) & latencies<=baseline(2));

        EEG=pop_loadset('/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/allSVZ/S4.set');
        switch typeCond
            case 'brainBrain'
                EEG.CAT.Conn.dDTF08=[];
                EEG.CAT.Conn.dDTF08=connStruct(1:8,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
                % EEG.CAT.Conn.dDTF08=EEG.CAT.Conn.dDTF08-repmat(baseVals,[1 1 1 length(latencies)]);
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(1:8,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.S(:,:,:,baseidx),4);
                % EEG.CAT.Conn.S=EEG.CAT.Conn.S-repmat(baseVals,[1 1 1 length(latencies)]);
                maxVal4Plot=.0015;
            case 'EMGEMG'
                EEG.CAT.Conn.dDTF08=[];
                EEG.CAT.Conn.dDTF08=connStruct(9:16,9:16,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
                % EEG.CAT.Conn.dDTF08=EEG.CAT.Conn.dDTF08-repmat(baseVals,[1 1 1 length(latencies)]);
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(9:16,9:16,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.S(:,:,:,baseidx),4);
                % EEG.CAT.Conn.S=EEG.CAT.Conn.S-repmat(baseVals,[1 1 1 length(latencies)]);
                maxVal4Plot=.0025;
            case 'brainEMG'
                EEG.CAT.Conn.dDTF08=[];
                EEG.CAT.Conn.dDTF08=connStruct(9:16,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
                % EEG.CAT.Conn.dDTF08=EEG.CAT.Conn.dDTF08-repmat(baseVals,[1 1 1 length(latencies)]);
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(9:16,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.S(:,:,:,baseidx),4);
                % EEG.CAT.Conn.S=EEG.CAT.Conn.S-repmat(baseVals,[1 1 1 length(latencies)]);
                maxVal4Plot=.0005;
            case 'EMGbrain'
                EEG.CAT.Conn.dDTF08=[];
                EEG.CAT.Conn.dDTF08=connStruct(1:8,9:16,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
                % EEG.CAT.Conn.dDTF08=EEG.CAT.Conn.dDTF08-repmat(baseVals,[1 1 1 length(latencies)]);
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(1:8,9:16,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
                baseVals=median(EEG.CAT.Conn.S(:,:,:,baseidx),4);
                % EEG.CAT.Conn.S=EEG.CAT.Conn.S-repmat(baseVals,[1 1 1 length(latencies)]);
                maxVal4Plot=.0005;
        end
        % EEG.CAT.Conn.dDTF08=[];
        % EEG.CAT.Conn.dDTF08=connStruct(9:16,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
        % baseVals=median(EEG.CAT.Conn.dDTF08(:,:,:,baseidx),4);
        % % EEG.CAT.Conn.dDTF08=EEG.CAT.Conn.dDTF08-repmat(baseVals,[1 1 1 length(latencies)]);
        % EEG.CAT.Conn.S=[];
        % EEG.CAT.Conn.S=connStructS(9:16,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
        % baseVals=median(EEG.CAT.Conn.S(:,:,:,baseidx),4);
        % % EEG.CAT.Conn.S=EEG.CAT.Conn.S-repmat(baseVals,[1 1 1 length(latencies)]);
        EEG.CAT.curComps=1:size(EEG.CAT.Conn.dDTF08,1);
        EEG.CAT.nbchan=size(EEG.CAT.Conn.dDTF08,1);
        EEG.CAT.curComponentNames(1:8)={'B1','B2','B3','B4','B5','B6','B7','B8'};
        EEGtemp=EEG;
        EEGtemp.CAT.Conn=rmfield(EEG.CAT.Conn,{'Coh','iCoh','RPDC','GPDC','ffDTF','GGC','DTF','PDC','S'});
        EEGtemp = pop_vis_TimeFreqGrid(EEGtemp,'nogui');
        close(gcf);
        EEG.CAT.configs.vis_TimeFreqGrid=EEGtemp.CAT.configs.vis_TimeFreqGrid;
        if any(strcmpi(typeCond,{'brainBrain','EMGEMG'}))
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout=[];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_direct=0;
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.triu='dDTF08';
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.ut_clim=[-maxVal4Plot maxVal4Plot];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.tril='dDTF08';
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.lt_clim=[-maxVal4Plot maxVal4Plot];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.diag='S';
            maxValMeasS=3; %1.5;
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.d_clim=[-maxValMeasS maxValMeasS];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.clim=[];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_selection='Partial';
        else
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout=[];
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_direct=0;
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.estimator='dDTF08';
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.arg_selection='Full';
            EEG.CAT.configs.vis_TimeFreqGrid.MatrixLayout.clim=[-maxVal4Plot maxVal4Plot];
        end


        EEG.CAT.configs.vis_TimeFreqGrid.freqscale='log';%'linear';%'log'; % %Note: log scale doesn't seem to alter plotting result
        EEG.CAT.configs.vis_TimeFreqGrid.foilines=[8 13 30];%log([4 8 13 30]); %
        EEG.CAT.configs.vis_TimeFreqGrid.foilinecolor=[0 0 0]; %[1 0 0; 1 0 0; 0 0 1; 0 0 1];
        origFreqs=EEGtemp.CAT.configs.vis_TimeFreqGrid.freqValues;
        EEG.CAT.configs.vis_TimeFreqGrid.freqValues=origFreqs(origFreqs>=4);
        EEG.CAT.configs.vis_TimeFreqGrid.events{1}={0 'k' ':' 2};
        if strcmp(cond(1),'a')
            EEG.CAT.configs.vis_TimeFreqGrid.events{2}={0.5 'k' ':' 2};
        else
            EEG.CAT.configs.vis_TimeFreqGrid.events{2}={1 'k' ':' 2};
        end
        EEG.CAT.configs.vis_TimeFreqGrid.timeRange=[-0.5 1.5];
        EEG = pop_vis_TimeFreqGrid(EEG,'nogui');

        %Save the plot
        set(gcf,'Position',[0.034 0 0.966 0.934]); %[0.348 .309 .547 .588]);
        saveas(gcf,['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/plots/Connplot_' cond '_' typeCond '_' num2str(maxVal4Plot) '.fig']);
%         export_fig(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/plots/Connplot_' cond '_' typeCond '_' num2str(maxVal4Plot) '.png']);
        close all;
    end
end