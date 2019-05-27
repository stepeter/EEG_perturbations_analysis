%Permuting subset of subjects to test the min, average, and max samples per connection.
%This helps estimate how increases in subject sample size affect the sample size at each connection.

subjectInds=[1:7 9:13 15:17 19:20 22:33];
numConns_boot=cell(16,16);
cond='allSVZ'; group = 'all';
RootData='/usr/local/VR_connectivity/Data/';
goodClusts=3:10;
load('/usr/local/VR_connectivity/Data/STUDYtopo.mat');
load('/usr/local/VR_connectivity/Data/cluster_1222018.mat'); %_pruned.mat');
clustDat=cluster; allSubjs=[1:13 15:17 19:20 22:33];
for i=subjectInds
    subjcode = ['WMISM_' num2str(i)];
    subj_basepath = [RootData subjcode filesep];
    EEGsets_outpath = ['/usr/local/VR_connectivity/Data/regularGroupConn/' cond]; %allSVZ'; %[subj_basepath 'EEG_sets'];
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
    for j=icNums
        for k=icNums
            numConns_boot{j,k}=[numConns_boot{j,k} i];
        end
    end
end

%% From a range of 1-29 subjects, do 200 subject permutations and average the result
numPerms=200; 
finalPermMean=zeros(1,29); finalPermStd=finalPermMean;
finalPermMean_min=finalPermMean; finalPermStd_min=finalPermMean;
finalPermMean_max=finalPermMean; finalPermStd_max=finalPermMean;
for numSubjs=1:29
    permVals=zeros(1,numPerms);
    permVals_min=zeros(1,numPerms);
    permVals_max=zeros(1,numPerms);
    disp(numSubjs)
    for i=1:numPerms
        tmpInds=randperm(length(subjectInds),numSubjs);
        newSubjInds=subjectInds(tmpInds);
        
        %Check the number of connections for each off-diagonal and average
        %results
        numConns=[];
        for j=1:8
            for k=1:8
                if j~=k
                    tmpSum=0;
                    for p=newSubjInds
                        tmpSum=tmpSum+sum(p==numConns_boot{j,k});
                    end
                    numConns=[numConns tmpSum];
                end
            end
        end
        permVals(1,i)=mean(numConns);
        permVals_min(1,i)=min(numConns);
        permVals_max(1,i)=max(numConns);
    end
    finalPermMean(1,numSubjs)=mean(permVals); finalPermStd(1,numSubjs)=std(permVals);
    finalPermMean_min(1,numSubjs)=mean(permVals_min); finalPermStd_min(1,numSubjs)=std(permVals_min);
    finalPermMean_max(1,numSubjs)=mean(permVals_max); finalPermStd_max(1,numSubjs)=std(permVals_max);
end

%%
figure; hold on;
errorbar(1:29,finalPermMean,finalPermStd,'g');
errorbar(1:29,finalPermMean_min,finalPermStd_min,'r');
errorbar(1:29,finalPermMean_max,finalPermStd_max,'b');
hold off;
xlim([0 30]); ylim([-1 20]);

Fit_mean = polyfit(1:29,finalPermMean,2);
Fit_min = polyfit(1:29,finalPermMean_min,2);
Fit_max = polyfit(1:29,finalPermMean_max,2);
finalPerm.meanVal=finalPermMean;
finalPerm.maxVal=finalPermMean_max;
finalPerm.minVal=finalPermMean_min;
finalPerm.meanValstd=finalPermStd;
finalPerm.maxValstd=finalPermStd_max;
finalPerm.minValstd=finalPermStd_min;
%% Simulation of IC number vs. number of connections (should show quadratic behavior)
numConns=zeros(1,29);
for i=1:29
    for j=1:i
        for k=1:i
            if j~=k
                numConns(i)=numConns(i)+1;
            end
        end
    end
end
figure; plot(1:29,numConns);
Fit_ICnum = polyfit(1:29,numConns,2);
finalPerm.simFitICnum=numConns;
finalPerm.x=1:29;

save('/usr/local/VR_connectivity/Data/regularGroupConn/plots/numSubjsNumConns_fit.mat','finalPerm');