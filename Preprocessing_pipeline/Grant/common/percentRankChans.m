function [percentRank]=percentRankChans(EEG,weights)
    %Takes in EEG channels and ranks them for channel rejection measures
    %depending on weights used.
    %Weights is 4 column vector with weighting for each measure:
    	%correlation, standard deviation, kurtosis, range
    
    scores=zeros(length(EEG.chanlocs),4); %nCol=1;
    
    %Rank correlations (higher is worse)
    disp('Computing correlations...');
    [badchanindts,badchanlogts,c_cellts,IND,maxcts]=detectnoisychannelts(EEG);
    percentSecondsBad = sum(badchanlogts') / size(badchanlogts,2);
    percentSecondsBad = abs(percentSecondsBad); %make sure all positive
    corRank=tiedrank(percentSecondsBad)/length(percentSecondsBad);
    scores(:,1)=corRank;
    
    %Rank standard deviation (higher is worse)
    disp('Computing standard deviation...');
    try
        stdev = std(EEG.data,0,2);
    catch
        %catch memory error and compute std one channel at a time
        stdev = [];
        for i = 1:EEG.nbchan
            disp(['Computing standard deviation for channel ' num2str(i) '...']);
            stdev(i) = std(EEG.data(i,:));
        end
    end
    stdev=abs(stdev); %make sure all positive
    stdRank=tiedrank(stdev)/length(stdev);
    scores(:,2)=stdRank;
    
    %Rank Kurtosis (higher is worse)
    %disp('Computing kurtosis...');
    [REJ_EEG rej_chan_idx measure] = pop_rejchan(EEG,'elec',...
        [1:EEG.nbchan],'threshold',5,'norm',...
        'on','measure','kurt');
    kurtRank=tiedrank(measure)/length(measure);
    scores(:,3)=kurtRank;
    
    %Rank Range (higher is worse)
    disp('Computing range...');
    ranges=range(EEG.data,2);
    rangeRank=tiedrank(ranges)/length(ranges);
    scores(:,4)=rangeRank;
    
    %Weighted average across percent ranks for each channel
    %(If a weighting is zero, column is removed before averaging)
    percentRank=mean(scores(:,weights~=0).*(ones(size(scores,1),1)*weights(weights~=0)),2);
end