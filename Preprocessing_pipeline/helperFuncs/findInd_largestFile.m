function [maxBytesInd]=findInd_largestFile(rawEEG_inpath)

%Find .set file with most memory space (will be training trial)
listSets = dir([rawEEG_inpath filesep '*.xdf']);
maxBytesInd=-1; maxBytesVal=-1;
for i=1:length(listSets)
    if listSets(i,1).bytes>maxBytesVal
        maxBytesVal=listSets(i,1).bytes;
        maxBytesInd=i;
    end
end