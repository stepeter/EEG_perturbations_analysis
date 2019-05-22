function [EEGout,rankEEG]=runCCA_multConds(EEGin,segmentBounds)
%Breaks EEG into segments to run CCA on individual conditions.
%Inputs:      EEGin - input EEG set
%             segmentBounds - can specify start and end of segments in
%             numerical array [startSeg1 endSeg1; startSeg2 endSeg2; ...]
%             or can use argument 'boundary' to break up by condition if
%             initially merged sets (will use boundary events)


dataRank=rank(double(EEGin.data'));
if dataRank<size(EEGin.data,1)
    error('Data is not full rank! Remove interpolated channels and run again.');
end

if ischar(segmentBounds) && strcmp(segmentBounds,'boundary')
    %Get boundary events
    boundary_latencies=[];
    for i=1:length(EEGin.event)
        if strcmp(EEGin.event(i).type,'boundary')
            boundary_latencies=[boundary_latencies EEGin.event(i).latency];
        end
    end
    %Special condition for subject data
    if strcmp(EEGin.subject,'WMISM_1')
        boundary_latencies(end-1)=[];
    end
    boundary_latencies2=[1 boundary_latencies EEGin.pnts];
else
    boundary_latencies2=segmentBounds;
end

EEG_orig=EEGin;
%Loop through segments
for i=1:(length(boundary_latencies2)-1)
    %Split into segment
    EEG=pop_select(EEG_orig,'point',boundary_latencies2([i i+1]));
    %Run CCA
    disp(['Segment # ' num2str(i)]);
%     %Remove motion artifact 1st
%     disp('Remove motion artifact');
    [EEG,Wcca,rej_ccaComps]=CCA_cleanChans_smp(EEG,i);
    %Substitute the new data into EEG_orig
    intBounds=int64(boundary_latencies2([i i+1]));
    EEG_orig.data(:,intBounds(1,1):intBounds(1,2))=EEG.data;
    EEG_orig.Wcca{i}=Wcca;
    EEG_orig.rej_ccaComps{i}=rej_ccaComps;
%     
%     %Remove muscular/ocular artifact
%     disp('Remove muscular/ocular artifact');
%     [EEG,Wcca,rej_ccaComps]=CCA_cleanChans_smp(EEG);
%     %Substitute the new data into EEG_orig
%     intBounds=int64(boundary_latencies2([i i+1]));
%     EEG_orig.data(:,intBounds(1,1):intBounds(1,2))=EEG.data;
%     EEG_orig.Wcca2{i}=Wcca;
%     EEG_orig.rej_ccaComps2{i}=rej_ccaComps;
    close all;
end
EEGout=EEG_orig;

%Calculate final rank of data
rankEEG=rank(double(EEGout.data'));
end