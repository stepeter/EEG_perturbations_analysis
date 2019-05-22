%Cut up EEG data based on sync signal
function [EEG]=sync_cut_v2(EEG,badEpchs)
   if nargin==1
        badEpchs=-1;
   end
   type={EEG.event(1:end).type};
   isRightType=zeros(1,length(type));
   for j=1:length(type)
       if (type{1,j}==6) || (type{1,j}==32774)
           isRightType(1,j)=1;
       end
   end
   
   %Pick out correct latencies
   evs=(cell2mat({EEG.event(1:end).latency})-1)/EEG.srate;
   lats=evs(isRightType==1);
   dLats=diff(lats);
   btwnEpochs=find(dLats>1.5);
   
   %Create array of times to exclude
   exclude=zeros(length(btwnEpochs),2);
   for i=1:length(btwnEpochs)
       exclude(i,1)=lats(btwnEpochs(i));
       exclude(i,2)=lats(btwnEpochs(i)+1);
   end
   
   remEpchs=[];
    if badEpchs~=-1
        if max(badEpchs>length(exclude))
            error('Bad epoch value is out of range!');
        end
        for i=1:length(badEpchs)
            if badEpchs(i)==1
                exclude(1,1)=0;
            else
                exclude(badEpchs(i),1)=exclude(badEpchs(i)-1,1);
                remEpchs=[remEpchs badEpchs(i)-1];
            end
        end
        exclude(remEpchs,:)=[];
    end
    
    %Edge of first and last epochs
    exclude=[0 lats(1); exclude; lats(end) EEG.times(end)/1000];
   
   disp(['Excluding: ' num2str(sum(exclude(:,2)-exclude(:,1))/60) ' minutes'])
   %Remove EEG data during those times
   EEG=pop_select(EEG,'notime',exclude);
end
