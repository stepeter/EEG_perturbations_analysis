%Cut up EEG data based on sync signal
function [EEG]=sync_cut_v2(EEG,badEpchs)
   if nargin==1
        badEpchs=-1;
   end
   type={EEG.event(1:end).type};
   isRightType=zeros(1,length(type));
   for j=1:length(type)
       if strcmp(type{1,j},'6') || strcmp(type{1,j},'32774') || floor(mean(type{1,j}==6)) || floor(mean((type{1,j}==32774)))
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
       exclude(i,1)=lats(btwnEpochs(i))+2/EEG.srate;
       exclude(i,2)=lats(btwnEpochs(i)+1)-2/EEG.srate;
   end
   
   %Create array of starting events
   startTimes=[lats(1) lats(btwnEpochs(1:(end))+1)]; %(end-1))+1)]; %changed 4/17/16
   
   %Add in starting of pass events
   for index = 1 : length(startTimes)
       EEG.event(end+1) = EEG.event(1);    % Add new event to end of event list
       EEG.event(end).latency = startTimes(index)*EEG.srate;    % latency of event (one time pt less than event bc of conversion in line 15)
       EEG.event(end).type = 'startPass';    % Name of new events
   end

   EEG = eeg_checkset(EEG, 'eventconsistency');     % Check all events for consistency
   EEG = pop_editeventvals( EEG, 'sort', {'latency', [0] } );     % Re-sort events
   
   remEpchs=[];
   if badEpchs~=-1
        if max(badEpchs>length(exclude)+1)
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
    if exclude(size(exclude,1),2)==0
        exclude=[0 lats(1)-2/EEG.srate; exclude];
        exclude(size(exclude,1),2)=EEG.times(end)/1000;
    else
        exclude=[0 lats(1)-2/EEG.srate; exclude; lats(end) EEG.times(end)/1000];
    end
   
   disp(['Excluding: ' num2str(sum(exclude(:,2)-exclude(:,1))/60) ' minutes'])
   %Remove EEG data during those times
   EEG=pop_select(EEG,'notime',exclude);
end
