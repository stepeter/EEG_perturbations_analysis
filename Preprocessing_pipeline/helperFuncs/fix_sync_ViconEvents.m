%% fixing the vicon and biosemi sync issues  %%%%%%%%%
%  J. Lukos 05/22/2013, Army Research Laboratory %%%%%
%  Modified by S. Peterson 02/16/2017

function [new_eeg_trig_latencies,new_csv_trig_latencies,EEG,new_biomech] = fix_sync_ViconEvents(EEG, csv_data, csv_sync_ch_idx, csv_sample_rate, syncFreq)


%Inputs:
%       - EEG: EEG struct
%       - csv_data: sync signal (Matlab variable with sync signal for each condition)
%       - group: name of group using
%       - csv_sample_rate: srate for Vicon (should be 1000 Hz)
%       - syncFreq: frequency of sync signal

eeg_trig_latencies = []; eeg_trig_event= [];

%Find rising sync latencies
for i=1:length({EEG.event(1:end).type})
    syncEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['Sync Rising' '*']),'once'); %[trigPrefix 'CCW_*']),'once');
    if ~isempty(syncEv)
        eeg_trig_latencies=[eeg_trig_latencies EEG.event(1,i).latency];
        eeg_trig_event=[eeg_trig_event EEG.event(1,i).urevent];
    end
end

if eeg_trig_latencies(1) == 1
   eeg_trig_latencies(1) = []; %this means the sync was on when the recording started. We don't want this, we want rising edges
end
%offset sync so it is zero based
% csv_data{2,csv_sync_ch_idx} = csv_data{2,csv_sync_ch_idx}-min(csv_data{2,csv_sync_ch_idx});

%define csv trigger latencies
csv_trig_latencies = [];

%if sync channel returned to zero for one frame consider this to be
%noise and use only the first onset event
% csv_trig_latencies(find(diff(csv_trig_latencies) <= 2)+1) = [];
for i = 1:length(csv_data{2,csv_sync_ch_idx})-1
    if csv_data{2,csv_sync_ch_idx}(i) < 2 && csv_data{2,csv_sync_ch_idx}(i+1) > 2
        csv_trig_latencies(length(csv_trig_latencies)+1) = i+1;
    end
end



%% fixing sync
% csv file

new_biomech=csv_data{2,csv_sync_ch_idx}(1:csv_trig_latencies(1)-1)';
for loop=1:length(csv_trig_latencies)-1
    if csv_trig_latencies(loop+1)-csv_trig_latencies(loop)~=csv_sample_rate*syncFreq
        offset=(csv_trig_latencies(loop+1)-csv_trig_latencies(loop)-1)/(csv_sample_rate*syncFreq-1);
        tmp_time=csv_trig_latencies(loop):offset:csv_trig_latencies(loop+1)-1;
        new_data=[];
        %             for j=1:size(csv_data{2,k},1)
        new_data=interp1([csv_trig_latencies(loop):1:csv_trig_latencies(loop+1)-1],csv_data{2,k}(csv_trig_latencies(loop):csv_trig_latencies(loop+1)-1),tmp_time,'spline');
        %             end
        new_biomech=horzcat(new_biomech,new_data);
    else
        old_data=[csv_data{2,csv_sync_ch_idx}(csv_trig_latencies(loop):csv_trig_latencies(loop+1)-1)];
        new_biomech=horzcat(new_biomech,old_data);
    end
end
new_biomech=horzcat(new_biomech,csv_data{2,csv_sync_ch_idx}(csv_trig_latencies(end):end)');

new_csv_trig_latencies = [];
for i = 1:size(new_biomech,2)-1
    if new_biomech(i) < 2 && new_biomech(i+1) > 2
        new_csv_trig_latencies(length(new_csv_trig_latencies)+1) = i+1;
    end
end

% eeg file
new_eeg=EEG.data(:,1:eeg_trig_latencies(1)-1);
for loop=1:length(eeg_trig_latencies)-1
    if eeg_trig_latencies(loop+1)-eeg_trig_latencies(loop)~=EEG.srate*syncFreq
        offset=(eeg_trig_latencies(loop+1)-eeg_trig_latencies(loop)-1)/(EEG.srate*syncFreq-1);
        tmp_time=eeg_trig_latencies(loop):offset:eeg_trig_latencies(loop+1)-1;
        for j=1:size(EEG.data,1)
            eeg_data(j,:)=interp1([eeg_trig_latencies(loop):1:eeg_trig_latencies(loop+1)-1],EEG.data(j,eeg_trig_latencies(loop):eeg_trig_latencies(loop+1)-1),tmp_time);
        end
        new_eeg=horzcat(new_eeg,eeg_data);
        EEG.event(eeg_trig_event(loop+1)).latency=length(new_eeg);
    else
        old_eeg=[EEG.data(:,eeg_trig_latencies(loop):eeg_trig_latencies(loop+1)-1)];     
        new_eeg=horzcat(new_eeg,old_eeg);
        EEG.event(eeg_trig_event(loop+1)).latency=length(new_eeg);
    end
end
new_eeg=horzcat(new_eeg,EEG.data(:,eeg_trig_latencies(end):end));
EEG.data=new_eeg;
EEG.pnts=length(EEG.data);

new_eeg_trig_latencies = [];
%Find rising sync latencies
for i=1:length({EEG.event(1:end).type})
    syncEv=regexp(EEG.event(1,i).type, regexptranslate('wildcard',['Sync Rising' '*']),'once'); %[trigPrefix 'CCW_*']),'once');
    if ~isempty(syncEv)
        new_eeg_trig_latencies=[new_eeg_trig_latencies EEG.event(1,i).latency];
    end
end

%remove trailing sync events so eeg and csv have the same number
if length(eeg_trig_latencies) > length(csv_trig_latencies)
    eeg_trig_latencies = eeg_trig_latencies(1:length(csv_trig_latencies));
    new_eeg_trig_latencies = new_eeg_trig_latencies(1:length(new_csv_trig_latencies));
else
    csv_trig_latencies = csv_trig_latencies(1:length(eeg_trig_latencies));
    new_csv_trig_latencies = new_csv_trig_latencies(1:length(new_eeg_trig_latencies));
end

%%
   
%compare csv and EEG trig latencies
old_eeg_trig_time = eeg_trig_latencies/EEG.srate;
new_eeg_trig_time = new_eeg_trig_latencies/EEG.srate;
old_csv_trig_time = csv_trig_latencies/csv_sample_rate;
new_csv_trig_time = new_csv_trig_latencies/csv_sample_rate;
time_shift=[];
length(new_csv_trig_time)
length(new_eeg_trig_time)
time_shift(1,:) = old_csv_trig_time - old_eeg_trig_time;
time_shift(2,:) = new_csv_trig_time - new_eeg_trig_time;
y=figure; plot(diff(old_csv_trig_time)); hold on;  plot(diff(old_eeg_trig_time),'r'); plot(diff(new_csv_trig_time),'g.'); plot(diff(new_eeg_trig_time),'r.');
title('sync data','Interpreter','none'); legend('biomech old','EEG','biomech new'); 
% saveas(y,['plots/' csv_filename(firstletter:lastletter) '_synctime']);
x=figure; plot(old_csv_trig_time,time_shift(1,:),'ro'); hold on; plot(new_csv_trig_time,time_shift(2,:),'bo');
title('time shift','Interpreter','none'); legend('old','new'); 
% saveas(x,['plots/' csv_filename(firstletter:lastletter) '_timeshift']);
close(y); close(x);
