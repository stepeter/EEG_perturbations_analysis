%Run CCA on data
delay=1; %can set different delay for CCA if want
[W,r] = bsscca(EEG.data,delay);

%Take CCA sources and put them in EEG.icaact
EEG.icaweights=real(W);
EEG.icasphere=eye(EEG.nbchan);
EEG=eeg_checkset(EEG,'ica');

%Use channel rejection framework to remove bad CCA components
tmpInds=1:EEG.nbchan;
clear tmp_comp_struct;
for i = 1:length(tmpInds)
    tmp_comp_struct(i).labels = num2str(tmpInds(i));
end
for i = 1:length(EEG.chanlocs)
    str{i} = EEG.chanlocs(i).labels;
end
mean_ica_act = EEG.icaact(tmpInds,:)-repmat(mean(EEG.icaact(tmpInds,:),2),1,EEG.pnts);

EEG2=EEG;
EEG2.data=EEG2.icaact; rej_chan_idxAll=[]; done=0;
while ~done
    rej_chan_temp = eegchan_listdlg(EEG2,'PromptString','Select Chans to Reject',...
        'SelectionMode','multiple','ListString',str);
    rej_chan_idx=unique([rej_chan_idxAll rej_chan_temp]);
    data2 = nan(size(EEG2.data));
    data2(rej_chan_idx,:)=mean_ica_act(rej_chan_idx,:);
    eegplot(mean_ica_act,'submean','off','srate',EEG.srate,...
    'eloc_file',tmp_comp_struct,'events',EEG.event,'winlength',50,'data2',data2);
    h2 = gcf; waitfor(h2)
    rej_chan_idxAll=rej_chan_idx;
    prompt = sprintf('%s\n%s','Reject components?',...
        [num2str(EEG.nbchan-length(rej_chan_idx)) ' components would remain']);
    answer = questdlg(prompt,...
        'Reject CCA components?','YES (again)','YES (done)','CLEAR','YES (again)');
    if strcmp(answer, '')
        answer = 'YES (again)';
    end
    switch answer
        case 'YES (done)'
            done=1;
        case 'CLEAR'
            rej_chan_idxAll=[];
    end
            
end

%Remove bad components and subtract from data
EEG=pop_subcomp(EEG,rej_chan_idxAll);