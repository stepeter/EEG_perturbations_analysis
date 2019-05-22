function [EEGout,Wcca,rej_ccaComps]=CCA_cleanChans_smp(EEG,num)
%%Takes in EEG data, runs Canonical correlation analysis, and then allows
%%user to reject CCA components based on raw trace
done_CCA=0;
while ~done_CCA
    %Run CCA on data
    delay=1; %can set different delay for CCA if want
    
    %Added in for removing really bad frames just for CCA
%     EEG_temp=EEG;
%     if num==8
%         EEG_temp=pop_select(EEG,'notime',[569 574]);
%     end
%     switch num
%         case 2
%             EEG_temp=pop_select(EEG,'notime',[43 47; 51 55; 76 81; 604 608]);
%         case 4
%             EEG_temp=pop_select(EEG,'notime',[487 492]);
%         case 6
%             EEG_temp=pop_select(EEG,'notime',[389 393; 540 545]);
%         case 7
%             EEG_temp=pop_select(EEG,'notime',[268 272]);
%         case 8
%             EEG_temp=pop_select(EEG,'notime',[216 220]);
%     end
%     if num==2
%         EEG_temp=pop_select(EEG,'notime',[70 77; 78 81; 87 90; 96 98; 126 129; 134 137; 142 145; 159 162; 180 183; 191 192; 230 233; 242 245; 265 268;...
%             269 273; 285 288; 317 320; 322 325; 355 358; 374 377; 387 388; 413 417; 427 430; 436 438; 452 454; 458 461; 466 470; 491 493; 500 502;...
%             512 514; 527 529; 541 543; 558 559; 562 564; 570 573; 590 592; 602 604; 626 630]);
%     end
%     [W,r] = bsscca_local(EEG_temp.data,delay);
    
    [W,r] = bsscca_local(EEG.data,delay);
        
    %Take CCA sources and put them in EEG.icaact (clear first)
    EEG.icaweights=[]; EEG.icasphere=[]; EEG.icawinv=[]; EEG.icaact=[];
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
%         figure;
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

    answer = questdlg('Rerun CCA after removed components?',...
            'Rerun CCA?','YES','NO','NO');
    if strcmp(answer,'NO')
        done_CCA=1;
    end
end
EEGout=EEG;
Wcca=W; rej_ccaComps=rej_chan_idxAll;
end

function [W,r] = bsscca_local(X,delay)
% bsscca() - Blind Source Separation through Canonical Correlation Analysis
% From AAR toolbox (EEGLAB extension)
%
% Usage:
%   >> [W,r] = bsscca(X,delay)
%
% Inputs:
%   X     - data matrix (dxN, data channels are rowwise)
%   delay - delay at which the autocorrelation of the sources will be
%           maximized (def: 1)
%
% Output:
%   W     - separation matrix
%   r     - autocorrelation of the estimated sources at the given delay
%
% See also:
%   BSSCCA_IFC, AUTOBSS

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 2, delay = 1; end
if nargin < 1, 
    help cca;
    return;
end


[d,T] = size(X);

% correlation matrices
Y = X(:,delay+1:end);
X = X(:,1:end-delay);
Cyy = (1/T)*Y*Y';
Cxx = (1/T)*X*X';
Cxy = (1/T)*X*Y';
Cyx = Cxy';
invCyy = pinv(Cyy);

% calculate W
[W,r] = eig(pinv(Cxx)*Cxy*invCyy*Cyx);
r = sqrt(abs(real(r)));
[r,I] = sort(diag(r),'descend');
W = W(:,I)';
end