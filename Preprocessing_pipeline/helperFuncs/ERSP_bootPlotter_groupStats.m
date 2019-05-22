function ERSP_bootPlotter_groupStats(STUDY,ALLEEG,tlim,flim,clim,clusters,cond_labels,savePlots,savePath,basetime,evPlotLines,erspNums)
%%Plots ERSPs using bootstrap significance masking
%%Inputs:
%       tlim - vector of [min max] times to plot (msec)
%       flim - vector of [min max] freqs to plot (Hz)
%       clim - vector of [min max] color scale (dB)
%       clusters - vector of cluster indices to make plots of
%       cond_labels - cell array of trial names plotted
%       savePlots - 1: save plots; 0: don't save plots
%       savePath - where to save plots if saving data
%       basetime - vector of [min max] times for baseline in msec (if nan,
%                  no signficance masking)
%       evPlotLines - vecotr of event times to mark events on ERSP plots

for j=clusters

    [STUDY, erspdata, alltimes, allfreqs] = std_readersp(STUDY, ALLEEG,'clusters',j);
    % Bootstrapping
    alpha = 0.05;%Set a significance value
    
    groupErspdata=cell(1,1);
    for i=1:length(erspdata)
        groupErspdata{1}=cat(3,groupErspdata{1},erspdata{i});
    end
    baseidxT = find(alltimes>=basetime(1) & alltimes<=basetime(end)); %Take times that are in your baseline
    baselineT = mean(groupErspdata{1}(:,baseidxT,:),2);
    %Find group significance mask across all conditions
    temp_ersp = groupErspdata{1}(:,:,:)-repmat(baselineT,1,length(alltimes));
    if ~isnan(alpha)
        pboot = bootstat(temp_ersp,'mean(arg1,3);','boottype','shuffle',...
            'label','ERSP','bootside','both','naccu',200,...
            'basevect',baseidxT,'alpha',alpha,'dimaccu',2);
    end
    
    sigMask=ones(size(erspdata{1}(:,:,1)));
    % Set basetime to NaN if you don't want to significant mask
    % Otherwise set basetime to the limits of your cycle, ex. stride
    baseidx = find(alltimes>=basetime(1) & alltimes<=basetime(end)); %Take times that are in your baseline
    baselines=[];
    for i=1:length(erspdata)
        baselines(:,i,:)=mean(erspdata{i}(:,baseidx,:),2);
    end
    baseline=mean(baselines(:,erspNums,:),2);
    for c = erspNums %1:length(erspdata)
        %Calculate and subtract baseline
%         baseline = mean(erspdata{c}(:,baseidx,:),2);
        curr_ersp = erspdata{c}(:,:,:)-repmat(baseline,1,length(alltimes));
        
        %Bootstrap and significance mask
        if ~isnan(alpha)
%             pboot = bootstat(curr_ersp,'mean(arg1,3);','boottype','shuffle',...
%                 'label','ERSP','bootside','both','naccu',200,...
%                 'basevect',baseidx,'alpha',alpha,'dimaccu',2);         
            curr_ersp = mean(curr_ersp,3);
            curr_maskedersp = curr_ersp;
            curr_maskedersp(curr_ersp > repmat(pboot(:,1),[1 size(curr_ersp,2)]) & curr_ersp < repmat(pboot(:,2),[1 size(curr_ersp,2)])) = 0;
        else
            curr_ersp = mean(curr_ersp,3);
            curr_maskedersp = curr_ersp;
        end
        maskedersp(c) = {curr_maskedersp};
        ersp(c) = {curr_ersp};
        sigMask=sigMask.*(curr_maskedersp~=0);
    end
    
    for c = erspNums
        maskedersp{c}=maskedersp{c}.*sigMask;
    end
    
    %%%Now plot the ersps
    figure('name',['Cls ' num2str(j)]);
    set(gcf,'Position',[1 1 1102 317]);
    % plots
    for k=erspNums %1:length(cond_labels)
%         if k==3
%             figure;
%         end
        h=subplot(1,length(erspNums),k); %2,2-mod(k,2)); %
        tftopo(maskedersp{k},alltimes,allfreqs,'vert',[evPlotLines],... 
            'limits',[tlim flim clim],'logfreq','native');
        title(cond_labels{k},'FontSize',16);
        ylimits = ylim; %(gca);
        set(gca,'YTick',log([4.01,8,13,30,50,99.4843]));%100-.05])); %4.009255933
        if k==1
            set(gca,'YTickLabel',{'4','8','13','30','50','100'},'Fontsize',12);
            ylabel(sprintf('Frequency\n(Hz)'),'Rotation',0,'fontsize',16,'fontweight','bold');
        else
            set(gca,'YTickLabel',{'','','','','','',''});
            ylabel('');
        end
%         set(gca,'XTick',[-500 0 500 990]);
%         set(gca,'XTickLabel',{'-500','0','500','1000'});
        if k~=2
            xlabel('');
        else
            xlabel('Time (msec)','Fontsize',16,'fontweight','bold');
        end
        ylhand = get(gca,'ylabel');
        %set(ylhand,'Rotation',0,'fontsize',16,'fontweight','bold');
        set(gca,'Fontsize',14,'fontweight','bold');
        set(gca,'box','off');
        set(gca,'XTick',[]);
%         subPos=get(h,'pos');
%         set(h,'pos',subPos-(k-1)*[.08 0 0 0]);
%         % Add lines for events
%         for e = 1:length(event_times)
%             plot(event_times(e)*[1 1],ylimits, 'k:', 'linewidth',2)
%         end
    end
    if savePlots==1
        saveas(gcf,[savePath 'clus' num2str(j) '_ersp.fig']);
        saveas(gcf,[savePath 'clus' num2str(j) '_ersp.jpg'])
    end
end

disp('ERSPs Plotted!');