function Spec_bootPlotter_designs(STUDY,ALLEEG,tlim,flim,clim,clusters,cond_labels,savePlots,savePath,basetime,evPlotLines,designs)
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
    if length(designs)>1
        for c=1:length(designs)
            [STUDY] = std_selectdesign(STUDY, ALLEEG, designs(c));
            [STUDY, specTemp, allfreqs] = std_readspec(STUDY, ALLEEG,'clusters',j);
            %         [STUDY, erspTemp, alltimes, allfreqs] = std_readersp(STUDY, ALLEEG,'clusters',j);
            specdata{c}=specTemp{1};
        end
    else
        [STUDY, specdata, allfreqs] = std_readspec(STUDY, ALLEEG,'clusters',j);
    end
        
    %%%Now plot the ersps
    figure('name',['Cls ' num2str(j)]);
    
    % plots
    hold all
    for k=1:length(cond_labels)
        inds=find(allfreqs>=4 & allfreqs<=100);
        plot(allfreqs(inds),specdata{k}(inds),'Linewidth',2);
%         h=subplot(1,length(cond_labels),k);
%         tftopo(maskedersp{k},alltimes,allfreqs,'vert',[evPlotLines],... 
%             'limits',[tlim flim clim],'logfreq','native');
        title('Power Spectrum','FontSize',16); %cond_labels{k},'FontSize',16);
%         ylimits = ylim; %(gca);
%         set(gca,'YTick',log([4.01,8,13,30,50,80,100-.05])); %4.009255933
%         if k==1
%             set(gca,'YTickLabel',{'4','8','13','30','50','80','100'},'Fontsize',12);
%         else
%             set(gca,'YTickLabel',{'','','','','','',''});
%             ylabel('');
%         end
%         set(gca,'XTick',[-500 0 500 990]);
%         set(gca,'XTickLabel',{'-500','0','500','1000'});
%         if k~=2
%             xlabel('');
%         else
%             xlabel('Time (msec)','Fontsize',14);
%         end
%         ylhand = get(gca,'ylabel');
%         set(ylhand,'fontsize',14);
%         subPos=get(h,'pos');
%         set(h,'pos',subPos-(k-1)*[.08 0 0 0]);
%         % Add lines for events
%         for e = 1:length(event_times)
%             plot(event_times(e)*[1 1],ylimits, 'k:', 'linewidth',2)
%         end
    end
    ylim([40 80]);
    xlim([4 100]);
    hold off;
    legend(cond_labels{1},cond_labels{2},cond_labels{3});
    if savePlots==1
        saveas(gcf,[savePath 'clus' num2str(j) '_spec.fig']);
        saveas(gcf,[savePath 'clus' num2str(j) '_spec.jpg'])
    end
end

disp('ERSPs Plotted!');