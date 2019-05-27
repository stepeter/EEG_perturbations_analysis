condS={'allSVZ','allWVZ','pullsStn','pullsWalk'};
typeCondS={'brainBrain','EMGEMG'};
tlim = [-0.5 1.5];
flim = [4 50];
clim = [-3 3];
savePlots = 1;
brainLabels = {'LO','RO','LS','ACC','RS','PP','SMA','AP'};
emgLabels = {'LTA','LSO','LMG','LPL','RTA','RSO','RMG','RPL'};
savePath = '/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/erspPlots/';
for condInds=1:4
    for typeCondInds=1:2
        cond=condS{condInds}; %'pullsStn'; %close all;
        typeCond=typeCondS{typeCondInds}; %'brainBrain'; %'brainBrain','EMGEMG','brainEMG','EMGbrain' (from -> to)

        load(['/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/' cond '/median/connStructS_baseSub.mat']);

        EEG=pop_loadset('/usr/local/VR_connectivity/Data/regularGroupConn_hiModelOrder/allSVZ/S4.set');
        switch typeCond
            case 'brainBrain'
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(1:8,1:8,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
            case 'EMGEMG'
                EEG.CAT.Conn.S=[];
                EEG.CAT.Conn.S=connStructS(9:16,9:16,:,:); %(9:16,9:16,:,:); %(1:8,1:8,:,:); %
        end
        
        %Plot for each cluster
        for frodo = 1:8
            if strcmp(typeCond,'brainBrain')
                labelName = brainLabels{frodo};
            else
                labelName = emgLabels{frodo};
            end
            figure('name',[labelName]);
            set(gcf,'Position',[284 222 412 324]);
            maskedersp = squeeze(EEG.CAT.Conn.S(frodo,frodo,:,:));
            alltimes = EEG.CAT.Conn.erWinCenterTimes;
            allfreqs = EEG.CAT.Conn.freqs;
            if any(strcmp(cond,{'allSVZ','allWVZ'}))
                evPlotLines = [0 0.5];
            else
                evPlotLines = [0 1];
            end
            tftopo(maskedersp,alltimes,allfreqs,'vert',[],...
                'limits',[tlim flim clim],'logfreq','native');
            title('','FontSize',16);
            hold on;
            plot3([evPlotLines(1) evPlotLines(1)],log([4 50]),[4 4],'--k','linewidth',3);
            plot3([evPlotLines(2) evPlotLines(2)],log([4 50]),[4 4],'--k','linewidth',3);
            
            set(gca,'YTick',log([4.01,8,13,30,50]));
            set(gca,'YTickLabel',{'4','8','13','30','50'},'Fontsize',12);
            set(gca,'Fontsize',14,'fontweight','bold');
            set(gca,'box','off','linewidth',3);
            ylabel(''); xlabel('');
            
            %Save the plot
            if savePlots==1
                saveas(gcf,[savePath cond '_' labelName '.fig']);
                export_fig([savePath cond '_' labelName '.png']);
%                 saveas(gcf,[savePath cond '_' labelName '.jpg']);
            end
            close all;
        end
    end
end