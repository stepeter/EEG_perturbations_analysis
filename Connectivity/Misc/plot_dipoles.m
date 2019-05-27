%Plot dipoles and centroids
close all;
STUDY.cluster = clustDat_new;
colors2use={[1 0 0], [0 1 0], [0 0 1], [1 0 1], [0 1 1], [.8549 .6471 .1255], [1 .4902 0], ... %[1 1 0], [1 .4902 0], ...
    [1 .0784 .5765], [.8 0 1], [.6 0 0], [0 0 0]};
centrDipsTogether=1; %1 - plot centroids and dipoles together; 0 - only dipoles
clusters_to_plot=10:-1:3; %3:10;
savePath='/usr/local/VR_connectivity/Data/neuroimage_revisions2/dipole_plots/';

diplotfig(STUDY,ALLEEG,clusters_to_plot, colors2use(1:length(clusters_to_plot)),savePath,centrDipsTogether);