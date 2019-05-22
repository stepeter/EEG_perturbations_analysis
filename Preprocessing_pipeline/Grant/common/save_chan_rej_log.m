function save_chan_rej_log(savefilename, chan_rej_log)

fid = fopen([savefilename '.txt'], 'w');

[r, c] = size(chan_rej_log);

for ic = 1:c
    fprintf(fid, '%s\t%s\n', 'Method: ', chan_rej_log{1,ic});
    fprintf(fid, '%s\t%d\n', 'Num channels rejected: ', chan_rej_log{2,ic});
    fprintf(fid, '%s\t', 'Rejected channels: '); 
    nchans = chan_rej_log{2,ic};
    for k = 1:(nchans-1)
        fprintf(fid, '%s\t', char(chan_rej_log{3,ic}(k)));
    end
    fprintf(fid, '%s\n', char(chan_rej_log{3,ic}(k+1)));
    fprintf(fid, '%s\t%d\n', 'Total number of channels rejected: ', chan_rej_log{4,ic});
    fprintf(fid, '%s\t%d\n\n', 'Remaining channels: ', chan_rej_log{5,ic});    
end

fclose(fid);

% save chan_rej_log variable to .mat
save([savefilename '.mat'], 'chan_rej_log');