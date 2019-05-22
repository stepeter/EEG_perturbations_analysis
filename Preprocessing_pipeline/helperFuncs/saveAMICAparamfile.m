function []=saveAMICAparamfile(ICA_path,pFlux,filename,num_chans,num_frames,doPCA,numPCs,maxthreads,group)
%Creates param file for AMICA given inputs
%doPCA: 1 - do PCA down to numPCs before ICA; 0 - no PCA
fid = fopen([ICA_path filesep group '_amica.param'],'w');
fprintf(fid,['files ' pFlux filename '\n']);
fprintf(fid,['outdir ' pFlux 'files_AMICA' group '/\n']);
fprintf(fid,'num_models 1\n');
fprintf(fid,'num_mix_comps 3\n');
fprintf(fid,'pdftype 0\n');
fprintf(fid,'block_size 128\n');
fprintf(fid,'max_iter 1500\n');
fprintf(fid,'num_samples 1\n');
fprintf(fid,'data_dim %d\n',num_chans);
fprintf(fid,'field_dim %d\n',num_frames);
fprintf(fid,'field_blocksize 1\n');
fprintf(fid,'share_comps 0\n');
fprintf(fid,'share_start 100\n');
fprintf(fid,'comp_thresh 0.990000\n');
fprintf(fid,'share_iter 100\n');
fprintf(fid,'lrate 0.100000\n');
fprintf(fid,'minlrate 1.000000e-08\n');
fprintf(fid,'lratefact 0.500000\n');
fprintf(fid,'rholrate 0.050000\n');
fprintf(fid,'rho0 1.500000\n');
fprintf(fid,'minrho 1.000000\n');
fprintf(fid,'maxrho 2.000000\n');
fprintf(fid,'rholratefact 0.500000\n');
fprintf(fid,'kurt_start 3\n');
fprintf(fid,'num_kurt 5\n');
fprintf(fid,'kurt_int 1\n');
fprintf(fid,'do_newton 1\n');
fprintf(fid,'newt_start 50\n');
fprintf(fid,'newt_ramp 10\n');
fprintf(fid,'newtrate 1.000000\n');
fprintf(fid,'do_reject 1\n');
fprintf(fid,'numrej 15\n'); %number of iterations to perform rejection during ICA
fprintf(fid,'rejsig 3.000000\n');
fprintf(fid,'rejstart 1\n');
fprintf(fid,'rejint 1\n');
fprintf(fid,'max_threads %d\n',maxthreads);
fprintf(fid,'writestep 10\n');
fprintf(fid,'write_nd 1\n');
fprintf(fid,'write_LLt 1\n');
fprintf(fid,'decwindow 1\n');
fprintf(fid,'max_decs 3\n');
fprintf(fid,'update_A 1\nupdate_c 1\nupdate_gm 1\nupdate_alpha 1\nupdate_mu 1\nupdate_beta 1\n');
fprintf(fid,'invsigmax 100.000000\n');
fprintf(fid,'invsigmin 0.000000\n');
fprintf(fid,'do_rho 1\n');
fprintf(fid,'load_rej 0\nload_W 0\nload_c 0\nload_gm 0\nload_alpha 0\nload_mu 0\nload_beta 0\nload_rho 0\nload_comp_list 0\n');
fprintf(fid,'do_mean 1\ndo_sphere 1\n');
fprintf(fid,'doPCA %d\n',doPCA);
fprintf(fid,'pcakeep %d\n',numPCs);
fprintf(fid,'pcadb 30.000000\n');
fprintf(fid,'byte_size 4\ndoscaling 1\nscalestep 1\n');

fclose(fid);
