function []=saveCUDAICAparamfile(ICA_path,pGran,filename,num_chans,num_frames,group)
fid = fopen([ICA_path filesep 'DataMerge_' group '_cudaica.sc'],'w');
naFile = [pGran filename];
fprintf(fid,['DataFile ' naFile '\n']);
fprintf(fid,['chans %d\n'],num_chans);
fprintf(fid,['frames %d\n'],num_frames);
fprintf(fid,'epochs 1\n');
WtsOut = [pGran 'files_CUDAICA' group '/DataMerge_' group '_cudaica.wts'];
SphOut = [pGran 'files_CUDAICA' group '/DataMerge_' group '_cudaica.sph'];
fprintf(fid,['WeightsOutFile ' [WtsOut '\n']]);
fprintf(fid,['SphereFile ' [SphOut '\n']]);
fprintf(fid,'sphering on\nbias on\nextended 0\npca 0\nlrate 1.0e-4\nblocksize 0\nstop 1.0e-7\nmaxsteps 1500\nposact off\nannealstep 0.98\nannealdeg 60\nmomentum 0\nverbose on');
fclose(fid);
