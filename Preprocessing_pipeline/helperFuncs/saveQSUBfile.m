function []=saveQSUBfile(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group)
%Creates qsub file to run AMICA on Flux

amicaBinaryPath='/nfs/hnldata/stepeter/eeglab13_5_4b/plugins/AMICA_15/amica15ub';

fid = fopen([ICA_path filesep group 'qsub.pbs'],'w');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'####  PBS preamble\n\n');
fprintf(fid,'#PBS -S /bin/sh\n');
fprintf(fid,['#PBS -N AMICA_15_' filename(1:end-4) '\n']);
fprintf(fid,['#PBS -M ' email '\n']);
fprintf(fid,'#PBS -m abe\n\n');
fprintf(fid,['#PBS -A ' account '\n']);
fprintf(fid,'#PBS -l qos=flux\n');
fprintf(fid,'#PBS -q flux\n\n');
fprintf(fid,'#PBS -l procs=%d,mem=%s,walltime=%s \n',procs,mem,walltime); %resources to allocate to job
fprintf(fid,'#PBS -j oe\n');
fprintf(fid,'#PBS -V\n\n');
fprintf(fid,'####  End PBS preamble\n\n');
fprintf(fid,'#  Put your job commands after this line\n');

%Actual job commands
fprintf(fid,['/nfs/hnldata/stepeter/mpich-install/bin/mpiexec -np ' num2str(procs) ' ' amicaBinaryPath ' ' pFlux 'DataMerge_' group '_amica.param']);


fclose(fid);

