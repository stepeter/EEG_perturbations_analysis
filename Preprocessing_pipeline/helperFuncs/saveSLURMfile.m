function []=saveSLURMfile(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group,maxthreads)
%Creates qsub file to run AMICA on Flux

amicaBinaryPath='/ufrc/dferris/s.peterson/test/AMICA_15/amica15ub';

fid = fopen([ICA_path filesep group 'qsub.sh'],'w');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,['#SBATCH --job-name=AMICA_15_' filename(1:end-4) ' # Job name\n']);
fprintf(fid,'#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n');
fprintf(fid,['#SBATCH --mail-user=' email '  # Where to send mail\n']);
fprintf(fid,['#SBATCH --ntasks=' num2str(maxthreads) '                   # Run a single task\n']);
fprintf(fid,'#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n',procs);
fprintf(fid,'#SBATCH --nodes=2            # Number of CPU cores per task\n');
fprintf(fid,['#SBATCH --mem=' mem '                  # Total memory limit\n']);
fprintf(fid,['#SBATCH --time=' walltime '             # Time limit hrs:min:sec\n']);
fprintf(fid,'#SBATCH --output=amicaRun_%%j.out     # Standard output and error log\n');
fprintf(fid,['#SBATCH --account=' account ' 	     # Account name\n']); %dferris
fprintf(fid,['#SBATCH --qos=' account '-b		     # Quality of service name\n\n']); %dferris
fprintf(fid,'# Run your program with correct path and command line options\n');

%Actual job commands
fprintf(fid,'module load ufrc\n');
fprintf(fid,'module load intel openmpi\n');
fprintf(fid,['mpirun ' amicaBinaryPath ' ' pFlux group '_amica.param']);

fclose(fid);