function []=saveSLURMfile_jobArray(ICA_path,pFlux,filename,email,account,procs,mem,walltime,group,arrayNums)
%Creates qsub file to run AMICA on Flux

amicaBinaryPath='/ufrc/dferris/s.peterson/eeglab13_5_4b/plugins/AMICA_15/amica15ub';

fid = fopen([ICA_path filesep group 'qsub.sh'],'w');
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,['#SBATCH --job-name=AMICA_15_' filename(1:end-4) ' # Job name\n']);
fprintf(fid,'#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n');
fprintf(fid,['#SBATCH --mail-user=' email '  # Where to send mail\n']);
fprintf(fid,'#SBATCH --ntasks=1                   # Run a single task\n');
fprintf(fid,'#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n',procs);
fprintf(fid,['#SBATCH --mem=' mem '                  # Total memory limit\n']);
fprintf(fid,['#SBATCH --time=' walltime '             # Time limit hrs:min:sec\n']);
fprintf(fid,'#SBATCH --output=amicaRun_%%j.out     # Standard output and error log\n');
fprintf(fid,['#SBATCH --account=' account ' 	     # Account name\n']); %dferris
fprintf(fid,['#SBATCH --qos=' account '		     # Quality of service name\n']); %dferris
fprintf(fid,['#SBATCH --array=' arrayNums ' 	     # Account name\n\n']); %dferris
fprintf(fid,'# Run your program with correct path and command line options\n');

%Actual job commands
fprintf(fid,'module load ufrc\n');
fprintf(fid,'module load intel\n');
fprintf(fid,['/ufrc/dferris/s.peterson/mpich-install/bin/mpiexec -np ' num2str(procs) ' ' amicaBinaryPath ' ' pFlux group  '_amica${SLURM_ARRAY_TASK_ID}.param']);

fclose(fid);