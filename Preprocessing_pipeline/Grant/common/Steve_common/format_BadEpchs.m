function badEpchs=format_BadEpchs(badEpchs)

%Function for manually creating badEpchs cell array (num subjs x num trials).
%Used to reduce clutter in main script
%Num trials are in alphabetical order

badEpchs{4,2}=1; %NS
badEpchs{6,3}=[1,29]; %S
badEpchs{2,3}=4; %S
badEpchs{7,2}=[1,44]; %NS
badEpchs{7,3}=[15,19]; %S
badEpchs{5,3}=[3]; %S
badEpchs{10,3}=[3]; %S
badEpchs{12,1}=[23,26]; %C
badEpchs{12,2}=[1,2,6,7,9,13,14,16,20,22,23,27,29,31,32];%24,26,28,30,31]; %NS
badEpchs{12,3}=[2:5,8,10,14,17,20,24,25]; %S
badEpchs{14,3}=1:8; %S
badEpchs{16,2}=[17]; %NS
badEpchs{18,2}=[1]; %NS
badEpchs{20,3}=[3,4,5,8,14,15]; %S

%Subjects that went NS-C-S: 1,2,4,5,10,11,12,15,18,19

%Subjects that went S-C-NS: 3,6,7,8,9,13,14,16,17,20
end
