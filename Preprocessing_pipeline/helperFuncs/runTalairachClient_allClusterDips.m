%Spits out an file to feed the talairach client to determine where each
%cluster dipole is in the brain (assumes same dipoles across conditions)

clusterNums=[3 4 6:12]; %select clusters to look at (can usually do up to 6 at a time)
outpath=[STUDY.filepath filesep 'talairachDips']; %folder to put output in

for clusternum=clusterNums
    disp(clusternum);
    numDips=length(STUDY.cluster(1,clusternum).comps);
    tal_coords=zeros(numDips+1,3);
    for i=1:numDips
        setnum=STUDY.cluster(1,clusternum).sets(1,i);
        dipnum=STUDY.cluster(1,clusternum).comps(1,i);
        mni_coords=ALLEEG(1,setnum).dipfit.model(1,dipnum).posxyz;
        tal_coords(i,:)=mni2tal(mni_coords); %round(mni2tal(mni_coords))
    end
    tal_coords(end,:)=mni2tal(STUDY.cluster(1,clusternum).dipole.posxyz); %add in mean dipole

    %Now write them out to a file
    fid=fopen([outpath filesep STUDY.filename(1:end-6) '_cluster_' num2str(clusternum) '.txt'],'w');
    for i=1:(numDips+1)
        fprintf(fid,'%g\t %g\t %g\n',tal_coords(i,1),tal_coords(i,2),tal_coords(i,3));
    end
    fclose(fid);
end

%% Read the file that Talairach client spits out and determine Brodmann areas of cluster
clusts={'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen','twenty'};
outpath='/home/stepeter/Desktop/MPT_060716_cortney';
for clusternum=clusterNums
    disp(clusternum)
    BA_data=importdata([outpath filesep STUDY.filename(1:end-6) '_cluster_' num2str(clusternum) '.td.txt']);
    BA_nums.(clusts{1,clusternum})=zeros(size(BA_data.textdata,1),1);
    for i=2:size(BA_data.textdata,1)
        temp=str2double(BA_data.textdata{i,9});
        if ~isempty(temp) && ~all(isnan(temp))
            BA_nums.(clusts{1,clusternum})(i-1,1)=temp(1,3);
        end
    end
end


%% Find Std Dev for each cluster (calculated as r=sqrt(x^2+y^2+z^2))
stds=zeros(length(clusterNums),1); q=0;
for clusternum=clusterNums
    q=q+1;
    disp(clusternum);
    numDips=length(STUDY.cluster(1,clusternum).comps);
    tal_coords=zeros(numDips+1,3);
    for i=1:numDips
        setnum=STUDY.cluster(1,clusternum).sets(1,i);
        dipnum=STUDY.cluster(1,clusternum).comps(1,i);
        mni_coords=ALLEEG(1,setnum).dipfit.model(1,dipnum).posxyz;
        tal_coords(i,:)=mni2tal(mni_coords); %round(mni2tal(mni_coords))
    end
    tal_coords(end,:)=mni2tal(STUDY.cluster(1,clusternum).dipole.posxyz); %add in mean dipole
    temp_std=std(tal_coords);
    stds(q,1)=sqrt(temp_std(:,1).^2+temp_std(:,2).^2+temp_std(:,3).^2);
end
