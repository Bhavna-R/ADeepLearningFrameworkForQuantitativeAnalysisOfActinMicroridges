%% A deep learning framework for quantitative analysis of actin microridges

%% Rajasekaran Bhavna1,2*, Mahendra Sonawane1

%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

%% 2 Current Address: Department of Data Science and Engineering, Indian Institute of Science Education and Research, Bhopal, Madhya Pradesh- 462066 

%% *Corresponding author email: bhavnarajasekaran@yahoo.com


function [PeridermCellsWithMicroridges_TracksLong, ALL_PeridermCellsWithMicroridges]= TrackPeridermCellCentroids(Cell_centroids,ExtractCellsWithMicroridges,VoxelSizeX_microns,TotalTime)

%% required outputs from: Periderm_wt_MicroridgesCellSegmentation.m

%% PeridermCellTracks_long (row-by-column  cell array)

%% Each row is a single tracked cell (grayscale periderm cell patterned with microridges)

%% PeridermCellsWithMicroridges_TracksLong: Trackslengths atleast 1/2(totaltime imaged) are considered

%% ALL_PeridermCellsWithMicroridges: All TrackLengths

pixel_distance=round(2.5/VoxelSizeX_microns);% 2.5 microns
distance_threshold=pixel_distance;% in pixel uints almost

tim=find(cellfun(@isempty,Cell_centroids));
if isempty(tim), time_length=(length(Cell_centroids)-1);
else
    time_length=(tim(1)-2);
end

for g=1:time_length
    Distbwtimes{1,g}=pdist2(Cell_centroids{g},Cell_centroids{g+1},'euclidean');
    Distbwtimes{g}(isnan(Distbwtimes{g}))=5000;
end

for g=1:time_length
    [r,c]= size(Distbwtimes{g});
    for x=1:r
        if (min(Distbwtimes{g}(x,:))<distance_threshold)
            [q w]=min(Distbwtimes{g}(x,:));
            gp{g}(x,1)=x;
            gp{g}(x,2)=w;
        else
            gp{g}(x,:)=[0 0];
        end
    end
end

for h=1:length(gp)
    gp{h}(~any(gp{h},2),:)=[];
end

getpairs=gp;
clear assemble
assemble=getpairs{1,1};

for  m=1:(length(getpairs)-2)
    clear c K f a b match* miss*
    [a,match1,match2]=intersect(getpairs{m+1}(:,1),assemble(:,m+1));
    [b,miss1,miss2]=setxor(getpairs{m+1}(:,1),assemble(:,m+1));

    if ~isempty(a)
        for d=1:length(a)
            assemble(match2(d),m+2)=getpairs{m+1}(match1(d),2);
        end
    else
        c=find(assemble(:,m+1));
        for s=1:length(c)
            assemble(c(s),m+2)=0;
        end
        if isempty(c)
            assemble(:,m+2)=0;
        end
    end

    if exist('b')
        for r=1:length(b)
            f=find(getpairs{m}(:,2)==b(r));
            if ~isempty(f)
                K=find(assemble(:,m)==getpairs{m}(f,1));
                assemble(K,m+1)=b(r);
            else
                assemble(end+1,m+1)=b(r);
            end
        end
    end
end

if exist('b')
    m=m+1;
    for r=1:length(b)
        f=find(getpairs{m}(:,1)==b(r));
        if ~isempty(f)
            K=find(assemble(:,m)==getpairs{m}(f,1));
            assemble(K,m+1)=getpairs{m}(f,2);
        end
    end
end

tracklinks=assemble;
tracklinks(all(~tracklinks,2),:)=[];
[tr_r tr_c]=size(tracklinks);

for r=1:tr_r
    trkind=tracklinks(r,:);
    for l=1:length(getpairs)
        if (trkind(l)~=0)
            ALL_PeridermCellsWithMicroridges{r,l}=ExtractCellsWithMicroridges{l}{trkind(l)};
        end
    end
end
CellTracklengths=(~cell2mat(cellfun(@isempty,ALL_PeridermCellsWithMicroridges,'UniformOutput',false)));
Considertracksind=(sum(CellTracklengths,2)>round(TotalTime/2));
PeridermCellsWithMicroridges_TracksLong=ALL_PeridermCellsWithMicroridges(Considertracksind,:);
