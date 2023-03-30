%% A deep learning framework for quantitative analysis of actin microridges

%% Rajasekaran Bhavna1,2*, Mahendra Sonawane1

%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

%% 2 Current Address: Department of Data Science and Engineering, Indian Institute of Science Education and Research, Bhopal, Madhya Pradesh- 462066 

%% *Corresponding author email: bhavnarajasekaran@yahoo.com



function [Cell_centroids, ExtractCellsWithMicroridges,ConvImg, VoxelSizeX_microns,TotalTime]=Periderm_wt_MicroridgesCellSegmentation(infile)
%% input: LSM file

%% infile=yolkcentre52hpf.lsm'; input microscopy file example

%% required addtional supporting files: LSMfileread.m, freqfilter_usecode.m, Correct_celldimensions.m

filepart=regexp(infile,'\.','split');
inputfile=LSMfileread(infile);
path=inputfile(1,1).filename;
path_parts=regexp(path,'/','split');
filepart=regexp(path_parts{end},'\.','split');
getpath = strfind(path, '/');
getdir= path(1:getpath(end));
image_bits=inputfile(1,1).bits;
Xpixels=inputfile(1,1).lsm.DimensionX;
Ypixels=inputfile(1,1).lsm.DimensionY;
ZSlicesinStack=inputfile(1,1).lsm.DimensionZ;
TotalTime=inputfile(1,1).lsm.DimensionTime
VoxelSizeX_microns=(inputfile(1,1).lsm.VoxelSizeX)*10^6;
VoxelSizeY_microns=(inputfile(1,1).lsm.VoxelSizeY)*10^6;
VoxelSizeZ_microns=(inputfile(1,1).lsm.VoxelSizeZ)*10^6;
Time_interval_seconds=inputfile(1,1).lsm.TimeInterval
xscale=VoxelSizeX_microns;
yscale=VoxelSizeY_microns;
zscale=VoxelSizeZ_microns;
duration_imaged=[num2str((TotalTime*Time_interval_seconds)/60),' mins']

%% set params
entl=3;
entu=7;
tmpent=1;
Sol=0.70;
Igauparam=3;
amp5param=2.8;
clbw=11;
dlbw=11;
setparams = struct('Sol',Sol,'entl',entl,'entu',entu,'tmpent',tmpent,'Igauparam',Igauparam,'amp5param',amp5param,'clbw',clbw,'dlbw',dlbw);

%cell_area_lowlim=3400;% flank cells
cell_area_lowlim=4500;% head/yolk cells
cell_area_upperlim=20000;
%%

for time=1:TotalTime

    for zslice=1:ZSlicesinStack
        Im(:,:,zslice)=inputfile(1,((time-1))*ZSlicesinStack+zslice).data;
        I{time}(:,:,zslice)=Im(:,:,zslice);
        Ientropy_overall(time,zslice)=entropy(Im(:,:,zslice));

        Ientropy(:,:,zslice)=entropyfilt(Im(:,:,zslice),true(3));

        if (Ientropy_overall(time,zslice)>entl && Ientropy_overall(time,zslice)<entu)

            tmp_entropy=Ientropy(:,:,zslice);
            tmp_entropy=tmp_entropy>tmpent;
            filteredentropy{time}(:,:,zslice)=uint16(tmp_entropy).* I{time}(:,:,zslice);
        else
            filteredentropy{time}(:,:,zslice)= zeros(Xpixels,Ypixels);
        end
        Ientrpmean{time}(:,:,:)= mean(filteredentropy{time}(:,:,:),3);
    end
end

for time=1:TotalTime

    Igaus=imgaussfilt(uint16(Ientrpmean{time}),Igauparam);

    freqfiltimg = freqfilter_usecode(Igaus,3,1);%default parameters
    freqfiltimg=uint16(freqfiltimg);

    amplitudeImage5=imgaussfilt(freqfiltimg,amp5param);
    imgl=imbinarize(uint16(amplitudeImage5));

    se = strel('disk',clbw);
    closeBW = imclose(imgl,se);

    se2 = strel('disk',dlbw);

    dilabw=imdilate(closeBW,se2);

    BWprp = bwpropfilt(dilabw,'Perimeter',1);
    masky=imcomplement(BWprp);
    L = bwlabel(masky);
    Bgaus = labeloverlay(Igaus,L);
    stats= regionprops(L,'Area','BoundingBox','ConvexImage','Centroid','Solidity');

    [~,index] = sortrows([stats.Area].'); stats = stats(index(end:-1:1)); clear index
    for s=1:length(stats), areastats(time,s)=stats(s).Area; end

    cells=sum([stats.Area]>cell_area_lowlim & [stats.Area]<cell_area_upperlim);

    if cells>0
        for l=1:length(stats)

            if (stats(l).Solidity > 0.82 && stats(l).Area>cell_area_lowlim && stats(l).Area<cell_area_upperlim)

                extractcells{time}{l}=imcrop(uint16(Ientrpmean{time}),[stats(l).BoundingBox]-1);
                isolatecell=extractcells{time}{l};
                ConvImg=stats(l).ConvexImage;
                ExtractCellsWithMicroridges{time}{l}= Correct_celldimensions(isolatecell,ConvImg);

                %figure;imshow(extractridgeswithincell{time}{l},[]);
                %title([num2str(stats(l).Area),' , ', num2str(stats(l).Solidity)]);
                Cell_centroids{time}(l,1:2)=stats(l).Centroid;
            else
                Cell_centroids{time}(l,1:2)=[NaN NaN];
            end
        end
    else
        continue;
    end
    disp(time);
    close all;
    clear Igaus freqfiltimg amplitudeImage5 imgl se se2 dilabw;
    clear BWprp masky L Bgaus stats;
end
