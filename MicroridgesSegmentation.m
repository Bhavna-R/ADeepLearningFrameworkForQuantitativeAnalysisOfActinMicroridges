%% A deep learning framework for quantitative analysis of actin microridges

%% Rajasekaran Bhavna1,2*, Mahendra Sonawane1

%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

%% 2 Current Address: Department of Data Science and Engineering, Indian Institute of Science Education and Research, Bhopal, Madhya Pradesh- 462066 

%% *Corresponding author email: bhavnarajasekaran@yahoo.com

s
function [Microridges_Binary, Microridges_Binary_Skeleton]=MicroridgesSegmentation(PeridermCellsWithMicroridges_TracksLong,ALL_PeridermCellsWithMicroridges)

%% required outputs from: TrackPeridermCellCentroids.m

%% option 1
%% [cell_nos,timepnts]=size(ALL_PeridermCellsWithMicroridges);
%% option 2
 [cell_nos,timepnts]=size(PeridermCellsWithMicroridges_TracksLong);
 sigmagradient=0.7;% parameter

for f=1:cell_nos
    for g=1:timepnts
        A=ALL_PeridermCellsWithMicroridges{f,g};
        % A=PeridermCellsWithMicroridges_TracksLong; %option2

        if ~isempty(A)
            J= imgaussfilt(double(A),0.7);
            fim=mat2gray(J,[0 max(max(J))]);
            gim=fim;
            gim(gim==0)=2;

            [imx,imy]=gaussgradient(gim,sigmagradient);

            [L_imxx L_imxy] = gaussgradient(imx,sigmagradient);
            [L_imyx L_imyy] = gaussgradient(imy,sigmagradient);

            Trace_lap=L_imxx+L_imyy;
            G_curv=Trace_lap;
            G_curv=G_curv.*(G_curv<0);
            S=1./(1+exp(-abs(sqrt(abs(G_curv)))));
            bin_S=imbinarize(S);
            M=J.*(bin_S);
            Microridges_Binary{f,g}= imbinarize(M);
            Microridges_Binary_Skeleton{f,g} = bwskel(Microridges_Binary{f,g});
        end
        clear A J M fim gim imx imy  L_imxx L_imxy L_imyx L_imyy Trace_lap G_curv S bin_S M
    end
end
