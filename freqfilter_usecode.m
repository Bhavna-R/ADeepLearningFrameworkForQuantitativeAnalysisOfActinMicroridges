%% https://www.mathworks.com/matlabcentral/fileexchange/40579-frequency-domain-filtering-for-grayscale-images
function freqfilterimg = freqfilter_usecode(image2d,thresh,n)
cim=double(image2d);
[r,c]=size(cim);

r1=2*r;
c1=2*c;

pim=zeros((r1),(c1));
kim=zeros((r1),(c1));

%padding
for i=1:r
    for j=1:c
        pim(i,j)=cim(i,j);
    end
end

%center the transform
for i=1:r
    for j=1:c
        kim(i,j)=pim(i,j)*((-1)^(i+j));
    end
end

fim=fft2(kim);
%%

%%

%n=1; %order for butterworth filter
%thresh=7;

% % function call for low pass filters
 %him=glp(fim,thresh); % gaussian low pass filter
% him=blpf(fim,thresh,n); % butterworth low pass filter

% % function calls for high pass filters
%him=ghp(fim,thresh); % gaussian high pass filter
him=bhp(fim,thresh,n);  %butterworth high pass filter

% % function call for high boost filtering
%him=hbg(fim,thresh);  % using gaussian high pass filter
% him=hbb(fim,thresh,n);  % using butterworth high pass filter

%inverse 2D fft
ifim=ifft2(him);

for i=1:r1
    for j=1:c1
        ifim(i,j)=ifim(i,j)*((-1)^(i+j));
    end
end
% removing the padding
for i=1:r
    for j=1:c
        rim(i,j)=ifim(i,j);
    end
end
% retaining the real parts of the matrix
rim=real(rim);
freqfilterimg=rim;
end