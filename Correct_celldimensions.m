%% A deep learning framework for quantitative analysis of actin microridges

%% Rajasekaran Bhavna1,2*, Mahendra Sonawane1

%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

%% 2 Current Address: Department of Data Science and Engineering, Indian Institute of Science Education and Research, Bhopal, Madhya Pradesh- 462066 

%% *Corresponding author email: bhavnarajasekaran@yahoo.com


function [extractridgeswithincell]= Correct_celldimensions(extractcells,ConvImg)

%% helper function required during periderm cell extraction

[ca,cb]=size(extractcells);
[la,lb]=size(ConvImg);

clear new_extractcells
if (ca==la && cb==lb)
    extractridgeswithincell =(extractcells).*uint16(ConvImg);
elseif (cb~=lb && ca==la)
    if lb>cb
        get_diff=lb-cb;
        new_extractcells=padarray(extractcells,[0,get_diff],'post');
    else
        get_diff=cb-lb;
        new_extractcells=padarray(ConvImg,[0,get_diff],'post');
    end
    extractridgeswithincell =(new_extractcells).*uint16(ConvImg);

elseif (cb==lb && ca~=la)
    if la>ca
        get_diff=la-ca;
        new_extractcells=padarray(extractcells,[get_diff,0],'post');
    else
        get_diff=ca-la;
        new_extractcells=padarray(ConvImg,[get_diff,0],'post');
    end
    extractridgeswithincell =(new_extractcells).*uint16(ConvImg);
elseif (cb~=lb && ca~=la)

    if (cb>lb && ca>la)
        get_diff1=cb-lb;
        get_diff2=ca-la;
        new_extractcells=padarray(extractcells,[get_diff1,get_diff2],'post');
    elseif (cb>lb && ca<la)
        get_diff1=cb-lb;
        get_diff2=la-ca;
        new_extractcells=padarray(extractcells,[get_diff1,get_diff2],'post');
    elseif (cb<lb && ca>la)
        get_diff1=lb-cb;
        get_diff2=ca-la;
        new_extractcells=padarray(extractcells,[get_diff1,get_diff2],'post');
    elseif (cb<lb && ca<la)
        get_diff1=lb-cb;
        get_diff2=la-ca;
        new_extractcells=padarray(extractcells,[get_diff1,get_diff2],'post');
    end
    extractridgeswithincell =(new_extractcells).*uint16(ConvImg);

end
