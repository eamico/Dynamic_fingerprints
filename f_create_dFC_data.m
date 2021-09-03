function dFCw_2D = f_create_dFC_data(configs)
%% This function computes dynamic FC profiles and stores them (upper triangular) in a 2D matrix dFCw_2D

wSize = configs.wSize; 
wjump = configs.wjump; 
numTP = configs.numTP;
winit_last = numTP-wSize+1;
numLayers = length(1:wjump:winit_last);
j=1;

dFCw_2D = nan(configs.numSubj,nnz(configs.mask_ut),numLayers);
ts_3D = importdata(configs.fMRI_file);
ts_3D = double(ts_3D);
for i=1:configs.numSubj
    ts = ts_3D(:,:,i);
    %% run sliding window dFC computation
    for l=1:numLayers 
        winit = ((l-1)*wjump)+1;
        wend = winit + wSize-1;
        ts_segment = ts(winit:wend,:);
        fc_segment = corr(ts_segment);                   
        fc_segment = (fc_segment + fc_segment')./2;
        fc_segment(logical(eye(size(fc_segment)))) = 0;
        fc_segment(isnan(fc_segment)) = 0;           
        dFCw_2D(j,:,l) = fc_segment(configs.mask_ut);            
    end
    j=j+1;   
end