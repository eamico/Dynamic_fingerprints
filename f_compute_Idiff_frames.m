function  [dID,k_Frames_test,k_Frames_retest] = f_compute_Idiff_frames(dFCw_2D_Test,dFCw_2D_Retest,configs)
%% This function evaluates dynamic brain fingerprints across temporal scales 
%% through dlself, dlothers and dldiff scores, when ranking dFC frames based 
%% on individual dIself, in descending order (see also Van De Ville et al., Science Advances 2021)

Ident_mat = cell(1,configs.numSubj);
Row_index = cell(1,configs.numSubj);
Col_index = cell(1,configs.numSubj);

disp('Frame selection, iterative')
for i= 1:configs.numSubj
    % extract test and retest frames for subject i
    aux_test = squeeze(dFCw_2D_Test(i,:,:));
    aux_retest = squeeze(dFCw_2D_Retest(i,:,:));
    % Iself "Identifiability sub-block" for subject i
    Ident_mat{i} = corr(aux_test,aux_retest); 
    % select best frames iteratively
    [Row_index{i},Col_index{i}] = f_select_frames_iter(Ident_mat{i});
end

dIdent = f_create_dIdent_tensor(dFCw_2D_Test,dFCw_2D_Retest); % create dynamic identifiability matrix
aux = dIdent;
nFrames = size(aux,1)./configs.numSubj;
k_Frames_test = zeros(nnz(configs.mask_ut),configs.numSubj);
k_Frames_retest = zeros(nnz(configs.mask_ut),configs.numSubj);
Ident_mean = zeros(configs.numSubj,configs.numSubj);
k_frames = configs.k_frames;
dID.Idiff = zeros(1,length(k_frames));
dID.Iself = zeros(1,length(k_frames));
dID.Iothers = zeros(1,length(k_frames));
dID.Iself_std = zeros(1,length(k_frames));
dID.Iothers_std = zeros(1,length(k_frames));
dID.Idiff_std = zeros(1,length(k_frames));
for k=1:length(k_frames)
    for s1=1:configs.numSubj
        for s2=1:configs.numSubj
            IndexRow = ((s1-1)*nFrames)+1:(s1*nFrames);
            IndexCol = ((s2-1)*nFrames)+1:(s2*nFrames);
            tmp = dIdent(IndexRow,IndexCol);
            tmp_fltd = tmp(Row_index{s1}(1:k_frames(k)),Col_index{s2}(1:k_frames(k)));
            Ident_mean(s1,s2) = nanmean(tmp_fltd(:));
        end
        if k==1 % retain Top best matching dynamic frame (test and retest)
            k_Frames_test(:,s1) = dFCw_2D_Test(s1,:,Row_index{s1}(1:k_frames(k)));
            k_Frames_retest(:,s1) = dFCw_2D_Retest(s1,:,Col_index{s1}(1:k_frames(k)));
        end
    end
    tmp = Ident_mean(:,:);
    mask_diag = logical(eye(configs.numSubj));
    Iself_vec = tmp(mask_diag);
    tmp(mask_diag) = nan;
    Iothers_vec = 0.5 .* (nanmean(tmp,2) + nanmean(tmp,1)');
    Idiff_vec = Iself_vec - Iothers_vec;
    dID.Iself(k) = nanmean(Iself_vec);
    dID.Iothers(k) = nanmean(Iothers_vec);
    dID.Idiff(k) = dID.Iself(k) - dID.Iothers(k);
    dID.Iself_std(k) = std(Iself_vec);
    dID.Iothers_std(k) = std(Iothers_vec); 
    dID.Idiff_std(k) = std(Idiff_vec);
end
