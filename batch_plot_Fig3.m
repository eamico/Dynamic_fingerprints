%% This script replicates the iterative frame selection result 
%% at different window lenghts for dynamic functional connectomes
%% as proposed in (Van De Ville et al., Science Advances 2021).
%% It also replicates Figure 3 of the paper, on 10 sample HCP subjects.
%% 
%% Enrico Amico, EPFL
%% version 1.2 October, 2021
%
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
%% Dimitri Van De Ville, Younes Farouj, Maria Giulia Preti Raphael Liegeois and Enrico Amico. 
%% When makes you unique: temporality of the human brain fingerprint. Science Advances, 2021. 

%% Plot Iterative Frame selection results
% EA, October 2021
%% initialize environment
clearvars;
close all;
clc
%% Configuration 
configs.wjump = 10; % Window step (in time points)
configs.numTP = 1200; % Number of time points
configs.TR = 720; % HCP TR
configs.parc = 'Schaefer400'; % Schaefer + Subcortical regions
    configs.Nparc = 419; % Number of brain regions 
    configs.mask_ut = triu(true(configs.Nparc,configs.Nparc),1); % Upper triangular mask 
    load('yeo_RS7_Schaefer400S.mat'); % Yeo Ordering of the matrices (Left and Right hemispheres together)  
configs.wSizeRange = [10 50 100 200 400 800];  % dFC window lenght (number of time points)
configs.numSubj = 10; % Number of sampled subjects      
configs.k_frames = 1:2:41; % Iterative frame selection range
n_wsize = length(configs.wSizeRange); % number of windows explored
Time_scales_str = {'7.2s' '36s' '72s' '144s' '288s' '576s'};
flags.FrameSelection_iter = 1; % execute iterative best frame selection
flags.plotIdiff_k=1; % Plot Fig 3A of Van de Ville et al, Science Advances
flags.plotStdFrames = 1; % Plot Fig 3B-C of Van de Ville et al, Science Advances

%% Iterative Frame Selection on Dynamic identifiability matrix (Top 1, 5, 10 etc.)
if flags.FrameSelection_iter==1
    n_wsize = length(configs.wSizeRange);
    Ident_mat = cell(n_wsize,configs.numSubj);
    dIdent = cell(1,n_wsize); 
    Row_index = cell(n_wsize,configs.numSubj);
    Col_index = cell(n_wsize,configs.numSubj);
    Top_Frame_test = zeros(nnz(configs.mask_ut),configs.numSubj,n_wsize);
    Top_Frame_retest = zeros(nnz(configs.mask_ut),configs.numSubj,n_wsize);
    dID_LR = cell(1,n_wsize);
    dID_RL = cell(1,n_wsize);
    disp('Computing dynamic Identifiability...')
    for t=1:n_wsize
        disp(Time_scales_str{t});
        configs.wSize = configs.wSizeRange(t);
        configs.fMRI_file = 'FCs_10S_R1_LR.mat';
        % create dynamic functional connectomes at specific time scale t
        dFCw_2D_Test_LR = f_create_dFC_data(configs);
        configs.fMRI_file = 'FCs_10S_R2_LR.mat';       
        dFCw_2D_Retest_LR = f_create_dFC_data(configs);
        configs.fMRI_file = 'FCs_10S_R1_RL.mat';
        dFCw_2D_Test_RL = f_create_dFC_data(configs);
        configs.fMRI_file = 'FCs_10S_R2_RL.mat';       
        dFCw_2D_Retest_RL = f_create_dFC_data(configs);
        % extract best k matching frames
        [dID_LR{t},k_Frames_testLR,k_Frames_retestLR] = f_compute_Idiff_frames(dFCw_2D_Test_LR,dFCw_2D_Retest_LR,configs);
        [dID_RL{t},k_Frames_testRL,k_Frames_retestRL] = f_compute_Idiff_frames(dFCw_2D_Test_RL,dFCw_2D_Retest_RL,configs);
        Top_Frame_test(:,:,t) = (k_Frames_testLR + k_Frames_testRL)./2; 
        Top_Frame_retest(:,:,t) = (k_Frames_retestLR + k_Frames_retestRL)./2;
    end
end

%% Plot dIdiff dIself dIothers across top k dFC frames (Fig. 3A Van De Ville et al.)
if flags.plotIdiff_k==1
    dIdiff_avg = zeros(n_wsize,length(configs.k_frames));
    dIself_avg = zeros(n_wsize,length(configs.k_frames));
    dIothers_avg = zeros(n_wsize,length(configs.k_frames));
    dIdiff_std = zeros(n_wsize,length(configs.k_frames));
    dIself_std = zeros(n_wsize,length(configs.k_frames));
    dIothers_std = zeros(n_wsize,length(configs.k_frames));
    for t=1:n_wsize
        dIdiff_avg(t,:) = 0.5.*(dID_LR{t}.Idiff + dID_RL{t}.Idiff);
        dIself_avg(t,:) = 0.5.*(dID_LR{t}.Iself + dID_RL{t}.Iself);
        dIothers_avg(t,:) = 0.5.*(dID_LR{t}.Iothers + dID_RL{t}.Iothers);
        dIdiff_std(t,:) =  0.5.*(dID_LR{t}.Idiff_std + dID_RL{t}.Idiff_std);
        dIself_std(t,:) =  0.5.*(dID_LR{t}.Iself_std + dID_RL{t}.Iself_std);
        dIothers_std(t,:) =  0.5.*(dID_LR{t}.Iothers_std + dID_RL{t}.Iothers_std);
    end
    
    figure
    bar_colors = {'k' 'c' 'm' 'r' 'g' 'b'};
    for t=1:n_wsize
        subplot(1,3,1);
        shadedErrorBar(configs.k_frames,dIself_avg(t,:).*100,(dIself_std(t,:)./sqrt(configs.numSubj)).*100,'lineprops', {[bar_colors{t} '-o'],'markerfacecolor',bar_colors{t}}); hold on;      
        ylabel('dIdiff score (%)'); xlabel('Top K FC_t frames'); axis square
        subplot(1,3,2);
        shadedErrorBar(configs.k_frames,dIothers_avg(t,:).*100,(dIothers_std(t,:)./sqrt(configs.numSubj)).*100,'lineprops', {[bar_colors{t} '-o'],'markerfacecolor',bar_colors{t}}); hold on;               
        ylabel('dIdiff score (%)'); xlabel('Top K FC_t frames'); axis square
        subplot(1,3,3);
        shadedErrorBar(configs.k_frames,dIdiff_avg(t,:).*100,(dIdiff_std(t,:)./sqrt(configs.numSubj)).*100,'lineprops', {[bar_colors{t} '-o'],'markerfacecolor',bar_colors{t}}); hold on;
        ylabel('dIdiff score (%)'); xlabel('Top K FC_t frames'); axis square
    end
end
%% Plot std and median of Top frames across subjects (Fig. 3B, 3C, Van De Ville et al.)
if flags.plotStdFrames==1
    
    SD_3D = zeros(configs.Nparc,configs.Nparc,n_wsize);
    figure, 
    for t=1:n_wsize
        aux_test = Top_Frame_test(:,:,t);
        aux_retest = Top_Frame_retest(:,:,t);
        SD_test = squeeze(std(aux_test,[],2));
        SD_retest = squeeze(std(aux_retest,[],2));
        SD_mat = zeros(configs.Nparc,configs.Nparc);
        SD_mat(configs.mask_ut) = (SD_test + SD_retest)./2;
        SD_mat = SD_mat + SD_mat';
        SD_3D(:,:,t) = SD_mat;
        ul = prctile(SD_mat(configs.mask_ut),95); ll = prctile(SD_mat(configs.mask_ut),5);
        subplot(2,3,t); imagesc(SD_mat(yeoOrder,yeoOrder),[ll ul]); 
        title({'Edgewise SD',['Top frame: ' Time_scales_str{t}]}); axis square; colorbar
    end
    
    Med_3D = zeros(configs.Nparc,configs.Nparc,n_wsize);
    figure, 
    for t=1:n_wsize
        aux_test = Top_Frame_test(:,:,t);
        aux_retest = Top_Frame_retest(:,:,t);
        Med_test = squeeze(median(aux_test,2));
        Med_retest = squeeze(median(aux_retest,2));
        Med_mat = zeros(configs.Nparc,configs.Nparc);
        Med_mat(configs.mask_ut) = (Med_test + Med_retest)./2;
        Med_mat = Med_mat + Med_mat';
        Med_3D(:,:,t) = Med_mat;
        ul = prctile(Med_mat(configs.mask_ut),95); ll = prctile(Med_mat(configs.mask_ut),5);
        subplot(2,3,t); imagesc(Med_mat(yeoOrder,yeoOrder),[-.7 .7]); axis square; colorbar
        title({'Edgewise Median',['Top frame: ' Time_scales_str{t}]}); axis square; colorbar
    end

end   
