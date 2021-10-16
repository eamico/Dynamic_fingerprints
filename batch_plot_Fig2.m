%% This script creates the dynamic identifiability matrices 
%% at different window lenghts of dynamic functional connectomes
%% as proposed in (Van De Ville et al., Science Advances 2021).
%% It also replicates Figure 2 of the paper, on 10 sample HCP subjects.
%% 
%% Enrico Amico, EPFL
%% version 1.2 October, 2021
%
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
%% Dimitri Van De Ville, Younes Farouj, Maria Giulia Preti Raphael Liegeois and Enrico Amico. 
%% When makes you unique: temporality of the human brain fingerprint. Science Advances, 2021. 

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
configs.wSizeRange = [10 50 100 200 400 800];  % dFC window lenght (number of time points)
configs.numSubj = 10; % Number of sampled subjects      
Time_scales_str = {'7.2s' '36s' '72s' '144s' '288s' '576s'};
n_wsize = length(configs.wSizeRange); % number of windows explored


flags.ComputeIdent = 1;
%% 1.0 Dynamic Identifiability matrix at different window lengths 
if flags.ComputeIdent==1
    Ident_mat = cell(1,n_wsize);
    dIdent_LR = cell(1,n_wsize);
    dIdent_RL = cell(1,n_wsize);   
    disp('Computing dynamic Identifiability (might take a while)...')
    for t=1:n_wsize
        disp(Time_scales_str{t});
        configs.wSize = configs.wSizeRange(t);
        configs.fMRI_file = 'FCs_10S_R1_LR.mat';
        dFCw_2D_Test_LR = f_create_dFC_data(configs);
        configs.fMRI_file = 'FCs_10S_R2_LR.mat';       
        dFCw_2D_Retest_LR = f_create_dFC_data(configs);
        dIdent_LR{t} = f_create_dIdent_tensor(dFCw_2D_Test_LR,dFCw_2D_Retest_LR);   
        configs.fMRI_file = 'FCs_10S_R1_RL.mat';
        dFCw_2D_Test_RL = f_create_dFC_data(configs);
        configs.fMRI_file = 'FCs_10S_R2_RL.mat';       
        dFCw_2D_Retest_RL = f_create_dFC_data(configs);
        dIdent_RL{t} = f_create_dIdent_tensor(dFCw_2D_Test_RL,dFCw_2D_Retest_RL);
    end
end
%% 2.0 Extract dynamic Iself, Iothers and Idiff

Ident_mean_LR = zeros(configs.numSubj,configs.numSubj,n_wsize);
Ident_mean_RL = zeros(configs.numSubj,configs.numSubj,n_wsize);

mask_diag = logical(eye(configs.numSubj));
dIdiff = nan(1,n_wsize);
dIself = nan(1,n_wsize);
dIothers = nan(1,n_wsize);


for t=1:n_wsize
    aux = dIdent_LR{t};
    nFrames = size(aux,1)./configs.numSubj;
    for s1=1:configs.numSubj
        for s2=1:configs.numSubj
            IndexRow = ((s1-1)*nFrames)+1:(s1*nFrames);
            IndexCol = ((s2-1)*nFrames)+1:(s2*nFrames);
            tmp = dIdent_LR{t}(IndexRow,IndexCol);
            Ident_mean_LR(s1,s2,t) = nanmean(tmp(:));
            tmp = dIdent_RL{t}(IndexRow,IndexCol);
            Ident_mean_RL(s1,s2,t) = nanmean(tmp(:));
        end          
    end
    Ident_mean = 0.5.*(Ident_mean_LR+Ident_mean_RL); 
    tmp = Ident_mean(:,:,t);
    dIself(t) = nanmean(tmp(mask_diag));
    dIothers(t) = nanmean(tmp(~mask_diag));
    dIdiff(t) = dIself(t) - dIothers(t);  
end
%% 3.0 Plot Dynamic Identifiability matrix
% Note that the dID matrices should not be averaged across LR and RL sessions
% (since they represent different resting-state time frames). However, the
% average across LR and RL session can be performed at the summary statistic stage 
% (i.e. dIself, dIothers, dIdiff computation, see lines 80-84)
figure, 
for t=1:n_wsize
    subplot(2,3,t); imagesc(dIdent_LR{t},[.1 .6]); colorbar; title([int2str(configs.wSizeRange(t).*configs.TR./1000) 's']);
    xlabel('TFrames x Subj retest'); ylabel('TFrames x Subj test');axis square; 
end
suptitle('Dynamic Identifiability matrix');
