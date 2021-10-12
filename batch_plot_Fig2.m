%% This script creates the dynamic identifiability matrices 
%% at different window lenghts of dynamic functional connectomes
%% as proposed in (Van De Ville et al., Science Advances 2021).
%% It also replicates Figure 2 of the paper, on 10 sample HCP subjects.
%% 
%% Enrico Amico, EPFL
%% version 1.1 September, 2021
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

n_wsize = length(configs.wSizeRange); % number of windows explored


flags.ComputeIdent = 1;
%% 1.0 Dynamic Identifiability matrix at different window lengths 
if flags.ComputeIdent==1
    Ident_mat = cell(1,n_wsize);
    % now average the dy
    dIdent_LR = cell(1,n_wsize);
    dIdent_RL = cell(1,n_wsize);   
    disp('Computing dynamic Identifiability (might take a while)...')
    for t=1:n_wsize
        disp(t);
        configs.wSize = configs.wSizeRange(t);
        configs.fMRI_file = 'FCs_10S_R1_LR.mat';
        dFCw_2D_Test = f_create_dFC_data(configs);
        configs.fMRI_folder = 'FCs_10S_R2_LR.mat';       
        dFCw_2D_Retest = f_create_dFC_data(configs);
        dIdent_LR{t} = f_create_dIdent_tensor(dFCw_2D_Test,dFCw_2D_Retest);   
        configs.fMRI_folder = 'FCs_10S_R1_RL.mat';
        dFCw_2D_Test = f_create_dFC_data(configs);
        configs.fMRI_folder = 'FCs_10S_R2_RL.mat';       
        dFCw_2D_Retest = f_create_dFC_data(configs);
        dIdent_RL{t} = f_create_dIdent_tensor(dFCw_2D_Test,dFCw_2D_Retest);
    end
end
%% 2.0 Extract dynamic Iself, Iothers and Idiff
% Note that these scores (as well as success rates) might fluctuate in small sample sizes

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
%% 3.0 Plot Dynamic Identifiability matrix matrix
% Note that the dID matrices should not be averaged across LR and RL sessions
% (since they represent different resting-state time frames). However, the
% average across LR and RL session can be performed at the summary statistic stage 
% (i.e. dIself, dIothers, dIdiff computation, see lines 83-87)
figure, 
for t=1:n_wsize
    subplot(2,3,t); imagesc(dIdent_LR{t},[.1 .6]); colorbar; title([int2str(configs.wSizeRange(t).*configs.TR./1000) 's']);
    xlabel('TFrames x Subj retest'); ylabel('TFrames x Subj test');axis square; 
end
suptitle('Dynamic Identifiability matrix');
