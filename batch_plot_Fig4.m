%% This script extracts the "time scales of brain fingerprints"  
%% as proposed in (Van De Ville et al., Science Advances 2021).
%% It also replicates Figure 4 of the paper.
%% 
%% Enrico Amico, EPFL
%% version 1.2 October, 2021
%
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
%% Dimitri Van De Ville, Younes Farouj, Maria Giulia Preti Raphael Liegeois and Enrico Amico. 
%% When makes you unique: temporality of the human brain fingerprint. Science Advances, 2021. 

%% Initialize the environment
clearvars;
close all;
%clc
%% Configuration
configs.TR = 720; % HCP TR
load('yeo_RS7_Schaefer400S.mat'); % Yeo Ordering of the matrices (Left and Right hemispheres together)  
configs.Nparc = 419; % Number of brain regions
configs.mask_ut = triu(true(configs.Nparc,configs.Nparc),1); % Upper triangular mask 
configs.numSubj = 100; % total number of subjects        
configs.wSizeRange = [10 50 100 200 400 800]; % dFC window lenght (number of time points)
    n_wsize = length(configs.wSizeRange); % number of windows explored
    
configs.FileName = 'Fig4_material.mat'; % load relevant files
% This mat files contains the edgewise intraclass correlation matrices for
% the different time scales (see Van De Ville et al., Science Advances for details)

load(configs.FileName);
    

    
flags.TimeScale1 = 1; % Replicate Panel 4A, Van De Ville et al., Science Advances 
flags.TimeScale2 = 1; % Replicate Panel 4B, Van De Ville et al., Science Advances
 
%% Time scales of brain fingerprints
if flags.TimeScale1==1
    ICC_avg = zeros(configs.Nparc,n_wsize);
    Good_ICC_nodal = zeros(configs.Nparc,n_wsize);
    Max_ICC_nodal = zeros(configs.Nparc,n_wsize);    
    % Plot ICC across time scales
    ICC_3D = zeros(configs.Nparc,configs.Nparc,n_wsize);
    figure,
    for t=1:n_wsize
        ICC_mat = zeros(configs.Nparc); ICC_mat(configs.mask_ut) = 0.5.*(ICC_2D_LR{t} + ICC_2D_RL{t}); ICC_mat = ICC_mat + ICC_mat';
        ICC_avg(:,t) = sum(ICC_mat)./(configs.Nparc-1);
        Max_ICC_nodal(:,t) = max(ICC_mat,[],2);     
        Good_ICC = ICC_mat.*(ICC_mat>.4); % according to Van De Ville et al., 2021
        ICC_3D(:,:,t) = ICC_mat>.4;
        Good_ICC_nodal(:,t) = sum(Good_ICC); 
        subplot(2,3,t); imagesc(Good_ICC(yeoOrder,yeoOrder),[-.1 .7]); colorbar; axis square; title([int2str(configs.wSizeRange(t).*(configs.TR./1000)) ' s']);
    end 
    yeo_nets = max(yeoROIs);
    yeo_names = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'SUBC'};
    yeo_Color = jet(yeo_nets);
    % Compute ICC value within and between Yeo network
    ICC_yeoMat = nan(yeo_nets,yeo_nets,n_wsize);
    figure,
    for t=1:n_wsize
        ICC_tmp = ICC_3D(:,:,t);              
        for yeo1=1:yeo_nets
            for yeo2=1:yeo_nets
                if yeo1==yeo2
                    mask_wyeo = triu(true(nnz(yeoROIs==yeo1)),1);
                    yeo_subset = ICC_tmp(yeoROIs==yeo1,yeoROIs==yeo2);
                    ICC_yeoMat(yeo1,yeo2,t) = (nnz(yeo_subset(mask_wyeo))./nnz(mask_wyeo)).*100;
                else
                    denominator = nnz(yeoROIs==yeo1).*nnz(yeoROIs==yeo2);
                    yeo_subset = ICC_tmp(yeoROIs==yeo1,yeoROIs==yeo2);                    
                    ICC_yeoMat(yeo1,yeo2,t) = (nnz(yeo_subset(:))./denominator).*100;
                end
            end
        end
        subplot(2,3,t), imagesc(ICC_yeoMat(:,:,t)); colorbar; hold on
        title([int2str(configs.wSizeRange(t).*configs.TR./1000) ' s']);
        xlabel('YeoICC'); ylabel('YeoICC');axis square; 
    end
     
end
%% Time scale 2
if flags.TimeScale2==1
    ICC_3D = zeros(configs.Nparc,configs.Nparc,n_wsize);
    for t=1:n_wsize
        ICC_mat = zeros(configs.Nparc); ICC_mat(configs.mask_ut) = 0.5.*(ICC_2D_LR{t} + ICC_2D_RL{t}); ICC_mat = ICC_mat + ICC_mat';
        ICC_avg(:,t) = sum(ICC_mat)./(configs.Nparc-1);
        Max_ICC_nodal(:,t) = max(ICC_mat,[],2);     
        Good_ICC = ICC_mat.*(ICC_mat>.4); % according to Van De Ville et al., 2021
        ICC_3D(:,:,t) = ICC_mat>.4;
        Good_ICC_nodal(:,t) = sum(Good_ICC>0)./(configs.Nparc-1); 
    end     
    % ICC trajectory per YeoNetwork
    yeo_nets = max(yeoROIs);
    yeo_names = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'SUBC'};
    yeo_Color = jet(yeo_nets);
    figure
    for yy=1:yeo_nets
        subplot(1,2,1);
        aux =mean(ICC_avg(yeoROIs==yy,:),1);
        plot(configs.wSizeRange*(configs.TR./1000),aux,'-o','color',yeo_Color(yy,:),'MarkerFaceColor',yeo_Color(yy,:)); hold on;
        axis square;  
        if yy==yeo_nets
            ylabel('ICC'); xlabel('Time (s)'); legend(yeo_names); title('Avg ICC');
        end
        subplot(1,2,2);
        aux =mean(Good_ICC_nodal(yeoROIs==yy,:),1);
        plot(configs.wSizeRange*(configs.TR./1000),aux,'-o','color',yeo_Color(yy,:),'MarkerFaceColor',yeo_Color(yy,:)); hold on;
        axis square; ylabel('ICC'); xlabel('Time (s)'); title('ICC>0.4')  
    end
end



