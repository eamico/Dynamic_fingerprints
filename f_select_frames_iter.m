function [Row_index,Col_index] = f_select_frames_iter(Ident_mat)
%% This function select best dynamic functional connectome frames for brain fingerprints 
%% across temporal scales, when ranking dFC frames based on the individual dIself score
%% (see also Van De Ville et al., Science Advances 2021)
% note that, for the frame selection based on dIself
% other methods other than the mean can be explored (i.e. max,median, etc.)
Ident_Rows= nanmean(Ident_mat,2);
Ident_Cols = nanmean(Ident_mat);
[~,Col_index] = sort(Ident_Cols,'descend');
[~,Row_index] = sort(Ident_Rows,'descend');

            
