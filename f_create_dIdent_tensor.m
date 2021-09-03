function dIdent = f_create_dIdent_tensor(FC_test_3D,FC_retest_3D)
%% Create dIdent tensor (Subjects test, Subjects retest, time scales), and store it as a supramarginal matrix (nsubj x nwin)
nwin = size(FC_test_3D,3);
nsubj = size(FC_test_3D,1);
FC_test_2D = zeros(size(FC_test_3D,2),nwin.*nsubj);
FC_retest_2D = zeros(size(FC_test_3D,2),nwin.*nsubj);
count = 1;
for s=1:nsubj
    for t=1:nwin
        aux_test = squeeze(FC_test_3D(s,:,t));
        aux_retest = squeeze(FC_retest_3D(s,:,t));
        FC_test_2D(:,count) = aux_test;
        FC_retest_2D(:,count) = aux_retest; 
        count = count +1;
    end
end
dIdent = corr(FC_test_2D,FC_retest_2D);