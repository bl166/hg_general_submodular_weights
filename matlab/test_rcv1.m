function [] = test_rcv1(para_p, delta, random_seed)

addpath('subfunc');
load('data/rcv1/data.mat'); % R: |E|x|V|, y: 1x|V|
% load('data/rcv1/data_nw+50_sb1+0.020000_sb2+0.000400.mat')
% load('data/rcv1/data_nw+30_sb1+0.020000_sb2+0.001000.mat')
load('data/rcv1/data_nw+30_sb1+0.010000_sb2+0.001000.mat')
% random_seed = 0;
rng(random_seed);

%%
mode_str = 'ct';
mode_num = 3;
dec_outloop_ = 1e-3;
err_inloop_ = 1e-3;
% dec_outloop = 1e-3;
% err_inloop = 1e-6;
maxiter_outer_ = 100;
maxiter_inner_ = 10000;

folder = ['data/rcv1/fin526_rng', num2str(random_seed), '_ww_out', num2str(-log10(dec_outloop_)),...
        '_in', num2str(-log10(err_inloop_)), '_maxin_', num2str(maxiter_inner_), '_maxout_', num2str(maxiter_outer_),'/'];


file = strcat(folder, mode_str, '_p', string(para_p), '_delta', string(delta), '/');
mkdir(file)
file_log = strcat(file, 'log.txt');
file_log_id = fopen(file_log, 'w');
fclose(file_log_id);
diary(file_log);

[n_e, n_v, card, kappa, R, incidence_list, parameter_list, mu] = hg_para(R, para_p, mode_str, delta, 'w', 'text');

%% rw-based clique
f_rwc = hg_rw_laplacian('c');
tic
[~, eigvec_rwc, eigval_rwc] = f_rwc(R, 'std', 2);
toc
vmin_rwc = eigvec_rwc(:, 2)';
fprintf('%f %f\n', eigval_rwc);

[labels_rwc, NCut_rwc] = general_optthreshold(incidence_list, parameter_list, mu, n_v, n_e, mode_num, delta, vmin_rwc);
err_rwc = comp_err(labels_rwc, y);

%% submodular
tic

maxiter_inner = maxiter_inner_;
maxiter_outer = maxiter_outer_;
dec_outloop = dec_outloop_;
err_inloop = err_inloop_;
warmstart = vmin_rwc;

[labels, NCut, vmin] = submodular_hypergraph_partition(incidence_list, parameter_list, mu, ...
    n_v, n_e, mode_num, delta, dec_outloop, err_inloop, warmstart, maxiter_inner, maxiter_outer);
toc

err = comp_err(labels, y);

fprintf('baseline %f %f\n', NCut_rwc, err_rwc);
fprintf('proposed %f %f\n', NCut, err);

diary off
%% save results
save(strcat(file, 'vmin_rwc.mat'), 'vmin_rwc');
save(strcat(file, 'labels_rwc.mat'), 'labels_rwc');
save(strcat(file, 'vmin.mat'), 'vmin');
save(strcat(file, 'labels.mat'), 'labels');

save(strcat(file, 'NCut_rwc.mat'), 'NCut_rwc');
save(strcat(file, 'err_rwc.mat'), 'err_rwc');
save(strcat(file, 'NCut.mat'), 'NCut');
save(strcat(file, 'err.mat'), 'err');

end
