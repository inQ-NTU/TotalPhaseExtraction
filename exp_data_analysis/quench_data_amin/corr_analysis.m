addpath('../../classes/')
addpath('../../plotting_func/')
load('ext_com_phase_data_uncoupled_batch2.mat')

corr_suite = class_1d_correlation(ext_com_phase_all_shots);
fourier_corr = diag(corr_suite.fourier_correlation());
g1_corr = corr_suite.g1_corr();
cov_matrix = corr_suite.covariance_matrix();