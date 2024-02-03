load('many_shots_10nk_1000samples_n20.mat')
fidelity_stats_10nk = mitigated_fidelity;
mean_fid_10nk = mean(fidelity_stats_10nk);
std_fid_10nk = std(fidelity_stats_10nk);
cov_in_10nk = cov_in;
cov_out_10nk = cov_out;
ext_com_phase_10nk = mitigated_ext_com_phase;
com_phase_10nk = com_phase; 

load('many_shots_20nk_1000samples_n20.mat')
fidelity_stats_20nk = mitigated_fidelity;
mean_fid_20nk = mean(fidelity_stats_20nk);
std_fid_20nk = std(fidelity_stats_20nk);
cov_in_20nk = cov_in;
cov_out_20nk = cov_out;
ext_com_phase_20nk = mitigated_ext_com_phase;
com_phase_20nk = com_phase; 

load('many_shots_30nk_1000samples_n20.mat')
fidelity_stats_30nk = mitigated_fidelity;
mean_fid_30nk = mean(fidelity_stats_30nk);
std_fid_30nk = std(fidelity_stats_30nk);
cov_in_30nk = cov_in;
cov_out_30nk = cov_out;
ext_com_phase_30nk = mitigated_ext_com_phase;
com_phase_30nk = com_phase; 

load('many_shots_40nk_1000samples_n20.mat')
fidelity_stats_40nk = mitigated_fidelity;
mean_fid_40nk = mean(fidelity_stats_40nk);
std_fid_40nk = std(fidelity_stats_40nk);
cov_in_40nk = cov_in;
cov_out_40nk = cov_out;
ext_com_phase_40nk = mitigated_ext_com_phase;
com_phase_40nk = com_phase; 

load('many_shots_50nk_1000samples_n20.mat')
fidelity_stats_50nk = mitigated_fidelity;
mean_fid_50nk = mean(fidelity_stats_50nk);
std_fid_50nk = std(fidelity_stats_50nk);
cov_in_50nk = cov_in;
cov_out_50nk = cov_out;
ext_com_phase_50nk = mitigated_ext_com_phase;
com_phase_50nk = com_phase; 

%average and standard deviation of fidelity statistics
temperature = [10,20,30,40,50];
mean_fidelities = [mean_fid_10nk, mean_fid_20nk, mean_fid_30nk, mean_fid_40nk, mean_fid_50nk];
std_fidelities = [std_fid_10nk, std_fid_20nk, std_fid_30nk, std_fid_40nk, std_fid_50nk]; 

%define z grid
z_grid = linspace(-50,50,100);

%Plotting
figure
g = tight_subplot(2,2,[.12 .12],[.15 .1],[.1 .1]);

axes(g(1))
%random = floor(rand()*1000);
plot(z_grid, com_phase_20nk(614,:),'Color','Black', 'LineStyle', '-.')
hold on
plot(z_grid, ext_com_phase_20nk(614,:), 'Color', 'red')
ylabel('$\phi_+ (z)$', 'Interpreter','latex')
xlabel('$z\; (\mu m)$', 'Interpreter','latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);


axes(g(2))
histogram(fidelity_stats_30nk, 'Normalization','probability','FaceColor', [0.8500 0.3250 0.0980])
ylabel('$P(F)$','Interpreter','latex')
xlabel('$F$','Interpreter','latex')
xlim([0.8,1])
xticks([0.8,0.9,1])
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);

axes(g(3))
errorbar(temperature, mean_fidelities, std_fidelities, 'or--')
ylabel('$\langle F \rangle$', 'Interpreter','latex')
xlabel('$T_+ \; (\rm nK)$','Interpreter','latex')
ylim([0.8,1.02])
yticks([0.8,0.9,1])
xlim([8,52])
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);


axes(g(4))
plot(z_grid, cov_in_50nk(50,:),'Color','black','LineStyle','-.')
hold on
plot(z_grid, cov_out_50nk(50,:), 'r-')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);
ylb = ylabel('$g_2(z,0)$','Interpreter','latex');    


set(g, 'FontName', 'Times', 'FontSize', 16)

figure
f = tight_subplot(1,2,[.12 .12],[.2 .15],[.1 .15]);
axes(f(1))
imagesc(z_grid, z_grid, cov_in_50nk)
xlabel('$z\; (\mu m)$', 'Interpreter','latex')
ylabel('$z^\prime\; (\mu m)$', 'Interpreter','latex')
clim([-0.5,1])
axes(f(2))
imagesc(z_grid, z_grid, cov_out_50nk)
cb = colorbar;
cb.Position = cb.Position + [0.1,0,0,0];
xlabel('$z\; (\mu m)$', 'Interpreter','latex')
clim([-0.5,1])
set(f, 'FontName', 'Times', 'FontSize', 16)
colormap(gge_colormap)




