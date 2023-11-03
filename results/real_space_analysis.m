load('real_space_reconstruction_10nk_500samples.mat')
raw_fidelity_10nk = raw_fidelity;
mitigated_fidelity_10nk = mitigated_fidelity;
cov_in_10nk = cov_in;
cov_out_10nk = cov_out;

load('real_space_reconstruction_50nk_500samples.mat')
raw_fidelity_50nk = raw_fidelity;
mitigated_fidelity_50nk = mitigated_fidelity;
cov_in_50nk = cov_in;
cov_out_50nk = cov_out;

load('real_space_reconstruction_20nk_500samples.mat')
raw_fidelity_20nk = raw_fidelity;
mitigated_fidelity_20nk = mitigated_fidelity;

load('real_space_reconstruction_30nk_500samples.mat')
raw_fidelity_30nk = raw_fidelity;
mitigated_fidelity_30nk = mitigated_fidelity;

load('real_space_reconstruction_40nk_500samples.mat')
raw_fidelity_40nk = raw_fidelity;
mitigated_fidelity_40nk = mitigated_fidelity;

temperature_list = [10,20,30,40,50];
mean_mitigated_fidelity = [mean(mitigated_fidelity_10nk), mean(mitigated_fidelity_20nk),...
    mean(mitigated_fidelity_30nk), mean(mitigated_fidelity_40nk), mean(mitigated_fidelity_50nk)];
mean_raw_fidelity = [mean(raw_fidelity_10nk), mean(raw_fidelity_20nk), mean(raw_fidelity_30nk),...
    mean(raw_fidelity_40nk), mean(raw_fidelity_50nk)];
std_mitigated_fidelity = [std(mitigated_fidelity_10nk), std(mitigated_fidelity_20nk),...
    std(mitigated_fidelity_30nk), std(mitigated_fidelity_40nk), std(mitigated_fidelity_50nk)];
std_raw_fidelity = [std(raw_fidelity_10nk), std(raw_fidelity_20nk), std(raw_fidelity_30nk),...
    std(raw_fidelity_40nk), std(raw_fidelity_50nk)];


%%%%Plotting%%%%
figure
g = tight_subplot(2,2,[.12 .1],[.1 .05],[.1 .05]);

axes(g(1))
histogram(raw_fidelity_10nk, 'Normalization','probability')
hold on
histogram(mitigated_fidelity_10nk, 'Normalization','probability')
ylabel('$P(F)$','Interpreter','latex')
xlabel('$F$','Interpreter','latex')
xlim([0.6,1])
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(g(2))
errorbar(temperature_list, mean_raw_fidelity, std_raw_fidelity, '--o')
hold on
errorbar(temperature_list, mean_mitigated_fidelity, std_mitigated_fidelity, '.-', 'MarkerSize',20);
ylabel('$\langle F \rangle$', 'Interpreter','latex')
xlabel('$T_+ \; (\rm nK)$','Interpreter','latex')
ylim([0.4,1.05])
xlim([8,52])
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);


z_grid = linspace(-50,50,100);
axes(g(3))
plot(z_grid, cov_in_10nk(50,:), '-.','Color','black')
hold on
plot(z_grid, cov_out_10nk(50,:),'r-')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylb = ylabel('$g^{(2)}$','Interpreter','latex');    
ylb.Position = ylb.Position - [0.1,0,0];
ylim([-0.25,0.5])
yticks([-0.25,0,0.25,0.5])
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(g(4))
plot(z_grid, cov_in_50nk(50,:),'-.','Color','black')
hold on
plot(z_grid, cov_out_50nk(50,:), 'r-')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);
ylb = ylabel('$g^{(2)}$','Interpreter','latex');    
ylb.Position = ylb.Position - [0.1,0,0];


set(g, 'FontName', 'Times', 'FontSize', 16)
