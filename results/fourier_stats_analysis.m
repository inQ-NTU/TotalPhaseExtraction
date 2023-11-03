load('fourier_stats_20nK_500samples.mat')
mean_abs_fourier_coeffs_in_20nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_20nk = mean_abs_fourier_coeffs_out;

load('fourier_stats_30nK_500samples.mat')
mean_abs_fourier_coeffs_in_30nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_30nk = mean_abs_fourier_coeffs_out;

load('fourier_stats_40nK_500samples.mat')
mean_abs_fourier_coeffs_in_40nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_40nk = mean_abs_fourier_coeffs_out;

load('fourier_stats_50nK_500samples.mat')
mean_abs_fourier_coeffs_in_50nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_50nk = mean_abs_fourier_coeffs_out;

 load('fourier_stats_60nK_500samples.mat')
mean_abs_fourier_coeffs_in_60nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_60nk = mean_abs_fourier_coeffs_out;

load('fourier_stats_70nK_500samples.mat')
mean_abs_fourier_coeffs_in_70nk = mean_abs_fourier_coeffs_in;
mean_abs_fourier_coeffs_out_70nk = mean_abs_fourier_coeffs_out;

l = 100;
n = 1:15;

fitfun = fittype( @(a,x) a./x.^2);

fitted_curve_in_20nk = fit(n',mean_abs_fourier_coeffs_in_20nk',fitfun);
fitted_curve_out_20nk = fit(n',mean_abs_fourier_coeffs_out_20nk',fitfun);
fitted_curve_in_30nk = fit(n',mean_abs_fourier_coeffs_in_30nk',fitfun);
fitted_curve_out_30nk = fit(n',mean_abs_fourier_coeffs_out_30nk',fitfun);

fitted_curve_in_40nk = fit(n',mean_abs_fourier_coeffs_in_40nk',fitfun);
fitted_curve_out_40nk = fit(n',mean_abs_fourier_coeffs_out_40nk',fitfun);

fitted_curve_in_50nk = fit(n',mean_abs_fourier_coeffs_in_50nk',fitfun);
fitted_curve_out_50nk = fit(n',mean_abs_fourier_coeffs_out_50nk',fitfun);

fitted_curve_in_60nk = fit(n',mean_abs_fourier_coeffs_in_60nk',fitfun);
fitted_curve_out_60nk = fit(n',mean_abs_fourier_coeffs_out_60nk',fitfun);

fitted_curve_in_70nk = fit(n',mean_abs_fourier_coeffs_in_70nk',fitfun);
fitted_curve_out_70nk = fit(n',mean_abs_fourier_coeffs_out_70nk',fitfun);


input_coeffvals = [coeffvalues(fitted_curve_in_20nk), coeffvalues(fitted_curve_in_30nk),...
    coeffvalues(fitted_curve_in_40nk), coeffvalues(fitted_curve_in_50nk), coeffvalues(fitted_curve_in_60nk),...
    coeffvalues(fitted_curve_in_70nk)];

output_coeffvals = [coeffvalues(fitted_curve_out_20nk), coeffvalues(fitted_curve_out_30nk),...
    coeffvalues(fitted_curve_out_40nk), coeffvalues(fitted_curve_out_50nk), coeffvalues(fitted_curve_out_60nk),...
    coeffvalues(fitted_curve_out_70nk)];

temperature_list = [20,30,40,50,60,70];

linear_fitfun = fittype(@(a,b,x) a+b*x);
fitted_linear_in = fit(temperature_list', input_coeffvals',linear_fitfun);
fitted_linear_out = fit(temperature_list', output_coeffvals', linear_fitfun);


n = 1:15;
g = tight_subplot(1,2,[.12 .12],[.2 .1],[.1 .08]);

axes(g(1))
plot(n, fitted_curve_in_20nk(n), '-.', 'Color','black')
hold on
plot(n, fitted_curve_in_20nk(n), 'x','Color','black')
plot(mean_abs_fourier_coeffs_out_20nk, 'o','Color','red')
plot(n, fitted_curve_out_20nk(n), 'Color', 'red')
xlabel('$\rm Mode \; index$', 'Interpreter','latex')
ylabel('$\langle|\Phi_n|^2\rangle$','Interpreter','latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);


f(1) = axes('Position',[.25 .5 .15 .3]);
box on
plot(n(5:15), fitted_curve_in_20nk(n(5:15)), 'Color','black','LineStyle','-.')
hold on
plot(n(5:15), mean_abs_fourier_coeffs_out_20nk(5:15), 'ro')
ax = gca;
ax.YAxis.Exponent = -2;



axes(g(2))
plot(temperature_list, fitted_linear_in(temperature_list), '-.','Color','black')
hold on
plot(temperature_list, input_coeffvals, 'x','Color','Black')
plot(temperature_list, fitted_linear_out(temperature_list), 'Color','red')
plot(temperature_list, output_coeffvals, 'o', 'Color','red')
xlabel('$T_+ \; (\rm nK)$', 'Interpreter', 'latex')
ylabel('$\alpha(T)$', 'Interpreter','latex')
xlim([20,70])
ylim([0.5,2.6])
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

set(g,'FontName','Times', 'FontSize', 16)
set(f,'FontName','Times', 'FontSize', 14)