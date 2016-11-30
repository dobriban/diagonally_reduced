%Plot shrinkage and MSE for in-sample vs out-of-sample denoising
%% compute shrinkage
gamma = 1/2;
delta = 1/2;
ells = linspace(0,10*sqrt(gamma),100)'; %spikes, BBP at sqrt(gamma)
[lambda,cos2] = standard_spiked_forward(ells*delta,gamma);
blp_shr = ells./(1+ells.*delta);
beta = 1+gamma./(delta.*ells);
eblp_shr = ells.*cos2.*beta./lambda;
eblp_shr_oos = ells.*cos2./(1+delta*ells.*cos2);
%% plot
a = {'-','--','-.',':'};    savefigs =1;
rng(2);
figure, hold on
h1 = plot(ells,blp_shr ,'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ells,eblp_shr ,'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{3});
h3 = plot(ells,eblp_shr_oos ,'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{2});
legend([h1,h2,h3],{'BLP','Opt EBLP','Opt EBLP-OOS'},'location','Best')
set(gca,'fontsize',20)
xlim([min(ells),max(ells)])
xlabel('Pop Spike')
ylabel('\eta^*')
%%
if savefigs==1
    filename = sprintf( './shr_blp_eblp_oos_gamma=%.2f_delta=%.2f.png',gamma,delta);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
end

%% compute MSE
blp_mse = ells./(1+delta*ells);
eblp_mse = ells-delta*(ells.*cos2.*beta).^2./lambda;
eblp_mse_oos = ells.*(1+delta*ells.*cos2.*(1-cos2))./(1+delta*ells.*cos2);
%% plot
a = {'-','--','-.',':'};    savefigs =1;
rng(2);
figure, hold on
h1 = plot(ells,blp_mse ,'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ells,eblp_mse ,'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{3});
h3 = plot(ells,eblp_mse_oos ,'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{2});
legend([h1,h2,h3],{'BLP','Opt EBLP','Opt EBLP-OOS'},'location','Best')
set(gca,'fontsize',20)
xlim([min(ells),max(ells)])
xlabel('Pop Spike')
ylabel('MSE')
%%
if savefigs==1
    filename = sprintf( './mse_blp_eblp_oos_gamma=%.2f_delta=%.2f.png',gamma,delta);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
end