%Test formulas for denoising MSE in a MC simulation
% Unreduced noise model
%% parameters
gamma = 1/2;
ells = linspace(0,20*sqrt(gamma),40)'; %spikes
n=1e3;
p = floor(n*gamma);
nMonte = 50;
delta = 1/2;
eta = 1;
mse_blp_mc =zeros(length(ells),1);
mse_eblp_mc = zeros(length(ells),1);
mse_opt_eblp_mc= zeros(length(ells),1);
mse_eta = zeros(length(ells),1);

%% finite sample simulation
x = tic;
print_iter=1;
for i=1:length(ells)
    if print_iter==1
        str = sprintf('spike %d out of %d. time: %.2f sec \n',i,length(ells), toc(x));
        fprintf(str);
    end
    blp =  zeros(nMonte,1); eblp =  zeros(nMonte,1);
    opt_eblp =  zeros(nMonte,1); eta_err =  zeros(nMonte,1);
    for l=1:nMonte
        Z = randn(n,1); Z = Z/norm(Z);
        u = randn(p,1); u = u/norm(u);
        X_0=randn(n,p); %nxp 
        Signal = ells(i)^(1/2)*Z*u';
        D = binornd(1,delta,[n,p]);
        X = D.*Signal +n^(-1/2)*X_0;
        [Z_hat, S, U_hat] = svd(X);
                
        Xhat_blp = ells(i)*delta/(ells(i)*delta^2+1)*(u*u')*(X(1,:)');
        Xhat_eblp = eta*(U_hat(:,1)*U_hat(:,1)')*(X(1,:)');
        
        [~,c2] = standard_spiked_forward(ells(i)*delta,gamma);
        opt_blp_shrink = delta.*ells(i).*c2./(delta^2.*ells(i)+1);
        Xhat_opt_eblp = opt_blp_shrink*(U_hat(:,1)*U_hat(:,1)')*(X(1,:)');
        
        blp(l) = n*norm(Xhat_blp -Signal(1,:)')^2;
        eblp(l) = n*norm(Xhat_eblp -Signal(1,:)')^2;
        opt_eblp(l) = n*norm(Xhat_opt_eblp -Signal(1,:)')^2;
    end
    mse_blp_mc(i) = mean(blp);
    mse_eblp_mc(i) = mean(eblp);
    mse_opt_eblp_mc(i) = mean(opt_eblp);
end

%% compute MSE theoretically
[mse_blp, ~, mse_opt_eblp,mse_eta_eblp] =  compute_mse_denoising_unred_noise(ells,gamma,delta,eta);

% plot MSE
a = {'-','--','-.',':'};    savefigs =1;
rng(2);
figure, hold on
h1 = plot(ells,mse_blp,'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ells,mse_blp_mc,'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(ells,mse_eta_eblp,'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(ells,mse_eblp_mc,'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
h5 = plot(ells,mse_opt_eblp,'linewidth',4,'color',rand(1,3));
set(h5,'LineStyle',a{1});
h6 = plot(ells,mse_opt_eblp_mc,'linewidth',4,'color',rand(1,3));
set(h6,'LineStyle',a{2});
xlabel('Pop Spike')
ylabel('MSE')
legend([h1,h2,h3,h4,h5,h6],{'BLP','BLP MC','\eta BLP','\eta BLP MC','OEBLP','OEBLP MC'},'location','Best')
%legend([h1,h2,h3,h4],{'BLP','BLP MC','EBLP','EBLP MC'},'location','Best')
set(gca,'fontsize',20)
xlim([min(ells),max(ells)])


%%
if savefigs==1
    filename = sprintf( './denoising_MSE_spiked_unred_noise_delta=%.2f_gamma= %.2f_nMonte=%d.png',delta,gamma,nMonte);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
end