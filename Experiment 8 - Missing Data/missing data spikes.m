%Test formulas for the spike and the angles
%When there is missing data
%% Toeplitz
gamma = 1/2;
n=2*1e2;
p = floor(n*gamma);
rho = 0.5;
r = rho.^(0:1:p-1);
Sigma = toeplitz(r);
t = eig(Sigma);
w = ones(p,1)/p;
ell_arr = linspace(0,15*sqrt(gamma),50); %spikes
delta_arr = [1/3;2/3;1];

lambda_num = zeros(length(ell_arr),length(delta_arr));  lambda_MC = lambda_num;
c_right_num = zeros(length(ell_arr),length(delta_arr)); c_right_MC = c_right_num;
c_left_num = zeros(length(ell_arr),length(delta_arr)); c_left_MC = c_left_num;
nMonte = 10;

%% compute and finite sample simulation
x = tic;
print_iter=1;
for i=1:length(ell_arr)
    if print_iter==1
        str = sprintf('spike %d out of %d. time: %.2f sec \n',i,length(ell_arr), toc(x));
        fprintf(str);
    end
    for j=1:length(delta_arr)
        [lambda_num(i,j),c_right_num(i,j),c_left_num(i,j)] = general_spiked_forward(delta_arr(j)^2*ell_arr(i), delta_arr(j)*t, w, gamma);
        lambda =  zeros(nMonte,1); c_right = lambda; c_left  = lambda;
        for l=1:nMonte
            u = randn(n,1); u = u/norm(u);
            v = randn(p,1); v = v/norm(v);
            X_0=randn(n,p);
            D = binornd(1,delta_arr(j),[n,p]);
            X = ell_arr(i)^(1/2)*u*v'+n^(-1/2)*X_0*diag(t.^(1/2));
            X = D.*X;
            [U, S, V] = svd(X);
            lambda(l) = S(1)^2;
            c_left(l) = (U(:,1)'*u)^2;
            c_right(l) = (V(:,1)'*v)^2;
        end
        lambda_MC(i,j) = mean(lambda);
        c_right_MC(i,j) = mean(c_right);
        c_left_MC(i,j) = mean(c_left);
    end
end

%% plot comparison: spike
a = {'-','--','-.',':'};       plotfigs=1;  savefigs =1;
rng(2);
figure, hold on
h1 = plot(ell_arr,lambda_num(:,3),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ell_arr,lambda_MC(:,3),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(ell_arr,lambda_num(:,2),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(ell_arr,lambda_MC(:,2),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
h5 = plot(ell_arr,lambda_num(:,1),'linewidth',4,'color',rand(1,3));
set(h5,'LineStyle',a{1});
h6 = plot(ell_arr,lambda_MC(:,1),'linewidth',4,'color',rand(1,3));
set(h6,'LineStyle',a{2});
xlabel('$\ell$','Interpreter','latex')
ylabel('$\lambda$','Interpreter','latex')
legend([h1,h2,h3,h4,h5,h6],{'\delta=1','\delta=1 MC','\delta=2/3','\delta=2/3 MC','\delta=1/3','\delta=1/3 MC'},'location','Best')
xlim([min(ell_arr),max(ell_arr)])
set(gca,'fontsize',20)

%%
if savefigs==1
    filename = sprintf( './projected_spike_toeplitz_gamma=%.2f_rho=%.2f_n=%d_nMC=%d.png',gamma,rho,n,nMonte);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end

%% angles
rng(2);
figure, hold on
h1 = plot(ell_arr,c_right_num(:,3),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ell_arr,c_right_MC(:,3),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(ell_arr,c_right_num(:,2),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(ell_arr,c_right_MC(:,2),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
h5 = plot(ell_arr,c_right_num(:,1),'linewidth',4,'color',rand(1,3));
set(h5,'LineStyle',a{1});
h6 = plot(ell_arr,c_right_MC(:,1),'linewidth',4,'color',rand(1,3));
set(h6,'LineStyle',a{2});
xlabel('$\ell$','Interpreter','latex')
ylabel('$c^2$','Interpreter','latex')
legend([h1,h2,h3,h4,h5,h6],{'\delta=1','\delta=1 MC','\delta=2/3','\delta=2/3 MC','\delta=1/3','\delta=1/3 MC'},'location','Best')
xlim([min(ell_arr),max(ell_arr)])
set(gca,'fontsize',20)
%%
if savefigs==1
    filename = sprintf( './empirical_cosine_comparison_toeplitz_gamma= %.2f_rho=%.2f_n=%d_nMC=%d.png',gamma,rho,n,nMonte);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end
