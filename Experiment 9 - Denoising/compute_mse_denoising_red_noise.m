function [mse_blp, mse_eblp, mse_opt_eblp,mse_eta_eblp] =  compute_mse_denoising_red_noise(ells,gamma,delta,eta)

if ~exist('delta','var')
    delta = 1;
end

if ~exist('eta','var')
    eta = 0;
end
%cos2= @(ell) max(0,(1-gamma./ell.^2)./(1+gamma./ell));
err_blp = @(ell) ell./(delta.*ell+1);
[lambda,cos2] = standard_spiked_forward(delta*ells,gamma);
beta = 1 + gamma./(delta.*ells);
beta(ells==0)=1;

err_eblp = @(ell,eta) ell+ eta.^2.*delta.*lambda ...
    -2*delta*eta.*ell.*cos2.*beta;

%BLP with true evecs
mse_blp = err_blp(ells); 

%optimal population BLP shrinkage coeffs
blp_shrink = ells./(delta.*ells+1);

%MSE using mis-specified pop shr coef for emp evecs
mse_eblp = err_eblp(ells,blp_shrink);
 
%optimal emp BLP shrinkage coeffs
opt_blp_shrink = ells.*cos2.*beta./lambda;

%MSE using emp shr coef for emp evecs
mse_opt_eblp = err_eblp(ells,opt_blp_shrink);

%MSE using eta as a shr coef for emp evecs
mse_eta_eblp =  err_eblp(ells,eta);