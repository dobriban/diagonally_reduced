function [mse_blp, mse_eblp, mse_opt_eblp,mse_eta_eblp] =  compute_mse_denoising_unred_noise(ells,gamma,delta,eta)

if ~exist('delta','var')
    delta = 1;
end

if ~exist('eta','var')
    eta = 0;
end
%cos2= @(ell) max(0,(1-gamma./ell.^2)./(1+gamma./ell));
err_blp = @(ell) ell./(delta^2.*ell+1);
%[lambda,cos2] = standard_spiked_forward(delta*ells,gamma);
[lambda,cos2] = standard_spiked_forward(delta^2*ells,gamma);

err_eblp = @(ell,eta) ell+ eta.^2.*lambda ...
    -2*eta.*delta.*cos2.*(ell+gamma/delta^2);

%BLP with true evecs
mse_blp = err_blp(ells);  

%optimal population BLP shrinkage coeffs
blp_shrink = delta*ells./(delta^2.*ells+1);

%MSE using mis-specified pop shr coef for emp evecs
mse_eblp = err_eblp(ells,blp_shrink);
 
%optimal emp BLP shrinkage coeffs
%opt_blp_shrink = cos2.*(delta.*ells+gamma)./lambda;
opt_blp_shrink = delta.*ells.*cos2./(delta^2.*ells+1);


%MSE using emp shr coef for emp evecs
mse_opt_eblp = err_eblp(ells,opt_blp_shrink);

%MSE using eta as a shr coef for emp evecs
mse_eta_eblp =  err_eblp(ells,eta);
