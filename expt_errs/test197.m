        function main
        startnow;
        rand(22,116);
        rand(22,116);

tic
        expt_errs
toc
        stopnow
        end
%
%
%
%
%
        function expt_errs
%
        ndeltas = 10
        ngams = 10
        nmonte = 200;
%
        gams = linspace(.1,1,ngams);
        gams = fliplr(gams)
%
        deltas = linspace(.1,1,ndeltas)
%
        m=1200;

        errs_hat_fr_table = zeros(ngams,ndeltas);
        errs_pred_fr_table = zeros(ngams,ndeltas);
%
        errs_hat_op_table = zeros(ngams,ndeltas);
        errs_pred_op_table = zeros(ngams,ndeltas);
%
        sdev_op_table = zeros(ngams,ndeltas);
        sdev_fr_table = zeros(ngams,ndeltas);

        errs_fr_cube = zeros(ngams,ndeltas,nmonte);
        errs_hat_fr_cube = zeros(ngams,ndeltas,nmonte);
%
        errs_op_cube = zeros(ngams,ndeltas,nmonte);
        errs_hat_op_cube = zeros(ngams,ndeltas,nmonte);
%
        for igam = 1:ngams
%
        gam = gams(igam);
        n=floor(m/gam);

        for idelta = 1:ndeltas
%
        delta = deltas(idelta);
%
        [errs_fr,errs_hat_fr,err_pred_fr,errs_op,errs_hat_op,...
           err_pred_op] = expt_cov_err(gam,m,n,delta,nmonte);

        errs_fr_cube(igam,idelta,1:nmonte) = errs_fr;
        errs_hat_fr_cube(igam,idelta,1:nmonte) = errs_hat_fr;
        errs_pred_fr_table(igam,idelta) = err_pred_fr;

        errs_op_cube(igam,idelta,1:nmonte) = errs_op;
        errs_hat_op_cube(igam,idelta,1:nmonte) = errs_hat_op;
        errs_pred_op_table(igam,idelta) = err_pred_op;

        save data11.mat

    end

    end

%
%   compute average errors
%

        errs_fr_table = mean(errs_fr_cube,3)
        errs_hat_fr_table = mean(errs_hat_fr_cube,3)
        errs_pred_fr_table
%
        errs_op_table = mean(errs_op_cube,3)
        errs_hat_op_table = mean(errs_hat_op_cube,3)
        errs_pred_op_table
%
%   compute standard deviations of empirical errors
%
        for igam=1:ngams
%
        for idelta=1:ndeltas
%
        var_op = var(errs_op_cube(igam,idelta,1:nmonte));
        sdev_op_table(igam,idelta) = sqrt(var_op)

        var_fr = var(errs_fr_cube(igam,idelta,1:nmonte));
        sdev_fr_table(igam,idelta) = sqrt(var_fr)
    end

    end

        save data11.mat

        end
%
%
%
%
%
        function [errs_fr,errs_hat_fr,err_pred_fr,errs_op,errs_hat_op,...
           err_pred_op] = expt_cov_err(gam,m,n,delta,nmonte)
%
        errs_op = zeros(1,nmonte); errs_hat_op = zeros(1,nmonte);
        errs_fr = zeros(1,nmonte); errs_hat_fr = zeros(1,nmonte);
%
        k=1;
%
        evmin = 1 + (sqrt(gam) + 2)/delta;
        evs_noisy = evmin + [k-1:-1:0];
        ells = evs_noisy;
%
        [err_pred_op,ierr_op] = cov_err_oper(ells,k,delta,gam);
        [err_pred_fr,ierr_fr] = cov_err_frob(ells,k,delta,gam);

        for ijk = 1:nmonte
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);
%
        [xhat_op,sxhat_op,err_hat_op,xhat_fr,sxhat_fr,...
           err_hat_fr] = cov_shrink_mat2(ys,inds,m,n,delta,k,'d');
%
%   operator norm loss
%
        errs_op(ijk) = norm(xhat_op - rcov);
        errs_hat_op(ijk) = err_hat_op;
%
%   frobenius norm loss
%
        errs_fr(ijk) = norm(xhat_fr - rcov,'fro')^2;
        errs_hat_fr(ijk) = err_hat_fr;

    end
        end
%
%
%
%
%
        function [xhat_op,sxhat_op,err_pred_op,xhat_fr,sxhat_fr,...
           err_pred_fr] = cov_shrink_mat2(ys,inds,m,n,delta,k_est,syst)
        sxhat_op = zeros(1,k_est); ells=zeros(1,k_est); 
        sxhat_fr = zeros(1,k_est);
%
%
%                             description:
%
%   This function computes a shrinkage estimator of the covariance of data
%   with missing values, where each coordinate is sampled with probability
%   delta (delta is given by user, and hence may be an estimate). The 
%   estimate is obtained by taking the sample covariance of the 
%   missing data ys, applying a linear operator to get an unbiased
%   estimator of the covariance of the xs, and then adjusting the spectrum 
%   of this unbiased estimator.
%
%   The shrinker is supposed to be asymptotically optimal for the specified 
%   loss. It estimates the signal covariance, NOT the signal + noise 
%   covariance (which can be estimated by adding the identity matrix).
%
%   The mean of the distribution is assumed to be zero (the code does not 
%   perform mean subtraction). The noise is assumed to be white with
%   variance 1.
%
%
%                           input parameters:
%
%   ys - the m-by-n matrix of data with zeros in the missing entries
%
%   inds - the m-by-n matrix with 1's in observed entries, 0's elsewhere
%
%   m - the dimensionality of each data vector
%
%   n - the number of samples
%
%   delta - the (perhaps estimated) sampling probability
%
%   k_est - an a prior estimate of the rank (i.e. the maximum allowed rank
%      for the rank)
%
%   loss - a character that specifies which loss function to minimize
%      'o' - operator norm loss
%      'f' - frobenius norm loss
%
%   syst - a character that specifies which unbiased estimator of x to use
%      (that is, which linear operator to apply to the covariance of ys)
%      'd' - the delta system
%      'n' - the nijs system (i.e. available case estimator)
%
%
%                           output parameters:
%
%   xhat2 - the m-by-m estimated covariance matrix, obtained by shrinking
%      the unbiased estimator using the optimal operator norm shrinker
%
%   uxhat - the m-by-m matrix whose columns are the eigenvectors of xhat2
%
%   sxhat2 - the m-dimensional vector of eigenvalues of xhat2
%

%
%   . . . take sample covariance of the Y vectors
%
        yhat = ys*ys' / n;
%
%   invert system to get estimate of xs covariance
%
        if (syst == 'd')
            xhat = yhat2xhat_delta(yhat,delta,m);
        end
%
        if (syst == 'n')
            nijs = inds*inds';
            xhat = yhat2xhat_nijs(yhat,nijs,m,n);
        end

        [uxhat,sxhat] = eigsmart(xhat,m,k_est);
%
%   apply optimal non-linearity to each eigenvalue
%
        gam=m/n;

%
%   operator norm loss
%
        for i=1:k_est
%
        sxhat_op(i) = cov_emp2oper(sxhat(i),delta,gam);
        ells(i) = sxhat_op(i);
        sxhat_op(i) = sxhat_op(i) - 1;
    end

        [err_pred_op,ierr_op]=cov_err_oper(ells,k_est,delta,gam);
        xhat_op = uxhat * diag(sxhat_op) * uxhat';
%
%
%   frobenius norm loss
%

        for i=1:k_est
%
        [sxhat_fr(i),ells(i)] = cov_emp2frob(sxhat(i),delta,gam);
        sxhat_fr(i) = sxhat_fr(i) - 1;
    end
        [err_pred_fr,ierr_fr]= cov_err_frob(ells,k_est,delta,gam);
        xhat_fr = uxhat * diag(sxhat_fr) * uxhat';

        end
%
%
%
%
%
        function [err_pred,ierr] = cov_err_oper(ells,k,delta,gam)
%
        ell = max(ells);

        ierr=0
        [sx1,sx_bulk,cc,sy1,sy_bulk] = cov_pop2emp2(ell,delta,gam);

        if (ell < 1+sqrt(gam)/delta)
           cc = 0;
           ierr=1;
        end

        err_pred =  (ell-1) * sqrt(1 - cc^2)

        end
%
%
%
%
%
        function [err_pred,ierr] = cov_err_frob(ells,k,delta,gam)
%
%   ells - the eigenvalue of signal plus noise
%

        err_pred = 0;
        for i=1:k
%
        ierr=0;
        [sx1,sx_bulk,cc,sy1,sy_bulk] = cov_pop2emp2(ells(i),delta,gam);
     
        if (ells(i) < 1 + sqrt(gam)/delta)
           cc = 0;
           ierr=1;
        end

        err_pred = err_pred + (1-cc^4)*(ells(i)-1)^2;
    end

        end
%
%
%
%
%
        function test_finite_angles()
%
        gam = .8
        m=500
        n=floor(m/gam)
        k=1
        delta = 1
%
        evmin = 1 + (sqrt(gam) + 1)/delta + 5
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,ell] = ...
           prep_gauss(m,n,k,delta,evs_noisy);
%
        [ys_hat,uy,vy,sy_sh,sy] = svd_shrink_mat(ys,m,n,delta,k,'f');

        [uy,sy,vy] = svdsmart(ys,m,n,m);
        chk0 = norm(uy*diag(sy)*vy' - ys,'fro')


        end
%
%
%
%
%
        function [x_filt,eta,err_hat] = blp_cheat_out(ys,m,n,gam,delta,...
           ell,cosl,uy,y_out)
%
        [eta,err_hat] = eta_fmla_out(delta,ell,cosl,gam)
        prod = sum(uy .* y_out)
        coeff = eta*prod
%
        x_filt = coeff*uy;

        end
%
%
%
%
%
        function [eta,err_hat] = eta_fmla_out(delta,ell,cosl,gam)
%
        sinl = sqrt(1-cosl^2);
%%%        chk0 = cosl^2 + sinl^2 - 1

%
        eta = ell*cosl^2 / (delta*ell*cosl^2 + 1)

        err_hat = ell * (delta*ell*cosl^2*sinl^2 + 1) ...
           / (delta*ell*cosl^2 + 1);

        end
%
%
%
%
%
        function [x_filt,err_hat,eta] = blp_cheat_in(ys,m,n,gam,delta,...
           ell,cosl,uy)
%
        x_filt = zeros(m,n);
%
        [eta,err_hat] = eta_fmla_in(delta,ell,cosl,gam);

        for i=1:n
%
        y_i = ys(1:m,i);
        prod = sum(y_i .* uy);
        coeff = prod * eta;

        x_filt(1:m,i) = coeff * uy;

    end


        end
%
%
%
%
%
        function [eta,err_hat] = eta_fmla_in(delta,ell,cosl,gam)
%
        s_pop = sqrt(ell);
        [s_emp,cosl,cosr] = svd_pop2emp_delta(s_pop,gam,delta);
        tsq = s_emp^2;

        chk0 = (1+delta*ell)*(1 + gam/(delta*ell)) - tsq

        bet = 1 + gam/(delta*ell);
        eta = ell * cosl^2 / (delta*ell + 1);
        err_hat = ell - delta*(ell*cosl^2*bet)^2 / tsq;


%%%        err_hat = ell - (ell * cosl^2 * bet)^2 * delta / tsq
        end
%
%
%
%
%
        function [x_filt,err_hat,eta] = blp_in(ys,m,n,gam,delta)
%
        x_filt = zeros(m,n);
%
        [uy,sy,vy] = svdsmart(ys,m,n,1);
        sy = sy / sqrt(n * delta)
        s_pop = svd_emp2pop_new(sy,gam)
        s_pop = s_pop / sqrt(delta)
%
        [s_emp,cosl,cosr] = svd_pop2emp_delta(s_pop,gam,delta)
        chk0 = sy-s_emp
%
        ell = s_pop^2;

        [x_filt,err_hat,eta] = blp_cheat_in(ys,m,n,gam,delta,ell,cosl,uy);


        end
%
%
%
%
%
        function [xfilt,err_hat,eta] = blp_oracle(ys,m,n,gam,delta,u,ell)
%
        xfilt = zeros(m,n);
        svpop = sqrt(ell);
%
        [svemp2,cosl,cosr] = svd_pop2emp_delta(svpop,gam,delta);
%
        [uy,svemp,vy] = svdsmart(ys,m,n,1);
        svemp=svemp/sqrt(n*delta);
%
        cosl2 = abs(sum(u .* uy));

        eta = ell / (delta*ell + 1);

        for i=1:n
%
        y_i = ys(1:m,i);
        prod = sum(y_i .* u);
        coeff = prod * eta;

        xfilt(1:m,i) = coeff * u;
    end

        err_hat = eta;

        end
%
%
%
%
%
        function x_hat = js_het2_many(ys,inds,m,n,sig)
%
%   ys - subsampled Gaussian data; each sampled entry has standard
%      deviation sig
%
        [x,ds,kcounts] = many2one(ys,inds,m,n,sig);
        [x_hat,a_hats,bs,dstars] = js_het2(x,ds,m)


        end
%
%
%
%
%
        function x_hat = js_het3_many(ys,inds,m,n,sig)
%
%   ys - subsampled Gaussian data; each sampled entry has standard
%      deviation sig
%
        [x,ds,kcounts] = many2one(ys,inds,m,n,sig);
        [x_hat,a_hat,err_hats,bs] = js_het3(x,ds,m)


        end
%
%
%
%
%
        function [x,ds,kcounts] = many2one(ys,inds,m,n,sig)
%
%   ys - subsampled Gaussian data; each sampled entry has standard
%      deviation sig
%
%
%   code returns the sample mean and effective variance D of each entry
%
        ds = zeros(m,1);
        kcounts = zeros(m,1);
        x = mean_estim(ys,inds,m,n);

        for i=1:m
%
        kcounts(i) = sum(inds(i,1:n));
        ds(i) = sig^2 / kcounts(i);
    end



        end
%
%
%
%
%
        function [x_hat,a_hat,err_hats,bs] = js_het3(x,ds,m)
%
%   estimates A once, then plugs it into the oracle
%
        a_hats = zeros(m,1);a_inits = zeros(m,1);es = zeros(m,1);
        fvals = zeros(m,1);dstars = zeros(m,1);es_all = zeros(m,m);
        bs = zeros(m,1);x_hat = zeros(m,1);

        ss = sum(x.^2,2);

%
%   compute the values E_j; each one is an unbiased estimator of A
%
        for j=1:m
%
        es(j) = ss(j) - ds(j);
    end

        a_init = mean(es);

        nds = ones(m,1);
        [a_hat,ierr] = root_newton(es,ss,ds,nds,m,a_init);
%
        [x_hat,err_hats,bs] = js_het_oracle(x,ds,m,a_hat)

        end
%
%
%
%
%
        function [x_hat,err_hats,bs] = js_het_oracle(x,ds,m,aa)
%
        x_hat = zeros(m,1);bs = zeros(m,1);err_hats = zeros(m,1);

        for ii=1:m
%
        bs(ii) = ds(ii) / (aa + ds(ii));

        bs(ii) = max(bs(ii),0);
        bs(ii) = min(bs(ii),1);

        err_hats(ii) = ds(ii)*(1 - bs(ii));
%
        x_hat(ii) = (1 - bs(ii)) * x(ii);
    end

        
        end
%
%
%
%
%
        function [x_hat,a_hats,bs,dstars] = js_het2(x,ds,m)
%
        a_hats = zeros(m,1);a_inits = zeros(m,1);es = zeros(m,1);
        fvals = zeros(m,1);dstars = zeros(m,1);es_all = zeros(m,m);
        bs = zeros(m,1);x_hat = zeros(m,1);

        ss = sum(x.^2,2);

        for ii=1:m

        ii

        nds = ones(m,1);
        nds(ii) = 3;
%
%   compute the values E_j; each one is an unbiased estimator of A
%
        es = 0*es;
        for j=1:m
%
        es(j) = (ss(j) - nds(j)*ds(j)) / nds(j);
    end

        es_all(1:m,ii) = es;
%
%   use newton to find the MLE for A, starting at the mean of the E_j's
%
        a_inits(ii) = mean(es);
        [a_hats(ii),ierrs(ii)] = root_newton(es,ss,ds,nds,m,a_inits(ii));
%
%   find the Fisher information at estimated A_hat and compute d^*
%
        for j=1:m
%
        [fvals(j),der1,der2] = eval_fisher(a_hats(ii),ds(j),nds(j));
    end

        dstars(ii) = 2*(a_hats(ii) + ds(ii))^2*sum(fvals);

        xx = (dstars(ii) - 4) / dstars(ii);
        yy = ds(ii) / (a_hats(ii) + ds(ii));
        bs(ii) = xx * yy;
%
        x_hat(ii) = (1 - bs(ii)) * x(ii);
    end

        end
%
%
%
%
%
        function [root,ierr] = root_newton(es,ss,ds,n_is,m,aa0)
%
%   want to find root of the objective function 
%
        nmax=100;
        aas = zeros(1,nmax);
        vals = zeros(1,nmax);
%
        thresh = 1d-13;

        aas(1) = aa0;

        for ijk = 1:nmax

        root = aas(ijk);
        [val,der1,der2] = eval_obje2(aas(ijk),es,ds,n_is,m);

        if (abs(val) <= thresh)
           root = aas(ijk);
           break;
        end

        ierr=0;
        if (ijk == nmax)
           ierr=1
        end

%%%        val = aas(ijk)^2;
%%%        der1 = 2*aas(ijk);

        aas(ijk+1) = aas(ijk) - val/der1;

    end

        return

        errs = aas(1:ijk) - aas(ijk)
        chk0 = eval_obje2(aas(ijk),es,ds,n_is,m)
        end
%
%
%
%
%
        function [val,der1,der2] = eval_obje2(aa,es,ds,n_is,m)
%
        fvals = zeros(m,1);
        fders1 = zeros(m,1);
        fders2 = zeros(m,1);
%
%   compute Fisher informations at each point
%
        for i=1:m
%
        [fvals(i),fders1(i),fders2(i)] = eval_fisher(aa,ds(i),n_is(i));
    end
%
%   weights are normalized Fisher information values
%
        whts = fvals / sum(fvals);
        zz = sum(whts .* es);
%
        yy = sum(es.*fders1)/sum(fvals);
        yy_der = sum(es.*fders2)/sum(fvals) ...
           - sum(es.*fders1)*sum(fders1)/sum(fvals)^2;
%
        xx = sum(fders1)/sum(fvals);
        xx_der = sum(fders2)/sum(fvals) - sum(fders1)^2/sum(fvals)^2;

        val = aa - zz;
        der1 = 1 - (yy - zz*xx);
        der2 = - (yy_der - zz*xx_der + (der1-1)*xx);
        end
%
%
%
%
%
        function [fval,der1,der2] = eval_fisher(a_hat,d,n)
%
        xx = 2 * (a_hat + d)^2;
        fval = n / xx;
        
        der1 = -n / (a_hat + d)^3;
        der2 = 3*n / (a_hat + d)^4;
        end
%
%
%
%
%
        function err_pred = js_err_fmla(aa,sig)
%
        dvar = sig^2;
        err_pred = dvar * aa/(dvar+aa);

        end
%
%
%
%
%
        function [rmu_jst,sig_n] = js_vector(xs,m,n,sig,rmu)
%
%   James-Stein estimator for full data with n samples; so the 
%   standard deviation gets scaled down by a factor of sqrt(n), and then
%   the usual JS estimator is applied, and then data is scaled up again
%
        xmean = mean(xs,2);
        sig_n = sig / sqrt(n);
%
%   normalize each coordinate to have standard deviation 1
%
        xmean = xmean / sig_n;
%
%   shrink the normalized vector, then scale back up by standard deviation
%
        [rmu_jst,err_hat,b_hat] = js_simple(xmean,m);
        rmu_jst  = rmu_jst * sig_n;

        end
%
%
%
%
%
        function [rmu_jst,err_hat,b_hat] = js_simple(x,m)
%
%   compute classical James-Stein shrinkage estimator on vector x in R^m,
%   assuming standard deviation 1
%
        rmu_mle = mean(x,2);
%
        ss = sum(x.^2);
        b_hat = (m-2) / ss ;
%
        shr = 1 - b_hat;
        rmu_jst = shr * rmu_mle;

        err_hat = 1 - (m-2)^2 / m / ss;

        end
%
%
%
%
%
        function rmu_shr = js_cheat(ys,m,n,inds,aa,sig)
%
%
%   computes shrinkage estimator rmu_shr of the mean, knowing sig and aa
%   (aa is the variance of the gaussian prior on the mean, which is assumed
%   to be mean zero)
%
%   the counts per row can be different
%
%   sig - the standard deviation of the noise
%
%   aa - the variance in the prior
%

        vars = zeros(m,1);
        sdevs = zeros(m,1);
        isums = zeros(m,1);
        bs = zeros(m,1);
        shrs = zeros(m,1);
        rmu_shr = zeros(m,1);

        rmu_mle = mean_estim(ys,inds,m,n);
        for i=1:m
%
        isums(i) = sum(inds(i,1:n));
        vars(i) = sig^2/isums(i);
        sdevs(i) = sqrt(vars(i));
%
        bs(i) = vars(i) / (vars(i) + aa);
        shrs(i) = 1 - bs(i);
        rmu_shr(i) = shrs(i) * rmu_mle(i);
    end

        end
%
%
%
%
%
        function test_angles
%
        gam = .8
        m=200
        n=floor(m/gam)
        k=1
%
        delta = 1
        evmin = 1 + (sqrt(gam) + 2)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);


        nsets=n;
        u_all = zeros(m,nsets);
        c_all = zeros(1,nsets);

        ll = n-1;

        iset_j = zeros(1,ll);
        for j=1:nsets
%
        j
        iset_j = randperm(n);
        iset_j = iset_j(1:ll);
%%%        iset_j(1:j-1) = 1:j-1;
%%%        iset_j(j:n-1) = j+1:n;

%
        xs_j = xs(1:m,iset_j);
        [u_j,s_j,v_j] = svdsmart(xs_j,m,n,1);

        if (u_j(1) < 0)
           u_j = -u_j;
        end

        c_j = abs(sum(u_j .* u))

        u_all(1:m,j) = u_j;
        c_all(j) = c_j;
    end

        u2 = sum(u_all,2);
        u2 = median(u_all,2);
        u2 = u2/norm(u2);


        [uu,ss,vv] = svdsmart(xs,m,n,1);

        max(c_all)
        cos_new = abs(sum(u2 .* u))
        cc = abs(sum(uu .* u))

        end
%
%
%
%
%
        function test_uneven()
%
        gam = .8
        m=1000
        n=floor(m/gam)
        k=10
%
        delta = .2
        evmin = 1 + (sqrt(gam) + 10)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);


%
%   add mean vector rmu to xs_clean and xs
%
        rmu = ones(m,1);
        xs_clean = xs_clean + repmat(rmu,1,n);
        xs = xs_clean + rnoise;
    
        kk=100
        pvec = [kk+1:m+kk]';
        pvec=pvec.^2;
        pvec = pvec/max(pvec)
%%%        plot(pvec,'*')

        ps = repmat(pvec,1,n);
        inds = rand_inds_fast2(m,n,ps);
%
        ys = xs .* inds;


        ymu_hat = mean_estim(ys,inds,m,n);

        ss = sum(inds,2);


        

        ymu_hat2 = ymu_hat.*sqrt(ss)/n


figure; plot(ymu_hat2,'*')

        stopnow

        yhat = ys*ys'/n;

        nijs=inds*inds';
        xhat = yhat2xhat_nijs(yhat,nijs,m,n);
        
        norm(xhat,'fro')
        norm(rcov2,'fro')


        res = xhat-rcov2;
        imagesc(abs(res))

stopnow
        plot(res(2,1:m),'*')

        end
%
%
%
%
%
        function test_whiten()
%
        gam = .8
        m=800
        n=floor(m/gam)
        k=10
%
        delta = .3
        evmin = 1 + (sqrt(gam) + 10)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);

%%%        [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,rnoise,...
%%%           inds,ys,probs] = prep_discr(m,n,k,delta,evmin);

        k_est=k+20


%
%   make diagonal covariance for noise
%
        dw_vals = sqrt([1:m]);
%%%        dw_vals = [1:m].^2;
%%%        dw_vals=ones(1,m);
        dw_vals = 2*dw_vals / max(dw_vals)
        dw = diag(dw_vals);
        rnoise2 = dw*rnoise;
%
        xs2 = xs_clean + rnoise2;
        ys2 = xs2 .* inds;
%
        est_mean = mean_estim(ys2,inds,m,n);
        ys3 = (ys2 - repmat(est_mean,1,n)).*inds;
%
        ycov = ys2*ys2'/(n*delta);
        vars = diag(ycov);

%%%        figure; plot(sqrt(vars),'*')
%%%        hold on; plot(dw_vals,'*')


        ys_wh = diag(1./sqrt(vars)) * ys2;


        est_mean = mean_estim(ys_wh,inds,m,n);
        ys_wh = (ys_wh - repmat(est_mean,1,n)).*inds;


%%%figure; imagesc(ys2)
%%%figure; imagesc(ys_wh)

        norm(ys_wh(1,1:n))^2 / (n*delta)
        norm(ys_wh(10,1:n))^2 / (n*delta)
        norm(ys_wh(100,1:n))^2 / (n*delta)
%

        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys3,m,n,k_est,delta);
        ys_raj = ys_raj + repmat(est_mean,1,n);
%
        [ys_raj2,s_hat,u,v] = raj_shrink_mat(ys_wh,m,n,k_est,delta);
        ys_raj2 = diag(sqrt(vars)) * ys_raj2;
        ys_raj2 = ys_raj2 + repmat(est_mean,1,n);
%
        [xs3,xhat3,sig_hat3] = denoise_matr(ys_wh,inds,m,n,delta,k_est,...
           'f','d');
        xs3 = diag(sqrt(vars)) * xs3;
%
        [xs4,xhat4,sig_hat4] = denoise_matr(ys_wh,inds,m,n,delta,k_est,...
           'o','d');
        xs4 = diag(sqrt(vars)) * xs4;
%
        err_raj = norm(ys_raj - xs_clean,'fro')
        err_raj2 = norm(ys_raj2 - xs_clean,'fro')
        err_xs3 = norm(xs3 - xs_clean,'fro')
        err_xs4 = norm(xs4 - xs_clean,'fro')
%
        norm(xs_clean,'fro')


        end
%
%
%
%
%
        function test_frey()
%
        m=560
        gam=.8
        n=floor( m/ gam)
%
        xs_clean = load_frey(n);
%
        sig=30
        delta = .2
%
%   add noise and subsample
%

        rnoise = sig*randn(m,n);

        dd = diag ( sqrt([1:m]));
%%%        rnoise = dd*rnoise;


%
        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        xs = xs_clean + rnoise;


%%%        display_frey(xs(1:m,1),1); colormap(gray)
%%%stopnow

        ys = xs .* inds;
        est_mean = mean_estim(ys,inds,m,n);
        ys2 = (ys - repmat(est_mean,1,n)).*inds;
%

%%%        k_est = 200;
        [uys,sys,vys] = svdsmart(ys,m,n,m);
        ell=sys.^2/n;
        [k_est,sigma_hat] = KN_rankEst(ell,n,1,.95,m)


        [xs2,xhat2,sig_hat2] = denoise_matr17(ys,inds,m,n,delta,k_est);

        [xs3,xhat3,sig_hat3] = denoise_matr(ys,inds,m,n,delta,k_est,...
           'f','n');


        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);
        [ys_raj2,s_hat,u,v] = raj_shrink_mat(ys2,m,n,k_est,delta);
        ys_raj2 = ys_raj2 + repmat(est_mean,1,n);

        err2 = norm(xs2 - xs_clean,'fro')
        err3 = norm(xs3 - xs_clean,'fro')
        err_raj = norm(ys_raj - xs_clean,'fro')
        err_raj2 = norm(ys_raj2 - xs_clean,'fro')
        norm(xs_clean,'fro')

        stopnow
        end
%
%
%
%
%
        function data = load_frey(n)
%
        load ../frey_rawface.mat

        data2=double(ff);
        clear ff;

        [m,n2]=size(data2);

        inds = randperm(n2);
        inds = inds(1:n);

        data = data2(1:m,inds(1:n));

%%%        ifig=1
%%%        display_frey(data(1:m,n),ifig)

        end
%
%
%
%
%
        function display_frey(ff,ifig)
%
        ff2 = reshape(ff,20,28)';
        figure(ifig);
        imagesc(ff2);
        end
%
%
%
%
%
        function test_mnist()
%
        m=784
        gam=.8
        n=floor( m/ gam)

        xs_clean = load_mnist(n);

        sig=2.5
        delta = .5
%
%   add noise and subsample
%

        rnoise = sig*randn(m,n);

        dd = diag ( sqrt([1:m]));
        rnoise = dd*rnoise;
%
        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        xs = xs_clean + rnoise;

%%%        aa = display_mnist(xs(1:m,10),1); colormap(gray)
%%%        aa = display_mnist(xs_clean(1:m,10),2); colormap(gray)
%%%stopnow
        ys = xs .* inds;

        est_mean = mean_estim(ys,inds,m,n);
        ys2 = (ys - repmat(est_mean,1,n)).*inds;

%

%%%        k_est = 200;
        [uys,sys,vys] = svdsmart(ys,m,n,m);
        ell=sys.^2/n;
        [k_est,sigma_hat] = KN_rankEst(ell,n,1,.95,m)



        [xs2,xhat2,sig_hat2] = denoise_matr17(ys,inds,m,n,delta,k_est);

        [xs3,xhat3,sig_hat3] = denoise_matr(ys,inds,m,n,delta,k_est,...
           'f','n');


        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);
        [ys_raj2,s_hat,u,v] = raj_shrink_mat(ys2,m,n,k_est,delta);
        ys_raj2 = ys_raj2 + repmat(est_mean,1,n);

        err2 = norm(xs2 - xs_clean,'fro')
        err3 = norm(xs3 - xs_clean,'fro')
        err_raj = norm(ys_raj - xs_clean,'fro')
        err_raj2 = norm(ys_raj2 - xs_clean,'fro')
        norm(xs_clean,'fro')
stopnow

%%%        norm(ys_raj,'fro')
%%%        norm(xs2,'fro')

save data199.mat

        ii=11
        aa = display_mnist(ys_raj(1:m,ii),1);
        aa = display_mnist(ys_raj2(1:m,ii),2);
        aa = display_mnist(xs2(1:m,ii),3);
        aa = display_mnist(xs_clean(1:m,ii),4);

        end
%
%
%
%
%
        function [data,labels] = load_mnist(n)
%
        load ../mnist_all.mat

        data2 = zeros(n*10,784);
        labels2 = zeros(1,n*10);

        kcount=0;

%%%        figure(1); colormap('gray')

        for i=0:9
%
        ichar = num2str(i);
        vname2 = strcat('train',ichar);
        arr = eval(vname2);
        [n1,m] = size(arr);

        inds = randperm(n1);
        inds=inds(1:n);

        data2(kcount+1:kcount+n,1:m) = arr(inds,1:m);
        labels2(kcount+1:kcount+n) = i;
        kcount=kcount+n;



%%%        subplot(2,5,i+1); imagesc(reshape(data2(kcount,1:m),28,28)')
    end

        iperm = randperm(10*n);
        inds = iperm(1:n);
        data = data2(inds,1:m)';
        labels = labels2(inds);


        return

        jj=99
        imagesc(reshape(data(1:m,jj),28,28)')
        labels(jj)

        end
%
%
%
%
%
        function aa = display_mnist(vec,ifig)
%
        aa = reshape(vec,28,28)';
        figure(ifig);
        imagesc(aa)

        end
%
%
%
%
%
        function aas = make_mat17(ys,inds,m,n,rcov)
%
        aas = zeros(m,m,n);

        for i=1:n

        i
%
%   i^th projection matrix
%
        kdim = sum(inds(1:m,i));
        proj_i = zeros(kdim,m);


        ind = 0;
        for j=1:m
%
        if (inds(j,i) == 1)
           ind = ind+1;
           proj_i(ind,j) = 1;
        end
    end


        bb = proj_i * rcov * proj_i' + eye(kdim);

        tol = 1e-10;
        bb = .5*(bb + bb');
        [b_inv,err1,err2] = pseudo_inv(bb,kdim,tol);

%%%        err1
%%%        err2

        aas(1:m,1:m,i) = rcov - rcov * proj_i' * b_inv * proj_i * rcov;
    end

        end
%
%
%
%
%
        function filt = make_filter17(rcov,proj,m,rnu)
%
        bb = proj * rcov * proj + rnu^2 * eye(m);
        b_inv2 = pinv(bb);

        filt = rcov * proj * b_inv2;


        end
%
%
%
%
%
        function test_usps()
%
        rand(1,10);

        delta = .6;
        gam = .5;
        sig = 2

%
%   load matrix of digits
%
        [data,npix,nusps] = usps_data;
        data = 10*data/max(data(:));
        m=npix
        n = floor( m/gam )


        [uu,rr]=qr(randn(m,m));
        uu=eye(m);
        data_old = data;
        data = uu*data;

%
%   select n random columns
%
        ichos = randperm(nusps);
        ichos = ichos(1:n);
        xs_clean = data(1:npix,ichos);
%%%        imagesc(xs_clean)

%
%   add noise and subsample
%
        rnoise = sig*randn(m,n);
%      
        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        xs = xs_clean + rnoise;


%%%figure; imagesc(reshape(xs(:,1),16,16))
%%%figure; imagesc(reshape(xs_clean(:,1),16,16))
%%%stopnow

        ys = xs .* inds;

        [uys,sys,vys] = svdsmart(ys,m,n,m);
        ell=sys.^2/n;
        [k_est,sigma_hat] = KN_rankEst(ell,n,1,.95,m)

        est_mean = mean_estim(ys,inds,m,n);
%%%norm(est_mean)
%%%figure; imagesc(reshape(est_mean,16,16))
%%%stopnow
        ys2 = (ys - repmat(est_mean,1,n)).*inds;

%
%%%        k_est = 30;




        [xs2,xhat2,sig_hat2] = denoise_matr17(ys,inds,m,n,delta,k_est);

        [xs3,xhat3,sig_hat3] = denoise_matr(ys,inds,m,n,delta,k_est,...
           'f','n');

        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);
        [ys_raj2,s_hat,u,v] = raj_shrink_mat(ys2,m,n,k_est,delta);
        ys_raj2 = ys_raj2 + repmat(est_mean,1,n);

        err2 = norm(xs2 - xs_clean,'fro')
        err3 = norm(xs3 - xs_clean,'fro')
        err_raj = norm(ys_raj - xs_clean,'fro')
        err_raj2 = norm(ys_raj2 - xs_clean,'fro')
        norm(xs_clean,'fro')

%%%        stopnow
        figure;
        ii=11
        subplot(2,2,1); imagesc(reshape(uu'*xs_clean(1:m,ii),16,16))
        subplot(2,2,2); imagesc(reshape(uu'*xs(1:m,ii),16,16))
        subplot(2,2,3); imagesc(reshape(uu'*xs2(1:m,ii),16,16))
        subplot(2,2,4); imagesc(reshape(uu'*ys_raj2(1:m,ii),16,16))

        set(figure(1),'Position',[500,500,440,400])
        colormap('gray')

stopnow





        [xs1,xhat1,sig_hat1] = denoise_matr17(ys,inds,m,n,delta,k_est);



%%%        [xs1,xhat1] = denoise_matre17(ys,inds,m,n,delta,k_est);
%%%        [xs2,xhat2] = denoise_matre(ys,inds,m,n,delta,k_est,...
%%%           'f','d');
%
        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);

        norm(xs1-xs2)
        norm(xhat1-xhat2)
        err1 = norm(xs1 - xs_clean,'fro')

        err_raj = norm(ys_raj - xs_clean,'fro')











        res1 = xs1 - xs_clean;
        figure; histogram(res1(:))


        figure; histogram(rnoise(:))

stopnoe
        res_raj = ys_raj - xs_clean;
        figure; histogram(res_raj(:))

stopnow
%
        norm(xs_clean,'fro')
        norm(xs1,'fro')
        norm(xs2,'fro')

        rank(xs2)
        rank(xs1)



%%%stopnow
        figure;
        ii=25
        subplot(2,2,1); imagesc(reshape(xs_clean(1:m,ii),16,16))
        subplot(2,2,2); imagesc(reshape(xs(1:m,ii),16,16))
        subplot(2,2,3); imagesc(reshape(xs1(1:m,ii),16,16))
        subplot(2,2,4); imagesc(reshape(ys_raj(1:m,ii),16,16))

        set(figure(1),'Position',[200,500,440,400])
        colormap('gray')

        end
%
%
%
%
%
        function [data,npix,nusps] = usps_data
%
        load ../mat_files/usps_all.mat
        data = real(data);
%
        load ../mat_files/usps_all.mat;
        data=double(data);
        data = reshape(data,256,[]);
%
        [npix,nusps] = size(data)
% 

        end
%
%
%
%
%
        function test_filt()
%
        gam = .8
        m=100
        n=floor(m/gam)
        k=10
%
        delta = .5
        evmin = 1 + (sqrt(gam) + 15)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);

%%%        [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,rnoise,...
%%%           inds,ys,probs] = prep_discr(m,n,k,delta,evmin);

        k_est=k+20



        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);

        [ys_raj2,relmse_hat,mse_hat] = optshrink(ys,k_est);

        chk0 = norm(ys_raj2/delta-ys_raj)

%%%stopnow


        [xs2,xhat2,sig_hat2] = denoise_matr(ys,inds,m,n,delta,k_est,...
           'f','d');


        [xs1,xhat1,sig_hat1] = denoise_matr17(ys,inds,m,n,delta,k_est);
%%%        [xs1,xhat1] = denoise_matre17(ys,inds,m,n,delta,k_est);
%%%        [xs2,xhat2] = denoise_matre(ys,inds,m,n,delta,k_est,...
%%%           'f','d');
%
        [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta);

        norm(xs1-xs2)
        norm(xhat1-xhat2)
        err1 = norm(xs1 - xs_clean,'fro')
        err2 = norm(xs2 - xs_clean,'fro')
        err3 = norm(ys_raj - xs_clean,'fro')
%
        norm(xs_clean,'fro')
        norm(xs1,'fro')
        norm(xs2,'fro')

        rank(xs2)
        rank(xs1)


        stopnow

        figure; plot(svds(xs1,100),'*')
        hold on; plot(svds(xs2,100),'*')


%
        sig_hat1
        sig_hat2


        xs200 = wfilter(ys,inds,m,n,xhat2,1);
        xs100 = wfilter(ys,inds,m,n,xhat1,1);

 xs200-xs100
%%%        xs2-xs200

        end
%
%
%
%
%
        function [xs_filt,xhat2,sig_hat] = denoise_matr17(ys,...
           inds,m,n,delta,k_est)
%
%   estimates covariance matrix and applies Wiener filtering; 
%   starts by estimating standard deviation of noise and rescaling signal
%
%
%
%   . . . estimate mean and subtract it off from the data
%
        est_mean = mean_estim(ys,inds,m,n);
        for i=1:n
%
        ys(1:m,i) = ys(1:m,i) - est_mean.*inds(1:m,i);
    end
%

%
%   estimate noise level and rescale signal so noise has unit variance
%
        sig_hat = estim_sig(ys,m,n,delta)
%%%        sig_hat=1;
        ys = ys/sig_hat;

        [xs_filt,xhat2] = denoise_matre17(ys,inds,m,n,delta,k_est);
%
%   scale estimator back up by sig_hat
%
        xs_filt = sig_hat*xs_filt;

%
%   add back the estimated mean
%
        for i=1:n
%
        xs_filt(1:m,i) = xs_filt(1:m,i) + est_mean;
    end



        end
%
%
%
%
%
        function [xs_filt,xhat2] = denoise_matre17(ys,inds,m,n,delta,...
           k_est)
%
%   estimates covariance matrix and applies Wiener filtering; supposes data
%   has already been rescaled by the standard deviation, so effective noise
%   has unit variance
%
%
%   uses wiener filtering with D-transform estimator of covariance of xs
%
        [xhat2,uxhat,sxhat2] = rov_shrink_mat(ys,m,n,k_est,delta);
        xhat2=.5*(xhat2+xhat2');

        rnu=1
        xs_filt = wfilter(ys,inds,m,n,xhat2,rnu);

        end
%
%
%
%
%
        function test_new()
%
        gam = .8
        m=1000
        n=floor(m/gam)
        k=10
%
        delta = .1
        evmin = 1 + (sqrt(gam) + 5)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
%%%        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
%%%           prep_gauss(m,n,k,delta,evs_noisy);

        [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,rnoise,...
           inds,ys,probs] = prep_discr(m,n,k,delta,evmin);

        k_est=k+50
        [xhat2,uxhat2,sxhat2] = cov_shrink_mat17(ys,inds,m,n,...
           delta,k_est,'f','d');


        [xhat3,uxhat3,sxhat3] = cov_shrink_mat(ys,inds,m,n,...
           delta,k_est,'f','d');

        err2 = norm(xhat2 - rcov,'fro')
        err3 = norm(xhat3 - rcov,'fro')
        norm(rcov,'fro')


        norm(xhat2,'fro')
        norm(xhat3,'fro')
        end
%
%
%
%
%
        function [xhat2,uxhat,sxhat2] = cov_shrink_mat17(ys,inds,m,n,...
           delta,k_est,loss,syst)
        sxhat2 = zeros(1,k_est);

%
%   . . . take sample covariance of the Y vectors
%
        yhat = ys*ys' / n;
%
%   invert system to get estimate of xs covariance
%
        if (syst == 'd')
            xhat = yhat2xhat_delta(yhat,delta,m);
        end
%
        if (syst == 'n')
            nijs = inds*inds';
            xhat = yhat2xhat_nijs(yhat,nijs,m,n);
        end
%
        [uxhat,sxhat17] = eigsmart(xhat,m,k_est);
        [uxhat17,sxhat] = svdsmart(ys,m,n,k_est);
        sxhat=sxhat/sqrt(n*delta);
%
%   apply optimal non-linearity to each eigenvalue
%
        gam=m/n;

        if (loss == 'o')
        for i=1:k_est
%
        sxhat2(i) = cov_emp2oper(sxhat(i),delta,gam);
        sxhat2(i) = sxhat2(i) - 1;
    end

        xhat2 = uxhat * diag(sxhat2) * uxhat';
    end
%
        if (loss == 'f')
        for i=1:k_est
%   
        sxhat2(i) = cov_emp2frob17(sxhat(i)^2,delta,gam);
        sxhat2(i) = sxhat2(i) - 1;
    end

        xhat2 = uxhat * diag(sxhat2) * uxhat';
    end

        end
%
%
%
%
%
        function rshrunk = cov_emp2frob17(rlam,delta,gam)
%
%   applies the asymptotically optimal nonlinearity for frobenius norm loss
%   to rlam in case of data missing with probability delta for m/n = gam
%
         rlam2 = rlam;
%%%%        rlam2 = delta*rlam + 1 - delta;
%
%   ell - the true population eigenvalue of x's (assuming rlam is the 
%   limiting empirical eigenvalue for x's
%
%
        [ell,cee,ess] = cov_emp2pop(rlam2,gam);
        ell = (ell - 1 + delta) / delta;
        rshrunk = ell*cee^2 + ess^2;

        end
%
%
%
%
%

        function [cov_hat,uu,s_hat] = rov_shrink_mat(ys,m,n,k_est,delta)
%
%
%                           description:
%
%   This function returns the estimator of a clean matrix X from the noisy
%   subsampled matrix ys by following the OptShrink algorithm of Raj Rao 
%   Nadakuditi.
%
%
%                         input parameters:
%
%   ys - the m-by-n matrix of data with zeros in the missing entries
%
%   m - the dimensionality of each data vector
%
%   n - the number of samples
%
%   delta - the (perhaps estimated) sampling probability
%
%   k_est - an a prior estimate of the rank (i.e. the maximum allowed value
%      for the rank)
%
%
%                             output parameters:
%
%   ys_raj - the denoised matrix
%
%   u - m-by-k_est matrix whose columns are the top left singular vectors
%      of ys_hat (and ys)
%
%   v - n-by-k_est matrix whose columns are the top right singular vectors
%      of ys_hat (and ys)
%
%   s_hat - the top k_est singular values of the estimated matrix ys_raj
%
%
%
        s_hat = zeros(1,k_est);
%
%   . . . rescale so variance of each noise entry is 1/n
%
        yhat = ys*ys'/n;
        xhat=yhat2xhat_delta(yhat,delta,m);
        [uu,ss] = eigsmart(xhat,m,k_est);

        ys = ys/sqrt(n*delta);
%
%   apply RRN shrinkage to each singular value
%
        


        [u,s,v] = svdsmart(ys,m,n,m);




        [uuu,ss] = eigsmart(xhat,m,m);
        s = sqrt(delta*ss + 1 - delta);



%%%        chk0 = norm(u*diag(s)*v' - ys,'fro')

        for i=1:k_est
%
        s_hat(i) = rov_emp2frob(delta,s,m,n,k_est,i);
    end
%
%   reconstruct the covariance matrix
%
        u=u(1:m,1:k_est);
%%%        cov_hat = u * diag(s_hat) * u';
        cov_hat = uu * diag(s_hat) * uu';

        end
%
%
%
%
%
        function eval = rov_emp2frob(delta,syw,m,n,k_est,ii)
%
%   applies the asymptotically optimal nonlinearity for frobenius norm loss
%   to rlam in case of data missing with probability delta for m/n = gam
%
        gam=m/n;
        bulk_edge = 1+sqrt(gam);
        if (syw(ii) <= bulk_edge)
           eval = 0;
           return;
        end

        [s_hat,dtop,dbot,t1,t2,t3,t4] = raj_emp2frob(syw,m,n,k_est,ii);
        eval = -2*t2/dbot/delta;
        
        end
%
%
%
%
%
        function compare_raj_svd()
%
        gam = .5
        m=600
        n=floor(m/gam)
        k=2
%
        delta = .15
        evmin = 1 + (sqrt(gam) + 15)/delta
%%%        evmin = 1 + sqrt(gam)/delta + 150
        evs_noisy = evmin + [k-1:-1:0]
%
        [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,rcov2,evs_clean] = ...
           prep_gauss(m,n,k,delta,evs_noisy);

%%%        [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,rnoise,...
%%%           inds,ys,probs] = prep_discr(m,n,k,delta,evmin);

        [utr,strue,vtr] = svdsmart(xs_clean,m,n,m);
        strue = strue/sqrt(n);

        yw = ys / sqrt(n*delta);
        [uyw,syw,vyw] = svdsmart(yw,m,n,m);
%
        k_est=k+40
        ii=k_est
%
%   frobenius norm shrinkage
%
        s_raj = raj_emp2frob(syw,m,n,k_est,ii);
        s_raj = s_raj/sqrt(delta);
%
        s_gav = svd_emp2frob(syw(ii),gam);
        s_gav = s_gav/sqrt(delta);
%
%   left and right cosines:
%
        ci_left = sum(uyw(1:m,ii) .* utr(1:m,ii));
        ci_rght = sum(vyw(1:n,ii) .* vtr(1:n,ii));
%
        s_raj
        s_gav
        strue(ii)*ci_left*ci_rght



%
%   operator norm shrinkage
%

        s_raj = raj_emp2oper(syw,m,n,k_est,ii);
        s_raj = s_raj/sqrt(delta);
%
        s_gav = svd_emp2oper(syw(ii),gam);
        s_gav = s_gav/sqrt(delta);
%
        s_raj
        s_gav
        strue(ii)


        end
%
%
%
%
%
        function s_hat = raj_emp2oper(s,m,n,k,ii)
%
%
%   estimate i^th singular value using RRN's formula
%
%
%   estimate i^th singular value using RRN's formula
%
        s_tail = s(k+1:m);
        si = s(ii);
%
        n1=n-k;
        m1=m-k;

        t1 = sum(1./(si^2 - s_tail.^2))*si/n1 + (n1-m1)/si/n1;
        t2 = sum(1./(si^2 - s_tail.^2))*si/m1;
%
        dd = t1*t2;
%
        s_hat = 1/sqrt(dd);


        end
%
%
%
%
%
        function [ys_raj,s_hat,u,v] = raj_shrink_mat(ys,m,n,k_est,delta)
%
%
%                           description:
%
%   This function returns the estimator of a clean matrix X from the noisy
%   subsampled matrix ys by following the OptShrink algorithm of Raj Rao 
%   Nadakuditi.
%
%
%                         input parameters:
%
%   ys - the m-by-n matrix of data with zeros in the missing entries
%
%   m - the dimensionality of each data vector
%
%   n - the number of samples
%
%   delta - the (perhaps estimated) sampling probability
%
%   k_est - an a prior estimate of the rank (i.e. the maximum allowed rank
%      for the rank)
%
%
%                             output parameters:
%
%   ys_raj - the denoised matrix
%
%   u - m-by-k_est matrix whose columns are the top left singular vectors
%      of ys_hat (and ys)
%
%   v - n-by-k_est matrix whose columns are the top right singular vectors
%      of ys_hat (and ys)
%
%   s_hat - the top k_est singular values of the estimated matrix ys_raj
%
%
%
        s_hat = zeros(1,k_est);
%
%   . . . rescale so variance of each noise entry is 1/n
%
        ys = ys/sqrt(n);
%
%   apply RRN shrinkage to each singular value
%
        [u,s,v] = svd(ys,'econ');
        s = diag(s);
%%%        chk0 = norm(u*diag(s)*v' - ys,'fro')
%%%        k_est

        for i=1:k_est
%
        s_hat(i) = raj_emp2frob(s,m,n,k_est,i);
    end
%
%   scale up again and reconstruct matrix
%
        s_hat = s_hat*sqrt(n);
        u=u(1:m,1:k_est);
        v=v(1:n,1:k_est);
        ys_raj = u * diag(s_hat) * v';


%
%   finally, rescale estimator to account for subsampling
%
        s_hat = s_hat/delta;
        ys_raj = ys_raj/delta;

        end
%
%
%
%
%
        function [s_hat,dtop,dbot,t1,t2,t3,t4] = raj_emp2frob(s,m,n,k,i)
%%%        function s_hat = raj_emp2frob(s,m,n,k,i)
%
%   estimate i^th singular value using RRN's formula
%
%
%   estimate i^th singular value using RRN's formula
%
        s_tail = s(k+1:m);
        si = s(i);
%
        n1=n-k;
        m1=m-k;

        t1 = sum(1./(si^2 - s_tail.^2))*si/n1 + (n1-m1)/si/n1;
        t2 = sum(1./(si^2 - s_tail.^2))*si/m1;
%
        t3 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/m1;
        t4 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/n1;
        t4 = t4 - 2*(n1-m1)/si^2/n1 + (n1-m1)/si^2/n1;
%
        dbot = t1*t3 + t2*t4;
        dtop = t1*t2;
%
        s_hat = -2*dtop/dbot;



        end
%
%
%
%
%
        function [s_hat,dtop,dbot,t1,t2,t3,t4] = raj_emp2frob_slow(s,...
           m,n,k,ii)
%
        si = s(ii);
        n1=n-k
        m1=m-k
        s_tail = s(k+1:m);

        xmat = zeros(m1,n1);

        for i=1:m1
%
        xmat(i,i) = s_tail(i);
    end

        xxt = xmat*xmat';
        xtx = xmat'*xmat;

        eye_n = eye(n1);
        eye_m = eye(m1);

        aa = si^2*eye_m - xxt;
        aa_inv = zeros(m1,m1);

        for i=1:m1
%
        aa_inv(i,i) = 1/aa(i,i);
    end

        bb = si^2*eye_n - xtx;
        bb_inv = zeros(n1,n1);

        for i=1:n1
%
        bb_inv(i,i) = 1/bb(i,i);
    end

        chk0 = norm(aa_inv*aa - eye(m1))
        chk0 = norm(bb_inv*bb - eye(n1))
%
        t2 = trace(aa_inv)*si/m1
        t1 = trace(bb_inv)*si/n1
%
        gg = si^2*eye_n - xtx;
        gg = gg.^2;
        gg_inv = zeros(n1,n1);
        for i=1:n1
%
        gg_inv(i,i) = 1/gg(i,i);
    end

        hh = si^2*eye_m - xxt;
        hh=hh.^2;
        hh_inv = zeros(m1,m1);
        for i=1:m1
%
        hh_inv(i,i) = 1/hh(i,i);
    end
%

        chk0 = norm(gg_inv*gg - eye(n1))
        chk0 = norm(hh_inv*hh - eye(m1))
%

        t4 = -2*trace(gg_inv)*si^2/n1 + t1/si
        t3 = -2*trace(hh_inv)*si^2/m1 + t2/si

        dbot = t1*t3 + t2*t4;
        dtop = t1*t2;
%
        s_hat = -2*dtop/dbot;


        end
%
%
%
%
%
        function si_hat = svd_emp2pop_raj(s,m,n,k,i)
%
%
%   estimate i^th singular value using the D transform
%
        s_tail = s(k+1:m);
        si = s(i);
%
        n1=n-k;
        m1=m-k;

        t1 = sum(1./(si^2 - s_tail.^2))*si/n1;
        t2 = sum(1./(si^2 - s_tail.^2))*si/m1;
%
        dhat = t1*t2;
        si_hat = 1/sqrt(dhat);
%

        end
%
%
%
%
%

        function [ys_hat,uy,vy,sy2,sy] = svd_shrink_mat(ys,m,n,delta,...
           k_est,loss)
%
%
%                             description:
%
%   This function computes a singular value shrinkage estimator of the 
%   a data matrix with missing values, where each coordinate is sampled 
%   with probability delta (delta is given by user, and hence may be an 
%   estimate). The estimate is obtained by adjusting the asymptotically
%   optimal shrinkers of Donoho-Gavish to the case of missing data.
%
%   The noise is assumed to be white with every entry having variance 1 
%   (note that this is NOT how the noise is scaled in the Donoho-Gavish 
%   paper or the Nadakuditi paper).
%
%
%                           input parameters:
%
%   ys - the m-by-n matrix of data with zeros in the missing entries
%
%   m - the dimensionality of each data vector
%
%   n - the number of samples
%
%   delta - the (perhaps estimated) sampling probability
%
%   k_est - an a prior estimate of the rank (i.e. the maximum allowed rank
%      for the rank)
%
%   loss - a character that specifies which loss function to minimize
%      'o' - operator norm loss
%      'f' - frobenius norm loss
%
%
%                             output parameters:
%
%   ys_hat - the denoised matrix
%
%   uy - m-by-k_est matrix whose columns are the top left singular vectors
%      of ys_hat (and ys)
%
%   vy - n-by-k_est matrix whose columns are the top right singular vectors
%      of ys_hat (and ys)
%
%   sy2 - the top k_est singular values of the estimated matrix ys_hat
%
%   sy - the top k_est singular values of the input matrix ys
%
%
        sy2 = zeros(1,k_est);
        gam=m/n;
%
%   ys2 - data matrix normalized so it is effectively signal plus noise 
%      with variance 1, divided by sqrt(n)
%
        [uy,sy,vy] = svdsmart(ys,m,n,k_est);
        sy_w = sy / sqrt(n * delta);

%%%        chk0 = norm(uy*diag(sy)*vy' - ys2,'fro')
%%%        sy

%
%   apply shrinkage estimator to spectrum of ys2
%
        if (loss == 'f')
        for i=1:k_est
%
        sy2(i) = svd_emp2frob(sy_w(i),gam);
    end
%
    end
%
        if (loss == 'o')
        for i=1:k_est
%
        sy2(i) = svd_emp2oper(sy_w(i),gam);
    end
%
    end

%
%   reconstruct the estimator
%
        sy2 = sy2 * sqrt(n/delta);
        ys_hat = uy * diag(sy2) * vy';

        end
%
%
%
%
%
        function [xhat2,uxhat,sxhat2,err_pred] = cov_shrink_mat(ys,inds,m,n,...
           delta,k_est,loss,syst)
        sxhat2 = zeros(1,k_est); ells=zeros(1,k_est);
%
%
%                             description:
%
%   This function computes a shrinkage estimator of the covariance of data
%   with missing values, where each coordinate is sampled with probability
%   delta (delta is given by user, and hence may be an estimate). The 
%   estimate is obtained by taking the sample covariance of the 
%   missing data ys, applying a linear operator to get an unbiased
%   estimator of the covariance of the xs, and then adjusting the spectrum 
%   of this unbiased estimator.
%
%   The shrinker is supposed to be asymptotically optimal for the specified 
%   loss. It estimates the signal covariance, NOT the signal + noise 
%   covariance (which can be estimated by adding the identity matrix).
%
%   The mean of the distribution is assumed to be zero (the code does not 
%   perform mean subtraction). The noise is assumed to be white with
%   variance 1.
%
%
%                           input parameters:
%
%   ys - the m-by-n matrix of data with zeros in the missing entries
%
%   inds - the m-by-n matrix with 1's in observed entries, 0's elsewhere
%
%   m - the dimensionality of each data vector
%
%   n - the number of samples
%
%   delta - the (perhaps estimated) sampling probability
%
%   k_est - an a prior estimate of the rank (i.e. the maximum allowed rank
%      for the rank)
%
%   loss - a character that specifies which loss function to minimize
%      'o' - operator norm loss
%      'f' - frobenius norm loss
%
%   syst - a character that specifies which unbiased estimator of x to use
%      (that is, which linear operator to apply to the covariance of ys)
%      'd' - the delta system
%      'n' - the nijs system (i.e. available case estimator)
%
%
%                           output parameters:
%
%   xhat2 - the m-by-m estimated covariance matrix, obtained by shrinking
%      the unbiased estimator using the optimal operator norm shrinker
%
%   uxhat - the m-by-m matrix whose columns are the eigenvectors of xhat2
%
%   sxhat2 - the m-dimensional vector of eigenvalues of xhat2
%

%
%   . . . take sample covariance of the Y vectors
%
        yhat = ys*ys' / n;
%
%   invert system to get estimate of xs covariance
%
        if (syst == 'd')
            xhat = yhat2xhat_delta(yhat,delta,m);
        end
%
        if (syst == 'n')
            nijs = inds*inds';
            xhat = yhat2xhat_nijs(yhat,nijs,m,n);
        end

%
%%%        [uxhat,sxhat] = eigsort(xhat,m);
%%%        uxhat=uxhat(1:m,1:k_est);
        [uxhat,sxhat] = eigsmart(xhat,m,k_est);
%%%        chk0 = norm(xhat*uxhat - uxhat*diag(sxhat))
%
%   apply optimal non-linearity to each eigenvalue
%
        gam=m/n;

        if (loss == 'o')
        for i=1:k_est
%
        sxhat2(i) = cov_emp2oper(sxhat(i),delta,gam);
        ells(i) = sxhat2(i);
        sxhat2(i) = sxhat2(i) - 1;
    end

        [err_pred,ierr]=cov_err_oper(ells,k_est,delta,gam);
        xhat2 = uxhat * diag(sxhat2) * uxhat';
    end
%
        if (loss == 'f')
        for i=1:k_est
%
        [sxhat2(i),ells(i)] = cov_emp2frob(sxhat(i),delta,gam);
        sxhat2(i) = sxhat2(i) - 1;
    end
        [err_pred,ierr]= cov_err_frob(ells,k_est,delta,gam);
        xhat2 = uxhat * diag(sxhat2) * uxhat';
    end

        end
%
%
%
%
%
        function [u,s] = eigsmart(a,m,k)
%
%   returns the top k eigenvectors/values of the m-by-m symmetric matrix a.
%   Note that s is the vector of eigenvalues, not a diagonal matrix.
%
        if (m/k > 10)
%
        [u,s] = eigs(a,k);
        s = sort(diag(s),'descend');
    end

        if (m/k <= 10)
%
        [u,s] = eig(a);
        [s,ii] = sort(diag(s),'descend');
        s=s(1:k);
        u = u(1:m,ii(1:k));
    end


        end
%
%
%
%
%
        function [xs_clean,u,s,v,rnoise,xs,ys,inds,rcov,...
           rcov2,evs_clean] = prep_gauss(m,n,k,delta,evs_noisy)
%
%   prepares covariance, full and subsampled gaussian data with specified
%   spectrum and sampling probability delta
%
        evs_clean = evs_noisy - 1;
        sig=1;
        [rcov,u] = rand_cov3(m,k,evs_clean);
%%%        [rcov,u] = rand_cov_hada(m,k,evs_clean);
%
        [xs,xs_clean,rnoise] = draw_gaussian(u,evs_clean,m,k,n,sig);
        [u,s,v] = svdsmart(xs_clean,m,n,k);
%
        rcov2= rcov + eye(m);
%
        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        ys = xs.*inds;



        end
%
%
%
%
%
        function [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,...
           rnoise,inds,ys,probs,ichosen] = prep_discr2(m,n,k,delta,evmin)
%
%   prepares covariance, full and subsampled data from a discrete
%   disctribution (k+1 vectors sampled at random) with specified
%   minimum eigenvalue of covariance and sampling probability delta
%

        evmin_clean = evmin - 1;
        probs = rand(1,k+1);
        probs = probs/sum(probs);
        [rcov,vecs,u,evs_clean] = cov_discrete(evmin_clean,k,m,probs);
        evs_noisy = evs_clean + 1;
%
        rcov2 = rcov + eye(m);
        sig=1;
        [xs,xs_clean,rnoise,ichosen] = draw_discrete(vecs,probs,k,m,n,sig);
        [u,s,v] = svdsmart(xs_clean,m,n,k);

        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        ys = xs.*inds;



        end
%
%
%
%
%
        function [rcov,rcov2,vecs,u,s,v,evs_clean,xs,xs_clean,...
           rnoise,inds,ys,probs] = prep_discr(m,n,k,delta,evmin)
%
%   prepares covariance, full and subsampled data from a discrete
%   disctribution (k+1 vectors sampled at random) with specified
%   minimum eigenvalue of covariance and sampling probability delta
%

        evmin_clean = evmin - 1;
        probs = rand(1,k+1);
        probs = probs/sum(probs);
        [rcov,vecs,u,evs_clean] = cov_discrete(evmin_clean,k,m,probs);
        evs_noisy = evs_clean + 1;
%
        rcov2 = rcov + eye(m);
        sig=1;
        [xs,xs_clean,rnoise,ichosen] = draw_discrete(vecs,probs,k,m,n,sig);
        [u,s,v] = svdsmart(xs_clean,m,n,k);

        ps = delta*ones(m,n);
        inds = rand_inds_fast2(m,n,ps);
%
        ys = xs.*inds;



        end
%
%
%
%
%
        function sig_hat = estim_sig(ys,m,n,delta)
%
%   uses the Donoho-Gavish estimator of the noise level
%
        gam = m/n;
        svs = svd(ys);
        svs = svs/sqrt(n*delta);

        xmed = mp_median(gam);
        sig_hat = median(svs) / sqrt(xmed);


        end
%
%
%
%
%
        function xmed = mp_median(gam)
%
%   returns median value of MP density with parameter gam in (0,1]
%
        a = (1-sqrt(gam))^2;
        b = (1+sqrt(gam))^2;
        x = (a+b)/2;

        for ijk=1:500
%
        cdfx = mp_cdf(x,gam);

        if (abs(cdfx - .5) < 10*eps)
           break;
        end
%
        if (cdfx > .5)
           b_old=b;
           b=x;
           x = (a+b)/2;
        end
%
        if (cdfx <= .5)
           a_old=a;
           a=x;
           x = (a+b)/2;
        end

    end
        val = mp_cdf(x,gam);
        xmed = x;
        end
%
%
%
%
%
        function rint = mp_cdf(x,gam)
%
%   integrate MP density from (1-sqrt(gam))^2 up to x
%
        a=(1-sqrt(gam))^2;
        b=(1+sqrt(gam))^2;
%
        rr = (x - a) * (b - x);
        y1 = asin((2*x-a-b)/(b-a)) - asin(-1);
        y2 = asin( ((a+b)*x - 2*a*b) / (x*(b-a)) ) - asin(-1);
%
        rint = sqrt(rr) + (a+b)*y1/2 - sqrt(a*b)*y2;
        rint = rint/(2*pi*gam);

        end
%
%
%
%
%
        function val = mp_eval(t,gam)
%
        x0 = (1-sqrt(gam))^2;
        x1 = (1+sqrt(gam))^2;
        val = (x1 - t) * (t - x0);
        val = sqrt(val) / (2*pi*t*gam);

        end
%
%
%
%
%
        function [u,s,v] = svdsmart(a,m,n,k)
%
        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
        function s_hat = svd_emp2frob(s,gam)
%
%   for an empirical singular value s, returns the shrunken singular value
%   s_hat that is optimal for frobenius loss; data is missing uniformly
%   at random with probabiliy delta, gam is the aspect ratio, and we assume
%   the data matrix has been normalized by 1/sqrt(n*delta) already
%
        bulk_edge = 1 + sqrt(gam);
        if (s <= bulk_edge)
           s_hat = 0;
           return;
        end
%
        s_pop = svd_emp2pop_new(s,gam);
        [sy17,cosl,cosr,bulk_edge] = svd_pop2emp(s_pop,gam);
        s_hat = s_pop * cosl * cosr;

        end
%
%
%
%
%
        function [sv_emp,cosl,cosr,bulk_edge] = svd_pop2emp(rlam,gam)
%
%   returns limiting empirical values of eigenvalue and left and right
%   angles between empirical singular vectors and population singular
%   vectors if population singular value is rlam
%
%   sv_pop - the population singular value WITHOUT NOISE and completely
%      unnormalized (it has not been divided by sqrt(n) or sqrt(delta))
%
%   sv_emp - the returned value; equals the predicted singular value of
%      the data matrix (with missing entries), divided by sqrt(n*delta)
%
%   note that sv_emp is the inverse function to svd_emp2pop
%
        sv_emp = sqrt( (rlam + 1/rlam) * (rlam + gam/rlam) );
%
        cosl = sqrt( (1-gam/rlam^4) / (1 + gam/rlam^2) );
        cosr = sqrt( (rlam^4 - gam) / (rlam^4 + rlam^2));

        bulk_edge = 1 + sqrt(gam);


        end
%
%
%
%
%
        function s_hat = svd_emp2oper(s,gam)
%
%   s is an empirical  singular value of a data matrix of the form:
%
%                        X + Z/sqrt(n)
%
%   where the singular values of X do not depend on n, and Z is noise
%   with unit variance. X is m-by-n, and gam = m/n. With this scaling, the
%   shrinker does not depend on n.
%
%
%   . . . if s is below bulk_edge, it's set to 0
%
        bulk_edge = 1 + sqrt(gam);
        if (s <= bulk_edge)
           s_hat = 0;
           return;
        end
%
%   . . . otherwise, it's returned to the population value
%
        s_hat = svd_emp2pop_new(s,gam);

        end
%
%
%
%
%
        function x = svd_emp2pop_new(y,gam)
%
%   given empirical singular value y of m-by-n matrix with gam=m/n,
%   returns the population singular value WITHOUT NOISE and not
%   normalized (not divided by sqrt(n) or sqrt(delta)
%
%   note that this is the inverse function to the svd_pop2emp
%
        t1 = y^2 - gam - 1;
        x = t1 + sqrt(t1^2 - 4*gam);
        x = sqrt(x/2);



        return
%
%   perform sanity check:
%
        t2 = (1 + x^2)*(1+gam/x^2);
        chk0 = t2 - y^2

        end
%
%
%
%
%
        function [sv_emp,cosl,cosr] = svd_pop2emp_delta(sv_pop,gam,...
           delta)
%
%   returns limiting empirical values of eigenvalue and left and right
%   angles between empirical singular vectors and population singular
%   vectors if population singular value is rlam
%
%   sv_pop - the population singular value WITHOUT NOISE and completely
%      unnormalized (it has not been divided by sqrt(n) or sqrt(delta));
%      in other words, it is the standard deviation of the data 
%      along one of the population principal components
%
%   sv_emp - the returned value; equals the predicted singular value of
%      the data matrix (with missing entries), divided by sqrt(n*delta)
%
%   note that sv_emp is the inverse function to svd_emp2pop
%
        rlam = sv_pop*sqrt(delta);
        [sv_emp,cosl,cosr,bulk_edge] = svd_pop2emp(rlam,gam);

        end
%
%
%
%
%
        function xs_filt = denoise_dumb(ys,inds,m,n,sig,delta)
%
%   uses wiener filtering with unbiased estimator of covariance of xs
%
        yhat = ys * ys' / n;
        xhat = yhat2xhat_delta(yhat,delta,m);
        xhat = xhat - eye(m);
        rnu=1
        xs_filt = wfilter(ys,inds,m,n,xhat,rnu);

        end
%
%
%
%
%
        function [xs_filt,xhat2,sig_hat] = denoise_matr(ys,...
           inds,m,n,delta,k_est,loss,syst)
%
%   estimates covariance matrix and applies Wiener filtering; 
%   starts by estimating standard deviation of noise and rescaling signal
%
%
%
%   . . . estimate mean and subtract it off from the data
%
        est_mean = mean_estim(ys,inds,m,n);
        for i=1:n
%
        ys(1:m,i) = ys(1:m,i) - est_mean.*inds(1:m,i);
    end
%

%
%   estimate noise level and rescale signal so noise has unit variance
%
        sig_hat = estim_sig(ys,m,n,delta);
%%%        sig_hat=1;
        ys = ys/sig_hat;

        [xs_filt,xhat2] = denoise_matre(ys,inds,m,n,delta,k_est,loss,syst);

%
%   scale estimator back up by sig_hat
%
        xs_filt = sig_hat*xs_filt;

%
%   add back the estimated mean
%
        for i=1:n
%
        xs_filt(1:m,i) = xs_filt(1:m,i) + est_mean;
    end



        end
%
%
%
%
%
        function [xs_filt,xhat2] = denoise_matre(ys,inds,m,n,delta,...
           k_est,loss,syst)
%
%   estimates covariance matrix and applies Wiener filtering; supposes data
%   has already been rescaled by the standard deviation, so effective noise
%   has unit variance


%   uses wiener filtering with shrinkage estimator of covariance of xs,
%   with specified loss function and linear system
%
        [xhat2,uxhat,sxhat2] = cov_shrink_mat(ys,inds,m,n,delta,k_est,...
           loss,syst);
        xhat2=.5*(xhat2+xhat2');

        rnu=1
        xs_filt = wfilter(ys,inds,m,n,xhat2,rnu);

        end
%
%
%
%
%
        function est_mean = mean_estim(ys,inds,m,n)

        est_mean = zeros(m,1);

%%%        chk0 = norm(ys.*inds - ys,'fro')
%%%        chk0 = norm(ys.*(1-inds),'fro')

        for i = 1:m
%
        icount=0;
        for j=1:n
%
        est_mean(i) = est_mean(i) + ys(i,j);
        icount = icount+inds(i,j);
    end
        est_mean(i) = est_mean(i)/icount;
    end

        est_mean2 = sum(ys,2) ./ sum(inds,2);


        return
        chk0 = norm(est_mean - est_mean2)

        end
%
%
%
%
%
        function xs_filt = wfilter(ys,inds,m,n,rcov,rnu)
%
        xs_filt = zeros(m,n);
%
        for i=1:n
%
%   i^th projection matrix
%
        kdim = sum(inds(1:m,i));
        ii = find(inds(1:m,i) == 1);
        bb =  rcov(ii,ii) + rnu^2 * eye(kdim);

        x1 = ys(ii,i);
        x2 = bb \ x1;
%%%        x2 = mldivide(bb, x1);
        xs_filt(1:m,i) = rcov(1:m,ii) * x2;

    end
        
%%%        xs_filt 

        end
%
%
%
%
%
        function [xs_filt,xhat2] = denoise_matre_slow(ys,inds,m,n,delta,...
           k_est,loss,syst)
%
%   estimates covariance matrix and applies Wiener filtering; supposes data
%   has already been rescaled by the standard deviation, so effective noise
%   has unit variance


%   uses wiener filtering with shrinkage estimator of covariance of xs,
%   with specified loss function and linear system
%
        [xhat2,uxhat,sxhat2] = cov_shrink_mat(ys,inds,m,n,delta,k_est,...
           loss,syst);
        xhat2=.5*(xhat2+xhat2');

        xs_filt = wfilter_slow(ys,inds,m,n,xhat2);

        end
%
%
%
%
%
        function xs_filt = wfilter_slow(ys,inds,m,n,rcov)
%
        xs_filt = zeros(m,n);
%
        for i=1:n
%
%   i^th projection matrix
%
        kdim = sum(inds(1:m,i));
        proj_i = zeros(kdim,m);

        ind = 0;
        for j=1:m
%
        if (inds(j,i) == 1)
           ind = ind+1;
           proj_i(ind,j) = 1;
        end

    end

        bb = proj_i * rcov * proj_i' + eye(kdim);

        tol = 1e-10;
        bb = .5*(bb + bb');
        [b_inv,err1,err2] = pseudo_inv(bb,kdim,tol);

%%%        err1
%%%        err2

        aa = rcov * proj_i' * b_inv;
        xs_filt(1:m,i) = aa * proj_i * ys(1:m,i);

    end
        

        end
%
%
%
%
%
        function [a_inv,err1,err2] = pseudo_inv(a,m,tol)
%
        a = .5*(a + a');
        [u,s] = eigsort(a,m);
%%%        chk0 = norm(a - u*diag(s)*u','fro')

        s_inv = 0*s;
        for i=1:m
%
        if (s(i) <= tol)
          break;
        end

        s_inv(i) = 1/s(i);
    end


        a_inv = u * diag(s_inv) * u';
        a_inv = .5*(a_inv + a_inv');

        err1 = norm(a * a_inv * a - a,'fro');
        err2 = norm(a_inv * a * a_inv - a_inv,'fro');


        end
%
%
%
%
%
        function [xs_opt,u_opt,s_opt,v_opt] = call_optspace(ys,inds,m,n)
%

        ys_sparse = sparse(ys);
        niter=200
        tol=1e-6
        [u_opt,s_opt,v_opt,dist] = OptSpace(ys_sparse,[],niter,tol);

        s_opt

        xs_opt = u_opt * s_opt * v_opt';

        end
%
%
%
%
%
        function ell = compute_ell2(rlam,delta,gam)
%
%   rlam - asymptotic empirical value of unbiased estimator of covariance 
%   of the x's
%
%   ell - the true population eigenvalue of x's covariance
%
        rlam2 = (rlam-1 + 1/delta)*delta;
        ell = compute_ell(rlam2,gam);
        ell = ( ell - 1 + delta ) / delta;

        end
%
%
%
%
%
        function [rshrunk,ell] = cov_emp2frob(rlam,delta,gam)
%
%   applies the asymptotically optimal nonlinearity for frobenius norm loss
%   to rlam in case of data missing with probability delta for m/n = gam
%
        rlam2 = delta*rlam + 1 - delta;

        bulk_edge = (1+sqrt(gam))^2;
        if (rlam2 <= bulk_edge)
%
        rshrunk = 1;
        ell=1;
        return
    end

%
%   ell - the true population eigenvalue of x's (assuming rlam is the 
%   limiting empirical eigenvalue for x's
%
%
        [ell,cee,ess] = cov_emp2pop(rlam2,gam);
        ell = (ell - 1 + delta) / delta;
        rshrunk = ell*cee^2 + ess^2;

        end
%
%
%
%
%
        function [ell,cee,ess] = cov_emp2pop(rlam,gam)
%
%   if rlam is the asymptotic empirical eigenvalue, this returns the 
%   population eigenvalue and cosine of the angle between empirical
%   eigenvector and top population eigenvector
%
%   These are three values that are needed to compute any of the optimal
%   shrinkers
%
%   . . . find the population eigenvalue
%
        ell = compute_ell(rlam,gam);
%
%   . . . find the angles
%
        [sx_spike,sx_bulk,cc] = cov_pop2emp(ell,gam);
%%%        chk0 = rlam - sx_spike
%
%   sy1 - the population value for the y's
%
%
%   cee, ess - cosine and sine of angles
%
        cee = compute_cee(ell,gam);
        ess = compute_ess(cee);
%%%        chk0 = cee^2 + ess^2 - 1


        end
%
%
%
%
%
        function [sx1,sx2,cc] = cov_pop2emp(rlam,gam)
%
%   returns predicted asymptotic values of first and second eigenvalues
%   as well as angles, when the top eigenvalue (of signal plus noise,
%   NOT the signal alone) is rlam, and gam is ratio between row and column 
%   dimensions
%
        sx1 = rlam*(1 + gam/(rlam - 1));
        sx2 = (1 + sqrt(gam))^2;
%
        r1 = 1 - gam/(rlam-1)^2;
        r2 = 1 + gam/(rlam-1);
        cc = sqrt(r1/r2);
%
        end
%
%
%
%
%
        function [sx1,sx_bulk,cc,sy1,sy_bulk] = cov_pop2emp2(rlam,...
           delta,gam)
%
%   returns predicted asymptotic values of first and second eigenvalues
%   as well as angles, when the top eigenvalue (of signal plus noise,
%   NOT the signal alone) is rlam, delta is uniform
%   sampling rate, and gam is ratio between row and column dimensions
%
        rlam2 = delta*rlam + 1 - delta;
        [sy1,sy_bulk,cc] = cov_pop2emp(rlam2,gam);
%
        sy_bulk = (1+sqrt(gam))^2;
%
        sx1 = sy1/delta + 1 - 1/delta;
        sx_bulk = (1 + sqrt(gam))^2/delta + 1 - 1/delta;
%

        end
%
%
%
%
%
        function rlam_new = shrink_oper_full(rlam,gam)
%
%   the optimal eigenvalue shrinker for operator norm loss with full data 
%   (no missing entries)
%
        thresh = (1 + sqrt(gam))^2;

        if (rlam <= thresh)
           rlam_new = 1;
           return
        end
%
        rlam_new = compute_ell(rlam,gam);

        end
%
%
%
%
%
        function [xs,xs_clean,rnoise] = draw_gaussian(u,spec,m,k,n,sig)
%
%   draws a rank k gaussian blob with covariance rcov = u * spec * u',
%   and then adds to it white gaussian noise of variance sig^2
%
        xs_clean = u * diag(sqrt(spec)) * randn(k,n);
        rnoise = randn(m,n) * sig;
        xs = xs_clean + rnoise;
%

        end
%
%
%
%
%
        function [xs,xs_clean,rnoise,ichosen] = draw_discrete(aa,...
           probs,k,m,n,sig)
%
%   randomly selects columns of aa; i^th column is selected with
%   probability probs(i), and there are k+1 columns. A total of n points
%   are selected. Then Gaussian noise with variance sig^2 is added.
%
        ichosen = zeros(1,n);
        xs_clean = zeros(m,n);
%
        for ijk=1:n

        r = rand();

        for i=1:k+1
%
        if (r <= sum(probs(1:i)))
           xs_clean(1:m,ijk) = aa(1:m,i);
           ichosen(ijk) = i;
           break;
        end
    end
    end

        rnoise = randn(m,n);
        xs = xs_clean + sig*rnoise;

        end
%
%
%
%
%
        function [rcov,aa,u,s] = cov_discrete(rlam,k,m,probs)
%
%   This function generates k+1 random vectors in R^m with probs-weighted
%   mean zero and there covariance rcov, as well as the eigendecomposition
%   of the covariance. The MIMUMUM eigenvalue of rcov is specified as rlam.
%
%                           input parameters:
%
%   rlam - the minumum eigenvalue of the covariance 
%
%   k - the rank of the covariance
%
%   m - the dimensionality of the vectors
%
%   probs - the 1-by-(k+1) dimensional vector of probabilities
%
%
%                          output parameters:
%
%   rcov - the m-by-m covariance matrix
%
%   aa - the m-by-(k+1) matrix of vectors
%
%   u - the m-by-k matrix of eigenvectors of rcov
%
%   s - the k-dimensional vector of corresponding eigenvalues of rcov
%
%
        rcov = zeros(m,m);
        rmu = zeros(m,1);
%
%   generate k+1 random vectors in R^m and center them (so ps-average is 0)
%
        aa = randn(m,k+1);

        for i=1:k+1
%
        rmu = rmu + probs(i) * aa(1:m,i);
    end

%%%        chk0 = norm(aa*probs' - rmu)

        for i=1:k+1
%
        aa(1:m,i) = aa(1:m,i) - rmu;
    end

%%%        chk0 = norm(sum(aa*probs'))

        for i=1:k+1
%
        rcov = rcov + probs(i) * aa(1:m,i) * aa(1:m,i)';
    end

        rcov = .5*(rcov + rcov');
%%%        [u,s] = eigsort(rcov,m);
        [u,s] = eigs(rcov,k);
        s = diag(s);
%
%   rescale so minimum non-zero eigenvalue is rlam
%
        s_min = s(k);
        aa = aa / sqrt(s(k)) * sqrt(rlam);
        rcov = rcov / s_min * rlam;
        rcov = .5*(rcov + rcov');
%%%        [u,s] = eigsort(rcov,m)
        s = s / s_min * rlam;


        end
%
%
%
%
%
        function [rcov,u] = rand_cov3(m,k,spec)
%
%   This function generates an m-by-m covariance matrix of rank k with 
%   user-specified spectrum.
%
%
%                         input parameters:
%
%   m - the dimensionality of the matrix
%
%   k - the rank of the matrix
%
%   spec - the vector of length k containing the spectrum 
%
%
%                        output parameters:
%
%   rcov - the m-by-m covariance matrix, with random eigenvectors
%
%   u - the m-by-k matrix of eigenvectors of rcov (random orthonormal set)
%

        u = orthoset_small(m,k);
        rcov = u * diag(spec) * u';
        rcov = .5 * (rcov + rcov');

        end
%
%
%
%
%
        function u = orthoset_small(m,k)
%
%   The function produces an m-by-k matrix u with orthonormal columns.
%
        u = zeros(m,k);
        g = randn(m,k);
%
        u(1:m,1) = g(1:m,1) / norm(g(1:m,1));

        for ijk = 2:k
%
        vec = g(1:m,ijk);

        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

%%%        chk0 = sum(vec .* vec2)

        u(1:m,ijk) = vec/norm(vec);

    end

    end


        end
%
%
%
%
%
        function [rcov,u,spec] = rand_cov_hada(m,k,spec)
%
%   The functions returns an m-dimensional covariance matrix of rank k
%   with specified spectrum whose eigenvectors are randomly selected
%   columns of the hadamard matrix. The point is that the diagonal elements
%   are all equal.
%
        u = hadam_sub(m,k);
        rcov = u * diag(spec) * u';
        rcov = .5 * (rcov + rcov');

        end
%
%
%
%
%
        function hsub = hadam_sub(m,ncols)
%
%   hsub has ncols random columns of the 2^k by 2^k hadamard matrix;
%   so its columns are orthonormal, and every coordinate has equal abs val
%
        hsub=zeros(m,ncols);
        iperm = randperm(m);
        inds = iperm(1:ncols);
%
        for i=1:m
%
        for j=1:ncols
%
        jj=iperm(j);
        ijdot = bitdot(i-1,jj-1);
        hsub(i,j) = (-1)^ijdot;
    end
    end

        hsub = hsub/sqrt(m);
        chk0 = norm(hsub'*hsub - eye(ncols),'fro')

        end
%
%
%
%
%
        function hmat = hadam_transf(m)
%
        hmat = zeros(m,m);
%
        for i=1:m
%
        for j=1:m
%
        ijdot = bitdot(i-1,j-1);
        hmat(i,j) = (-1)^ijdot;
    end
    end

        hmat = hmat/sqrt(m);

        end
%
%
%
%
%
        function kdot = bitdot(m,n)
%
        [kexpm,lm] = binexp(m);
        [kexpn,ln] = binexp(n);
        len = min([lm,ln]);
        kprod = kexpm(1:len) .* kexpn(1:len);
        kdot = sum(kprod);

        end
%
%
%
%
%
        function [kexp,len] = binexp(m)
%
        if (m==0)
           len=1;
           kexp=0;
           return;
        end
%
        if (m==1)
           len=1;
           kexp=1;
           return;
        end

        len = floor(log2(m))+1;
        kexp = zeros(1,len);
        m2 = m;
%
        for ijk = 1:len+10
%
        ii=floor(log2(m2));
        kexp(ii+1) = 1;
        m2 = m2 - 2^ii;

        if(m2 == 0)
          break;
        end

    end

        end
%
%
%
%
%
        function [xhat_delta,xhat_nijs] = cov_shrink_first(ys,m,n,...
           delta_hat,gam,nijs)
%
%   This function first shrinks the sample covariance of the projected
%   data, and then inverts the linear systems (both delta and nijs systems)
%   to obtain the final estimates of the covariance of the x's
%
%
%                         input parameters:
%
%   ys - the m-by-n matrix of observations, with zeroes in the missing
%      entries
%   m,n - the dimensionality and number of samples, respectively, of data
%
%   delta_hat - the estimated value of delta, the sampling probability
%
%   gam - the ratio of m/n
%
%   nijs - the m-by-m matrix of mutual counts 
%
%
%                       output parameters:
%
%   xhat_delta - the estimated m-by-m covariance obtained by inverting the
%      delta system
%   xhat_nijs - the estimated m-by-m covariance obtained by inverting the
%      nijs system
%
%
        yhat = ys*ys' / n;
        b = yhat / delta_hat;
        [vb,sb] = eigsort(b,m);
%%%        chk0 = norm(b*vb - vb*diag(sb),'fro')

        sb2 = 0*sb;
        for i=1:m

%%%        sb2(i) = shrink_spec(sb(i),gam);
%%%        sb2(i) = shrink_oper_full(sb(i),gam);
        sb2(i) = cov_emp2oper(sb(i),1,gam);
        sb2(i) = sb2(i)-1;
    end
%%%        figure;plot(sb,'*')
%%%        hold on;plot(sb2,'*')
%%%        stopnow
        b2 = vb*diag(sb2)*vb' * delta_hat;

%
%   invert systems to get estimators for x's covariance
%
        xhat_nijs = yhat2xhat_nijs(b2,nijs,m,n);
        xhat_delta = adjinv_mat(b2,delta_hat,m);

        end
%
%
%
%
%
        function xhat = yhat2xhat_nijs(yhat,nijs,m,n)
%
        xhat = zeros(m,m);
        iok = find(nijs > 0);
        xhat(iok) = n*yhat(iok) ./ nijs(iok);

        end
%
%
%
%
%
        function xhat = yhat2xhat_delta(yhat,delta,m)
%
        yhat = yhat/delta;
        ydiag = diag(yhat);
        xhat = yhat + (delta-1)*diag(ydiag);
        xhat = xhat/delta;

        end
%
%
%
%
%
        function x2 = softhresh99(x,thresh)
%
        if (x <= thresh)
            x2 = 0;
            return;
        end
%
        x2 = x - thresh - 1;
        
        end
%
%
%
%
%
        function [sxhat,sxhat2,xhat,xhat2] = shrinkmat_soft(...
           ys,m,n,delta,gam,tconst)       
        sxhat2 = zeros(1,m);
%
%   take sample covariance of the Y vectors and rescale so noise is white
%
        yhat = ys*ys' / (n*delta);
        ydiag = diag(yhat);
%
        xhat = yhat + (delta-1)*diag(ydiag);
        xhat = xhat/delta;
%
%%%        [uxhat,sxhat] = eigs(xhat,kvecs);
%%%        sxhat = diag(sxhat);
        [uxhat,sxhat] = eigsort(xhat,m);
        chk0 = norm(xhat*uxhat - uxhat*diag(sxhat));
%
%   apply non-linearity to each eigenvalue
%
        xnorm = norm(xhat);
        xtrace = sum(diag(xhat));
        thresh = tconst*sqrt(xtrace*xnorm*log(2*m)/n) / delta;
        for i=1:m
%
        sxhat2(i) = softhresh99(sxhat(i),thresh);
    end

        xhat2 = uxhat * diag(sxhat2) * uxhat';

        end
%
%
%
%
%
        function rshrunk = cov_emp2oper(rlam,delta,gam)
%
%   applies the asymptotically optimal nonlinearity for operator norm loss
%   to rlam in case of data missing with probability delta for m/n = gamma
%
        rlam2 = delta*rlam + 1 - delta;

        bulk_edge = (1+sqrt(gam))^2;
        if (rlam2 <= bulk_edge)
%
        rshrunk = 1;
        ell=1;
        return
    end

        ell = compute_ell(rlam2,gam);
        rshrunk = (ell - 1 + delta) / delta;

%%%        rshrunk = ell/delta;

        return
%
%   perform sanity check:
%
        aaa=compute_ell2(rlam,delta,gam);
        chk0 = aaa-rshrunk

        end
%
%
%
%
%
        function u = orthoset_mod_small(x,m,k)
%
%   produces an m-by-k matrix u with orthonormal columns
%
        u = zeros(m,k);
        u(1:m,1) = x(1:m,1) / norm(x(1:m,1));

        for ijk = 2:k
%
        vec = x(1:m,ijk);
        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);
%
        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;
        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;
%
        u(1:m,ijk) = vec/norm(vec);

    end

    end


        end
%
%
%
%
%
        function adjinv = adjinv_mat2(amat,deltas,m)
%
        adjinv = amat;
        for i=1:m
%
        adjinv(i,i) = adjinv(i,i) * deltas(i);
    end

        for j=1:m
%
        for i=1:m
%
        adjinv(i,j) = adjinv(i,j) / deltas(i) / deltas(j);
    end
%
    end

        end
%
%
%
%
%
        function adjmat = adjust_mat2(amat,deltas,m)
%
%   takes covariance matrix amat and returns covariance adjmat for
%   subsampled data, with subsampling probabilities given by deltas
%
        adjmat = amat;
        for j=1:m
%
        for i=1:m
%
        adjmat(i,j) = amat(i,j) * deltas(i) * deltas(j);
    end
%
    end

        for i=1:m
%
        adjmat(i,i) = adjmat(i,i)/deltas(i);
    end

        return
%
%   to do it in vectorized form:
%
        adjmat2 = diag(deltas) * amat * diag(deltas) + ...
           (diag(deltas - deltas.^2))*diag(diag(amat));

        end
%
%
%
%
%
        function adjmat = adjust_mat(amat,delta,m)
%
%   takes covariance matrix amat and returns covariance adjmat for
%   subsampled data, with subsampling probability given by delta; note
%   that each coordinate has the SAME probability delta of being sampled
%
        adjmat = delta^2 * amat;

        for i=1:m
%
        adjmat(i,i) = adjmat(i,i) / delta;
    end

        end
%
%
%
%
%
        function adjinv = adjinv_mat(amat,delta,m)
%
        adjinv = amat;
        for i=1:m
%
        adjinv(i,i) = delta * adjinv(i,i);
    end

        adjinv = adjinv / delta^2;

        end
%
%
%
%
%
        function u = orthoset_mod(x,m)
%
%   produces an m-by-k matrix u with orthonormal columns
%
        u = zeros(m,m);
        u(1:m,1) = x(1:m,1) / norm(x(1:m,1));

        for ijk = 2:m
%
        vec = x(1:m,ijk);
        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);
%
        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;
        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;
%
        u(1:m,ijk) = vec/norm(vec);

    end

    end


        end
%
%
%
%
%
        function [u,s] = eigsort(a,m)
%
        [u,s] = eig(a);
        s = diag(s);
        [s,ii] = sort(s,'descend');
        u = u(1:m,ii);

%%%        chk0 = norm(u*diag(s)*u' - a,'fro') / norm(a,'fro')
        end
%
%
%
%
%
        function ess = compute_ess(cee)
%
        ess = sqrt( 1 - cee^2);

        end
%
%
%
%
%
        function cee = compute_cee(ell,gam)
%
        rnum = 1 - gam/(ell - 1)^2;
        rden = 1 + gam/(ell - 1);
        cee = sqrt( rnum / rden );

        end
%
%
%
%
%
        function ell = compute_ell(rlam,gam)
%
        t1 = rlam + 1 - gam;
        t2 = sqrt( t1^2 - 4*rlam );

        ell = .5 * (t1 + t2);
        

        end
%
%
%
%
%
        function rlam = compute_lam(ell,gam)
%
        rlam = ell * (1 + gam / (ell - 1));

        end
%
%
%
%
%
        function inds = rand_inds_fast2(m,n,ps)
%
        ff = rand(m,n);
        inds = real(ff <= ps);

        end
%
%
%
%
%
        function [inds,ninds,indsc] = rand_inds_fast(m,n,ps)
%
        ff = rand(m,n);
        inds = real(ff <= ps);
        ninds = sum(inds(:));

        indsc = ones(m,n) - inds;

        end
%
%
%
%
%
        function [inds,ninds,indsc,ifbad] = rand_inds(m,n,ps)

        inds = zeros(m,n);
        irows = zeros(1,m);
        icols = zeros(1,n);
        
        ninds = 0;

        for j = 1:n
 
        for i = 1:m

        r = rand();

        if (r < ps(i,j))
           inds(i,j) = 1;
           ninds = ninds + 1;
           icols(j) = icols(j) + 1;
           irows(i) = irows(i) + 1;
        end
    end

    end


        indsc = ones(m,n) - inds;


        mcounts = inds * inds';

        ifbad = 0;

        chkmin = min(min(mcounts));
        if (chkmin <= 1)
           ifbad = 1
        end


mean(mean(inds))
mean(ps(:))

        end
%
%
%
%
%
        function startnow
%
        delete out13
        diary('out13')
        diary on
%
%%%        format shortE
        format longE


%
        rng('default');

        end
%
%
%
%
%
        function stopnow
%
        diary off
        stop

        end
