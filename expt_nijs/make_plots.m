        function make_plots()
%
%   plot the data returned by test_limits code
%
        load data200.mat
        set(groot,'defaultLineLineWidth','remove')
        set(groot,'defaultLineLineWidth',2)


        ndeltas
        ngams

        subplot(1,2,1)

        igam = 1;
        gam1 = gams(igam)
%
        evmins = 1 + (sqrt(gam1) + 1)./deltas + 5
        evs_noisy = zeros(2,ndeltas)

        for j=1:ndeltas
%
        evs_noisy(1:2,j) = evmins(j) + [k-1:-1:0]
    end
        ells = evs_noisy-1
        ss1 = sum(ells,1)

        vals1 = err_op_means(igam,:) ./ ss1;
        hold on; plot(deltas,vals1,'b-.','LineWidth',2)

        igam = 5;
        gam2 = gams(igam)
%
        evmins = 1 + (sqrt(gam1) + 1)./deltas + 5
        evs_noisy = zeros(2,ndeltas)

        for j=1:ndeltas
%
        evs_noisy(1:2,j) = evmins(j) + [k-1:-1:0]
    end
        ells = evs_noisy-1
        ss2 = sum(ells,1)

        vals2 = err_op_means(igam,:) ./ ss2;
        hold on; plot(deltas,vals2,'r--','LineWidth',2)



        igam = 10;
        gam3 = gams(igam)
%        evmins = 1 + (sqrt(gam1) + 1)./deltas + 5
        evs_noisy = zeros(2,ndeltas)

        for j=1:ndeltas
%
        evs_noisy(1:2,j) = evmins(j) + [k-1:-1:0]
    end
        ells = evs_noisy-1
        ss3 = sum(ells,1)

        vals3 = err_op_means(igam,:) ./ ss3;
        hold on; plot(deltas,vals3,'g','LineWidth',2)

        legend('\gamma=.1','\gamma=.45','\gamma=.9')

        xlim([.09,1])
        xlabel('\delta')
        ylabel('difference in solutions')

%
%
%

        subplot(1,2,2)

        igam = 1;
        gam1 = gams(igam)
%
        vals1 = err_fr_means(igam,:) ./ ss1;
        hold on; plot(deltas,vals1,'b-.','LineWidth',2)


        igam = 5;
        gam2 = gams(igam)
%
        vals2 = err_fr_means(igam,:) ./ ss2;
        hold on; plot(deltas,vals2,'r--','LineWidth',2)


        igam = 10;
        gam3 = gams(igam)
%
        vals3 = err_fr_means(igam,:) ./ ss3;
        hold on; plot(deltas,vals3,'g','LineWidth',2)

        legend('\gamma=.1','\gamma=.45','\gamma=.9')

        xlim([.09,1])
        xlabel('\delta')
        ylabel('difference in solutions')


        set(figure(1),'Position',[500,500,1075,400])
        savefig(figure(1),'fig200')
