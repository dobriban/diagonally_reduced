        function make_plots()
%
%   plot the data returned by test_limits code
%
        load data201.mat
        set(groot,'defaultLineLineWidth','remove')
        set(groot,'defaultLineLineWidth',2)


        ndels
        dels




        figure(1)

        subplot(1,2,1)
%
        idel = 1
        vals1 = evals(idel,:);

        del1 = dels(idel);

        hold on; plot(vals1,err_mean_shr(idel,:),'b-.','LineWidth',2)
        hold on; plot(vals1,err_mean_opt(idel,:),'r--','LineWidth',2)

        legend('EBLP error','OptSpace error')
        xlabel('$\ell$','Interpreter','latex')
        ylabel('average error')
        title('\delta = .1')
%

        subplot(1,2,2)
%
        idel = 7
        vals2 = evals(idel,:);

        del2 = dels(idel)


        hold on; plot(vals1,err_mean_shr(idel,:),'b-.','LineWidth',2)
        hold on; plot(vals1,err_mean_opt(idel,:),'r--','LineWidth',2)

        legend('EBLP error','OptSpace error')
        xlabel('$\ell$','Interpreter','latex')
        ylabel('average error')
        title('\delta = .7')

%


        del1
        del2

        set(figure(1),'Position',[500,500,1075,400])
        savefig(figure(1),'fig201b')



