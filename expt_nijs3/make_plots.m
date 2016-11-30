        function make_plots()
%
        load data203.mat
        set(groot,'defaultLineLineWidth','remove')
        set(groot,'defaultLineLineWidth',2)

        ndeltas
        n_ns

        figure(1);

        idel = 2;
        del1 = deltas(idel)

        vals1 = err_op_means(:,idel)
        hold on; plot(ns,vals1,'b-.','LineWidth',2)

%

        idel = 3;
        del2 = deltas(idel)

        vals2 = err_op_means(:,idel)
        hold on; plot(ns,vals2,'r--','LineWidth',2)

%

        idel = 4;
        del3 = deltas(idel)

        vals3 = err_op_means(:,idel)
        hold on; plot(ns,vals3,'g','LineWidth',2)



        legend('\delta=.2','\delta=.3','\delta=.4')

        xlabel('$n$','Interpreter','latex')
        ylabel('difference in solutions')

        set(figure(1),'Position',[500,500,850,700])
        savefig(figure(1),'fig203')

        del1
        del2


        end