choose_vars       = [12 7 2 3 4 5 9 10 8 11];
[data,vmed,vars] = est_irf_fev(nirf,bound,ndraws, nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% ==============================================================================================================================================================================================
%----------------------------------------------------------------------
% Plot the responses of KPSS and TFP adjusted. 
%----------------------------------------------------------------------
fonttype          ='Arial';
ftsizeaxis        = 11;
titlefontsize     = 10;
%data = data*100;
[nvars,temp]=size(vars);
N = nvars;
[nirf,temp,temp]=size(data);
R=round(N/2);

%% Figure 2 of the replication
a=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG5_CGV_RESTAT');   


figure(a)
for n=1:2 
        subplot(1,2,n);            
        plot(1:nirf,data(:,nvars+n,1), 'LineWidth',2);
        hold on                                        
        plot(1:nirf,data(:,nvars+n,2),'r--','LineWidth',2);
        hold on                            
        plot(1:nirf,data(:,nvars+n,3),'r--','LineWidth',2);
        hold on
        plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
        ylabel('percent','FontSize',12)
        xlabel('Horizon (quarters)','FontSize',12)
        title(num2str(vars(n,:)),'Interpreter','tex','FontSize',titlefontsize)
        
        set(gca,'XTick',[0;4;8;12;16;20])
        set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
        set(gca, 'FontName', fonttype);
        set(gca, 'FontSize', ftsizeaxis);
        set(gca,'Layer','top');

        box off 
        axis tight
        hold off;                                                                 
end

       
%% Figure 3 of the replication
b=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG6_CGV_RESTAT');   
figure(b);
for n=3:N
        subplot(R-1,2,n-2);            
        plot(1:nirf,data(:,nvars+n,1), 'LineWidth',2); hold on                                        
        plot(1:nirf,data(:,nvars+n,2),'r--','LineWidth',2);hold on                            
        plot(1:nirf,data(:,nvars+n,3),'r--','LineWidth',2);hold on
        plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
        ylabel('percent','FontSize',12)
        xlabel('Horizon (quarters)','FontSize',12)
        title(num2str(vars(n,:)),'Interpreter','tex','FontSize',titlefontsize) %'FontWeight','bold',

        set(gca,'XTick',[0;4;8;12;16;20])
        set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
        set(gca, 'FontName', fonttype);
        set(gca, 'FontSize', ftsizeaxis);
        set(gca,'Layer','top');
        box off 
        axis tight
        hold off;                               
end
