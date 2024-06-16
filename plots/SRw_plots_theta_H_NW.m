function SRw_plots_theta_H_NW()
    
    close all;

    sigmaH_angles = -10:1:25;
    
    magn_factors_idx = 1; %looking only for regional stress, k=0;
    
    mus = {};
    for j=1:6
        mus{j} = (j-1)*0.3;
    end
    mus{7} = 1e6;
    mus{8} = "rand";
    mu_idx = [1:3:16 17 18];
    
    mu_leg = ["\mu=0.0","\mu=0.3","\mu=0.6","\mu=0.9","\mu=1.2","\mu=1.5","\mu=inf","random \mu"];

    SRw = zeros(length(mus),length(sigmaH_angles));
    
    slip_models = ["JB2"];%;"RB2";"XB2"];
    
    % looking only at no uncertainty case for now
    realization = 1;
    randang = 0;
    shear_ext = "";
    dist_ext = "";
    
    for sm = 1:1
        slip_model_name = slip_models(sm);
        for sa = 1:length(sigmaH_angles)
            string_stress_angles = num2str(sigmaH_angles(sa),'%02.f');
            SRw(:,sa) = getSRw(string_stress_angles,mu_idx,magn_factors_idx,slip_model_name,realization,randang,shear_ext,dist_ext);
        end
    end
    
    xLimits = [min(sigmaH_angles) max(sigmaH_angles)];

    maxSRw = max(SRw,[],'all');
    minSRw = min(SRw,[],'all');
    rangeSRw = maxSRw-minSRw;
    yLimits = [0 40];%[minSRw-0.1*rangeSRw maxSRw+0.1*rangeSRw];
    yTicks = 0:10:40;
    
    xLabel = '\theta_{H,reg} [\circ]';
    yLabel = 'SR_w [\circ]';
    
    nb_plots = length(magn_factors_idx);
    t = tiledlayout(nb_plots,1);
    ax = gobjects(1,nb_plots);
    
    sgtitle('SRw of \theta_{r}-\theta_{H} as a function of regional stress orientation (any L)');


    for d = 1:nb_plots
        ax(d) = nexttile;
        
        xline(0,'LineStyle','--','Color','k', 'LineWidth', 1.5);
        hold on
        xline(15,'LineStyle','--','Color','k', 'LineWidth', 1.5);
        xline(14,'LineStyle',':','Color','k', 'LineWidth', 1.5);
        for m = 1:length(mus)

            y_plot = SRw(m,:);
            p = plot(ax(d), sigmaH_angles, y_plot, 'Color', colour(m), 'LineWidth', 1.5, DisplayName=mu_leg(m));
            p.LineStyle = '-';
            p.Marker = 'x';
            p.MarkerSize = 11;
            %p.MarkerFaceColor = colour(m);
            hold on
            
            if (isfloat(mus{m}))
                if (mus{m} == 0.6)
                    p.Marker = '.';
                    p.MarkerSize = 20;
                else
                    uistack(p, 'bottom');
                end
            elseif (strcmp(mus{m},"rand"))
                p.Marker = 'none';
                p.LineStyle = ':';
                p.Color = 'k';
                x_fill = [sigmaH_angles, fliplr(sigmaH_angles)];
                y_fill = [y_plot, ones(size(y_plot))*yLimits(2)];
                hf = fill(x_fill, y_fill, [0.7, 0.7, 0.7], 'FaceAlpha', 0.5,'LineStyle','none');
                uistack(hf, 'bottom');
            end
            
            xlim(xLimits);
            ylim(yLimits);
            yticks(yTicks);

            h = gca;
            h.YAxis.MinorTick = 'On';
            h.YAxis.MinorTickValues = 0:2:40;

        end

        set(legend,'Location', 'eastoutside');
        lh = findobj(gcf, 'Type', 'legend');
        lh.AutoUpdate = 'off';
        lh.PlotChildren(8:10) = [];
        lh.PlotChildren = lh.PlotChildren([1, 2, 3, 4, 5, 8, 6, 7, 9]);
        lh.PlotChildren(1) = [];
        lh.PlotChildren(end) = [];

        title('100m');
        ylabel(yLabel);
        xlabel(xLabel);
    
    end

    figurePosition = [1,121,960,960/2];
    set(gcf, 'Position', figurePosition);
    fontsize(gcf,12,'points');

    figname = "SRw_allL_reg_orientation_review";
    savefig(figname+".fig");
    print('-vector', '-depsc', figname+".eps");
    print('-dpng', figname+".png", '-r600');
            
end
