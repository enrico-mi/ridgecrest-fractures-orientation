function SRw_plots_magnfactor_NW()
    
    addpath( fullfile('.','computes') );
    addpath( fullfile('.','plots') );
    
    close all;

    sigmaH_angles = 14;

    mus = {};
    for j=1:6
        mus{j} = (j-1)*0.3;
    end
    mus{7} = 1e6;
    mus{8} = "rand";
    mu_idx = [1:3:16 17 18];

    mu_leg = ["\mu=0.0","\mu=0.3","\mu=0.6","\mu=0.9","\mu=1.2","\mu=1.5","\mu=inf","random \mu"];

    magn_factors = 0.:0.1:1;
    
    SRw = zeros(length(mus),length(sigmaH_angles),length(magn_factors));
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
            
            for mf = 1:length(magn_factors)
                SRw(:,sa,mf) = getSRw(string_stress_angles,mu_idx,mf,slip_model_name,realization,randang,shear_ext,dist_ext);
                
            end
        end
    end
    
    xLimits = [0 1];

    maxSRw = max(SRw,[],'all');
    minSRw = min(SRw,[],'all');
    rangeSRw = maxSRw-minSRw;
    yLimits = [0 40];%[minSRw-0.1*rangeSRw maxSRw+0.1*rangeSRw];
    yTicks = 0:10:40;
    
    xLabel = 'k [-]';
    yLabel = 'SR_w [\circ]';
    
    nb_plots = length(sigmaH_angles);
    t = tiledlayout(nb_plots,1);
    ax = gobjects(1,nb_plots);
    
    sgtitle('SRw of \theta_{r}-\theta_{H} with coseismic stresses intensity (L>100m, depth = 100m)');

    for sa = 1:nb_plots
        ax(sa) = nexttile;
        for m = 1:length(mus)
            
            y_plot = squeeze(SRw(m,sa,:))';
            p = plot(ax(sa), magn_factors, y_plot, 'Color', colour(m), 'LineWidth', 1.5, DisplayName=mu_leg(m));
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
                x_fill = [magn_factors, fliplr(magn_factors)];
                y_fill = [y_plot, ones(size(y_plot))*yLimits(2)];
                hf = fill(x_fill, y_fill, [0.7, 0.7, 0.7], 'FaceAlpha', 0.5,'LineStyle','none');
                uistack(hf, 'bottom');
            end            

        end
        xlim(xLimits);
        ylim(yLimits);
        yticks(yTicks);
            
        h = gca;
        h.YAxis.MinorTick = 'On';
        h.YAxis.MinorTickValues = 0:2:40;

        % displaying and re-ordering legend
        set(legend,'Location', 'eastoutside');
        lh = findobj(gcf, 'Type', 'legend');
        lh.AutoUpdate = 'off';
        lh.PlotChildren = lh.PlotChildren([1, 2, 3, 4, 5, 8, 6, 7, 9]);
        lh.PlotChildren(1) = [];
        lh.PlotChildren(end) = [];

        ylabel(yLabel);
        title('N14E');
        xlabel(xLabel);
    
    end

    figurePosition = [1,121,960,960/2];
    set(gcf, 'Position', figurePosition);
    fontsize(gcf,12,'points');

    figname = "SRw_100m_N14E_tau_magnitude_review";
    savefig(figname+".fig");
    print('-vector', '-depsc', figname+".eps");
    print('-dpng', figname+".png", '-r600');
            
end