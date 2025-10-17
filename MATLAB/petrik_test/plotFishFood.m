function [] = plotFishFood(fig, X, Y, Fish, Food, current, it, dt, titleString )
    % Plot food concentration
    %ax1 = subplot(1,2,1);
    %ax2 = subplot(1,2,2);
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    t = tiledlayout(2,2);

    % Fish
    ax1 = nexttile;
    img1 = pcolor(X, Y, Fish);
    img1.FaceColor='interp';
    img1.EdgeColor='none';
    cmap1 = cmocean('dens');
    colormap(ax1,cmap1); cb1 = colorbar;
    ylabel(cb1, '$gww/m^3$','FontSize',14,'interpreter','latex');
    set(cb1, 'TickLabelInterpreter','latex');
    axis tight
    pbaspect([2 1 1])
    xlabel(ax1,'x [km]', 'interpreter','latex');
    ylabel(ax1,'y [km]', 'interpreter','latex');
    title(ax1,'Fish Concentration in $gww/m^3$', 'interpreter','latex');
    caxis([0 0.1]);

    % Food
    ax2 = nexttile;
    img2 = pcolor(X,Y,Food);
    img2.FaceColor='interp';
    img2.EdgeColor='none';
    cmap2 = cmocean('algae');
    cb2 = colorbar;
    ylabel(cb2, '$gww/m^3$','FontSize',14,'interpreter','latex');
    set(cb2, 'TickLabelInterpreter','latex');
    axis tight
    pbaspect([2 1 1])
    xlabel(ax2,'x [km]', 'interpreter','latex');
    ylabel(ax2,'y [km]', 'interpreter','latex');
    title(ax2,'Food Concentration in $gww/m^3$', 'interpreter','latex');

    % Fish direction vectors
    ax3 = nexttile([1 2]);
    [n m o] = size(current);
    if o > 2

        u_holder = current(:,:,3) - current(:,:,4);
        v_holder= current(:,:,2) - current(:,:,1);
        current(:,:,1) = u_holder;
        current(:,:,2) = v_holder;
    end
    xvel = squeeze(current(:,:,1));
    yvel = squeeze(current(:,:,2));
    % Scaling the vectors
    % magnitude = sqrt(xvel.^2+yvel.^2);
    % maxMag = max(magnitude(:));
    % minMag = min(magnitude(:));
    % xvel_scaled = 0.1 + (xvel - minMag) / (maxMag - minMag) * (1 - 0.1);
    % yvel_scaled = 0.1 + (yvel - minMag) / (maxMag - minMag) * (1 - 0.1);    
    % xvel_scaled = 1000 .* Fish .* xvel;
    % yvel_scaled = 1000 .* Fish .* yvel;
    xvel(xvel > 0) = 1; xvel(xvel < 0) = -1;
    yvel(yvel > 0) = 1; yvel(yvel < 0) = -1;
    % Plotting
    img3 = quiver(X, Y, xvel, yvel, 1);
    axis tight
    pbaspect([5 1 1])
    xlabel(ax3,'x [km]', 'interpreter','latex');
    ylabel(ax3,'y [km]', 'interpreter','latex');
    title(ax3,'Overall Current (Vc+Vf) Direction Vectors', 'interpreter','latex');


    % Overall Title
    title(t,[ titleString, ', Day = ' num2str(it*dt, '%2.2f'), ...
    ', numfish=', num2str( sum(Fish(:)) ), ...
    ', numfood=', num2str( sum(Food(:))) ], 'interpreter','latex', ...
    'FontSize',20);

    % xlabel(t,'x [km]');
    % ylabel(t,'y [km]');
    % %zlabel('z');
    % title(t,[ titleString, ', Day = ' num2str(it*dt, '%4.2f'), ...
    % ', numfish=', num2str( sum(Fish(:)) ), ...
    % ', numfood=', num2str( sum(Food(:))) ]);
    %colorbar
    %clim([0 0.7]);
    %zlim([0 1.2]);
    drawnow;

    %fig = figure(); 
    % tcl = tiledlayout(fig,1,1); 
    % 
    % ax = nexttile(tcl); 
    % 
    % surf1Colormap = cmocean('dens'); 
    % tc = colormapToTruecolor(surf1Colormap,Fish);
    % hsurf1 = surf(ax,X,Y,Fish,tc,'FaceColor','interp','EdgeAlpha',0.3);
    % clim([0 1]);
    % % create second surface using colors from the autumn colormap
    % hold on
    % surf2Colormap = cmocean('algae'); 
    % tc = colormapToTruecolor(surf2Colormap,Food);
    % hsurf2 = surf(ax,X,Y,Food,tc,'FaceColor','interp','EdgeAlpha',0.3);
    % clim([0 1]);
    % %rotate(hsurf2,[0 1 0],45)
    % %axis equal
    % view([-3,7])
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')
    % title(['Fish concentration profile, iter = ' num2str(it, '%2.2f'), ...
    % ', numfish=', num2str( sum(Fish(:)) )]);
    % zlim([0 1.2]);
    % 
    % % Assign colormaps (which wont' affect the surfaces) and add colorbars
    % % A second, hidden axes is needed to host the second colormap.
    % axHidden = axes(tcl,'visible','off','HandleVisibility','off');
    % colormap(ax,surf1Colormap)
    % colormap(axHidden,surf2Colormap)
    % cb1 = colorbar(ax);
    % cb1.Layout.Tile = 'east';
    % cb1.Label.String = 'Fish Conc.';
    % cb2 = colorbar(axHidden);
    % cb2.Layout.Tile = 'east';
    % cb2.Label.String = 'Food Conc.';
    % drawnow; 

end

function tc = colormapToTruecolor(map,ZData)
% map is a n-by-3 colormap matrix
% ZData is a k-by-w matrix of ZData (or CData, I suppose)
% tc is a kxwx3 Truecolor array based on map and ZData values.
tcIdx = round(rescale(ZData,1,height(map)));
tc = reshape(map(tcIdx,:),[size(ZData),3]);
end