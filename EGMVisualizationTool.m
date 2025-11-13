% Generate a synthetic test matrix for threebody_view
t = linspace(0, 20, 500)'; % time from 0 to 20
n = length(t);

% Body 1: small circular orbit in XY plane
r1 = [cos(t), sin(t), 0.1*sin(2*t)];
v1 = [-sin(t), cos(t), 0.2*cos(2*t)];

% Body 2: larger, slower orbit tilted in XZ plane
r2 = [2*cos(0.5*t), 0.5*sin(0.5*t), 2*sin(0.5*t)];
v2 = [-sin(0.5*t), 0.5*cos(0.5*t), cos(0.5*t)];

% Body 3: elliptical orbit in XY plane
r3 = [1.5*cos(0.7*t), 0.8*sin(0.7*t), 0.3*cos(0.5*t)];
v3 = [-1.5*0.7*sin(0.7*t), 0.8*0.7*cos(0.7*t), -0.3*0.5*sin(0.5*t)];

% Combine into full matrix: [t, body1(6), body2(6), body3(6)]
M = [t, r1, v1, r2, v2, r3, v3];

% Test
threebody_view(M);

function threebody_view(M)
% THREEBODY_VIEW - Clean interactive 3-body visualizer with visible slider & state panels

    if nargin < 1 || size(M,2) ~= 19
        error('Input must be [N x 19].');
    end
    N = size(M,1);
    t  = M(:,1);
    r1 = M(:,2:4);  v1 = M(:,5:7);
    r2 = M(:,8:10); v2 = M(:,11:13);
    r3 = M(:,14:16);v3 = M(:,17:19);
    tailLen = max(10, round(0.05 * N));

    %% --- FIGURE SETUP ---
    f = uifigure('Name','3-Body Viewer','Color','k','Position',[100 100 1200 750]);

    % ---- 3D AXES ----
    ax = uiaxes(f);
    ax.Position = [380 120 780 580]; % fixed pixel placement
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    enableDefaultInteractivity(ax);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    ax.XColor=[0.7 0.7 0.7]; ax.YColor=[0.7 0.7 0.7]; ax.ZColor=[0.7 0.7 0.7];
    ax.GridColor=[0.25 0.25 0.25]; ax.GridAlpha=0.5;
    ax.MinorGridColor=[0.15 0.15 0.15]; ax.MinorGridAlpha=0.4;
    title(ax,'3-Body Simulation','Color',[0.95 0.95 0.95]);
    view(ax,35,20);

    % ---- AXIS LIMITS ----
    allR=[r1; r2; r3];
    mins=min(allR,[],1); maxs=max(allR,[],1);
    span=max(maxs-mins); if span<=0, span=1; end
    ctr=(mins+maxs)/2; mrg=0.6;
    xlim(ax,[ctr(1)-mrg*span,ctr(1)+mrg*span]);
    ylim(ax,[ctr(2)-mrg*span,ctr(2)+mrg*span]);
    zlim(ax,[ctr(3)-mrg*span,ctr(3)+mrg*span]);

    %% --- VISUAL STYLE ---
    C1=[1.0 0.45 0.2];
    C2=[0.25 0.7 1.0];
    C3=[0.4 1.0 0.45];

    hT1 = plot3(ax,nan,nan,nan,'-','Color',C1*0.8,'LineWidth',2);
    hT2 = plot3(ax,nan,nan,nan,'-','Color',C2*0.8,'LineWidth',2);
    hT3 = plot3(ax,nan,nan,nan,'-','Color',C3*0.8,'LineWidth',2);
    hP1 = scatter3(ax,r1(1,1),r1(1,2),r1(1,3),70,C1,'filled','MarkerEdgeColor','w');
    hP2 = scatter3(ax,r2(1,1),r2(1,2),r2(1,3),70,C2,'filled','MarkerEdgeColor','w');
    hP3 = scatter3(ax,r3(1,1),r3(1,2),r3(1,3),70,C3,'filled','MarkerEdgeColor','w');
    legend(ax,[hP1 hP2 hP3],{'Body 1','Body 2','Body 3'},...
        'TextColor',[0.9 0.9 0.9],'Location','northeast','Color','k');

    %% --- LABELS & PANELS ---
    uilabel(f,'Text','State Vectors','FontColor',[0.9 0.9 0.9],...
        'FontWeight','bold','BackgroundColor','k','Position',[30 675 300 25],'FontSize',13);

    txt1 = uitextarea(f,'Position',[30 580 320 85],'Editable','off',...
        'FontName','Consolas','FontSize',11,'BackgroundColor','k','FontColor',C1);
    txt2 = uitextarea(f,'Position',[30 470 320 85],'Editable','off',...
        'FontName','Consolas','FontSize',11,'BackgroundColor','k','FontColor',C2);
    txt3 = uitextarea(f,'Position',[30 360 320 85],'Editable','off',...
        'FontName','Consolas','FontSize',11,'BackgroundColor','k','FontColor',C3);

    % Time label centered below slider
    lbl = uilabel(f,'Position',[500 65 250 22],...
        'Text',sprintf('t = %.3g',t(1)),...
        'FontColor',[0.9 0.9 0.9],'BackgroundColor','k',...
        'HorizontalAlignment','center','FontSize',12);

    %% --- SLIDER ---
    sld = uislider(f);
    sld.Position = [380 50 780 3];
    sld.Limits = [1 N];
    sld.Value = 1;
    sld.MajorTicks = [];
    sld.ValueChangingFcn = @(~,evt) update(round(evt.Value));
    sld.ValueChangedFcn  = @(~,evt) update(round(evt.Value));

    %% --- INITIAL DRAW ---
    update(1);

    %% --- UPDATE FUNCTION ---
    function update(k)
        k = max(1,min(N,k));
        idx = max(1,k-tailLen+1):k;

        % trails
        set(hT1,'XData',r1(idx,1),'YData',r1(idx,2),'ZData',r1(idx,3));
        set(hT2,'XData',r2(idx,1),'YData',r2(idx,2),'ZData',r2(idx,3));
        set(hT3,'XData',r3(idx,1),'YData',r3(idx,2),'ZData',r3(idx,3));

        % points
        set(hP1,'XData',r1(k,1),'YData',r1(k,2),'ZData',r1(k,3));
        set(hP2,'XData',r2(k,1),'YData',r2(k,2),'ZData',r2(k,3));
        set(hP3,'XData',r3(k,1),'YData',r3(k,2),'ZData',r3(k,3));

        % labels
        lbl.Text = sprintf('t = %.4g   (sample %d/%d)',t(k),k,N);
        txt1.Value = {
            sprintf('Body 1:')
            sprintf('  x = %.4g   y = %.4g   z = %.4g',r1(k,1),r1(k,2),r1(k,3))
            sprintf('  vx = %.4g  vy = %.4g  vz = %.4g',v1(k,1),v1(k,2),v1(k,3))
        };
        txt2.Value = {
            sprintf('Body 2:')
            sprintf('  x = %.4g   y = %.4g   z = %.4g',r2(k,1),r2(k,2),r2(k,3))
            sprintf('  vx = %.4g  vy = %.4g  vz = %.4g',v2(k,1),v2(k,2),v2(k,3))
        };
        txt3.Value = {
            sprintf('Body 3:')
            sprintf('  x = %.4g   y = %.4g   z = %.4g',r3(k,1),r3(k,2),r3(k,3))
            sprintf('  vx = %.4g  vy = %.4g  vz = %.4g',v3(k,1),v3(k,2),v3(k,3))
        };
        drawnow limitrate;
    end
end
