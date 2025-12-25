clc; clear; close all;

%% ----------------  Data generation  ----------------
t_all = linspace(0, 1, 7278);
phi_all = zeros(size(t_all));
v_all   = zeros(size(t_all));
a_all   = zeros(size(t_all));
PI = pi;
t_stages = [0, 30, 70, 125, 165, 210, 250, 305, 345, 360] / 360;

for i = 1:length(t_all)
    t = t_all(i);
    %% 1 Dwell
    if t >= t_stages(1) && t < t_stages(2)
        phi = 0; v = 0; a = 0;
    %% 2 Modified sine
    elseif t >= t_stages(2) && t < t_stages(3)
        t2_start = t_stages(2); t2_end = t_stages(3);
        delta_phi = 30; delta_t = t2_end - t2_start; tn = (t - t2_start)/delta_t;
        if tn <= 1/8
            phi = 1/(4+PI)*(PI*tn - 1/4*sin(4*PI*tn));
            v   = PI/(4+PI)*(1 - cos(4*PI*tn));
            a   = 4*PI^2/(4+PI)*sin(4*PI*tn);
        elseif tn <= 7/8
            phi = 1/(4+PI)*(2 + PI*tn - 9/4*sin((4*PI*tn+PI)/3));
            v   = PI/(4+PI)*(1 - 3*cos((4*PI*tn+PI)/3));
            a   = 4*PI^2/(4+PI)*sin((4*PI*tn+PI)/3);
        else
            phi = 1/(4+PI)*(4 + PI*tn - 1/4*sin(4*PI*tn));
            v   = PI/(4+PI)*(1 - cos(4*PI*tn));
            a   = 4*PI^2/(4+PI)*sin(4*PI*tn);
        end
        phi = phi*delta_phi; v = v*(delta_phi/delta_t); a = a*(delta_phi/delta_t^2);
    %% 3 5th-order polynomial
    elseif t >= t_stages(3) && t < t_stages(4)
        tn = (t - t_stages(3))/(t_stages(4)-t_stages(3));
        phi = (10*tn^3 - 15*tn^4 + 6*tn^5)*30 + 30;
        v   = (30*tn^2 - 60*tn^3 + 30*tn^4)*(30/(t_stages(4)-t_stages(3)));
        a   = (60*tn - 180*tn^2 + 120*tn^3)*(30/(t_stages(4)-t_stages(3))^2);
    %% 4 Cosine acceleration (fall)
    elseif t >= t_stages(4) && t < t_stages(5)
        tn = (t - t_stages(4))/(t_stages(5)-t_stages(4));
        phi = 60 - 0.5*(1 - cos(PI*tn))*60;
        v   = -(PI/2)*sin(PI*tn)*(60/(t_stages(5)-t_stages(4)));
        a   = -(PI^2/2)*cos(PI*tn)*(60/(t_stages(5)-t_stages(4))^2);
    %% 5 Dwell
    elseif t >= t_stages(5) && t < t_stages(6)
        phi = 0; v = 0; a = 0;
    %% 6 7th-order polynomial
    elseif t >= t_stages(6) && t < t_stages(7)
        tn = (t - t_stages(6))/(t_stages(7)-t_stages(6));
        phi = (35*tn^4 - 84*tn^5 + 70*tn^6 - 20*tn^7)*60;
        v   = (140*tn^3 - 420*tn^4 + 420*tn^5 - 140*tn^6)*(60/(t_stages(7)-t_stages(6)));
        a   = (420*tn^2 - 1680*tn^3 + 2100*tn^4 - 840*tn^5)*(60/(t_stages(7)-t_stages(6))^2);
    %% 7 Dwell
    elseif t >= t_stages(7) && t < t_stages(8)
        phi = 60; v = 0; a = 0;
    %% 8 Sinusoidal acceleration (fall)
    elseif t >= t_stages(8) && t < t_stages(9)
        tn = (t - t_stages(8))/(t_stages(9)-t_stages(8));
        phi = 60 - (tn - 1/(2*PI)*sin(2*PI*tn))*60;
        v   = -(1 - cos(2*PI*tn))*(60/(t_stages(9)-t_stages(8)));
        a   = -2*PI*sin(2*PI*tn)*(60/(t_stages(9)-t_stages(8))^2);
    %% 9 Dwell
    else
        phi = 0; v = 0; a = 0;
    end
    phi_all(i) = phi; v_all(i) = v; a_all(i) = a;
end

%% ================= Three subplots + top legend =================
lw = 4;
axisSize   = 20;
labelSize  = 22;
panelMarkSize = 24;     % a/b/c font size
legendSize = 20;
clr = lines(9);         % 9 colors for 9 stages

stageName = {'Stage 1 - dwell', ...
             'Stage 2 - modified sine', ...
             'Stage 3 - 5th-order polynomial', ...
             'Stage 4 - cosine acceleration', ...
             'Stage 5 - dwell', ...
             'Stage 6 - 7th-order polynomial', ...
             'Stage 7 - dwell', ...
             'Stage 8 - sinusoidal acceleration', ...
             'Stage 9 - dwell'};
clr = lines(9);

fig = figure('Color','w','Units','pixels','Position',[80 80 1400 800]);
tlo = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
% Global title
try
    title(tlo, 'Motion design of the spatial cam mechanism', ...
          'FontName','Times New Roman','FontSize',26,'FontWeight','bold');
catch
    sgtitle('Motion design of the spatial cam mechanism', ...
            'FontName','Times New Roman','FontSize',26,'FontWeight','bold');
end

%% 1) Displacement (solid line)
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
hPhi = gobjects(9,1);
for k = 1:9
    idx = (t_all >= t_stages(k)) & (t_all < t_stages(k+1));
    if k==9, idx = (t_all >= t_stages(k)); end
    hPhi(k) = plot(ax1, t_all(idx), phi_all(idx), ...
        'Color', clr(k,:), 'LineWidth', lw, 'LineStyle','-');
end
set(ax1,'FontSize',axisSize,'LineWidth',1.2,...
    'FontName','Times New Roman','GridAlpha',0.15);
ylabel(ax1,'Angular displacement (°)',...
       'FontSize',labelSize,'FontName','Times New Roman');

% Force ticks to include 60 and 70
yl1 = ylim(ax1);
ylim(ax1, [yl1(1), max(70, max(phi_all))*1.05]);
yt = get(ax1,'YTick'); yt = unique([yt 60 70]); yt = sort(yt);
set(ax1,'YTick',yt);

% Panel label a
text(ax1,0.02,0.92,'\bf a','Units','normalized',...
     'FontSize',panelMarkSize,'FontName','Times New Roman');

%% 2) Velocity (dashed line)
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
hV = gobjects(9,1);
for k = 1:9
    idx = (t_all >= t_stages(k)) & (t_all < t_stages(k+1));
    if k==9, idx = (t_all >= t_stages(k)); end
    hV(k) = plot(ax2, t_all(idx), v_all(idx), ...
        'Color', clr(k,:), 'LineWidth', lw, 'LineStyle','--');
end
set(ax2,'FontSize',axisSize,'LineWidth',1.2,...
    'FontName','Times New Roman','GridAlpha',0.15);
ylabel(ax2,'Angular velocity v (°/s)',...
       'FontSize',labelSize,'FontName','Times New Roman');
ylim(ax2,[-1.1*max(abs(v_all)), 1.1*max(abs(v_all))]);

text(ax2,0.02,0.92,'\bf b','Units','normalized',...
     'FontSize',panelMarkSize,'FontName','Times New Roman');

%% 3) Acceleration (dotted line)
ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
hA = gobjects(9,1);
for k = 1:9
    idx = (t_all >= t_stages(k)) & (t_all < t_stages(k+1));
    if k==9, idx = (t_all >= t_stages(k)); end
    hA(k) = plot(ax3, t_all(idx), a_all(idx), ...
        'Color', clr(k,:), 'LineWidth', lw, 'LineStyle',':');
end

% Black dashed connections between stages (keep as in original)
idx3_end  = find(t_all >= t_stages(4),1,'first')-1; idx4_start = idx3_end+1;
plot(ax3, [t_all(idx3_end) t_all(idx4_start)], ...
           [a_all(idx3_end) a_all(idx4_start)], ...
           'k--','LineWidth',1.5);
idx4_end  = find(t_all >= t_stages(5),1,'first')-1; idx5_start = idx4_end+1;
plot(ax3, [t_all(idx4_end) t_all(idx5_start)], ...
           [a_all(idx4_end) a_all(idx5_start)], ...
           'k--','LineWidth',1.5);

set(ax3,'FontSize',axisSize,'LineWidth',1.2,...
    'FontName','Times New Roman','GridAlpha',0.15);
ylabel(ax3,'Angular acceleration a (°/s^2)',...
       'FontSize',labelSize,'FontName','Times New Roman');
xlabel(ax3,'Time t (s)',...
       'FontSize',labelSize,'FontName','Times New Roman');

text(ax3,0.02,0.92,'\bf c','Units','normalized',...
     'FontSize',panelMarkSize,'FontName','Times New Roman');

%% Top legend under the main title (two rows)
rows    = 2;
numCols = ceil(numel(stageName)/rows);

lgd = legend(ax1, hPhi, stageName, ...
    'Orientation','horizontal', ...
    'NumColumns', numCols, ...
    'FontSize', legendSize, ...
    'Box','on', 'Interpreter','tex');

try
    lgd.Layout.Tile = 'north';
catch
    lgd.Location = 'northoutside';
end

try, lgd.ItemTokenSize = [18 9]; end
set(lgd,'FontName','Times New Roman');

% Optional export
% set(fig,'Renderer','painters');
% exportgraphics(fig,'cam_motion_tripanel_english.pdf','ContentType','vector');
% exportgraphics(fig,'cam_motion_tripanel_english.png','Resolution',600);
