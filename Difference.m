clc; clear; close all;

%% ===== Fonts & sizes (SCI-style English) =====
font_main = 'Times New Roman';     % or 'Arial'
fs_label  = 24;                    % axis label
fs_tick   = 20;                    % tick numbers
fs_stats  = 24;                    % stats box text
fs_annot  = 20;                    % point annotations

%% ===== Plot padding (y-limits) =====
pad_ratio_phi = 0.12;  % displacement error
pad_ratio_v   = 0.12;  % velocity error
pad_ratio_a   = 0.14;  % acceleration error
shrink_phi    = 1.3;   % >1 expand away from center; =1 none
shrink_v      = 1.0;
shrink_a      = 1.0;

%% -------------------- Data I/O --------------------
file_phi = 'C:\Users\robot\Desktop\Adams11\1.xls';
file_v   = 'C:\Users\robot\Desktop\Adams11\new2.xls';
file_a   = 'C:\Users\robot\Desktop\Adams11\2.10800.xls';

% Displacement
try
    M = readmatrix(file_phi);
catch
    M = xlsread(file_phi);
end
t_data   = M(:,1);
phi_data = M(:,2);

% Velocity
sim_v   = readmatrix(file_v,'FileType','text','NumHeaderLines',3);
t_sim_v = sim_v(:,1);
v_sim   = -sim_v(:,2);

% Acceleration
sim_a   = readmatrix(file_a,'FileType','text','NumHeaderLines',3);
t_sim_a = sim_a(:,1);
a_sim   = -sim_a(:,2);

%% -------------------- Theoretical values --------------------
PI = pi;
t_stage = [0 30 70 125 165 210 250 305 345 360]/360;
N = numel(t_data);
phi_th = zeros(N,1);
v_th   = zeros(N,1);
a_th   = zeros(N,1);

for i = 1:N
    t = t_data(i);
    if (t >= t_stage(1)) && (t < t_stage(2))
        phi=0; v=0; a=0;
    elseif (t >= t_stage(2)) && (t < t_stage(3))
        delta=30; T=t_stage(3)-t_stage(2); tn=(t-t_stage(2))/T;
        if tn<=1/8
            phi=1/(4+PI)*(PI*tn - 1/4*sin(4*PI*tn));
            v=PI/(4+PI)*(1 - cos(4*PI*tn))*(delta/T);
            a=4*PI^2/(4+PI)*sin(4*PI*tn)*(delta/T^2);
        elseif tn<=7/8
            phi=1/(4+PI)*(2 + PI*tn - 9/4*sin((4*PI*tn+PI)/3));
            v=PI/(4+PI)*(1 - 3*cos((4*PI*tn+PI)/3))*(delta/T);
            a=4*PI^2/(4+PI)*sin((4*PI*tn+PI)/3)*(delta/T^2);
        else
            phi=1/(4+PI)*(4 + PI*tn - 1/4*sin(4*PI*tn));
            v=PI/(4+PI)*(1 - cos(4*PI*tn))*(delta/T);
            a=4*PI^2/(4+PI)*sin(4*PI*tn)*(delta/T^2);
        end
        phi=phi*delta;
    elseif (t >= t_stage(3)) && (t < t_stage(4))
        T=t_stage(4)-t_stage(3); tn=(t-t_stage(3))/T;
        phi=(10*tn^3 - 15*tn^4 + 6*tn^5)*30 + 30;
        v=(30*tn^2 - 60*tn^3 + 30*tn^4)*(30/T);
        a=(60*tn - 180*tn^2 + 120*tn^3)*(30/T^2);
    elseif (t >= t_stage(4)) && (t < t_stage(5))
        T=t_stage(5)-t_stage(4); tn=(t-t_stage(4))/T;
        phi=60 - 0.5*(1 - cos(PI*tn))*60;
        v=-(PI/2)*sin(PI*tn)*(60/T);
        a=-(PI^2/2)*cos(PI*tn)*(60/T^2);
    elseif (t >= t_stage(5)) && (t < t_stage(6))
        phi=0; v=0; a=0;
    elseif (t >= t_stage(6)) && (t < t_stage(7))
        T=t_stage(7)-t_stage(6); tn=(t-t_stage(6))/T;
        phi=(35*tn^4 - 84*tn^5 + 70*tn^6 - 20*tn^7)*60;
        v=(140*tn^3 - 420*tn^4 + 420*tn^5 - 140*tn^6)*(60/T);
        a=(420*tn^2 - 1680*tn^3 + 2100*tn^4 - 840*tn^5)*(60/T^2);
    elseif (t >= t_stage(7)) && (t < t_stage(8))
        phi=60; v=0; a=0;
    elseif (t >= t_stage(8)) && (t < t_stage(9))
        T=t_stage(9)-t_stage(8); tn=(t-t_stage(8))/T;
        phi=60 - (tn - 1/(2*PI)*sin(2*PI*tn))*60;
        v=-(1 - cos(2*PI*tn))*(60/T);
        a=-2*PI*sin(2*PI*tn)*(60/T^2);
    else
        phi=0; v=0; a=0;
    end
    phi_th(i)=phi; v_th(i)=v; a_th(i)=a;
end

%% -------------------- Errors --------------------
delta_phi     = phi_data - phi_th;
v_sim_interp  = interp1(t_sim_v, v_sim, t_data, 'linear', 'extrap');
a_sim_interp  = interp1(t_sim_a, a_sim, t_data, 'linear', 'extrap');
diff_v        = v_th - v_sim_interp;
diff_a        = a_th - a_sim_interp;

%% -------------------- Average deviation percentage (your formula) --------------------
% 100 * sum(|Δ|) / sum(|measurement|)
trunc2 = @(x) fix(x*100)/100; % truncate to 2 decimals (toward zero)

num_phi = sum(abs(delta_phi)); den_phi = sum(abs(phi_data));
num_v   = sum(abs(diff_v));    den_v   = sum(abs(v_sim_interp));
num_a   = sum(abs(diff_a));    den_a   = sum(abs(a_sim_interp));

if den_phi>0, fprintf('Avg deviation %% (displacement): %.2f%%\n', trunc2(100*num_phi/den_phi));
else,         fprintf('Avg deviation %% (displacement): N/A (zero denominator)\n'); end
if den_v>0,   fprintf('Avg deviation %% (velocity):     %.2f%%\n', trunc2(100*num_v/den_v));
else,         fprintf('Avg deviation %% (velocity):     N/A (zero denominator)\n'); end
if den_a>0,   fprintf('Avg deviation %% (acceleration): %.2f%%\n', trunc2(100*num_a/den_a));
else,         fprintf('Avg deviation %% (acceleration): N/A (zero denominator)\n'); end

%% -------------------- Colors --------------------
c_phi = [0 0.4470 0.7410];
c_v   = [0.8500 0.3250 0.0980];
c_a   = [0.10 0.60 0.60];
c_pos = [0.85 0.10 0.10];
c_neg = [0.10 0.40 0.90];

%% ========= Figure 1: Angular displacement error =========
fig1 = figure('Color','w','Visible','on');
try, set(fig1,'WindowState','maximized'); catch, set(fig1,'Units','normalized','OuterPosition',[0 0 1 1]); end
ax1 = axes(fig1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'FontName',font_main,'FontSize',fs_tick);

plot(ax1, t_data, delta_phi, '-', 'LineWidth', 2.2, 'Color', c_phi);
[max_pos1,pos_idx1] = max(delta_phi); [min_neg1,neg_idx1] = min(delta_phi);
plot(ax1, t_data(pos_idx1), max_pos1, 'o', 'MarkerSize', 9, 'MarkerFaceColor', c_pos, 'MarkerEdgeColor','k','LineWidth',1.1);
plot(ax1, t_data(neg_idx1), min_neg1, 's', 'MarkerSize', 9, 'MarkerFaceColor', c_neg, 'MarkerEdgeColor','k','LineWidth',1.1);
yline(ax1,0,'-','HandleVisibility','off');

xlim(ax1,[min(t_data) max(t_data)]);
set_nice_ylim(ax1, delta_phi, pad_ratio_phi, shrink_phi);

text_adjacent(ax1, t_data(pos_idx1), max_pos1, sprintf('Max positive deviation: %.2f°', trunc2(max_pos1)), c_pos, fs_annot, font_main);
text_adjacent(ax1, t_data(neg_idx1), min_neg1, sprintf('Max negative deviation: %.2f°', trunc2(min_neg1)), c_neg, fs_annot, font_main);

xlabel(ax1,'t (s)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);
ylabel(ax1,'Angular displacement error (°)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);

add_stats_multiline(ax1, '\epsilon_{\varphi}', ...
    max(abs(delta_phi)), mean(abs(delta_phi)), rms(delta_phi), ...
    '^{\circ}', fs_stats, font_main);

%% ========= Figure 2: Angular velocity error =========
fig2 = figure('Color','w','Visible','on');
try, set(fig2,'WindowState','maximized'); catch, set(fig2,'Units','normalized','OuterPosition',[0 0 1 1]); end
ax2 = axes(fig2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'FontName',font_main,'FontSize',fs_tick);

plot(ax2, t_data, diff_v, '-', 'LineWidth', 2.2, 'Color', c_v);
[max_pos2,pos_idx2] = max(diff_v); [min_neg2,neg_idx2] = min(diff_v);
plot(ax2, t_data(pos_idx2), max_pos2, 'o','MarkerSize',9,'MarkerFaceColor',c_pos,'MarkerEdgeColor','k','LineWidth',1.1);
plot(ax2, t_data(neg_idx2), min_neg2, 's','MarkerSize',9,'MarkerFaceColor',c_neg,'MarkerEdgeColor','k','LineWidth',1.1);
yline(ax2,0,'-','HandleVisibility','off');

xlim(ax2,[min(t_data) max(t_data)]);
set_nice_ylim(ax2, diff_v, pad_ratio_v, shrink_v);

text_adjacent(ax2, t_data(pos_idx2), max_pos2, sprintf('Max positive deviation: %.2f°/s', trunc2(max_pos2)), c_pos, fs_annot, font_main);
text_adjacent(ax2, t_data(neg_idx2), min_neg2, sprintf('Max negative deviation: %.2f°/s', trunc2(min_neg2)), c_neg, fs_annot, font_main);

xlabel(ax2,'t (s)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);
ylabel(ax2,'Angular velocity error (°/s)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);

add_stats_multiline(ax2, '\epsilon_{\dot{\varphi}}', ...
    max(abs(diff_v)), mean(abs(diff_v)), rms(diff_v), ...
    '^{\circ}/\mathrm{s}', fs_stats, font_main);

%% ========= Figure 3: Angular acceleration error =========
fig3 = figure('Color','w','Visible','on');
try, set(fig3,'WindowState','maximized'); catch, set(fig3,'Units','normalized','OuterPosition',[0 0 1 1]); end
ax3 = axes(fig3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
set(ax3,'FontName',font_main,'FontSize',fs_tick);

plot(ax3, t_data, diff_a, '-', 'LineWidth', 2.2, 'Color', c_a);
[max_pos3,pos_idx3] = max(diff_a); [min_neg3,neg_idx3] = min(diff_a);
plot(ax3, t_data(pos_idx3), max_pos3, 'o','MarkerSize',9,'MarkerFaceColor',c_pos,'MarkerEdgeColor','k','LineWidth',1.1);
plot(ax3, t_data(neg_idx3), min_neg3, 's','MarkerSize',9,'MarkerFaceColor',c_neg,'MarkerEdgeColor','k','LineWidth',1.1);
yline(ax3,0,'-','HandleVisibility','off');

xlim(ax3,[min(t_data) max(t_data)]);
set_nice_ylim(ax3, diff_a, pad_ratio_a, shrink_a);

text_adjacent(ax3, t_data(pos_idx3), max_pos3, sprintf('Max positive deviation: %.2f°/s^2', trunc2(max_pos3)), c_pos, fs_annot, font_main);
text_adjacent(ax3, t_data(neg_idx3), min_neg3, sprintf('Max negative deviation: %.2f°/s^2', trunc2(min_neg3)), c_neg, fs_annot, font_main);

xlabel(ax3,'t (s)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);
ylabel(ax3,'Angular acceleration error (°/s^2)','FontWeight','bold','FontName',font_main,'FontSize',fs_label);

add_stats_multiline(ax3, '\epsilon_{\ddot{\varphi}}', ...
    max(abs(diff_a)), mean(abs(diff_a)), rms(diff_a), ...
    '^{\circ}/\mathrm{s}^{2}', fs_stats, font_main);

%% ==================== Local functions ====================
function set_nice_ylim(ax, y, pad_ratio, shrink_factor)
% add vertical padding and optional symmetric expansion
if nargin<3 || isempty(pad_ratio),     pad_ratio = 0.10; end
if nargin<4 || isempty(shrink_factor), shrink_factor = 1.0; end
yl = [min(y) max(y)];
if yl(1)==yl(2), yl = yl + [-1 1]*eps; end
r  = range(yl);
yl = [yl(1)-pad_ratio*r, yl(2)+pad_ratio*r]; % base padding
c  = mean(yl);
yl = c + (yl - c) * shrink_factor;           % expand/contract
ylim(ax, yl);
end

function h = text_adjacent(ax, x0, y0, str, color, fs_annot, font_main)
% place text near a point with auto side selection + draggable
xr = xlim(ax); yr = ylim(ax);
dx = 0.02*(xr(2)-xr(1));   % horizontal offset (2% width)
dy = 0.03*(yr(2)-yr(1));   % vertical offset (3% height)
if x0 < mean(xr), x = x0 + dx; ha = 'left';  else, x = x0 - dx; ha = 'right'; end
if y0 < mean(yr), y = y0 + dy; va = 'bottom'; else, y = y0 - dy; va = 'top';    end
x = min(max(x, xr(1)+0.01*(xr(2)-xr(1))), xr(2)-0.01*(xr(2)-xr(1)));
y = min(max(y, yr(1)+0.01*(yr(2)-yr(1))), yr(2)-0.01*(yr(2)-yr(1)));
h = text(ax, x, y, str, ...
    'FontSize',fs_annot,'FontName',font_main,'FontWeight','bold', ...
    'Color', color,'HorizontalAlignment',ha,'VerticalAlignment',va, ...
    'BackgroundColor',[1 1 1 0.88],'EdgeColor','k','LineWidth',0.8, ...
    'Clipping','on','HitTest','on','PickableParts','all');
make_text_draggable(h);
end

function make_text_draggable(h)
set(h,'ButtonDownFcn',@startDrag);
    function startDrag(~,~)
        fig = ancestor(h,'figure'); ax = ancestor(h,'axes');
        set(fig,'WindowButtonMotionFcn',@dragging, ...
                'WindowButtonUpFcn',@stopDrag, ...
                'Pointer','hand');
        uistack(h,'top');
        function dragging(~,~)
            cp = ax.CurrentPoint(1,1:2);
            xr = xlim(ax); yr = ylim(ax);
            cp(1) = min(max(cp(1), xr(1)), xr(2));
            cp(2) = min(max(cp(2), yr(1)), yr(2));
            set(h,'Position',cp);
        end
        function stopDrag(~,~)
            set(fig,'WindowButtonMotionFcn','', 'WindowButtonUpFcn','', 'Pointer','arrow');
        end
    end
end

function h = add_stats_multiline(ax, eps_label_latex, maxabs, mae, rmse, unit_latex, fs_stats, font_main)
% Stats box (3 lines): max|ε_*|, MAE, RMSE (units; 2-decimal truncation)
xr = xlim(ax); yr = ylim(ax);
x_box = xr(1) + 0.88*(xr(2)-xr(1));
y_box = yr(1) + 0.08*(yr(2)-yr(1));
trunc2 = @(x) fix(x*100)/100;

maxabs_t = trunc2(maxabs);
mae_t    = trunc2(mae);
rmse_t   = trunc2(rmse);

str = sprintf([ ...
    '$\\boldmath \\max\\,\\left| %s \\right|$: %.2f$%s$' ...
    '\n$\\boldmath \\mathrm{MAE}$: %.2f$%s$' ...
    '\n$\\boldmath \\mathrm{RMSE}$: %.2f$%s$' ], ...
    eps_label_latex, maxabs_t, unit_latex, ...
    mae_t,           unit_latex, ...
    rmse_t,          unit_latex);

h = text(ax, x_box, y_box, str, ...
    'FontName', font_main, 'FontSize', fs_stats, 'FontWeight','bold', ...
    'Color',[0 0 0], 'EdgeColor',[0.55 0 0], 'BackgroundColor',[1 1 1 0.88], ...
    'LineWidth',1.2, 'Margin',0.5, ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
    'Interpreter','latex', 'Clipping','on','HitTest','on','PickableParts','all');
make_text_draggable(h);
end
