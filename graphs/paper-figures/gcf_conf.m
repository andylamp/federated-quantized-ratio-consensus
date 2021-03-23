% get the gcf relative to your monitor
% get(1,'Position')

% set the gcf
%set(gcf,'position',[1400, 500, 600, 420])

% better
% set(gcf,'position',[1400, 500, 850, 450])
set(gcf,'position',[1400, 500, 850, 500])

% set font for the current figure in all bits of it.
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)

ylabel('Node State Variables ($q_j^s[k]$)', 'interpreter', "latex", "fontsize", 18)
xlabel('Number of Iterations ($k$)', 'interpreter', "latex", "fontsize", 18)

ylabel('Iterations', 'interpreter', "latex", "fontsize", 18)
xlabel('Nodes', 'interpreter', "latex", "fontsize", 18)

% specific plot configs

% -- converge stats with vertical fix -- %
set(gcf,'position',[1400, 500, 850, 350])

% -- iterations + time with vertical fix -- %
set(gcf,'position',[1400, 500, 850, 300])