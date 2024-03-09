function cleanSubplot(varargin)
% cleanSubplot(fontsize, linewidth)

if nargin == 0
    fs = 20;
    lw = 2;
elseif nargin == 1
    fs = varargin{1};
    lw = 2;
elseif nargin == 2
    fs = varargin{1};
    lw = varargin{2};
end
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1];

set(gcf,'color',bgc);
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'FontName','Arial')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);

end