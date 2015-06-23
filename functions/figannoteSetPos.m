%   to annote a figure with a textbox

function figannoteSetPos(objfig,str,fontsize,where)
% example for where: [0.35 0.9 0.3 0.1]
    figure(objfig);
    an=annotation('textbox',where);
    set(an,'String', str);
    set(an,'VerticalAlignment','middle');
    set(an,'HorizontalAlignment','left');
    set(an,'FontSize',fontsize);
    set(an,'FitHeightToText','on');
    set(an,'Interpreter','latex');
    set(an,'LineStyle','none');