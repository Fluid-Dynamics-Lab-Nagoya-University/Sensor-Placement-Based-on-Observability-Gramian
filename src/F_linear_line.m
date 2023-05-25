function F_linear_line(hfig,hplot,order,strtmp)
    figure(hfig);
    xtmp    = hplot.XData;
    ytmp    = hplot.YData;
%     order   =
    clr     = hplot.Color;
    
    isodd   = (mod(size(xtmp,2),2) == 1);
    if isodd
        xc  = xtmp(fix(end/2+0.5));
        yc  = ytmp(fix(end/2+0.5));
    else
        xc  = (xtmp(end/2)+xtmp(end/2+1))/2;
        yc  = (ytmp(end/2)+ytmp(end/2+1))/2;
    end
    
    loglog(xtmp,((xtmp).^order)*yc/xc.^order,'DisplayName','','Color',clr);

    ifdisp  = true;
    if ifdisp
        txttmp = text(xc,yc,['  \propto ',strtmp],'Fontsize',15,...
            'FontName','Times New Roman','Interpreter','tex'); 
%         hfig.TextBox
    end

end