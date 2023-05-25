function [fig_handle] = F_ffvis_sensor(data,rmin,rmax,sens_col,sens_row)
    fig_handle = figure();
    imshow(data,[rmin rmax],...
        'Border','tight','InitialMagnification',500);
    daspect([1 1 1])
    set(gcf,'PaperPositionMode','auto')
    
    hold on; 
    [p,~]=size(sens_col);
    if p>0
%         s = scatter(sens_col,sens_row,[],[0 0 0]);
%         s.LineWidth = 1;
%         s.MarkerFaceColor = 'none';
        plot(sens_col,sens_row,'ko','MarkerSize',24,...
            'LineWidth', 2.5, 'MarkerFaceColor','none');
        a = [1:p]'; b = num2str(a); c = cellstr(b);
        dx = -7; dy = -1;
%         text(sens_col+dx, sens_row+dy, c,'FontSize',20);
    end
    axis off
end