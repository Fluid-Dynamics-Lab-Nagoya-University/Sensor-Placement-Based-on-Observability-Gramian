function [f1,p1,p2] = F_fig_SST_sensors_color_compare(sensi,mask, sensorsA,sensorsB)
f1=figure();
f1.Units='pixels';
movegui('center')
%   display_sensors Template for displaying image data
%   Displays enso sensors on white w/ black continents
snapshot = NaN*zeros(numel(mask),1);
x=sensi';
snapshot(mask==1) = x;
C = reshape(real(snapshot),size(mask))';
shading interp;
jetmod=jet(256);
jetmod(1,:)=0;
colormap(jetmod);
%%
%{
    c=colorbar;    %<========color bar
    if verLessThan('matlab','9.8')
    else
        f1.OuterPosition=[201 401 560 400];
        axis off
        c=colorbar();
        c.FontSize = 14;
%         c.Label.String = 'Standard deviation of Temperature [℃] ';
        c.Label.String = 'Surface Temperature [℃] ';
%         c.Label.String = 'RMSE of Temperature [K]';
        c.Label.FontName = 'TimesNewRoman';
        c.Label.FontSize = 16;
%         c.Location='northoutside'

        daspect([1.3 1 1])
    end
%}
daspect([1.3 1 1])
axis off
caxis([-5 30]);
%%
set(gca, 'FontName', 'Times','Color','white', 'FontSize', 20);

pivotA=sensorsA;
pivotB=sensorsB;

sensors_location = zeros(size(mask));
P = zeros(size(x));
P(pivotA)=1:length(pivotA);
sensors_location(mask==1) = P;
S = reshape(real(sensors_location)',numel(mask),1);
[~,ICA,~] = unique(S);

%! align Ilin with pivot somehow
[I1,J1] = ind2sub(size(sensors_location'),ICA(2:end));

Q = zeros(size(x));
Q(pivotB)=1:length(pivotB);
sensors_location(mask==1) = Q;
S = reshape(real(sensors_location)',numel(mask),1);
[~,ICB,~] = unique(S);
[I2,J2] = ind2sub(size(sensors_location'),ICB(2:end));
hold on;
[p,~]=size(sensorsA);

p1 = plot(J1,I1,'o','MarkerSize', 8,'LineWidth', 1.5, 'MarkerFaceColor', 'none','MarkerEdgeColor','k','MarkerFaceColor','none');
p2 = plot(J2,I2,'x','MarkerSize', 8,'LineWidth', 1.5, 'MarkerFaceColor', 'none','MarkerEdgeColor','r','MarkerFaceColor','none');
hold off;
%     axis off;
end
