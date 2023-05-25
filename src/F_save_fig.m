function F_save_fig(file_name)
%     cd_a=cd;
%     cd ../
%    date_now=datestr(now,'mmdd_HHMMSS')
%    file_name=(sprintf('fig_error_%s_%s',data_type,date_now));
gcfnow = gcf;
% gcfnow.InvertHardcopy = 'off';
ax = gca;
% ax.Color = 'none';
% gcfnow.Color = 'none';
    saveas(gcfnow, file_name,'png');
%     saveas(gcf, file_name,'pdf');
if verLessThan('matlab','9.8')
else
    exportgraphics(gcfnow,[file_name,'.pdf'],'BackgroundColor','none','ContentType','vector')
end
    saveas(gcfnow, file_name);
if ispc==1
%     saveas(gcfnow, file_name,'meta');
%     cd(cd_a)
% print -dmeta -vector
copygraphics(gcfnow,'BackgroundColor','none')
end
disp(['saved to --> ',file_name])
end