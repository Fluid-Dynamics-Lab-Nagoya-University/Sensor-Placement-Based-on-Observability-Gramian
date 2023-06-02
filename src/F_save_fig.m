function F_save_fig(file_name)
gcfnow = gcf;
saveas(gcfnow, file_name,'png');
if verLessThan('matlab','9.8')
else
    exportgraphics(gcfnow,[file_name,'.pdf'],'BackgroundColor','none','ContentType','vector')
end
saveas(gcfnow, file_name);
if ispc==1
    copygraphics(gcfnow,'BackgroundColor','none')
end
disp(['saved to --> ',file_name])
end