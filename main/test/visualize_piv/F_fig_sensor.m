function [fig_handle] = F_fig_sensor(back,data_size,rmin,rmax,loc_normal,sensors)
%% usage
% fig1 = F_fig_sensor(back,zeros(121,121),-10,20,loc.normal,sensors);
% cs = '18deg'; case_num = char(cs);
% F_wingdraw(case_num,wingx,wingy);
%%
[y,x]=size(data_size);
[data]=F_vec_to_pic(back(:,1),loc_normal,x,y);
[sens_row,sens_col]=F_sens_trim(data, loc_normal, sensors);
[fig_handle] = F_ffvis_sensor(data,rmin,rmax,sens_col,sens_row);
% F_colorbar_PIV
map = F_colormap_org('jet',255);
% map(round(end/2)-2:round(end/2)+2,:) = repmat([0.4 0.4 0.4],[5 1]);
colormap(map);
end