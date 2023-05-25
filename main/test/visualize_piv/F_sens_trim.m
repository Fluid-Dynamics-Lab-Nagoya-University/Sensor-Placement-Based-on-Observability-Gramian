function [sens_row,sens_col] = F_sens_trim(data, loc_normal, sensors)
    [y, x] = size(data);
    
    ind_normal=find(loc_normal==1);
    sens_row = []; sens_col = [];
    [p,~] = size(sensors);
    for i = 1:p
        sensor = sensors(i,1);
        ind_sensor=ind_normal(sensor);
        loc_sensor=false(y,x,1);
        loc_sensor(ind_sensor)=1;
        [a,b] = find(loc_sensor);
        sens_row = [sens_row; a];
        sens_col = [sens_col; b];
    end
end