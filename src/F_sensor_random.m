function [sensors]=F_sensor_random( U, p, s)


[n,~]=size(U);
if isempty(s)
    sensors = randperm(n,p)';
else
    sensors = randperm(round(n/s),p)';
end
sensors = [sensors(:)];

   
end