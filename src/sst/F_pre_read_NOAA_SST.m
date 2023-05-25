function [ Lat, Lon, time, mask, sst ] = F_pre_read_NOAA_SST( filename, maskname )
%read_data Read in ENSO data
%   Return meshgrid of lat,long,times and sea-surface temp

    sst  = ncread(filename, 'sst');
    time = ncread(filename, 'time'); %days
    % each element of time array is a new week, in units of days

    lat  = ncread(filename, 'lat');
    lon  = ncread(filename, 'lon');

    test = ncinfo(maskname);
    if find(ismember({test.Variables.Name},'mask'))
        mask = ncread(maskname, 'mask');
    elseif find(ismember({test.Variables.Name},'lsmask'))
        mask = ncread(maskname, 'lsmask');
    else
        error('cannot read mask from input mask data')
    end
    mask = mask(:,:,1);
    [Lat, Lon] = meshgrid(lat, lon);

    % % Do not plot continent mask
    % Lat(mask==0) = NaN;
    % Lon(mask==0) = NaN;
    
end
