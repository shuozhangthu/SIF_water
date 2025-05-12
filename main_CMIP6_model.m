function main_CMIP6_model()
    % MAIN FUNCTION: Process CMIP6 model data
    
    % 1. Initialize parameters
    model_params = initialize_parameters();
    
    % 2. Process each input file
    process_model_data(model_params);
end

function model_params = initialize_parameters()
    % Initialize model parameters and paths
    
    model_params = struct();
    
    % Model configuration
    model_params.model = 'ssp245'; % Can be changed to 'ssp585' 
    model_params.name = 'CESM2-WACCM';
    % 'ACCESS-CM2', 'ACCESS-ESM1-5','CanESM5','CMCC-CM2-SR5'，'CMCC-ESM2','IITM-ESM','MPI-ESM1-2-HR'  1850
    %  'AWI-CM-1-1-MR','CAMS-CSM1-0' ,'CAS-ESM2-0','FGOALS-g3','NESM3'  Start year:2015
    %  'CESM2-WACCM'  Start year: 0001
    model_params.variable = 'pr';  % Can be changed to 'tas' for temperature
    model_params.start_year = 1;   % Adjust based on model
    
    model_params.input_dir = fullfile('./input', model_params.model, ...
                                     model_params.name, model_params.variable);
    model_params.output_dir = fullfile('./output', model_params.model, ...
                                      model_params.name, 'processed_025deg', ...
                                      model_params.variable);
    
    if ~exist(model_params.output_dir, 'dir')
        mkdir(model_params.output_dir);
    end
    
    model_params.target_size = [720, 1440];  % 0.25 degree resolution
    model_params.latlim = [-90 90];
    model_params.lonlim = [-180 180];
end

function process_model_data(model_params)
    
    file_pattern = fullfile(model_params.input_dir, '*.nc');
    files = dir(file_pattern);
    
    if isempty(files)
        error('No input files found at: %s', model_params.input_dir);
    end
    
    for n = 1:length(files)
        process_single_file(files(n), model_params);
    end
end

function process_single_file(file_info, model_params)
    
    file_path = fullfile(model_params.input_dir, file_info.name);
    fprintf('Processing file: %s\n', file_info.name);
    
    try
        data = ncread(file_path, model_params.variable);
        time_data = ncread(file_path, 'time');
        
        base_date = datetime(model_params.start_year, 1, 1, 0, 0, 0);
        time_vector = base_date + days(time_data);
        
        year_num = length(time_vector) / 12;
        first_year = time_vector(1).Year;
        
        for year_idx = 1:year_num
            process_yearly_data(data, year_idx, first_year, model_params);
        end
        
    catch ME
        warning('File processing failed: %s\nError: %s', file_info.name, ME.message);
    end
end

function process_yearly_data(data, year_idx, first_year, model_params)

    start_idx = (year_idx-1)*12 + 1;
    end_idx = year_idx*12;
    yearly_data = data(:, :, start_idx:end_idx);
    
    switch model_params.variable
        case 'pr'  % Precipitation
            processed_data = nanmean(yearly_data, 3) * 3600 * 24 * 365;
            units = 'mm/year';
        case 'tas' % Temperature
            processed_data = nanmean(yearly_data, 3) - 273.15;
            units = 'degC';
    end
    
    rotated_data = rot90(processed_data);
    [~, cols] = size(rotated_data);
    
    % Longitude adjustment (0-360 to -180-180)
    half_cols = floor(cols/2);
    rearranged_data = [rotated_data(:, half_cols+1:cols), rotated_data(:, 1:half_cols)];
    
    % Resample to target resolution
    resized_data = imresize(rearranged_data, model_params.target_size, 'nearest');
    
    current_year = first_year + year_idx - 1;
    output_filename = sprintf('%s_%s_%s_%04d.nc', ...
                             model_params.model, ...
                             model_params.name, ...
                             model_params.variable, ...
                             current_year);
    output_path = fullfile(model_params.output_dir, output_filename);
    
    nccreate(output_path, model_params.variable, ...
             'Dimensions', {'lat', model_params.target_size(1), ...
                           'lon', model_params.target_size(2)}, ...
             'Format', 'classic');
    
    ncwrite(output_path, model_params.variable, resized_data);
    
    ncwriteatt(output_path, model_params.variable, 'units', units);
    ncwriteatt(output_path, model_params.variable, 'long_name', ...
              [model_params.variable ' annual mean']);
    

    lat = linspace(model_params.latlim(1), model_params.latlim(2), model_params.target_size(1));
    lon = linspace(model_params.lonlim(1), model_params.lonlim(2), model_params.target_size(2));
    
    nccreate(output_path, 'lat', 'Dimensions', {'lat', model_params.target_size(1)});
    nccreate(output_path, 'lon', 'Dimensions', {'lon', model_params.target_size(2)});
    
    ncwrite(output_path, 'lat', lat);
    ncwrite(output_path, 'lon', lon);

    ncwriteatt(output_path, 'lat', 'units', 'degrees_north');
    ncwriteatt(output_path, 'lon', 'units', 'degrees_east');
end