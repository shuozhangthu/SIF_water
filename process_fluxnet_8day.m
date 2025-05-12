function [composite_8day] = process_fluxnet_8day()
    % Main function to process FLUXNET2015 data into 8-day composites
    % Output: 
    %   composite_8day - 2D matrix [ (years*46) × sites ] with 8-day GPP values
    %                   NaN indicates failed quality control
    
    try
        % 1. Read site metadata
        [site_names, year_range] = read_site_info();
        
        % 2. Initialize output matrix
        composite_8day = initialize_composite_output(site_names, year_range);
        
        % 3. Process all sites
        composite_8day = process_all_sites(site_names, year_range, composite_8day);
        
    catch ME
        error('FLUXNET processing failed: %s', ME.message);
    end
end

function [site_names, year_range] = read_site_info()
    % Read site metadata from Excel file
    % Returns:
    %   site_names  - Cell array of site IDs
    %   year_range  - Matrix [start_year, end_year] per site
    
    info_file = './site_info_FLUXNET.xlsx';
    
    if ~exist(info_file, 'file')
        error('Site info file not found: %s', info_file);
    end
    
    % Read site names (column A) and year ranges
    [~, site_names] = xlsread(info_file, 'A2:A300');
    year_range = xlsread(info_file);
    
    % Validate data consistency
    if size(year_range,1) ~= length(site_names)
        error('Site/year range mismatch in input file');
    end
    year_range = year_range(:,1:2); 
end

function composite_8day = initialize_composite_output(site_names, year_range)
    % Initialize output matrix for 8-day composites
    % Matrix dimensions: [total_years*46 × num_sites]
    
    min_year = min(year_range(:,1));
    max_year = max(year_range(:,2));
    total_years = max_year - min_year + 1;
    num_sites = length(site_names);
    
    % Preallocate with NaN (46 periods/year × total years)
    composite_8day = nan(46 * total_years, num_sites); 
end

function composite_8day = process_all_sites(site_names, year_range, composite_8day)
    % Process data for all sites with progress tracking
    
    min_year = min(year_range(:,1));
    num_sites = length(site_names);
    
    fprintf('Processing %d sites...\n', num_sites);
    
    for j = 1:num_sites
        fprintf('Site %d/%d: %s... ', j, num_sites, site_names{j});
        tic;
        
        try
            % Build data file path
            data_file = sprintf('./FLUXSET_DD/FLX_%s_FLUXNET2015_FULLSET_DD.csv', site_names{j});
            
            if ~exist(data_file, 'file')
                warning('Missing data file: %s', data_file);
                continue;
            end
            
            % Read variable names from header
            [~, var_names] = xlsread(data_file);
            
            % Find target columns (GPP and QC)
            [target_col, qc_col] = find_target_columns(var_names);
            
            if isempty(target_col)
                warning('GPP columns not found in site %s', site_names{j});
                continue;
            end
            
            % Process site data
            composite_8day = process_site_data(data_file, site_names{j}, ...
                year_range(j,:), target_col, qc_col, min_year, composite_8day, j);
            
            fprintf('Done (%.1f sec)\n', toc);
        catch ME
            fprintf('ERROR processing site %s: %s\n', site_names{j}, ME.message);
        end
    end
end

function [target_col, qc_col] = find_target_columns(var_names)
    % Identify column indices for target variables
    % Returns:
    %   target_col - Column index for GPP_NT_VUT_REF and GPP_DT_VUT_REF
    %   qc_col     - Column index for corresponding QC flag
    
    target_col = find(strcmpi(var_names, 'GPP_NT_VUT_REF'), 1); % GPP_DT_VUT_REF
    qc_col = find(strcmpi(var_names, 'GPP_NT_VUT_REF_QC'), 1);
    
    % Validate column existence
    if isempty(target_col) || isempty(qc_col)
        warning('Required columns not found in data file');
    end
end

function composite_8day = process_site_data(data_file, site_name, site_years, ...
                          target_col, qc_col, min_year, composite_8day, site_idx)
    % Process individual site data and generate 8-day composites
    
    % Read numerical data
    data = xlsread(data_file);
    
    % Extract year from timestamp (YYYYMMDD format)
    data(:,1) = floor(data(:,1)./10000);
    
    % Process each year
    for year = site_years(1):site_years(2)
        year_data = data(data(:,1) == year, :);
        
        % Skip incomplete years
        if size(year_data, 1) < 365
            warning('Incomplete data for %s year %d (%d days)', ...
                   site_name, year, size(year_data,1));
            continue;
        end
        
        % Apply QC filters:
        % 1. Remove negative values
        % 2. Remove QC < 0.2 or missing values
        daily_values = year_data(1:365, target_col);
        daily_qc = year_data(1:365, qc_col);
        daily_values(daily_values < 0 | daily_qc < 0.2 | daily_values == -9999) = nan;
        
        % Generate 8-day composites
        yearly_composite = create_8day_composite(daily_values);
        
        % Store in output matrix
        year_offset = year - min_year;
        start_row = year_offset * 46 + 1;
        end_row = (year_offset + 1) * 46;
        composite_8day(start_row:end_row, site_idx) = yearly_composite;
    end
end

function composite_data = create_8day_composite(daily_values)
    % Generate 8-day composites from daily data with quality control
    % Rules:
    % Each composite requires ≥2 valid daily measurements
    
    composite_data = nan(46, 1); % Preallocate with NaN
    valid_threshold = 2;         % Minimum valid days per composite
    
    for period = 1:46
        start_day = (period-1)*8 + 1;
        end_day = min(period*8, 365);
        
        period_data = daily_values(start_day:end_day);
        valid_data = period_data(~isnan(period_data));
        
        % Only average if sufficient valid data exists
        if numel(valid_data) >= valid_threshold
            composite_data(period) = mean(valid_data, 'omitnan');
        end
    end
end