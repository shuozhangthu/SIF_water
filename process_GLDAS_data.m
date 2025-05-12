function process_GLDAS_data()
    % Main function to process GLDAS data: aggregate 3-hourly to 8-day data
    
    % 1. Configuration parameters
    inputDir = './input/';
    outputDir = './output/';
    variableNames = {'Rainf', 'Tair'};  % Variables to process (precipitation and temperature)
    startDate = datetime(2000,1,1);     % Start date
    endDate = datetime(2020,12,31);     % End date
    
    % 2. Create output directory if not exists
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    % 3. Process each variable
    for varIdx = 1:length(variableNames)
        currentVar = variableNames{varIdx};
        fprintf('Processing variable: %s\n', currentVar);
        
        fileList = dir(fullfile(inputDir, ['*', currentVar, '*.nc']));  % Assuming NetCDF format
        
        aggregateDataBy8Days(fileList, currentVar, outputDir, startDate, endDate);
    end
    
    fprintf('All processing completed!\n');
end

function aggregateDataBy8Days(fileList, varName, outputDir, startDate, endDate)
    
    allDates = startDate:endDate;
    eightDayIntervals = 1:8:length(allDates);
    
    for i = 1:length(eightDayIntervals)
        startIdx = eightDayIntervals(i);
        endIdx = min(startIdx+7, length(allDates));
        currentDates = allDates(startIdx:endIdx);
        
        fprintf('Processing 8-day period from %s to %s\n', ...
            datestr(currentDates(1)), datestr(currentDates(end)));
        
        periodFiles = {};
        for d = 1:length(currentDates)
            currentDate = currentDates(d);
            dateStr = datestr(currentDate, 'yyyymmdd');
            pattern = ['*', dateStr, '*.nc'];
            
            matchedFiles = dir(fullfile(fileList(1).folder, pattern));
            periodFiles = [periodFiles; fullfile({matchedFiles.folder}, {matchedFiles.name})'];
        end
        
        if ~isempty(periodFiles)
            aggregatedData = aggregateFilesData(periodFiles, varName);
            
            saveAggregatedData(aggregatedData, varName, outputDir, ...
                currentDates(1), currentDates(end));
        end
    end
end

function aggregatedData = aggregateFilesData(filePaths, varName)
    
    firstFile = filePaths{1};
    varInfo = ncinfo(firstFile, varName);
    dataSize = varInfo.Size;
    
    allData = zeros([dataSize, length(filePaths)], 'single');
    
    for i = 1:length(filePaths)
        allData(:,:,:,i) = ncread(filePaths{i}, varName);
    end
    
    switch varName
        case 'Rainf'  
            % Convert kg/m²/s to mm/3hour: 1 kg/m² = 1 mm, then multiply by seconds
            allData = allData * 3600 * 3;
            aggregatedData = sum(allData, 4, 'omitnan');
            
        case 'Tair'   
            % Convert Kelvin to Celsius
            allData = allData - 273.15;
            aggregatedData = mean(allData, 4, 'omitnan');
            
    end
end

function saveAggregatedData(data, varName, outputDir, startDate, endDate)
    % Save aggregated data to NetCDF file
    
    startStr = datestr(startDate, 'yyyymmdd');
    endStr = datestr(endDate, 'yyyymmdd');
    outputFile = fullfile(outputDir, [varName, '_8day_', startStr, '_', endStr, '.nc']);
    
    nccreate(outputFile, varName, 'Dimensions', {'lon', size(data,1), 'lat', size(data,2), 'time', 1});
    ncwrite(outputFile, varName, data);
    
    switch varName
        case 'Rainf'
            ncwriteatt(outputFile, varName, 'units', 'mm/8day');
            ncwriteatt(outputFile, varName, 'long_name', '8-day accumulated precipitation');
        case 'Tair'
            ncwriteatt(outputFile, varName, 'units', 'degC');
            ncwriteatt(outputFile, varName, 'long_name', '8-day average air temperature');
    end
    
    fprintf('Saved: %s\n', outputFile);
end