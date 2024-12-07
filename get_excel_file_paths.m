function [path1, path2, path3] = get_excel_file_paths()
    % This function returns the full paths of three specific Excel files.
    
    % Define the file names (inputs inside the function)
    file1 = 'D:\DataWMask.xlsx'; % Dataset file path
    file2 = 'C:\Users\acer\Desktop\breast image\'; % out put folder path
    file3 = 'data3.xlsx'; % Example file
    
    % If the files do not contain full paths, assume they are in the current directory
    if ~contains(file1, '\') && ~contains(file1, '/')
        path1 = fullfile(pwd, file1); % Add current directory path
    else
        path1 = file1; % If full path is already provided
    end

    if ~contains(file2, '\') && ~contains(file2, '/')
        path2 = fullfile(pwd, file2); % Add current directory path
    else
        path2 = file2; % If full path is already provided
    end

    if ~contains(file3, '\') && ~contains(file3, '/')
        path3 = fullfile(pwd, file3); % Add current directory path
    else
        path3 = file3; % If full path is already provided
    end
    
    % Display the paths (optional)
    disp('The full paths of the Excel files are:');
    disp(path1);
    disp(path2);
    disp(path3);
end
