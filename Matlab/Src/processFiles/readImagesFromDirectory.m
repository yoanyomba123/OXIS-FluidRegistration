function images = readImagesFromDirectory(directory)
    % obtain pattern to match files
    files = dir(directory);
    % read all images and store in a struct
    filesLen = length(files);
    
    % create cell struct
    imageDataStore = cell(1, filesLen);
    
    for i = 1: filesLen
        path = fullfile(files(i).folder + "/" + files(i).name); 
        imageDataStore{i} = im2double(dicomread(path)); 
    end
    images = imageDataStore;
end

