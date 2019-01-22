function [Template,Source] = loadImages(path)
% takes as input a path to the folder containting the images and loads them into the workspace. It is
% important to note that the provided file path must contain two files
% named Template.dcm and source.dcm
% Read In Images
templatePath = path + "/Template.dcm";
sourcePath = path + "/Source.dcm";

Template = im2double(dicomread(templatePath));
Source = im2double(dicomread(sourcePath));

Template =  imsharpen(Template,'Radius', 10, 'Amount', 10)
Source =  imsharpen(Source,'Radius', 10, 'Amount', 10)

end

