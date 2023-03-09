function PreProcData(matFilename)
% 输入图像数据，初始化并赋值data
global data;
global data2;
global data3;
global data4;

data = 0;
temp = load(matFilename{1});
fieldName = fieldnames(temp);
eval(['data = temp.' fieldName{1} ';']);
data = (2^12-1) * ((data - min(data,[],'all')) / (max(data,[],'all')-min(data,[],'all')));

if length(matFilename) > 1
    data2 = 0;
    temp = load(matFilename{2}, 'data');
    fieldName = fieldnames(temp);
    eval(['data2 = temp.' fieldName{1} ';']);
    data2 = (2^12-1) * ((data2 - min(data2,[],'all')) / (max(data2,[],'all')-min(data2,[],'all'))) .^ 1.2;
end

if length(matFilename) > 2
    data3 = 0;
    temp = load(matFilename{3}, 'data');
    fieldName = fieldnames(temp);
    eval(['data3 = temp.' fieldName{1} ';']);
    data3 = (2^12-1) * ((data3 - min(data3,[],'all')) / (max(data3,[],'all')-min(data3,[],'all'))) .^ 1.2;
end

if length(matFilename) > 3
    data4 = 0;
    temp = load(matFilename{4}, 'data');
    fieldName = fieldnames(temp);
    eval(['data4 = temp.' fieldName{1} ';']);
    data4 = (2^12-1) * ((data4 - min(data4,[],'all')) / (max(data4,[],'all')-min(data4,[],'all'))) .^ 1.2;
end
