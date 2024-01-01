function [ simA, Discr ] = SimulateA_OccluderEstimation(simuParams, downsamp_fact)
%GETLINEARFORWARDMODEL Simulates the forward operator (measurement matrix)
% for the monitor/camera experimental setup, using the set parameters
% provided in simulationParams (struct).
%
%   Usage:
%   [ simA, Discr ] = SimulateMonitorA4OccluderEst(simuParams, downsamp_fact)
%   Input:
%       * simuParams (struct obj):  Containing the parameters used for simulating the forward model A.
%                                   -simuParams.Occluder: occluder position in 3D-space
%                                   -simuParams.NumBlocks: occluder position in 3D-space
%                                   -simuParams.Ndiscr_mon: occluder position in 3D-space
%                                   -simuParams.numPixels: occluder position in 3D-space
%                                   -simuParams.D: occluder position in 3D-space
%                                   -simuParams.viewAngleCorrection: occluder position in 3D-space
%                                   -simuParams.IlluminationBlock_Size: occluder position in 3D-space
%                                   -simuParams.Mon_Offset: LCD position in (x,z)-plane
%       * downsamp_fact:            Amount of downsampling applied onto the 
%                                   camera measurement (this is reflected 
%                                   in the no. of rows of A).
%   Output:
%       * simA (N x Ndiscr_mon matrix): The forward model, transport matrix.
%       * Discr (matrix):               Camera FOV discretization (x,z)-plane on visible wall.

% 看返回值，最后返回的就是当前的这些参数，在FOV上形成的图案

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 11-Oct-2017 11:04: 
% v2.0 07-Nov-2018 16:34:07 (Clean-up and commented for sharing -JMB)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.

NumBlocks = simuParams.NumBlocks;
numPixels = simuParams.numPixels;
if nargin>1
    %numPixels = floor(numPixels * downsamp_fact + 1);
    numPixels = floor(numPixels/(2^downsamp_fact)); % ds_factor也就是用在这里计算pixel的数目
end
nDiscr = simuParams.Ndiscr_mon; % 对应前面的6吗？
D = simuParams.D;
FOV_cord = simuParams.FOV_cord;
Occluder = simuParams.Occluder; % 这个occluder不需要转换吗，没得到坐标轴上的实际位置啊
viewAngleCorrection = simuParams.viewAngleCorrection;
IlluminationBlock_size = simuParams.IlluminationBlock_Size;
Mon_Offset = simuParams.Mon_Offset;

[rowCount, colCount] = meshgrid(1:(NumBlocks(1)),1:(NumBlocks(2)));
rowCount = rowCount'; colCount = colCount';
% Note this not so natural ordering used below, due to the file naming
% convention for the monitor blocks used by Charlie's capture code.
ActiveBlock = [rowCount(:), colCount(:)]; % 先是row后是col，固定行，遍历列 ？

% Preallocate memory for array.
totalBlocks = NumBlocks(1)*NumBlocks(2);
simA = zeros(numPixels^2,totalBlocks);


% SIMULATE CAMERA MEASUREMENTS

parfor jj=1:totalBlocks %parfor 这里其实就是在构建对应的光传输矩阵
    Im = SimulateForwardModelPerBlock(NumBlocks, nDiscr, ActiveBlock(jj,:),...
        numPixels, D, FOV_cord, [Occluder],viewAngleCorrection,...
        IlluminationBlock_size, Mon_Offset);
    Im = flipud(Im);
    simA(:,jj) = Im(:); % 本质上干的事情一致
end
% for jj=1:totalBlocks%parfor 
%     Im = SimulateForwardModelPerBlock(NumBlocks, nDiscr, ActiveBlock(jj,:),...
%         numPixels, D, FOV_cord, [Occluder],viewAngleCorrection,...
%         IlluminationBlock_size, Mon_Offset);
%     Im = flipud(Im);
%     imshow(Im)
%     simA(:,jj) = Im(:); % 本质上干的事情一致
% end
% hidden = imread("pattern8.png")
% hidden = hidden(:,:,1)
% y = simA * double(hidden(:))
% y = reshape(y, [218 218])
% imshow(y, [])
% save('projection.mat', 'y')
Discr = [linspace(FOV_cord(1,1),FOV_cord(2,1),numPixels);
                    linspace(FOV_cord(1,2),FOV_cord(2,2),numPixels)]; % numPixels只会影响到这里，他其实根本没有去管FOV上的图案

end