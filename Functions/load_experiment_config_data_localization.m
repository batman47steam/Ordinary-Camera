%% This script details measurements configurations for the experimental tests
% performed.
%

filename = 'TestPosD11';
NumBlocks_sim = [29 36];    % Number of scene patches [hor, vert]
D = 1.03; % 相机和墙面之间的距离

% OCCLUDER DEFINITION flat障碍物的长，宽，高 宽当然就是0了，因为没什么厚度
Occ_size = [0.077 0 0.075]; % Size of the flat occluder (x,y,z)

% Please note here that the oocluder's y-position is stated as (D-p_y), where
% D = 1.03 (D: distance between hidden-scene and visible wall
% planes) and p_y is the distance between the visible wall plane and
% occluder plane.
if useEstimatedOccPos
    % Estimated occluder positions for the different scenes
    %             Occ_LLcorner = [0.4583    D-0.4892    0.2026];  % Mushroom
    %             Occ_LLcorner = [0.4569    D-0.4556    0.2080];  % Tommy
    %             Occ_LLcorner = [0.4733    D-0.4639    0.2072];  % BU red
    %             Occ_LLcorner = [0.4693    D-0.4674    0.2080];  % RGB
else
    Occ_LLcorner = [0.470 D-0.460 0.2040];  % True/measured occluder position 实际测量的障碍物的位置
end
% Define bottom-right and top-left edges of occluder
Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size]; % 上下顶点的坐标

%CAMERA CONFIG
FOV_size = [0.4372 0.4372];
FOV_LLCorner = [0.5128 0.0482];
FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size]; % 和上面的occluder一样，对应上下顶点的坐标

% MONITOR CONFIG
Mon_Offset = [0.0210 0.136];        % (x,z)-position of the lower-right edge of usuable portion of LCD screen 因为显示的图案不一定能占满整个显示屏

NumBlocks_col = NumBlocks_sim(2);   % Number of scene patches (horizontally)
NumBlocks_row = NumBlocks_sim(1);   % Number of scene patches (vertically)
ScreenSize = [0.408 0.3085];        % Size of the LCD display [m] 完整的显示屏显示大小
ScreenResolution = [1280 1024];     % [hor, ver] pixels
NumPixelsPerMonitorBlock = 35;      % Number of monitor pixels in each hidden-scene patch
PixelSize_m = (ScreenSize./ScreenResolution);   % Dimension of each scene patch [m] 每个像素对应为多少米

Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_sim(1)); % mod是除了以后再去取余，其实就是对实际的位置再进行一些修正
Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_sim(2));

IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock; % 就是一个patch块的大小



% Set parameters for simulating the forward model (i.e. A matrix)
simuParams.NumBlocks = NumBlocks_sim;
simuParams.Ndiscr_mon = Ndiscr_mon;
simuParams.numPixels = numPixels;
simuParams.D = D;
simuParams.FOV_cord = FOV_cord;
simuParams.Occluder = Occluder;
simuParams.viewAngleCorrection = viewAngleCorrection;
simuParams.IlluminationBlock_Size = IlluminationBlock_Size;
simuParams.Mon_Offset = Mon_Offset;


if ismac
    datafilepath = './Data/TestPosD11/';
elseif ispc
    datafilepath = '.\Data\TestPosD11\';
end
