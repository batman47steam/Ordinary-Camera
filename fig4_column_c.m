%% Scene reconstruction using differencing method (Fig4 column c)
%
% ------------------------------Pseudo-code------------------------------
% 1. Simulate light transport matrix A given scene geometry parameters.
%    
% 2. Solve TV regularized optimization problem to reconstruct scene
% (FISTA).
% -----------------------------------------------------------------------

% Last Modified by Charles Saunders at Boston University
% 09-Nov-2018 (Clean-up and commented for sharing)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.

%%
% Functions
addpath('Functions')

clear variables;

TestLetter = '10-25'; % Picks the test data to use
numPixels = 1280; % Number of pixels in camera measurement，像素数目由什么决定 (RGBG pattern转换以后就是1008*1008*3)

% Parameters
%Ndiscr_mon = 6; %Discretization of each scene patch,每个block里面又离散化为6个点
Ndiscr_mon = 6;
downsamp_factor = 0.1 % imresize的降采样率，直接乘上去的
%downsamp_factor = 3; %Downsampling of measurements 2^downsamp_factor
viewAngleCorrection = 0;  %True/False
useEstimatedOccPos = 0; %Use estimated occluder position or not

load_experiment_config_data_reconstruction % 这个是function里面的一个函数，还不知道干了什么，反正把东西都加载出来了

% Data path
if ismac
    calibParams.filepath = './Data/TestPosD11/';
elseif ispc
    calibParams.filepath = '.\Data\TestPosD11\';
end
%%%%%%% Setup %%%%%%%% (Nothing to manually set here)

% Wall/imaging plane wall point不就是FOV中心的位置
wall_point = [FOV_LLCorner(1) + FOV_size(1)/2,D,FOV_LLCorner(2)+FOV_size(2)/2];  %Point on plane D是相机和墙之间的距离，到底谁是原点
wall_vector_1 = [FOV_size(1)/2,0,0]; %Vector defining one direction of FOV (and extent)
wall_vector_2 = [0,0,FOV_size(2)/2]; %Vector defining the orthogonal direction (and extent)
wall_normal = cross(wall_vector_1,wall_vector_2); % 三维空间向量的叉积，确定法线
wall_normal = wall_normal./norm(wall_normal);

walln_points = floor(numPixels * downsamp_factor); % 主要目的是为了和加载进来的ROI图片的像素数目一致
%walln_points = floor(numPixels/(2^downsamp_factor)); %Number of points to render in each direction，这个是最后相机上拍了多少点吗

% Discretize imaging plane
f_imageplane = gpuArray((zeros(walln_points)));
wall_vec = (-1:2/(walln_points-1):1); % 以2/（walln_points-1) 为间隔，在（-1，1）之间生成walln_points-1个点
wall_matr(1,:) = gpuArray((wall_point(1) + wall_vec*wall_vector_1(1) + wall_vec*wall_vector_2(1))); % Fov_1 Fov_1 0
wall_matr(2,:) = gpuArray((wall_point(2) + wall_vec*wall_vector_1(2) + wall_vec*wall_vector_2(2))); % D     0     0
wall_matr(3,:) = gpuArray((wall_point(3) + wall_vec*wall_vector_1(3) + wall_vec*wall_vector_2(3))); % Fov_2 0     Fov_2

% Mon_offset(1)我理解的是屏幕有效显示区域low-right-corner的坐标, 这些block块谁先发光有影响吗
Monitor_xlim = [0 NumBlocks_col]*IlluminationBlock_Size(1) + Mon_Offset(1); %  IlluminationBlock对应每个35x35像素块的实际大小
Monitor_y = 0; % 相当于显示屏的位置是y轴的起点
Monitor_zlim = [0 NumBlocks_row]*IlluminationBlock_Size(2) + Mon_Offset(2);
Mon_xdiscr = (linspace(Monitor_xlim(1),Monitor_xlim(2),NumBlocks_col));
Mon_zdiscr = (linspace(Monitor_zlim(2),Monitor_zlim(1),NumBlocks_row)); % 指定了元素的个数, z这里颠倒了一下，所以应该是从最左上角开始的

wallparam.wall_matr = wall_matr;
wallparam.wall_point = wall_point;
wallparam.wall_vector_1 = wall_vector_1;
wallparam.wall_vector_2 = wall_vector_2;
wallparam.wall_normal = wall_normal;
wallparam.walln_points = walln_points;

numPixels = floor(numPixels * downsamp_factor+1)
%numPixels = floor(numPixels/(2^downsamp_factor));

%% Load data

% Available hidden test scenes data choose one [scene name: 'key']
%           * RGB bars scene:       'rgb'
%           * Text 'BU' scene:      'bu'
%           * Mushroom scene:       'mushroom'
%           * Tommy scene:          'tommy'

scene = 'T';

switch scene
    case 'mushroom' % 只有求逆的时候才会加载拍摄的图像，前面估计A的时候根本不需要拍摄的结果
        [test_image1,ground_truth1]=load_image1('image_test_mushroom20.mat',calibParams.filepath,downsamp_factor);   
        Occ_LLcorner = [0.4583 0.5408 0.2026]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e08 * [0.51    0.561    2.448]; %TV regularization parameter
        
        %Discretization dependent scaling of A matrix
        sr = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.5313e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 16250/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'tommy'
        [test_image1,ground_truth1]=load_image1('image_test_smilehat20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4569 0.5744 0.2080]; %Estimated occluder position from localization script 也是遮挡物最下角的坐标
        
        tv_reg_param = 1e07 * [5.72    6.76    5.72];  %TV regularization parameter
        
         %Discretization dependent scaling of A matrix
        sr = 1.1406e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.3594e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.9063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'bu'
        [test_image1,ground_truth1]=load_image1('image_test_bur20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4733   0.5661    0.2072]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e07*[5   25   25];  %TV regularization parameter
        
         %Discretization dependent scaling of A matrix
        sr = 12500/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 15625/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.7188e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'rgb'
        [test_image1,ground_truth1]=load_image1('image_test_colbar20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4693 0.5629 0.2080]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e06 * [52.5   50   47.5];  %TV regularization parameter
        
        %Discretization dependent scaling of A matrix
        sr = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'T'
        % 先随便的瞎读一下，反正simulateA的时候也不需要图片信息啊！
        %[test_image1,ground_truth1]=load_image1('image_test_colbar20.mat',calibParams.filepath,downsamp_factor);
        %Occ_LLcorner = [0.4693 0.5629 0.2080]; %Estimated occluder position from localization script
        %Occ_LLcorner = [1 0.8 0.912]
        Occ_LLcorner = [1 0.8 0.92]
        %Occ_LLcorner = [1.0147, 0.8207, 0.9181];
        % tv正则化先也不需要去管
        tv_reg_param = [1000   10   10];  %TV regularization parameter
        
        %Discretization dependent scaling of A matrix
        sr = 180.708/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 8.708/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sb = 8.708/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9;

        camera_capture = double(imread("pattern2.tif"));
        camera_capture = camera_capture(480:1568, 480:1568); % 这里的测量值也是人为的减去背景了
        camera_capture(camera_capture<0) = 0;
        camera_capture = fliplr(imresize(camera_capture, 0.2, 'bilinear')); % 坐标系统一，面向屏幕
        test_image1 = camera_capture;
        [r, c] = size(test_image1);
        rgbImage = zeros(r, c, 3);
        rgbImage(:, :, 1) = test_image1;
        rgbImage(:, :, 2) = test_image1;
        rgbImage(:, :, 3) = test_image1;
        test_image1 = rgbImage;
        %figure();
        %imshow(camera_capture, []);
        
end

%%
% Occluder Occ_size对应flat障碍物的长，宽，高 
occ_corner(1,:,1) = Occ_LLcorner; % 从这个来看，应该也是最下方处顶点的位置
occ_corner(2,:,1) = Occ_LLcorner + [Occ_size(1), 0, 0];
occ_corner(3,:,1) = Occ_LLcorner + [Occ_size(1), 0, Occ_size(3)];
occ_corner(4,:,1) = Occ_LLcorner + [0, 0, Occ_size(3)];

occ_corner(1,:,2) = Occ_LLcorner + [Occ_size(1)/2-0.0035, 0, 0]; % 这个是板子下面的腿
occ_corner(2,:,2) = Occ_LLcorner + [Occ_size(1)/2+0.0035, 0, 0];
occ_corner(3,:,2) = Occ_LLcorner + [Occ_size(1)/2-0.0035, 0, -Occ_LLcorner(3)]; % Occ_LLcorner(3) - Occ_LLcorner(3) = 0 所以z的起点就是地面啊
occ_corner(4,:,2) = Occ_LLcorner + [Occ_size(1)/2+0.0035, 0, -Occ_LLcorner(3)];
%%%%%%%%%%%%%% 上面是一个大的长方形，下面是一个小的长方形支撑着


%% Simulate Transport Matrix
disp('Simulating transport matrix...') % 关注下需要提供哪些参数，最后一个参数是moniter depth
%[simA] = simulate_A(wallparam, (occ_corner),simuParams, Mon_xdiscr,Mon_zdiscr, 0);
simA = load("nlos64.mat").simA;
%save('nlos128.mat', 'simA');

% 加载隐藏图片
hidden = imread("mushroom64.png");
imshow(hidden, [])
hidden = double(hidden)
hidden = hidden/255;
imshow(double(hidden));
% 分别对三个通道进行列化
r = hidden(:,:,1);
r = r(:);
g = hidden(:,:,2);
g = g(:);
b = hidden(:,:,3);
b = b(:);

% 读入手机拍摄的环境光下的background
background = imread("background128.jpg");

% h_r_mean = mean(r);
% h_g_mean = mean(g);
% h_b_mean = mean(b);

y1 = simA * double(r);
y2 = simA * double(g);
y3 = simA * double(b);

% proj_r_mean = mean(y1);
% proj_g_mean = mean(y2);
% proj_b_mean = mean(y3);

y1 = reshape(y1, [128 128]);
% figure();
% imshow(y1, []);
y2 = reshape(y2, [128 128]);
figure();
imshow(y2, []);
y3 = reshape(y3, [128 128]);
% figure();
% imshow(y3, []);
y(:,:,1) = y1;
y(:,:,2) = y2;
y(:,:,3) = y3;

% 先往measurement的结果上添加SBR
y = y + 10 * double(background);


figure();
imshow(y(:,:,2), []);


figure()
y = y/max(max(max(y)));
noisy_image = awgn(y, 50, 'measured'); % 加了measured，确保最后整体的psnr是50db
proj= noisy_image;
gt = hidden;
save("mushroom64_b.mat", 'simA', 'proj', 'gt');


% % 求伪逆测试结果
% A_inv = pinv(simA);
% y_r = noisy_image(:,:,1)
% y_g = noisy_image(:,:,2)
% y_b = noisy_image(:,:,3)
% 
% % 利用最小二乘去求初始结果
% r_init = A \ y_r(:)
% g_init = A \ y_g(:)
% b_init = A \ y_b(:)
% %x = A_inv * y1(:)
% x = simA \ y1(:);
% 
% %figure()
% %imshow(y, []); % 这里本身也就只是用了简单的Af, 还没有涉及到噪声的部分Af + b
% 
% % [r, c] = size(y);
% % rgbImage = zeros(r, c, 3);
% % rgbImage(:, :, 1) = y;
% % rgbImage(:, :, 2) = y;
% % rgbImage(:, :, 3) = y;
% % test_image1 = rgbImage;
% 
% pattern = imread("pattern2.tif");
% pattern = pattern(480:1568, 480:1568);
% pattern(pattern<0) = 0;
% pattern = fliplr(imresize(pattern, 0.2, 'bilinear')); % 就是说你的拍摄结果也要flip一下，改成是在墙后面看到的
% pattern = double(pattern);
% %pattern = imresize(pattern, 0.2, 'bilinear');
% % 加载由标定得到的A
% [simA] = load('A.mat');
% reconstruct = double(simA) \ double(pattern(:));
% restore = reshape(reconstruct, [8,8]);
% 
% save('nlos.mat', 'pattern', 'simA') % 直接保存成结构体倒是挺炫酷的
% 
% %imshow(fliplr(restore)); % 结果这里左右是否进行flip应该也无所谓
% figure()
% imshow(restore);



%% Reconstruction 
% sr, sg, sb 这些是进行重建的时候提供的，得到光传输矩阵以后，通过最优化的方法进行逆运算
%final_im1 = reconstruct_tv_it_cbg(simA,[sr,sg,sb],  test_image1, tv_reg_param, NumBlocks_sim, [0,0,0]); % test_image1换成restore看下
final_im1 = reconstruct_tv_it_cbg(simA,[sr,sg,sb],  proj, tv_reg_param, NumBlocks_sim, [0,0,0], gt); % test_image1换成restore看下
%% Plots
figure()

subplot(1,2,1)
%imshow(ground_truth1/255)
%title('Ground truths')
%subplot(1,2,2)
imshow(final_im1(:,:,1), [])
