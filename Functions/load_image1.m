function [image,ground_truth] = load_image1(file,path,downsampling)
% Loads an image, removes background and splits into colour channels


iml = load([path,file]); % 加载进来的是一个结构体
imload = double(iml.image); % imload是相机拍摄的结果
ground_truth = iml.ground_truth; % ground_truth是需要重建的目标
measure_shut = get_shutter_speed(iml); % 曝光相关的参数也加载进来了 3500



if ~exist('calshut','var')
    calshut = measure_shut;
end

imcol = get_color_image(imload); % 应该是把RGBG pattern的变为三通道的RGB图片

testimr = (calshut/measure_shut)*downsample_2(imcol(:,:,1),downsampling); % 降采样有点小奇怪，没细看, 16倍降采样
testimg = (calshut/measure_shut)*downsample_2(imcol(:,:,2),downsampling);
testimb = (calshut/measure_shut)*downsample_2(imcol(:,:,3),downsampling);

image(:,:,1) = testimr; % 分别取出每个通道对应RGB
image(:,:,2) = testimg;
image(:,:,3) = testimb;
end

function [sp] = get_shutter_speed(iml)

speed = iml.shutter_speed;

if isfield(iml,'averaged_iterations')
    sp = speed*iml.averaged_iterations;
else
    sp = speed*iml.iterations;
    
end
end

function [col_image] = get_color_image(im)
	red = im(1:2:end,1:2:end);
    green = (im(2:2:end,1:2:end)+im(1:2:end,2:2:end))/2;
    blue = im(2:2:end,2:2:end);
    
    col_image(:,:,1) = red;
    col_image(:,:,2) = green;
    col_image(:,:,3) = blue;
end

% 用这种方法来实现下采样，其实就是间隔一定的距离去取值
function [image] = downsample_2(im,iterations)
    %Downsampling 2^iterations times without averaging

    for i=1:iterations
        im1 = im(1:2:end,1:2:end);
        im2 = im(1:2:end,2:2:end);
        im3 = im(2:2:end,1:2:end);
        im4 = im(2:2:end,2:2:end);

        image = (im1+im2+im3+im4)/4;
        im = image;
    end
    
    if iterations == 0 
        image = im;
    end
end
