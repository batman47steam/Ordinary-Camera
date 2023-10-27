%% Simulate the A matrix using forward model
%% Charles Saunders and John Murray-Bruce at Boston University

function [A] = simulate_A(wallparam, occ_corner, simuParams, Mon_xdiscr, Mon_zdiscr, Monitor_depth)
%% Simulates the full A matrix. GPU accelerated

wall_matr = wallparam.wall_matr;
wall_point = wallparam.wall_point;
wall_vector_1 = wallparam.wall_vector_1;
wall_vector_2 = wallparam.wall_vector_2;
wall_normal = wallparam.wall_normal;
walln_points = wallparam.walln_points;

NumBlocks_row = simuParams.NumBlocks(1);
NumBlocks_col = simuParams.NumBlocks(2);
Ndiscr_mon = simuParams.Ndiscr_mon; % 每个大的patch又被离散化成了6x6=36个point

blockcount = 1;

A = (zeros(walln_points^2, NumBlocks_row*NumBlocks_col,'gpuArray')); % A是由墙面上离散化的点，和显示器的patch数目决定的

for mc = NumBlocks_col:-1:1 %For each scene patch 以block为主导的，也就是显示屏为主导, col从最大的index开始的 ？
    for mr = 1:NumBlocks_row % row从最小的index开始的
        lightposy = Monitor_depth;

        % View angle model
        if simuParams.viewAngleCorrection==1
            MM = ViewingAngleFactor([Mon_xdiscr(mc)-simuParams.IlluminationBlock_Size(1)/2, lightposy ,Mon_zdiscr(mr)-simuParams.IlluminationBlock_Size(2)/2],wall_matr(1,:), wall_matr(3,end:-1:1), simuParams.D);
        elseif simuParams.viewAngleCorrection==0
            MM = ones(length(wall_matr(1,:)),length(wall_matr(3,end:-1:1)));
        end

        % Discretize a monitor block 把每个block离散化为36个点
        [lightposx, lightposz] = meshgrid( Mon_xdiscr(mc):-simuParams.IlluminationBlock_Size(1)/Ndiscr_mon:Mon_xdiscr(mc)-simuParams.IlluminationBlock_Size(1)+(simuParams.IlluminationBlock_Size(1)/Ndiscr_mon)/10, Mon_zdiscr(mr):-simuParams.IlluminationBlock_Size(2)/Ndiscr_mon:Mon_zdiscr(mr)-simuParams.IlluminationBlock_Size(2)+(simuParams.IlluminationBlock_Size(2)/Ndiscr_mon)/10); 
        % Simulate the block，虽然一个block里面有36个点，但是更宏观的角度其实就等价于每次只亮一个像素点
        %image = simulate_block(wall_matr,wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, [lightposx(:),ones(Ndiscr_mon^2,1).*lightposy,lightposz(:)], occ_corner, repmat(MM,[1,1,Ndiscr_mon^2]));
        image = flipud(simulate_block(wall_matr,wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, [lightposx(:),ones(Ndiscr_mon^2,1).*lightposy,lightposz(:)], occ_corner, repmat(MM,[1,1,Ndiscr_mon^2])));
        %imshow(image, [])
        % Add to A matrix as column
        A(:,blockcount) = image(:); % 每次只亮一个像素点得到的光传输矩阵的每一列
        blockcount = blockcount + 1;
        
    end
end

A = gather(A);

end

        


function [M] = ViewingAngleFactor(MonitorPixel_xyz, FOV_xdiscr,FOV_zdiscr, D)
powexp = 18;%18;%20; %5;
ang = atan((MonitorPixel_xyz(3)-FOV_zdiscr)./(D));
Mz = cos(ang).^powexp;

Mx = cos(atan((MonitorPixel_xyz(1)-FOV_xdiscr)./(D))).^1;

M = (Mx'*Mz)';
end

