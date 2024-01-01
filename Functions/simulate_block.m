%% Simulates one monitor block. GPU accelerated.
%% Charles Saunders at Boston University
% wall_mat是墙上的所有点，然后light_source只是显示屏上那一小块上的光源
function [image] = simulate_block(wall_mat,wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos, occ_corner,MM)


    % Squares each element in input
    function [vecsqr] = calc(vec)
        vecsqr =  vec.^2;
    end

    % Vector from light position to wall position.
    % wall_mat([1,3],:)取出第一行和第三行, x的坐标和z的坐标，重复36次
    % light_source_pos(:,[1,3])取出第一列和第三列
    % 可以理解为从屏幕当前的block往墙面上的每一点去引射线，从而构成若干个向量
    vec = repelem(wall_mat([1,3],:),1,1,size(light_source_pos,1)) - repelem(reshape(light_source_pos(:,[1,3])',[2,1,length(light_source_pos)]),1,walln_points,1);
    % Distance squared
    vs = arrayfun(@calc, vec);
    
    % y component is constant, so: 这个简单，因为墙和监视器在y之间的距离就是D
    vs2 = ones(1,walln_points)*((wall_mat(2,1)-light_source_pos(1,2)).^2);

    % Calculate distance squared again (to normalize dot products and for
    % distance losses).
    tmp = arrayfun(@calc,repmat(vs(1,:,:) + vs2,[walln_points,1,1]) + repmat(permute(vs(2,:,:),[2,1,3]),[1,walln_points,1]));
    
    % Final intensity 最最简单的情况下mm相当于全是1，但是到墙上以后会根据距离有一个衰减
    intensity = MM.*repmat(vs2,[walln_points,1,length(light_source_pos)])./tmp;
    
    % Calculate occluder positions
    image = zeros(walln_points,walln_points,'gpuArray');
    if size(occ_corner)>0
    for i = 1:size(vec,3) % 遍历屏幕上那一个block处，所有的点光源？一个block离散为了36个点 occ_corner(:,:,1)是真正的遮挡物 occ_corner(:,:,2)是下面的底座
        occluder_image = simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos(i,:), occ_corner(:,:,1)); % 返回的是逻辑值，代表遮挡和未被遮挡
        for o = 2:size(occ_corner,3) % o 这个索引对应的不就是遮挡物的底座吗，相乘就是遮挡物和遮挡物的底座造成的最终结果
            occluder_image = occluder_image.*simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos(i,:), occ_corner(:,:,o));
        end
        
        image = image + intensity(:,:,i).*occluder_image; % 这样来看，肯定是1对应没有被遮挡，0对应着被遮挡，没被遮挡的就会作用在image上
    end
    else
        image = sum(intensity,3); % 一个block上36个点，所有点在wall上作用后的最终结果
    end

end



