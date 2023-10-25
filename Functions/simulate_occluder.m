%% Generates occluder mask
%% Charles Saunders at Boston University

function [ occluder_image ] = simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos, occ_corner)

    % Find projection of occluder corners onto wall
    occ_corner_proj = zeros(size(occ_corner,1),3);
    for i = 1:size(occ_corner,1) % 遍历遮挡物的每一个角点，start:显示屏光源 end:遮挡物角点 => 遮挡物的角点在墙面上的投影
        occ_corner_proj(i,:) = ray_plane_intersect_y(light_source_pos,occ_corner(i,:),wall_point,wall_normal); % 点光源, 障碍物的角点，FOV的中心点，墙面法线
    end

    % Coordinates in terms of the FOV
    % wall_point是FOV的中心，wall_vector_1其实就是FOV_size(1)/2
    % (wall_point-wall_vector) 其实就是FOV最下角的坐标吧
    % 乘上(walln_points/(wall_vector_1(1)*2)
    % 得到的是在126x126中对应的相对坐标，而不是实际以m为单位的坐标，因为感兴趣的只是相机拍到的那部分，126个点中的哪些个点
    occ_image_coords_1 = (occ_corner_proj(:,1) - (wall_point(1) - wall_vector_1(1))) * walln_points/(wall_vector_1(1)*2);
    occ_image_coords_2 = (occ_corner_proj(:,3) - (wall_point(3) - wall_vector_2(3))) * walln_points/(wall_vector_2(3)*2);
    
    % Generate mask
     k = convhull(occ_image_coords_1,occ_image_coords_2,'simplify', true);
     occluder_image = ~poly2mask(occ_image_coords_1(k),occ_image_coords_2(k),walln_points,walln_points);

    
end