%% This script hold the measurement configurations of the experimental tests
% performed. Given the experiment test letter, the correct dataset will be
% loaded for MATLAB processing.

% John Murray-Bruce at Boston University


switch TestLetter
    
    case 'D11'
        filename = 'TestPosD11';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;
        
        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            %% Estimated using occluder localization scripts
            % Projection for mushroom
            Occ_LLcorner = [0.4583 0.5408 0.2026];
            
            %Projection for Tommy
            %Occ_LLcorner = [0.4569 0.5744 0.2080];
            
            % Projection for BU red (bur)
            %Occ_LLcorner = [0.4733   0.5661    0.2072];
          
            % Projection for RGB bars (colbar)
            %Occ_LLcorner = [0.4693 0.5629 0.2080];
            
            %%
        else
            % True
            Occ_LLcorner = [0.475 D-0.460 0.214];
        end
        
        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
        
        %CAMERA CONFIG
        FOV_size = [0.4372 0.4372];
        FOV_LLCorner = [0.521 0.048];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];
        
        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];
        
        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024];
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);
        
        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));
        
        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
       
    case '10-25'
        filename = 'TestPosD11'; % file_name要改成对应的，而且要知道他加载的目的是什么
        subblocksperaxis = [1 1];
        NumBlocks_sim = [64 64].*subblocksperaxis;
        NumBlocks_cal = [64 64];
        D = 1.03;
        
        % OCCLUDER DEFINITION
        Occ_size = [0.077 0 0.077];
        if useEstimatedOccPos
            %% Estimated using occluder localization scripts
            % projection for H & T
            %Occ_LLcorner = [1 0.8 0.912]; % 这么严格，还会检查两个if语句是否一样

            % Projection for mushroom
            Occ_LLcorner = [0.4583 0.5408 0.2026];
            
            %Projection for Tommy
            %Occ_LLcorner = [0.4569 0.5744 0.2080];
            
            % Projection for BU red (bur)
            %Occ_LLcorner = [0.4733   0.5661    0.2072];
          
            % Projection for RGB bars (colbar)
            %Occ_LLcorner = [0.4693 0.5629 0.2080];
            
            %%
        else
            % True
            %Occ_LLcorner = [0.475 D-0.460 0.214];
            Occ_LLcorner = [1 0.8 0.912];
            %Occ_LLcorner = [0.9319 0.7112 0.9595]
            %Occ_LLcorner = [0.7, 0.5250, 0.8750];
        end
        
        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
        
        %CAMERA CONFIG
        FOV_size = [0.1914 0.1914];
        FOV_LLCorner = [1.115+0.08448 0.83+0.08448]; % 因为只在Fov中截取了一小块，所以计算时候的Fov往里面缩小一点
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];
        
        % MONITOR CONFIG, 这个取决于最后他到底是怎么得到每一块光源的
        Mon_Offset = [0.0820 0.756];
        
        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.345 0.195];
        ScreenResolution = [3840 2160];
        NumPixelsPerMonitorBlock = 33;
        PixelSize_m = (ScreenSize./ScreenResolution);
        
        % 这个应该是面对着屏幕的时候从右上角开始显示图片，所以要加上边缘的位置
        % Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        % Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        % 但是我其实是在面对屏幕时候的最左上角显示的，所以Mon_Offset直接就是起点

        % 不能够整除的时候
        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        %Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        % 和多少个像素组成一个block相关的，应该就是这里了
        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
        
    case 'D11Video'
        filename = 'TestPosD11Video';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;
        
        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.4750    0.5583    0.2000];
        else
            Occ_LLcorner = [0.477 D-0.466 0.21];
        end
        
        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
        
        %CAMERA CONFIG
        FOV_size = [0.4672 0.4672];
        FOV_LLCorner = [0.496 0.014];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];
        
        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];
        
        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);
        
        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));
        
        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
        
        
case {'Chair'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.503 D-0.436 0.244];
        else
            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];

    case {'ChairRealScene'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 ,1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
 
    
    case {'ChairAmb'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 ,1]; %Change to [2,2] to double resolution
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03-0.006;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
    
        
 case {'Chair3DScene'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 , 1]; %Change to [2,2] to double resolution
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.40 0.40];
        FOV_LLCorner = [0.557 0.062];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
 
        
    otherwise
        disp('Not available');
end


% Set simulation of forward model parameters!

simuParams.NumBlocks = NumBlocks_sim;
simuParams.Ndiscr_mon = Ndiscr_mon;
simuParams.numPixels = numPixels;
simuParams.D = D;
simuParams.FOV_cord = FOV_cord;
simuParams.Occluder = Occluder;
simuParams.viewAngleCorrection = viewAngleCorrection;

simuParams.IlluminationBlock_Size = IlluminationBlock_Size;
simuParams.Mon_Offset = Mon_Offset;
calibParams.scaling = 1;

if ismac
    calibParams.filepath = ['', filename,'/'];
elseif ispc
    calibParams.filepath = ['', filename,'\'];
end

calibParams.x_max = NumBlocks_cal(2);
calibParams.y_max = NumBlocks_cal(1);

