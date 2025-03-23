% 参数定义
Lx = 0.1; % 电池壳长度 (m)
Ly = 0.05; % 电池壳宽度 (m)
Lz = 0.02; % 电池壳高度 (m)
nx = 20; % x方向网格数
ny = 10; % y方向网格数
nz = 5; % z方向网格数
k = 200; % 热导率 (W/m·K)
T0 = 20; % 初始温度 (°C)
T_left = 100; % 左边界温度 (°C)
T_right = 50; % 右边界温度 (°C)
T_top = 80; % 上边界温度 (°C)
T_bottom = 30; % 下边界温度 (°C)

% 网格生成
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
z = linspace(0, Lz, nz);
[X, Y, Z] = meshgrid(x, y, z);

% 初始化温度场
T = T0 * ones(ny, nx, nz);

% 边界条件
T(:, 1, :) = T_left; % 左边界
T(:, end, :) = T_right; % 右边界
T(1, :, :) = T_top; % 上边界
T(end, :, :) = T_bottom; % 下边界

% 迭代求解
max_iter = 1000;
tol = 1e-6;
for iter = 1:max_iter
    T_old = T;
    for i = 2:ny-1
        for j = 2:nx-1
            for k = 2:nz-1
                T(i, j, k) = 0.1667 * (T_old(i+1, j, k) + T_old(i-1, j, k) + ...
                             T_old(i, j+1, k) + T_old(i, j-1, k) + ...
                             T_old(i, j, k+1) + T_old(i, j, k-1));
            end
        end
    end
    % 检查收敛
    if max(abs(T(:) - T_old(:))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end

% 结果可视化 (选择一个截面)
figure;
slice(X, Y, Z, T, Lx/2, Ly/2, Lz/2); % 在x, y, z的中点处切片
colorbar;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Temperature Distribution in Battery Shell');