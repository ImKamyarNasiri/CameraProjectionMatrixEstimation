image = imread('images/IMG_2216.jpg');

%%%same setting as Q2
world_points = [
    0, 6, 6;
    0, 6, 12;
    0, 12, 6;
    0, 12, 12;

    6, 6, 0;
    12, 12, 0;
    6, 12, 0;
    12, 12, 0;

    6, 0, 6;
    6, 0, 12;
    12, 0, 6;
    12, 0, 12;
];
%normalization
[num_rows, ~] = size(world_points);
avg_world = sum(world_points) / num_rows;
std_world = sqrt(sum(sum((world_points - repmat(avg_world, num_rows, 1)).^2)) / (3*num_rows));

world_points_n = (world_points - repmat(avg_world, num_rows, 1)) ./ std_world;

pixel_points = [
    2881, 1429;
    2602, 1248;
    2823, 1701;
    2560, 1510;

    3392, 1493;
    3612, 1336;
    3314, 1763;
    3543, 1586;

    3167, 1010;
    2869, 826;
    3402, 857;
    3105, 698;
];
%normalization
[num_rows, ~] = size(pixel_points);
avg_pixel = sum(pixel_points) / num_rows;
std_pixel = sqrt(sum(sum((pixel_points - repmat(avg_pixel, num_rows, 1)).^2)) / (2*num_rows));

pixel_points_n = (pixel_points - repmat(avg_pixel, num_rows, 1)) ./ std_pixel;

A = [];
%compute matrix A
for i = 1:size(world_points_n, 1)
    x_world = world_points_n(i, 1);
    y_world = world_points_n(i, 2);
    z_world = world_points_n(i, 3);
    
    x_pixel = pixel_points_n(i, 1);
    y_pixel = pixel_points_n(i, 2);   
    
    element_1 = [
        x_world, y_world, z_world, 1, ...
        0, 0, 0, 0, ...
        -x_pixel * x_world, -x_pixel * y_world, -x_pixel * z_world, -x_pixel
        ];

    element_2 = [
        0, 0, 0, 0, ...
        x_world, y_world, z_world, 1, ...
        -y_pixel * x_world, -y_pixel * y_world, -y_pixel * z_world, -y_pixel
        ];
     
    A = [A; element_1; element_2];
end
%svd of A
[U, S, V] = svd(A);
%pick the eigenvector corresponding to the smallest singular value
[~, index_of_min_singular_value] = min(diag(S));
P_n = V(:, index_of_min_singular_value);
P_n = reshape(P_n, 3, 4);

M1_1 = [
    1/std_pixel, 0, 0;
    0, 1/std_pixel, 0;
    0, 0, 1;
    ];
M1_2 = [
    1, 0, -avg_pixel(1);
    0, 1, -avg_pixel(2);
    0, 0, 1;
    ];

M1 = M1_1 * M1_2;

M2_1 = [
    1/std_world, 0, 0, 0;
    0, 1/std_world, 0, 0;
    0, 0, 1/std_world, 0;
    0, 0, 0, 1;
    ];
M2_2 = [
    1, 0, 0, -avg_world(1);
    0, 1, 0, -avg_world(2);
    0, 0, 1, -avg_world(3);
    0, 0, 0, 1;
    ];

M2 = M2_1 * M2_2;

%final version of P to work on unnormalized data
P = inv(M1) * P_n * M2;
% disp(P);
% disp(P(3, 1:3));
% disp(norm(P(3, 1:3)));
%%normalize w.r.t the last row of P
P = P / norm(P(3, 1:3));

%P_ is the 3x4 matrix, whild P is 3x3
P_ = P;
P = P(:, 1:3);
%compute R_z_theta, the first element of matrix R
theta = atan(P(3,1)/P(3,2));
R_z_theta = [
    cos(theta), sin(theta), 0;
    -sin(theta), cos(theta), 0;
    0, 0, 1;
];
%compute R_x_beta, the second element of matrix R
R_A = P * R_z_theta;
Beta = atan(R_A(3,2)/R_A(3,3));
R_x_beta = [
    1, 0, 0;
    0, cos(Beta), sin(Beta);
    0, -sin(Beta), cos(Beta);
];
%compute R_z_gamma, the third element of matrix R
R_B = R_A * R_x_beta;
gamma = atan(R_B(2,1)/R_B(2,2));
R_z_gamma = [
    cos(gamma), sin(gamma), 0;
    -sin(gamma), cos(gamma), 0;
    0, 0, 1;
];
% R is the transpose of multiplication of these three matrices
R = (R_z_theta * R_x_beta * R_z_gamma)';
%% inverse of R is equal with R transpose, det(R) = 1
% disp(det(R));
% disp(R);
% disp(inv(R));

% The diagonal of PQ should be positive, Q=R'
PQ = P * R_z_theta * R_x_beta * R_z_gamma;
% disp(PQ);

K = [
    3063.73307388903, 0, 2023.71937454517;
    0, 3050.01366910237, 1499.05930779402;
    0, 0, 1;
];

% s will give us: [I | -c]
s = (inv(P) * P_);
% disp(s);
%last column of s
T = - s(:, 4);
disp(norm(T));
% disp(T);

%another way of calculating rotation matrix, based on K
% r = inv(K) * P;
% disp(r);
% disp(det(r));



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Image 2
image_2 = imread('images/IMG_2217.jpg');

% imshow(image_2);
% dcm = datacursormode(gcf);
% set(dcm, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on', 'Enable', 'on');
% pause;
% cursor_info = getCursorInfo(dcm);
% for i = 1:length(cursor_info)
%     x = cursor_info(i).Position(1);
%     y = cursor_info(i).Position(2);
%     pixel_value = impixel(image_2, x, y);
%     disp(['Point ' num2str(i) ': (' num2str(x) ', ' num2str(y) ') - Pixel value: ' num2str(pixel_value)]);
% end

world_points_2 = [
    0, 6, 6;
    0, 6, 12;
    0, 12, 6;
    0, 12, 12;

    6, 6, 0;
    12, 12, 0;
    6, 12, 0;
    12, 12, 0;

    6, 0, 6;
    6, 0, 12;
    12, 0, 6;
    12, 0, 12;
];


[num_rows, ~] = size(world_points_2);
avg_world_2 = sum(world_points_2) / num_rows;
std_world_2 = sqrt(sum(sum((world_points_2 - repmat(avg_world_2, num_rows, 1)).^2)) / (3*num_rows));

world_points_n_2 = (world_points_2 - repmat(avg_world_2, num_rows, 1)) ./ std_world_2;

pixel_points_2 = [
    2684, 1426;
    2480, 1214;
    2652, 1706;
    2468, 1486;

    3240, 1542;
    3560, 1402;
    3184, 1806;
    3492, 1698;

    3048, 1002;
    2832, 794;
    3380, 894;
    3132, 682;
];
[num_rows, ~] = size(pixel_points_2);
avg_pixel_2 = sum(pixel_points_2) / num_rows;
std_pixel_2 = sqrt(sum(sum((pixel_points_2 - repmat(avg_pixel_2, num_rows, 1)).^2)) / (2*num_rows));

pixel_points_n_2 = (pixel_points_2 - repmat(avg_pixel_2, num_rows, 1)) ./ std_pixel_2;


A_2 = [];
for i = 1:size(world_points_n_2, 1)
    x_world = world_points_n_2(i, 1);
    y_world = world_points_n_2(i, 2);
    z_world = world_points_n_2(i, 3);

    x_pixel = pixel_points_n_2(i, 1);
    y_pixel = pixel_points_n_2(i, 2);   

    element_1 = [
        x_world, y_world, z_world, 1, ...
        0, 0, 0, 0, ...
        -x_pixel * x_world, -x_pixel * y_world, -x_pixel * z_world, -x_pixel
        ];

    element_2 = [
        0, 0, 0, 0, ...
        x_world, y_world, z_world, 1, ...
        -y_pixel * x_world, -y_pixel * y_world, -y_pixel * z_world, -y_pixel
        ];

    A_2 = [A_2; element_1; element_2];
end

[U_2, S_2, V_2] = svd(A_2);
[~, index_of_min_singular_value] = min(diag(S_2));
P_n_2 = V_2(:, index_of_min_singular_value);
P_n_2 = reshape(P_n_2, 3, 4);

M1_1 = [
    1/std_pixel, 0, 0;
    0, 1/std_pixel, 0;
    0, 0, 1;
    ];
M1_2 = [
    1, 0, -avg_pixel(1);
    0, 1, -avg_pixel(2);
    0, 0, 1;
    ];

M1 = M1_1 * M1_2;

M2_1 = [
    1/std_world, 0, 0, 0;
    0, 1/std_world, 0, 0;
    0, 0, 1/std_world, 0;
    0, 0, 0, 1;
    ];
M2_2 = [
    1, 0, 0, -avg_world(1);
    0, 1, 0, -avg_world(2);
    0, 0, 1, -avg_world(3);
    0, 0, 0, 1;
    ];

M2 = M2_1 * M2_2;
P_2 = inv(M1) * P_n_2 * M2;
% disp(P);
% disp(P(3, 1:3));
% disp(norm(P(3, 1:3)));
P_2 = P_2 / norm(P_2(3, 1:3));

P_2_ = P_2;
P_2 = P_2(:, 1:3);

theta = atan(P_2(3,1)/P_2(3,2));
R_z_theta_2 = [
    cos(theta), sin(theta), 0;
    -sin(theta), cos(theta), 0;
    0, 0, 1;
];

R_A_2 = P_2 * R_z_theta_2;
Beta = atan(R_A_2(3,2)/R_A_2(3,3));
R_x_beta_2 = [
    1, 0, 0;
    0, cos(Beta), sin(Beta);
    0, -sin(Beta), cos(Beta);
];

R_B_2 = R_A_2 * R_x_beta_2;
gamma = atan(R_B_2(2,1)/R_B_2(2,2));
R_z_gamma_2 = [
    cos(gamma), sin(gamma), 0;
    -sin(gamma), cos(gamma), 0;
    0, 0, 1;
];

R_2 = (R_z_theta_2 * R_x_beta_2 * R_z_gamma_2)';
% disp(R_2);
PQ_2 = P_2 * R_z_theta_2 * R_x_beta_2 * R_z_gamma_2;
% disp(Q);

K = [
    3063.73307388903, 0, 2023.71937454517;
    0, 3050.01366910237, 1499.05930779402;
    0, 0, 1;
];

s = (inv(P_2) * P_2_);
T_2 = - s(:, 4);
% disp(T_2);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image 3

image_3 = imread('images/IMG_2218.jpg');

% imshow(image_3);
% dcm = datacursormode(gcf);
% set(dcm, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on', 'Enable', 'on');
% pause;
% cursor_info = getCursorInfo(dcm);
% for i = 1:length(cursor_info)
%     x = cursor_info(i).Position(1);
%     y = cursor_info(i).Position(2);
%     pixel_value = impixel(image_3, x, y);
%     disp(['Point ' num2str(i) ': (' num2str(x) ', ' num2str(y) ') - Pixel value: ' num2str(pixel_value)]);
% end

world_points_3 = [
    0, 6, 6;
    0, 6, 12;
    0, 12, 6;
    0, 12, 12;

    6, 6, 0;
    12, 12, 0;
    6, 12, 0;
    12, 12, 0;

    6, 0, 6;
    6, 0, 12;
    12, 0, 6;
    12, 0, 12;
];

[num_rows, ~] = size(world_points_3);
avg_world_3 = sum(world_points_3) / num_rows;
std_world_3 = sqrt(sum(sum((world_points_3 - repmat(avg_world_3, num_rows, 1)).^2)) / (3*num_rows));

world_points_n_3 = (world_points_3 - repmat(avg_world_3, num_rows, 1)) ./ std_world_3;


pixel_points_3 = [
    2380, 1476;
    2284, 1239;
    2372, 1770;
    2272, 1530;

    2880, 1649;
    3268, 1587;
    2848, 1934;
    3212, 1868;

    2772, 1102;
    2644, 890;
    3148, 1038;
    3016, 830;
];

[num_rows, ~] = size(pixel_points_3);
avg_pixel_3 = sum(pixel_points_3) / num_rows;
std_pixel_3 = sqrt(sum(sum((pixel_points_3 - repmat(avg_pixel_3, num_rows, 1)).^2)) / (2*num_rows));

pixel_points_n_3 = (pixel_points_3 - repmat(avg_pixel_3, num_rows, 1)) ./ std_pixel_3;

A_3 = [];
for i = 1:size(world_points_n_3, 1)
    x_world = world_points_n_3(i, 1);
    y_world = world_points_n_3(i, 2);
    z_world = world_points_n_3(i, 3);

    x_pixel = pixel_points_n_3(i, 1);
    y_pixel = pixel_points_n_3(i, 2);   

    element_1 = [
        x_world, y_world, z_world, 1, ...
        0, 0, 0, 0, ...
        -x_pixel * x_world, -x_pixel * y_world, -x_pixel * z_world, -x_pixel
        ];

    element_2 = [
        0, 0, 0, 0, ...
        x_world, y_world, z_world, 1, ...
        -y_pixel * x_world, -y_pixel * y_world, -y_pixel * z_world, -y_pixel
        ];

    A_3 = [A_3; element_1; element_2];
end

[U_3, S_3, V_3] = svd(A_3);
[~, index_of_min_singular_value] = min(diag(S_3));
P_n_3 = V_3(:, index_of_min_singular_value);
P_n_3 = reshape(P_n_3, 3, 4);

M1_1 = [
    1/std_pixel, 0, 0;
    0, 1/std_pixel, 0;
    0, 0, 1;
    ];
M1_2 = [
    1, 0, -avg_pixel(1);
    0, 1, -avg_pixel(2);
    0, 0, 1;
    ];

M1 = M1_1 * M1_2;

M2_1 = [
    1/std_world, 0, 0, 0;
    0, 1/std_world, 0, 0;
    0, 0, 1/std_world, 0;
    0, 0, 0, 1;
    ];
M2_2 = [
    1, 0, 0, -avg_world(1);
    0, 1, 0, -avg_world(2);
    0, 0, 1, -avg_world(3);
    0, 0, 0, 1;
    ];

M2 = M2_1 * M2_2;

P_3 = inv(M1) * P_n_3 * M2;
P_3 = P_3 / norm(P_3(3, 1:3));

P_3_ = P_3;
P_3 = P_3(:, 1:3);

theta = atan(P_3(3,1)/P_3(3,2));
R_z_theta_3 = [
    cos(theta), sin(theta), 0;
    -sin(theta), cos(theta), 0;
    0, 0, 1;
];

R_A_3 = P_3 * R_z_theta_3;
Beta = atan(R_A_3(3,2)/R_A_3(3,3));
R_x_beta_3 = [
    1, 0, 0;
    0, cos(Beta), sin(Beta);
    0, -sin(Beta), cos(Beta);
];

R_B_3 = R_A_3 * R_x_beta_3;
gamma = atan(R_B_3(2,1)/R_B_3(2,2));
R_z_gamma_3 = [
    cos(gamma), sin(gamma), 0;
    -sin(gamma), cos(gamma), 0;
    0, 0, 1;
];

R_3 = (R_z_theta_3 * R_x_beta_3 * R_z_gamma_3)';
% disp(R_3);

PQ_3 = P_3 * R_z_theta_3 * R_x_beta_3 * R_z_gamma_3;
% disp(Q);

K = [
    3063.73307388903, 0, 2023.71937454517;
    0, 3050.01366910237, 1499.05930779402;
    0, 0, 1;
];

s = (inv(P_3) * P_3_);
T_3 = - s(:, 4);
% disp(T_3);


%define origin of world coordinates, and 3 axis
origin = [0, 0, 0];
x_axis = [1, 0, 0];
y_axis = [0, 1, 0];
z_axis = [0, 0, 1];

figure;
quiver3(origin(1), origin(2), origin(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 2);
hold on;

quiver3(origin(1), origin(2), origin(3), y_axis(1), y_axis(2), y_axis(3), 'g', 'LineWidth', 2);

quiver3(origin(1), origin(2), origin(3), z_axis(1), z_axis(2), z_axis(3), 'b', 'LineWidth', 2);

% First image camera coordinates origin is the negative of translation vector, and each row of
%rotation matrix is the axis of camera coordinates.
origin_camera = -T;
x_axis = R(1, :);
y_axis = R(2, :);
z_axis = R(3, :);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), y_axis(1), y_axis(2), y_axis(3), 'g', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), z_axis(1), z_axis(2), z_axis(3), 'b', 'LineWidth', 3);

% Second image camera coordinates origin is the negative of translation vector, and each row of
%rotation matrix is the axis of camera coordinates.
origin_camera = -T_2;
x_axis = R_2(1, :);
y_axis = R_2(2, :);
z_axis = R_2(3, :);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), y_axis(1), y_axis(2), y_axis(3), 'g', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), z_axis(1), z_axis(2), z_axis(3), 'b', 'LineWidth', 3);

% Third image camera coordinates origin is the negative of translation vector, and each row of
%rotation matrix is the axis of camera coordinates.
origin_camera = -T_3;
x_axis = R_3(1, :);
y_axis = R_3(2, :);
z_axis = R_3(3, :);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), x_axis(1), x_axis(2), x_axis(3), 'r', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), y_axis(1), y_axis(2), y_axis(3), 'g', 'LineWidth', 3);
quiver3(origin_camera(1), origin_camera(2), origin_camera(3), z_axis(1), z_axis(2), z_axis(3), 'b', 'LineWidth', 3);

axis equal;
grid on;

legend('X-axis', 'Y-axis', 'Z-axis');

hold off;
