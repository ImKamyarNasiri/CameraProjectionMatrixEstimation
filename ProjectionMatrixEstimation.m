image = imread('images/IMG_2216.jpg');

%define 12 world points
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
%normalize world coordinates
[num_rows, ~] = size(world_points);
avg_world = sum(world_points) / num_rows;
std_world = sqrt(sum(sum((world_points - repmat(avg_world, num_rows, 1)).^2)) / (3*num_rows));

world_points_n = (world_points - repmat(avg_world, num_rows, 1)) ./ std_world;

%%select points and find pixel coordinates correspond to that point
% imshow(image);
% dcm = datacursormode(gcf);
% set(dcm, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on', 'Enable', 'on');
% pause;
% cursor_info = getCursorInfo(dcm);
% for i = 1:length(cursor_info)
%     x = cursor_info(i).Position(1);
%     y = cursor_info(i).Position(2);
%     pixel_value = impixel(image, x, y);
%     disp(['Point ' num2str(i) ': (' num2str(x) ', ' num2str(y) ') - Pixel value: ' num2str(pixel_value)]);
% end

%define pixel cordinates 
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
%normalize pixel coordinates
[num_rows, ~] = size(pixel_points);
avg_pixel = sum(pixel_points) / num_rows;
std_pixel = sqrt(sum(sum((pixel_points - repmat(avg_pixel, num_rows, 1)).^2)) / (2*num_rows));

pixel_points_n = (pixel_points - repmat(avg_pixel, num_rows, 1)) ./ std_pixel;

A = [];
%%compute matrix A
for i = 1:size(world_points_n, 1)
    x_world = world_points_n(i, 1);
    y_world = world_points_n(i, 2);
    z_world = world_points_n(i, 3);
    
    x_pixel = pixel_points_n(i, 1);
    y_pixel = pixel_points_n(i, 2);   
    
    %first row for each sample
    element_1 = [
        x_world, y_world, z_world, 1, ...
        0, 0, 0, 0, ...
        -x_pixel * x_world, -x_pixel * y_world, -x_pixel * z_world, -x_pixel
        ];
    
    %second row for each sample
    element_2 = [
        0, 0, 0, 0, ...
        x_world, y_world, z_world, 1, ...
        -y_pixel * x_world, -y_pixel * y_world, -y_pixel * z_world, -y_pixel
        ];
     
    A = [A; element_1; element_2];
end
%compute SVD of A
[U, S, V] = svd(A);
%pick the eigenvector corresponding to the smallest singular value
[~, index_of_min_singular_value] = min(diag(S));
P_n = V(:, index_of_min_singular_value);
P_n = reshape(P_n, 3, 4);
%normalize w.r.t the last row of P
P_n = P_n / norm(P_n(:, 1:3));

%define normalization matrices of world and pixel coordinates
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

%%test phase
world_points_test = [
    0, 6, 8;
    8, 6, 0;
    3 , 6 , 0;
    0, 6, 3;
];

camera_points_test = [
    2720, 1332;
    3524, 1378;
    3271, 1579;
    3014, 1551;
];

%convert to homogenous coordinates
world_points_test_1 = [world_points_test'; ones(1, size(world_points_test', 2))];

estimated_pixels_test = P * world_points_test_1;
%obtain x, y by dividing by the third element
estimated_pixels_test = estimated_pixels_test(1:2, :) ./ estimated_pixels_test(3, :);
disp(estimated_pixels_test);

%error is abs of subtraction of true and predicted pixel coordinates
errors = abs(estimated_pixels_test - camera_points_test');
disp(errors);
mean_errors = mean(errors(:));
std_errors = std(errors(:));

disp(mean_errors);
disp(std_errors);
