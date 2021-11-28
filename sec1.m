close all
clear
clc
load('mats/coord_1080.mat');
load('mats/h_ref_1080.mat');
load('mats/h_1080.mat');

% load('mats/coord720.mat');
% load('mats/h_ref720.mat');
% load('mats/h720.mat');


% for a=1:6
%     img{a} = imread(sprintf('images/1080/img%d.jpeg',a));
%     img{a} = imread(sprintf('images/720/img%d.jpeg',a));
% end
% [coords img_p]= get_real_points_checkerboard_vmmc(9, 160, 1);
% [coords img_p]= get_real_points_checkerboard_vmmc(9, 110, 1);
% coords = transpose(coords);
% for i=1:6
%     xy{i} = get_user_points_vmmc(img{i});
%     h{i} = homography_solve_vmmc(coords, xy{i}); 
%     tform{i} = transpose(maketform('projective', (h{i})'));
%     t_img{i} = imtransform(img_p,tform{i});
%     [H{i}, rperr{i}] = homography_refine_vmmc(coords, xy{i}, h{i});
%     figure;
%     subplot(121);imshow(img{i});
%     subplot(122);imshow(t_img{i});
% end
% 
% save('mats/coord_1080.mat','xy');
% save('mats/h_ref_1080.mat', 'H');
% save('mats/h_1080.mat', 'h');
% save('mats/coord720.mat','xy');
% save('mats/h_ref720.mat', 'H');
% save('mats/h720.mat', 'h');


A6 = internal_parameters_solve_vmmc(h);
[R6 T6] = external_parameters_solve_vmmc(A6, h);

% refined

A_ref6 = internal_parameters_solve_vmmc(H);
[R_ref6 T_ref6] = external_parameters_solve_vmmc(A_ref6, H);
