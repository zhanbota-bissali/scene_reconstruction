%% SECTION 1
close all
clear
clc
load('mats/coord_1080.mat');
load('mats/h_ref_1080.mat');
load('mats/h_1080.mat');
% 
% load('mats/coord720.mat');
% load('mats/h_ref720.mat');
% load('mats/h720.mat');

A = internal_parameters_solve_vmmc(h);
[R T] = external_parameters_solve_vmmc(A, h);
% refined

A_ref = internal_parameters_solve_vmmc(H);
[R_ref T_ref] = external_parameters_solve_vmmc(A_ref, H);




%% SECTION 2

%%define dataset path
params.Directory    = fullfile('images/images_6/all');

%%Detector parameters
params.detector     =  'SIFT'; %'SIFT', 'SURF', 'DoH', 'KAZE'
params.nscales      =        10;
params.noctaves     =        3;
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      300;
params.th           =    0.001; % alternative to npoints.

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'DSP-SIFT'; % 'SIFT', 'SURF', 'DSP-SIFT'
params.desOnDecom   =    false; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false; % set to true to avoid orientation estimation.
% for DSP-SIFT
params.dsp.ns       =      15;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       7;% largest scale (relative to detection);    

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.8;
params.Metric       =  'SSD';
%% END OF PARAMETER DEFINITION %%

%% addpaths
addpath(genpath('./detectors/'));
addpath(genpath('./descriptors/'));
addpath(genpath('./toolbox/'));

%% preload dataset
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);

%% initialize (sort of)
ima{numImages}           = [];
points{numImages}        = [];
decomposition{numImages} = [];
features{numImages}      = [];

%% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end

%% detect & describe
for j = 1:numImages
%% Load and convert images %%
% scale = 0.8;
% ima{j}      =       imresize(readimage(params.Scene, j),scale);
ima{j}      =       readimage(params.Scene, j);
gima        =      im2double(rgb2gray(ima{j}));

%% PoI Detection %%
sprintf('Detecting for image: %d',j)
[points{j},decomposition{j}] =  myDetector(gima,params);

%% PoI Description %%
sprintf('Describing for image: %d',j)
[features{j},points{j}]      =  myDescriptor(points{j},decomposition{j},params);

%% show detections
figure(j)
imshow(ima{j}); hold on;
plot(points{j},'showOrientation',true);

end

%% SECTION 3

% Computation of the essential matrix and euclidean reconstruction

% include ACT_lite path
ACT_path = 'ACT_lite';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = 'extra_funs';
addpath(genpath(extra_funs_path));

vgg_path = 'vgg';
addpath(genpath(vgg_path));

warning off
disp('************************************* START')

MaxRatio = 0.8;
Metric  = params.Metric;

q_data = n_view_matching(points, features, ima, MaxRatio, Metric);
q_data = homogenize_coords(q_data);

n_points = size(q_data,2);
ncam = numImages;


% 2. Fundamental matrix using 1 and the last cameras
% ------------------------------------------------------------------------

q_2cams(:,:,1)=q_data(:,:,1); 
q_2cams(:,:,2)=q_data(:,:,ncam);

[F,P_2cams,Q_est,q_est] = MatFunProjectiveCalib(q_2cams); 

disp(['Mean re-projection error of initial reconstruction (2 cams)   = ' num2str( ErrorRetroproy(q_2cams,P_2cams,Q_est)/2 )]);
draw_reproj_error(q_2cams,P_2cams,Q_est);

%% 3.Projective Bundle Adjustment, projection matrices of of the cameras
% ------------------------------------------------------------------------
P = zeros(3,4,ncam);
P(:,:,1) = P_2cams(:,:,1);
P(:,:,ncam) = P_2cams(:,:,2);

for i = 1:ncam
    P(:,:,i) =  PDLT_NA(q_data(:,:,i), Q_est, 0,0);
end

disp(['Residual reprojection error resectioning   = ' num2str( ErrorRetroproy(q_data,P,Q_est)/2 )]);
draw_reproj_error(q_data,P,Q_est);

vp = ones(n_points,ncam);

[P_adj,Q_adj, q_adj] = BAProjectiveCalib(q_data,P,Q_est,vp);
disp(['Residual reprojection error Bundle Adjustment   = ' num2str( ErrorRetroproy(q_data,P_adj,Q_adj)/2 )]);
draw_reproj_error(q_data,P_adj,Q_adj);
q_adj2(:,:,1) = q_adj(:,:,1);
q_adj2(:,:,2) = q_adj(:,:,2);

%% 4. Re-computed Fundamental matrix between two of the cameras
% ------------------------------------------------------------------------

F_recomp = vgg_F_from_P(P(:,:,1), P(:,:,2));




%% 5. Obtain the essential matrix (E) from the fundamental matrix (F) and the
% intrinsic parameter matrices (K).
% ------------------------------------------------------------------------

K=zeros(3,3,2);
K(:,:,1)=A;
K(:,:,2)=A;
E =  normalize_matrix((K(:,:,1).')*F_recomp*K(:,:,2));

% ------------------------------------------------------------------------
% 4. Factorize the essential matrix with the 2 possible solutions for R. 
% Use the function factorize_E to obtain R_est(:,:,1) and R_est(:,:,2) and T_est.
% ------------------------------------------------------------------------
[R_est,T_est] = factorize_E(E);

% ------------------------------------------------------------------------
% Save the 4 solutions (R,t) in the structures Rcam(3,3,cam,sol), T(3,cam,sol),
% where cam indicates the camera number and sol indicates the solution number (1, 2, 3 or 4).
% ------------------------------------------------------------------------
Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);

for a=1:4
Rcam(:,:,1,a)=eye(3,3);
end
for a=1:4
    if a<3
        Rcam(:,:,2,a)=R_est(:,:,1);
    else
        Rcam(:,:,2,a)=R_est(:,:,2);
    end
end
Tcam(:,2,1)=T_est;
Tcam(:,2,2)=-T_est;
Tcam(:,2,3)=T_est;
Tcam(:,2,4)=-T_est;
% ------------------------------------------------------------------------
% 5. For each solution we obtain an Euclidean solution and we visualize it.
% ------------------------------------------------------------------------
npoints = size(q_adj2,2);
Q_euc = zeros(4,npoints,2); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;

for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc =  TriangEuc(Rcam(:,:,2,sol), Tcam(:,2,sol),K,q_adj2(:,:,:));
       
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc, K, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
    % Compute the projection matrices from K, Rcam, Tcam
    for k=1:2
        P_euc(:,:,k) = K(:,:,k)*[Rcam(:,:,k,sol) -Rcam(:,:,k,sol)*Tcam(:,k,sol)];
    end
    
    % Obtain the re-projected points q_rep
    for n=1:2
        q_rep(:,:,n) = P_euc(:,:,n)*Q_euc;
    end
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
      daspect([1, 1, 1]);
      pbaspect([1, 1, 1]);
      axis([-1000, 1000, -1000, 1000]);
    end
end

disp('************************************* END')




