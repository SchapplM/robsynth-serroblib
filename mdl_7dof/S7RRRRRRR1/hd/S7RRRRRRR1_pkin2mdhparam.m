% Convert vector of kinematic parameters to modified DH parameters of
% S7RRRRRRR1
%
% Input:
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [7x1]
%   Rotation around z
% b_mdh [7x1]
%   Translation along z
% alpha_mdh [7x1]
%   Rotation around x
% a_mdh [7x1]
%   Translation along x
% theta_mdh [7x1]
%   Rotation around z
% d_mdh [7x1]
%   Translation along z
% qoffset_mdh [7x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S7RRRRRRR1_pkin2mdhparam(pkin)


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t2 = -pi / 0.2e1;
t1 = pi / 0.2e1;
t3 = [0; t1; t2; t2; t1; t1; t2;];
alpha_mdh = t3;

% Aus parameters_mdh_a_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0;];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(1); 0; pkin(2); 0; pkin(3); 0; pkin(4);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0;];
qoffset_mdh = t1;
