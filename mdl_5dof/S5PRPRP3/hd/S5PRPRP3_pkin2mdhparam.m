% Convert vector of kinematic parameters to modified DH parameters of
% S5PRPRP3
%
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [5x1]
%   Rotation around z
% b_mdh [5x1]
%   Translation along z
% alpha_mdh [5x1]
%   Rotation around x
% a_mdh [5x1]
%   Translation along x
% theta_mdh [5x1]
%   Rotation around z
% d_mdh [5x1]
%   Translation along z
% qoffset_mdh [5x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S5PRPRP3_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(8,1)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_pkin2mdhparam: Kinematic parameters pkin have to be [8x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; t1; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(1); pkin(2); pkin(3); pkin(4);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [pkin(7); 0; pkin(8); 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [0; pkin(5); 0; pkin(6); 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; 0; 0; 0;];
qoffset_mdh = t1;
