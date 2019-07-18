% Convert vector of kinematic parameters to modified DH parameters of
% S4PRRR2
%
% Input:
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [4x1]
%   Rotation around z
% b_mdh [4x1]
%   Translation along z
% alpha_mdh [4x1]
%   Rotation around x
% a_mdh [4x1]
%   Translation along x
% theta_mdh [4x1]
%   Rotation around z
% d_mdh [4x1]
%   Translation along z
% qoffset_mdh [4x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S4PRRR2_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(2,1)}
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_pkin2mdhparam: Kinematic parameters pkin have to be [2x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [pi / 0.2e1; 0; 0; 0;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [0; 0; pkin(1); pkin(2);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [pi / 0.2e1; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; 0; 0;];
qoffset_mdh = t1;
