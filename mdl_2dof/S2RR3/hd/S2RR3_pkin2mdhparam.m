% Convert vector of kinematic parameters to modified DH parameters of
% S2RR3
%
% Input:
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [2x1]
%   Rotation around z
% b_mdh [2x1]
%   Translation along z
% alpha_mdh [2x1]
%   Rotation around x
% a_mdh [2x1]
%   Translation along x
% theta_mdh [2x1]
%   Rotation around z
% d_mdh [2x1]
%   Translation along z
% qoffset_mdh [2x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S2RR3_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(3,1)}
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_pkin2mdhparam: Kinematic parameters pkin have to be [3x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [0; 0;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(1);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(2); pkin(3);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0;];
qoffset_mdh = t1;
