% Convert vector of kinematic parameters to modified DH parameters of
% S1R1
%
% Input:
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [1x1]
%   Rotation around z
% b_mdh [1x1]
%   Translation along z
% alpha_mdh [1x1]
%   Rotation around x
% a_mdh [1x1]
%   Translation along x
% theta_mdh [1x1]
%   Rotation around z
% d_mdh [1x1]
%   Translation along z
% qoffset_mdh [1x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S1R1_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(1,1)}
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_pkin2mdhparam: Kinematic parameters pkin have to be [1x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [0;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [0;];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(1);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0;];
qoffset_mdh = t1;
