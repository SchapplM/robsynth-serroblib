% Convert vector of kinematic parameters to modified DH parameters of
% S2PP1
%
% Input:
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2]';
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
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S2PP1_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(1,1)}
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2PP1_pkin2mdhparam: Kinematic parameters pkin have to be [1x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [-pi / 0.2e1; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [-pi / 0.2e1; pi / 0.2e1;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(1);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [-pi / 0.2e1; pi / 0.2e1;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0;];
qoffset_mdh = t1;
