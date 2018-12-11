% Convert vector of kinematic parameters to modified DH parameters of
% S6RRPRRR14
%
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [6x1]
%   Rotation around z
% b_mdh [6x1]
%   Translation along z
% alpha_mdh [6x1]
%   Rotation around x
% a_mdh [6x1]
%   Translation along x
% theta_mdh [6x1]
%   Rotation around z
% d_mdh [6x1]
%   Translation along z
% qoffset_mdh [6x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = S6RRPRRR14_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(14,1)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_pkin2mdhparam: Kinematic parameters pkin have to be [14x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; pkin(6); pkin(7); pkin(8); t1; t1;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(1); pkin(2); pkin(3); pkin(4); pkin(5);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; pkin(14); 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(9); pkin(10); 0; pkin(11); pkin(12); pkin(13);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
qoffset_mdh = t1;
