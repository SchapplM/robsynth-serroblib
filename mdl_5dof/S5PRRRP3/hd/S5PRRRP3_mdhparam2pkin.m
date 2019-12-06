% Convert vector of modified DH parameters to kinematic parameter vector for
% S5PRRRP3
%
% Input:
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
%
% Output:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = S5PRRRP3_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); a_mdh(3); a_mdh(4); a_mdh(5); d_mdh(2); d_mdh(3); d_mdh(4); theta_mdh(1);];
pkin = t1;
