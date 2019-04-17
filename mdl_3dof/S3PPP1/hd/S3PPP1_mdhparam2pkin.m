% Convert vector of modified DH parameters to kinematic parameter vector for
% S3PPP1
%
% Input:
% beta_mdh [3x1]
%   Rotation around z
% b_mdh [3x1]
%   Translation along z
% alpha_mdh [3x1]
%   Rotation around x
% a_mdh [3x1]
%   Translation along x
% theta_mdh [3x1]
%   Rotation around z
% d_mdh [3x1]
%   Translation along z
% qoffset_mdh [3x1]
%   Offset on joint coordinate q
%
% Output:
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,theta1]';

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = S3PPP1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); a_mdh(3); theta_mdh(1);];
pkin  = t1;
