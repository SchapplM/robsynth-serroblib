% Convert vector of modified DH parameters to kinematic parameter vector for
% S7RRRRRRR1
%
% Input:
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
%
% Output:
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function pkin = S7RRRRRRR1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [d_mdh(1); d_mdh(3); d_mdh(5); d_mdh(7);];
pkin  = t1;
