% Convert vector of modified DH parameters to kinematic parameter vector for
% S2RR3
%
% Input:
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
%
% Output:
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = S2RR3_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); d_mdh(1); d_mdh(2);];
pkin = t1;
