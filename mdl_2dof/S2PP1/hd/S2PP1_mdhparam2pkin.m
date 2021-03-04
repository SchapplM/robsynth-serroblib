% Convert vector of modified DH parameters to kinematic parameter vector for
% S2PP1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2]';

% Quelle: HybrDyn-Toolbox
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = S2PP1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2);];
pkin = t1;
