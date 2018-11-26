% Convert vector of modified DH parameters to kinematic parameter vector for
% S2RR2
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
%   pkin=[d2]';

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function pkin = S2RR2_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [d_mdh(2);];
pkin  = t1;
