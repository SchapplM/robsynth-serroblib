% Convert vector of modified DH parameters to kinematic parameter vector for
% S3PRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function pkin = S3PRR1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); a_mdh(3); d_mdh(2); d_mdh(3);];
pkin  = t1;
