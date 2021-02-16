% Convert vector of modified DH parameters to kinematic parameter vector for
% S1P1
%
% Input:
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
%
% Output:
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = S1P1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [];
pkin = t1;
pkin = NaN;% Dummy-Wert, da pkin nicht leer sein soll
