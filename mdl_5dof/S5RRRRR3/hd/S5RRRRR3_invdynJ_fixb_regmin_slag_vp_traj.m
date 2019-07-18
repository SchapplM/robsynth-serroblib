% Calculate time series of minimal parameter regressor of inv. dyn. joint torques for
% S5RRRRR3
%
% Input:
% Q [NTx5]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx5]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx5]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% RV_Traj [NTx99]
%   time series of regressor matrices as vectors
%   see S5RRRRR3_invdynJ_fixb_regmin2vec.m

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV_Traj = S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,5]),
%$cgargs  coder.newtype('double',[inf,5]),
%$cgargs  coder.newtype('double',[inf,5]),
%$cgargs  zeros(3,1), zeros(5,1)}
assert(isreal(Q) && all(size(Q,2) == 5), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj: Q needs to be [NTx5] (double)');
assert(isreal(QD) && all(size(QD,2) == 5), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj: QD needs to be [NTx5] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 5), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj: QDD needs to be [NTx5] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_regmin_slag_vp_traj: Kinematic parameters pkin have to be [5x1] (double)');
  
%% Trajektorie der Regressor-Vektoren aufbauen
RV_Traj = NaN(size(Q,1), 99);
for ii = 1:size(Q,1)
  RV_Traj(ii,:) = S5RRRRR3_invdynJ_fixb_regmin2vec( ...
    S5RRRRR3_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end
