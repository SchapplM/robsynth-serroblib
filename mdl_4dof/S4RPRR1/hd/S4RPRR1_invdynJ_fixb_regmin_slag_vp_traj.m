% Calculate time series of minimal parameter regressor of inv. dyn. joint torques for
% S4RPRR1
%
% Input:
% Q [NTx4]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx4]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx4]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
%
% Output:
% RV_Traj [NTx20]
%   time series of regressor matrices as vectors
%   see S4RPRR1_invdynJ_fixb_regmin2vec.m

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV_Traj = S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,4]),
%$cgargs  coder.newtype('double',[inf,4]),
%$cgargs  coder.newtype('double',[inf,4]),
%$cgargs  zeros(3,1), zeros(7,1)}
assert(isreal(Q) && all(size(Q,2) == 4), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj: Q needs to be [NTx4] (double)');
assert(isreal(QD) && all(size(QD,2) == 4), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj: QD needs to be [NTx4] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 4), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj: QDD needs to be [NTx4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp_traj: Kinematic parameters pkin have to be [7x1] (double)');
  
%% Trajektorie der Regressor-Vektoren aufbauen
RV_Traj = NaN(size(Q,1), 20);
for ii = 1:size(Q,1)
  RV_Traj(ii,:) = S4RPRR1_invdynJ_fixb_regmin2vec( ...
    S4RPRR1_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end
