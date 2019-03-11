% Calculate time series of minimal parameter regressor of inv. dyn. joint torques for
% S3PRR1
%
% Input:
% Q [NTx3]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx3]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx3]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
%
% Output:
% RV_Traj [NTx14]
%   time series of regressor matrices as vectors
%   see S3PRR1_invdynJ_fixb_regmin2vec.m

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV_Traj = S3PRR1_invdynJ_fixb_regmin_slag_vp_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,3]),
%$cgargs  coder.newtype('double',[inf,3]),
%$cgargs  coder.newtype('double',[inf,3]),
%$cgargs  zeros(3,1), zeros(4,1)}
assert(isreal(Q) && all(size(Q,2) == 3), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp_traj: Q needs to be [NTx3] (double)');
assert(isreal(QD) && all(size(QD,2) == 3), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp_traj: QD needs to be [NTx3] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 3), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp_traj: QDD needs to be [NTx3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp_traj: Kinematic parameters pkin have to be [4x1] (double)');
  
%% Trajektorie der Regressor-Vektoren aufbauen
RV_Traj = NaN(size(Q,1), 14);
for ii = 1:size(Q,1)
  RV_Traj(ii,:) = S3PRR1_invdynJ_fixb_regmin2vec( ...
    S3PRR1_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end
