% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:24
% EndTime: 2019-12-05 17:36:24
% DurationCPUTime: 0.12s
% Computational Cost: add. (90->30), mult. (221->70), div. (0->0), fcn. (112->6), ass. (0->27)
t108 = sin(qJ(4));
t109 = cos(qJ(4));
t104 = sin(pkin(8));
t106 = cos(pkin(8));
t107 = cos(pkin(7));
t115 = -pkin(1) * t107 - pkin(2);
t93 = qJD(3) + (-pkin(3) * t106 - pkin(6) * t104 + t115) * qJD(1);
t105 = sin(pkin(7));
t99 = (pkin(1) * t105 + qJ(3)) * qJD(1);
t97 = t104 * qJD(2) + t106 * t99;
t120 = t108 * t93 + t109 * t97;
t110 = qJD(1) ^ 2;
t119 = t104 ^ 2 * t110;
t118 = qJD(1) * t104;
t117 = qJD(1) * t108;
t116 = t106 * qJD(1);
t114 = t104 * t117;
t113 = t109 * t118;
t112 = -t108 * t97 + t109 * t93;
t102 = t106 * qJD(2);
t100 = -qJD(4) + t116;
t98 = t115 * qJD(1) + qJD(3);
t95 = t104 * t99 - t102;
t90 = qJD(5) - t102 + (pkin(4) * t117 + t99) * t104;
t89 = -qJ(5) * t114 + t120;
t88 = -t100 * pkin(4) - qJ(5) * t113 + t112;
t1 = [t110 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t105 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t110, -t98 * t116, t98 * t118, (t104 * t95 + t106 * t97) * qJD(1), t97 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1, t109 ^ 2 * t119 / 0.2e1, -t109 * t108 * t119, -t100 * t113, t100 * t114, t100 ^ 2 / 0.2e1, -t112 * t100 + t95 * t114, t120 * t100 + t95 * t113, (-t108 * t89 - t109 * t88) * t118, t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t90 ^ 2 / 0.2e1;];
T_reg = t1;
