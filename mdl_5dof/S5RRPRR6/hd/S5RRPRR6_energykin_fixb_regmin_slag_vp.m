% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:10
% EndTime: 2019-12-05 18:36:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (175->32), mult. (261->77), div. (0->0), fcn. (146->8), ass. (0->36)
t139 = pkin(1) * qJD(1);
t116 = qJD(1) + qJD(2);
t117 = sin(pkin(9));
t118 = cos(pkin(9));
t130 = cos(qJ(2)) * t139;
t126 = qJD(3) - t130;
t101 = (-pkin(3) * t118 - pkin(7) * t117 - pkin(2)) * t116 + t126;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t131 = sin(qJ(2)) * t139;
t108 = t116 * qJ(3) + t131;
t136 = t108 * t118;
t138 = t120 * t101 + t123 * t136;
t137 = t108 * t116;
t113 = t116 ^ 2;
t114 = t117 ^ 2;
t135 = t113 * t114;
t134 = t116 * t117;
t133 = t116 * t120;
t132 = t118 * t116;
t129 = t114 * t137;
t128 = t117 * t133;
t127 = t123 * t134;
t111 = -qJD(4) + t132;
t125 = t123 * t101 - t120 * t136;
t122 = cos(qJ(5));
t119 = sin(qJ(5));
t115 = t118 ^ 2;
t109 = -qJD(5) + t111;
t106 = -t116 * pkin(2) + t126;
t104 = (-t119 * t120 + t122 * t123) * t134;
t103 = (t119 * t123 + t120 * t122) * t134;
t102 = (pkin(4) * t133 + t108) * t117;
t98 = -pkin(8) * t128 + t138;
t97 = -t111 * pkin(4) - pkin(8) * t127 + t125;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t113 / 0.2e1, t116 * t130, -t116 * t131, -t106 * t132, t106 * t134, (t114 + t115) * t137, t106 ^ 2 / 0.2e1 + (t115 / 0.2e1 + t114 / 0.2e1) * t108 ^ 2, t123 ^ 2 * t135 / 0.2e1, -t123 * t120 * t135, -t111 * t127, t111 * t128, t111 ^ 2 / 0.2e1, -t125 * t111 + t120 * t129, t138 * t111 + t123 * t129, t104 ^ 2 / 0.2e1, -t104 * t103, -t104 * t109, t103 * t109, t109 ^ 2 / 0.2e1, -(-t119 * t98 + t122 * t97) * t109 + t102 * t103, (t119 * t97 + t122 * t98) * t109 + t102 * t104;];
T_reg = t1;
