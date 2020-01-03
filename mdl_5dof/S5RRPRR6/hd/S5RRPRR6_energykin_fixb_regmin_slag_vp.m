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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:05:53
% EndTime: 2020-01-03 12:05:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (175->32), mult. (261->77), div. (0->0), fcn. (146->8), ass. (0->36)
t137 = pkin(1) * qJD(1);
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t114 = qJD(1) + qJD(2);
t129 = sin(qJ(2)) * t137;
t106 = t114 * qJ(3) + t129;
t116 = cos(pkin(9));
t134 = t106 * t116;
t115 = sin(pkin(9));
t128 = cos(qJ(2)) * t137;
t123 = qJD(3) - t128;
t99 = (-pkin(3) * t116 - pkin(7) * t115 - pkin(2)) * t114 + t123;
t136 = t118 * t99 + t121 * t134;
t135 = t106 * t114;
t111 = t114 ^ 2;
t112 = t115 ^ 2;
t133 = t111 * t112;
t132 = t114 * t115;
t131 = t114 * t118;
t130 = t116 * t114;
t127 = t112 * t135;
t126 = t115 * t131;
t125 = t121 * t132;
t109 = -qJD(4) + t130;
t124 = -t118 * t134 + t121 * t99;
t120 = cos(qJ(5));
t117 = sin(qJ(5));
t113 = t116 ^ 2;
t107 = -qJD(5) + t109;
t104 = -t114 * pkin(2) + t123;
t102 = (-t117 * t118 + t120 * t121) * t132;
t101 = (t117 * t121 + t118 * t120) * t132;
t100 = (pkin(4) * t131 + t106) * t115;
t96 = -pkin(8) * t126 + t136;
t95 = -t109 * pkin(4) - pkin(8) * t125 + t124;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t111 / 0.2e1, t114 * t128, -t114 * t129, -t104 * t130, t104 * t132, (t112 + t113) * t135, t104 ^ 2 / 0.2e1 + (t113 / 0.2e1 + t112 / 0.2e1) * t106 ^ 2, t121 ^ 2 * t133 / 0.2e1, -t121 * t118 * t133, -t109 * t125, t109 * t126, t109 ^ 2 / 0.2e1, -t124 * t109 + t118 * t127, t136 * t109 + t121 * t127, t102 ^ 2 / 0.2e1, -t102 * t101, -t102 * t107, t101 * t107, t107 ^ 2 / 0.2e1, -(-t117 * t96 + t120 * t95) * t107 + t100 * t101, (t117 * t95 + t120 * t96) * t107 + t100 * t102;];
T_reg = t1;
