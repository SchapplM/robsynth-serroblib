% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:39
% EndTime: 2019-12-31 19:01:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (131->35), mult. (313->83), div. (0->0), fcn. (195->8), ass. (0->33)
t126 = qJD(1) ^ 2;
t135 = t126 / 0.2e1;
t118 = sin(pkin(9));
t111 = (pkin(1) * t118 + pkin(6)) * qJD(1);
t125 = cos(qJ(3));
t116 = t125 * qJD(2);
t122 = sin(qJ(3));
t102 = qJD(3) * pkin(3) + t116 + (-pkin(7) * qJD(1) - t111) * t122;
t131 = qJD(1) * t125;
t133 = t122 * qJD(2) + t125 * t111;
t105 = pkin(7) * t131 + t133;
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t134 = t121 * t102 + t124 * t105;
t132 = qJD(1) * t122;
t130 = qJD(1) * qJD(3);
t119 = cos(pkin(9));
t129 = -pkin(1) * t119 - pkin(2);
t107 = t121 * t132 - t124 * t131;
t128 = t102 * t124 - t105 * t121;
t109 = (-pkin(3) * t125 + t129) * qJD(1);
t123 = cos(qJ(5));
t120 = sin(qJ(5));
t117 = qJD(3) + qJD(4);
t112 = t129 * qJD(1);
t108 = (t121 * t125 + t122 * t124) * qJD(1);
t106 = qJD(5) + t107;
t104 = t108 * t123 + t117 * t120;
t103 = t108 * t120 - t123 * t117;
t99 = pkin(4) * t107 - pkin(8) * t108 + t109;
t98 = pkin(8) * t117 + t134;
t97 = -pkin(4) * t117 - t128;
t1 = [t135, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t118 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t126, t122 ^ 2 * t135, t122 * t126 * t125, t122 * t130, t125 * t130, qJD(3) ^ 2 / 0.2e1, -t112 * t131 + (-t111 * t122 + t116) * qJD(3), -t133 * qJD(3) + t112 * t132, t108 ^ 2 / 0.2e1, -t108 * t107, t108 * t117, -t107 * t117, t117 ^ 2 / 0.2e1, t109 * t107 + t128 * t117, t109 * t108 - t134 * t117, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t106, -t103 * t106, t106 ^ 2 / 0.2e1, (-t120 * t98 + t123 * t99) * t106 + t97 * t103, -(t120 * t99 + t123 * t98) * t106 + t97 * t104;];
T_reg = t1;
