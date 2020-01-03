% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:07
% DurationCPUTime: 0.11s
% Computational Cost: add. (123->37), mult. (315->80), div. (0->0), fcn. (199->8), ass. (0->31)
t119 = sin(pkin(8));
t113 = (pkin(1) * t119 + qJ(3)) * qJD(1);
t120 = cos(pkin(9));
t116 = t120 * qJD(2);
t118 = sin(pkin(9));
t101 = t116 + (-pkin(6) * qJD(1) - t113) * t118;
t106 = t118 * qJD(2) + t120 * t113;
t130 = qJD(1) * t120;
t102 = pkin(6) * t130 + t106;
t123 = sin(qJ(4));
t125 = cos(qJ(4));
t132 = t123 * t101 + t125 * t102;
t131 = qJD(1) * t118;
t121 = cos(pkin(8));
t129 = -pkin(1) * t121 - pkin(2);
t109 = t123 * t131 - t125 * t130;
t128 = t125 * t101 - t123 * t102;
t108 = qJD(3) + (-pkin(3) * t120 + t129) * qJD(1);
t126 = qJD(1) ^ 2;
t124 = cos(qJ(5));
t122 = sin(qJ(5));
t112 = t129 * qJD(1) + qJD(3);
t110 = (t118 * t125 + t120 * t123) * qJD(1);
t107 = qJD(5) + t109;
t105 = -t118 * t113 + t116;
t104 = t122 * qJD(4) + t124 * t110;
t103 = -t124 * qJD(4) + t122 * t110;
t98 = t109 * pkin(4) - t110 * pkin(7) + t108;
t97 = qJD(4) * pkin(7) + t132;
t96 = -qJD(4) * pkin(4) - t128;
t1 = [t126 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t119 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t126, -t112 * t130, t112 * t131, (-t105 * t118 + t106 * t120) * qJD(1), t106 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1, t110 ^ 2 / 0.2e1, -t110 * t109, t110 * qJD(4), -t109 * qJD(4), qJD(4) ^ 2 / 0.2e1, t128 * qJD(4) + t108 * t109, -t132 * qJD(4) + t108 * t110, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t107, -t103 * t107, t107 ^ 2 / 0.2e1, (-t122 * t97 + t124 * t98) * t107 + t96 * t103, -(t122 * t98 + t124 * t97) * t107 + t96 * t104;];
T_reg = t1;
