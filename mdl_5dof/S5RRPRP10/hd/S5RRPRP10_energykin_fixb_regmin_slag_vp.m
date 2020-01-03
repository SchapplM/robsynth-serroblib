% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (110->34), mult. (256->74), div. (0->0), fcn. (122->4), ass. (0->30)
t136 = -pkin(2) - pkin(7);
t124 = qJD(1) ^ 2;
t135 = t124 / 0.2e1;
t123 = cos(qJ(2));
t134 = t123 * t124;
t121 = sin(qJ(2));
t126 = -qJ(3) * t121 - pkin(1);
t107 = (t136 * t123 + t126) * qJD(1);
t131 = t121 * qJD(1);
t130 = pkin(6) * t131 + qJD(3);
t108 = pkin(3) * t131 + t136 * qJD(2) + t130;
t120 = sin(qJ(4));
t122 = cos(qJ(4));
t133 = t122 * t107 + t120 * t108;
t132 = qJD(1) * t123;
t114 = -pkin(6) * t132 - qJD(2) * qJ(3);
t129 = qJD(1) * qJD(2);
t109 = pkin(3) * t132 - t114;
t128 = t121 * t129;
t127 = t123 * t129;
t125 = -t120 * t107 + t122 * t108;
t115 = qJD(4) + t131;
t113 = -qJD(2) * pkin(2) + t130;
t112 = t122 * qJD(2) - t120 * t132;
t111 = t120 * qJD(2) + t122 * t132;
t110 = (-pkin(2) * t123 + t126) * qJD(1);
t103 = t111 * pkin(4) + qJD(5) + t109;
t102 = -t111 * qJ(5) + t133;
t101 = t115 * pkin(4) - t112 * qJ(5) + t125;
t1 = [t135, 0, 0, t121 ^ 2 * t135, t121 * t134, t128, t127, qJD(2) ^ 2 / 0.2e1, pkin(1) * t134 - pkin(6) * t128, -t124 * pkin(1) * t121 - pkin(6) * t127, (t113 * t121 - t114 * t123) * qJD(1), t113 * qJD(2) + t110 * t132, -t114 * qJD(2) - t110 * t131, t110 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1, t112 ^ 2 / 0.2e1, -t112 * t111, t112 * t115, -t111 * t115, t115 ^ 2 / 0.2e1, t109 * t111 + t125 * t115, t109 * t112 - t133 * t115, -t101 * t112 - t102 * t111, t102 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1;];
T_reg = t1;
