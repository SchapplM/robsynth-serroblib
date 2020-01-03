% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:42
% EndTime: 2019-12-31 17:53:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->34), mult. (288->67), div. (0->0), fcn. (150->4), ass. (0->28)
t125 = qJD(1) ^ 2;
t133 = t125 / 0.2e1;
t132 = qJ(2) * t125;
t120 = sin(pkin(7));
t130 = qJD(1) * t120;
t110 = qJ(2) * t130 + qJD(3);
t108 = -pkin(6) * t130 + t110;
t121 = cos(pkin(7));
t129 = qJD(1) * t121;
t109 = (-pkin(6) + qJ(2)) * t129;
t123 = sin(qJ(4));
t124 = cos(qJ(4));
t131 = t123 * t108 + t124 * t109;
t117 = -qJD(1) * pkin(1) + qJD(2);
t128 = qJ(2) ^ 2 * t133;
t104 = -pkin(2) * t129 - qJ(3) * t130 + t117;
t102 = pkin(3) * t129 - t104;
t127 = t124 * t108 - t123 * t109;
t119 = t121 ^ 2;
t118 = t120 ^ 2;
t114 = t119 * t132;
t111 = t119 * t128;
t107 = (t120 * t124 - t121 * t123) * qJD(1);
t106 = (t120 * t123 + t121 * t124) * qJD(1);
t101 = qJD(4) * qJ(5) + t131;
t100 = -qJD(4) * pkin(4) + qJD(5) - t127;
t99 = t106 * pkin(4) - t107 * qJ(5) + t102;
t1 = [t133, 0, 0, -t117 * t129, t117 * t130, t118 * t132 + t114, t111 + t118 * t128 + t117 ^ 2 / 0.2e1, -t104 * t129, t110 * t130 + t114, -t104 * t130, t111 + t104 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1, t107 ^ 2 / 0.2e1, -t107 * t106, t107 * qJD(4), -t106 * qJD(4), qJD(4) ^ 2 / 0.2e1, t127 * qJD(4) + t102 * t106, -t131 * qJD(4) + t102 * t107, -t100 * qJD(4) + t99 * t106, t100 * t107 - t101 * t106, t101 * qJD(4) - t99 * t107, t101 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1;];
T_reg = t1;
