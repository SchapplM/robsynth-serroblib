% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:40
% EndTime: 2019-12-05 15:54:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (110->32), mult. (286->73), div. (0->0), fcn. (197->8), ass. (0->31)
t130 = cos(qJ(4));
t116 = sin(pkin(9));
t120 = sin(qJ(2));
t111 = qJD(2) * qJ(3) + t120 * qJD(1);
t125 = pkin(6) * qJD(2) + t111;
t104 = t125 * t116;
t117 = cos(pkin(9));
t105 = t125 * t117;
t119 = sin(qJ(4));
t129 = -t119 * t104 + t130 * t105;
t128 = qJD(2) * t116;
t127 = qJD(2) * t117;
t126 = qJD(1) * qJD(2);
t124 = -t130 * t104 - t119 * t105;
t122 = cos(qJ(2));
t123 = -t122 * qJD(1) + qJD(3);
t108 = (-pkin(3) * t117 - pkin(2)) * qJD(2) + t123;
t121 = cos(qJ(5));
t118 = sin(qJ(5));
t115 = qJD(4) + qJD(5);
t114 = t117 ^ 2;
t113 = t116 ^ 2;
t109 = -qJD(2) * pkin(2) + t123;
t107 = (t130 * t116 + t117 * t119) * qJD(2);
t106 = t119 * t128 - t130 * t127;
t99 = t106 * pkin(4) + t108;
t98 = -t118 * t106 + t121 * t107;
t97 = t121 * t106 + t118 * t107;
t96 = -t106 * pkin(7) + t129;
t95 = qJD(4) * pkin(4) - t107 * pkin(7) + t124;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t122 * t126, -t120 * t126, -t109 * t127, t109 * t128, (t113 + t114) * t111 * qJD(2), t109 ^ 2 / 0.2e1 + (t114 / 0.2e1 + t113 / 0.2e1) * t111 ^ 2, t107 ^ 2 / 0.2e1, -t107 * t106, t107 * qJD(4), -t106 * qJD(4), qJD(4) ^ 2 / 0.2e1, t124 * qJD(4) + t108 * t106, -t129 * qJD(4) + t108 * t107, t98 ^ 2 / 0.2e1, -t98 * t97, t98 * t115, -t97 * t115, t115 ^ 2 / 0.2e1, (-t118 * t96 + t121 * t95) * t115 + t99 * t97, -(t118 * t95 + t121 * t96) * t115 + t99 * t98;];
T_reg = t1;
