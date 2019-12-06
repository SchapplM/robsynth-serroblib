% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:05
% EndTime: 2019-12-05 17:29:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (134->35), mult. (330->83), div. (0->0), fcn. (191->8), ass. (0->29)
t114 = sin(pkin(8));
t117 = cos(pkin(8));
t118 = cos(pkin(7));
t124 = -pkin(1) * t118 - pkin(2);
t101 = qJD(3) + (-pkin(3) * t117 - qJ(4) * t114 + t124) * qJD(1);
t115 = sin(pkin(7));
t110 = (pkin(1) * t115 + qJ(3)) * qJD(1);
t107 = t114 * qJD(2) + t117 * t110;
t113 = sin(pkin(9));
t116 = cos(pkin(9));
t97 = t113 * t101 + t116 * t107;
t127 = t114 * t116;
t126 = qJD(1) * t114;
t125 = t117 * qJD(1);
t123 = t113 * t126;
t96 = t116 * t101 - t113 * t107;
t106 = t117 * qJD(2) - t114 * t110;
t105 = qJD(4) - t106;
t121 = qJD(1) ^ 2;
t120 = cos(qJ(5));
t119 = sin(qJ(5));
t111 = -qJD(5) + t125;
t109 = t124 * qJD(1) + qJD(3);
t104 = (-t113 * t119 + t116 * t120) * t126;
t103 = (t113 * t120 + t116 * t119) * t126;
t98 = pkin(4) * t123 + t105;
t95 = -pkin(6) * t123 + t97;
t94 = (-pkin(4) * t117 - pkin(6) * t127) * qJD(1) + t96;
t1 = [t121 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t115 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t121, -t109 * t125, t109 * t126, (-t106 * t114 + t107 * t117) * qJD(1), t107 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1, (t105 * t113 * t114 - t117 * t96) * qJD(1), (t105 * t127 + t117 * t97) * qJD(1), (-t113 * t97 - t116 * t96) * t126, t97 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1, t104 ^ 2 / 0.2e1, -t104 * t103, -t104 * t111, t103 * t111, t111 ^ 2 / 0.2e1, -(-t119 * t95 + t120 * t94) * t111 + t98 * t103, (t119 * t94 + t120 * t95) * t111 + t98 * t104;];
T_reg = t1;
