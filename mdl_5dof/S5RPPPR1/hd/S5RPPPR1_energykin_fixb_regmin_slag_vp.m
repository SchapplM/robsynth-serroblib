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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:31
% EndTime: 2020-01-03 11:20:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (134->35), mult. (330->83), div. (0->0), fcn. (191->8), ass. (0->29)
t113 = sin(pkin(7));
t108 = (pkin(1) * t113 + qJ(3)) * qJD(1);
t112 = sin(pkin(8));
t115 = cos(pkin(8));
t105 = t112 * qJD(2) + t115 * t108;
t111 = sin(pkin(9));
t114 = cos(pkin(9));
t116 = cos(pkin(7));
t122 = -pkin(1) * t116 - pkin(2);
t99 = qJD(3) + (-pkin(3) * t115 - qJ(4) * t112 + t122) * qJD(1);
t95 = t114 * t105 + t111 * t99;
t125 = t112 * t114;
t124 = qJD(1) * t112;
t123 = t115 * qJD(1);
t121 = t111 * t124;
t94 = -t111 * t105 + t114 * t99;
t104 = t115 * qJD(2) - t112 * t108;
t103 = qJD(4) - t104;
t119 = qJD(1) ^ 2;
t118 = cos(qJ(5));
t117 = sin(qJ(5));
t109 = -qJD(5) + t123;
t107 = t122 * qJD(1) + qJD(3);
t102 = (-t111 * t117 + t114 * t118) * t124;
t101 = (t111 * t118 + t114 * t117) * t124;
t96 = pkin(4) * t121 + t103;
t93 = -pkin(6) * t121 + t95;
t92 = (-pkin(4) * t115 - pkin(6) * t125) * qJD(1) + t94;
t1 = [t119 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t113 ^ 2 / 0.2e1 + t116 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t119, -t107 * t123, t107 * t124, (-t104 * t112 + t105 * t115) * qJD(1), t105 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1, (t103 * t111 * t112 - t115 * t94) * qJD(1), (t103 * t125 + t115 * t95) * qJD(1), (-t111 * t95 - t114 * t94) * t124, t95 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1, t102 ^ 2 / 0.2e1, -t102 * t101, -t102 * t109, t101 * t109, t109 ^ 2 / 0.2e1, -(-t117 * t93 + t118 * t92) * t109 + t96 * t101, (t117 * t92 + t118 * t93) * t109 + t96 * t102;];
T_reg = t1;
