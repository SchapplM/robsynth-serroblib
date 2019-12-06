% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:56
% EndTime: 2019-12-05 17:53:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (132->34), mult. (335->81), div. (0->0), fcn. (209->8), ass. (0->31)
t128 = qJD(1) ^ 2;
t134 = t128 / 0.2e1;
t121 = sin(pkin(8));
t115 = (pkin(1) * t121 + pkin(6)) * qJD(1);
t127 = cos(qJ(3));
t118 = t127 * qJD(2);
t125 = sin(qJ(3));
t108 = qJD(3) * pkin(3) + t118 + (-qJ(4) * qJD(1) - t115) * t125;
t132 = qJD(1) * t127;
t133 = t125 * qJD(2) + t127 * t115;
t109 = qJ(4) * t132 + t133;
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t101 = t120 * t108 + t122 * t109;
t131 = qJD(1) * qJD(3);
t123 = cos(pkin(8));
t130 = -pkin(1) * t123 - pkin(2);
t100 = t122 * t108 - t120 * t109;
t111 = qJD(4) + (-pkin(3) * t127 + t130) * qJD(1);
t126 = cos(qJ(5));
t124 = sin(qJ(5));
t119 = qJD(3) + qJD(5);
t116 = t130 * qJD(1);
t113 = (t120 * t127 + t122 * t125) * qJD(1);
t112 = (-t120 * t125 + t122 * t127) * qJD(1);
t104 = -t112 * pkin(4) + t111;
t103 = t124 * t112 + t126 * t113;
t102 = -t126 * t112 + t124 * t113;
t99 = t112 * pkin(7) + t101;
t98 = qJD(3) * pkin(4) - t113 * pkin(7) + t100;
t1 = [t134, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t121 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t128, t125 ^ 2 * t134, t125 * t128 * t127, t125 * t131, t127 * t131, qJD(3) ^ 2 / 0.2e1, -t116 * t132 + (-t125 * t115 + t118) * qJD(3), t116 * t125 * qJD(1) - t133 * qJD(3), -t100 * t113 + t101 * t112, t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1, t103 ^ 2 / 0.2e1, -t103 * t102, t103 * t119, -t102 * t119, t119 ^ 2 / 0.2e1, t104 * t102 + (-t124 * t99 + t126 * t98) * t119, t104 * t103 - (t124 * t98 + t126 * t99) * t119;];
T_reg = t1;
