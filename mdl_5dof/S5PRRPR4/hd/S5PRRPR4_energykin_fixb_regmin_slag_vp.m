% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:21
% EndTime: 2019-12-05 16:23:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (113->30), mult. (286->73), div. (0->0), fcn. (191->8), ass. (0->31)
t130 = qJD(2) ^ 2;
t137 = t130 / 0.2e1;
t125 = sin(qJ(3));
t126 = sin(qJ(2));
t119 = qJD(2) * pkin(6) + t126 * qJD(1);
t131 = qJ(4) * qJD(2) + t119;
t113 = qJD(3) * pkin(3) - t131 * t125;
t128 = cos(qJ(3));
t115 = t131 * t128;
t122 = sin(pkin(9));
t123 = cos(pkin(9));
t106 = t122 * t113 + t123 * t115;
t129 = cos(qJ(2));
t134 = t129 * qJD(1);
t136 = qJD(2) * (-qJD(2) * pkin(2) - t134);
t135 = qJD(3) * t119;
t133 = qJD(1) * qJD(2);
t132 = qJD(2) * qJD(3);
t105 = t123 * t113 - t122 * t115;
t118 = -t134 + qJD(4) + (-pkin(3) * t128 - pkin(2)) * qJD(2);
t127 = cos(qJ(5));
t124 = sin(qJ(5));
t121 = qJD(3) + qJD(5);
t117 = (t122 * t128 + t123 * t125) * qJD(2);
t116 = (-t122 * t125 + t123 * t128) * qJD(2);
t109 = -t116 * pkin(4) + t118;
t108 = t124 * t116 + t127 * t117;
t107 = -t127 * t116 + t124 * t117;
t104 = t116 * pkin(7) + t106;
t103 = qJD(3) * pkin(4) - t117 * pkin(7) + t105;
t1 = [qJD(1) ^ 2 / 0.2e1, t137, t129 * t133, -t126 * t133, t125 ^ 2 * t137, t125 * t130 * t128, t125 * t132, t128 * t132, qJD(3) ^ 2 / 0.2e1, -t125 * t135 - t128 * t136, t125 * t136 - t128 * t135, -t105 * t117 + t106 * t116, t106 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1, t108 ^ 2 / 0.2e1, -t108 * t107, t108 * t121, -t107 * t121, t121 ^ 2 / 0.2e1, (t127 * t103 - t124 * t104) * t121 + t109 * t107, -(t124 * t103 + t127 * t104) * t121 + t109 * t108;];
T_reg = t1;
