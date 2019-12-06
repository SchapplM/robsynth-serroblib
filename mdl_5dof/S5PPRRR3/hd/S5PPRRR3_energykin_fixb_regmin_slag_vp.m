% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:58
% EndTime: 2019-12-05 15:16:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->26), mult. (156->67), div. (0->0), fcn. (101->8), ass. (0->28)
t118 = qJD(3) ^ 2;
t129 = t118 / 0.2e1;
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t110 = sin(pkin(9));
t127 = qJD(1) * t110;
t128 = t114 * qJD(2) + t117 * t127;
t111 = cos(pkin(9));
t126 = qJD(1) * t111;
t116 = cos(qJ(4));
t125 = qJD(3) * t116;
t121 = t117 * qJD(2) - t114 * t127;
t124 = (-qJD(3) * pkin(3) - t121) * qJD(3);
t123 = qJD(3) * qJD(4);
t122 = t116 * t126;
t102 = qJD(3) * pkin(6) + t128;
t113 = sin(qJ(4));
t120 = t116 * t102 - t113 * t126;
t119 = qJD(1) ^ 2;
t115 = cos(qJ(5));
t112 = sin(qJ(5));
t109 = qJD(4) + qJD(5);
t104 = (t112 * t116 + t113 * t115) * qJD(3);
t103 = t112 * t113 * qJD(3) - t115 * t125;
t99 = (-pkin(4) * t116 - pkin(3)) * qJD(3) - t121;
t98 = pkin(7) * t125 + t120;
t97 = -t122 + qJD(4) * pkin(4) + (-pkin(7) * qJD(3) - t102) * t113;
t1 = [t119 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t110 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1) * t119, t129, t121 * qJD(3), -t128 * qJD(3), t113 ^ 2 * t129, t113 * t118 * t116, t113 * t123, t116 * t123, qJD(4) ^ 2 / 0.2e1, (-t113 * t102 - t122) * qJD(4) - t116 * t124, -t120 * qJD(4) + t113 * t124, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * t109, -t103 * t109, t109 ^ 2 / 0.2e1, (-t112 * t98 + t115 * t97) * t109 + t99 * t103, -(t112 * t97 + t115 * t98) * t109 + t99 * t104;];
T_reg = t1;
