% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:50
% EndTime: 2019-12-31 17:47:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (131->37), mult. (322->80), div. (0->0), fcn. (170->6), ass. (0->33)
t137 = qJD(1) ^ 2;
t145 = t137 / 0.2e1;
t144 = qJ(2) * t137;
t131 = sin(pkin(8));
t134 = cos(pkin(7));
t143 = t131 * t134;
t132 = sin(pkin(7));
t139 = -qJ(3) * t132 - pkin(1);
t114 = qJD(2) + ((-pkin(2) - qJ(4)) * t134 + t139) * qJD(1);
t142 = qJD(1) * t132;
t122 = qJ(2) * t142 + qJD(3);
t119 = pkin(3) * t142 + t122;
t133 = cos(pkin(8));
t111 = t133 * t114 + t131 * t119;
t141 = qJD(1) * t134;
t140 = qJ(2) ^ 2 * t145;
t120 = qJD(4) + (pkin(3) + qJ(2)) * t141;
t110 = -t131 * t114 + t133 * t119;
t136 = cos(qJ(5));
t135 = sin(qJ(5));
t130 = t134 ^ 2;
t129 = t132 ^ 2;
t128 = -qJD(1) * pkin(1) + qJD(2);
t125 = t130 * t144;
t123 = t130 * t140;
t121 = t133 * t141 + qJD(5);
t118 = qJD(2) + (-pkin(2) * t134 + t139) * qJD(1);
t116 = (t132 * t135 - t136 * t143) * qJD(1);
t115 = (t132 * t136 + t135 * t143) * qJD(1);
t112 = (pkin(4) * t133 + pkin(6) * t131) * t141 + t120;
t109 = pkin(6) * t142 + t111;
t108 = -pkin(4) * t142 - t110;
t1 = [t145, 0, 0, -t128 * t141, t128 * t142, t129 * t144 + t125, t123 + t129 * t140 + t128 ^ 2 / 0.2e1, t122 * t142 + t125, t118 * t141, -t118 * t142, t118 ^ 2 / 0.2e1 + t123 + t122 ^ 2 / 0.2e1, (t120 * t133 * t134 + t110 * t132) * qJD(1), (-t111 * t132 - t120 * t143) * qJD(1), (t110 * t131 - t111 * t133) * t141, t111 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1, t116 ^ 2 / 0.2e1, t116 * t115, t116 * t121, t115 * t121, t121 ^ 2 / 0.2e1, (-t135 * t109 + t136 * t112) * t121 - t108 * t115, -(t136 * t109 + t135 * t112) * t121 + t108 * t116;];
T_reg = t1;
