% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:08
% EndTime: 2019-12-31 19:06:12
% DurationCPUTime: 1.32s
% Computational Cost: add. (1424->205), mult. (2092->263), div. (0->0), fcn. (1360->10), ass. (0->143)
t97 = cos(qJ(4));
t157 = qJD(4) * t97;
t200 = -qJD(5) * t97 - t157;
t178 = sin(qJ(1));
t179 = cos(qJ(1));
t199 = g(1) * t178 - g(2) * t179;
t100 = qJD(4) ^ 2;
t87 = qJDD(1) - qJDD(3);
t98 = cos(qJ(3));
t171 = t87 * t98;
t150 = qJD(1) - qJD(3);
t197 = t150 ^ 2;
t95 = sin(qJ(3));
t198 = (t100 + t197) * t95 + t171;
t196 = -qJDD(2) + t199;
t153 = qJ(2) * qJD(1);
t99 = -pkin(1) - pkin(2);
t195 = qJD(3) * t153 - t99 * qJDD(1) - qJDD(2);
t94 = sin(qJ(4));
t96 = cos(qJ(5));
t166 = t94 * t96;
t93 = sin(qJ(5));
t45 = t93 * t97 + t166;
t88 = qJD(4) + qJD(5);
t190 = t88 * t45;
t194 = t150 * t190;
t116 = t200 * t96;
t43 = -t178 * t95 - t179 * t98;
t44 = -t178 * t98 + t179 * t95;
t134 = g(1) * t43 + g(2) * t44;
t192 = qJD(4) * t150;
t151 = qJD(1) * qJD(2);
t152 = qJ(2) * qJDD(1);
t64 = t99 * qJD(1) + qJD(2);
t191 = qJD(3) * t64 + t151 + t152;
t164 = t97 * t87;
t167 = t94 * t87;
t8 = t96 * t164 - t93 * t167 - t194;
t162 = t98 * qJ(2) + t95 * t99;
t141 = -t95 * qJ(2) + t98 * t99;
t160 = pkin(1) * qJDD(1);
t189 = t160 + t196;
t135 = g(1) * t44 - g(2) * t43;
t37 = -t95 * t153 + t64 * t98;
t183 = pkin(4) * t97;
t78 = -pkin(3) - t183;
t18 = -t150 * t78 - t37;
t165 = t96 * t97;
t169 = t93 * t94;
t42 = -t165 + t169;
t120 = t191 * t95 + t195 * t98;
t185 = pkin(3) * t87;
t11 = t120 + t185;
t158 = qJD(4) * t94;
t121 = -t150 * t158 + t164;
t6 = t121 * pkin(4) + t11;
t92 = qJ(4) + qJ(5);
t81 = cos(t92);
t188 = -t135 * t81 - t18 * t190 - t6 * t42;
t19 = t88 * t169 + t116;
t80 = sin(t92);
t187 = t135 * t80 + t18 * t19 - t6 * t45;
t186 = pkin(7) + pkin(8);
t184 = pkin(3) * t150;
t47 = -pkin(7) + t162;
t180 = pkin(8) - t47;
t38 = t98 * t153 + t64 * t95;
t24 = -pkin(7) * t150 + t38;
t143 = -pkin(8) * t150 + t24;
t17 = t143 * t97;
t177 = t17 * t96;
t149 = t150 * t169;
t27 = t150 * t165 - t149;
t29 = t45 * t150;
t176 = t29 * t27;
t36 = qJD(2) * t95 + t162 * qJD(3);
t175 = t36 * t150;
t174 = t38 * t150;
t86 = qJDD(4) + qJDD(5);
t173 = t42 * t86;
t172 = t45 * t86;
t170 = t150 * t94;
t90 = t94 ^ 2;
t161 = -t97 ^ 2 + t90;
t156 = qJD(5) * t93;
t155 = qJD(5) * t96;
t148 = pkin(4) * t158;
t146 = t150 * t157;
t142 = qJD(4) * t186;
t140 = qJD(4) * t180;
t46 = pkin(3) - t141;
t136 = -t38 + t148;
t16 = t143 * t94;
t7 = t116 * t150 + t88 * t149 - t93 * t164 - t87 * t166;
t131 = t19 * t29 + t45 * t7;
t15 = qJD(4) * pkin(4) - t16;
t130 = -t15 * t93 - t177;
t129 = t19 * t88 - t172;
t128 = t190 * t88 + t173;
t25 = t180 * t94;
t26 = t180 * t97;
t127 = t25 * t96 + t26 * t93;
t126 = t25 * t93 - t26 * t96;
t61 = t186 * t94;
t62 = t186 * t97;
t125 = -t61 * t96 - t62 * t93;
t124 = -t61 * t93 + t62 * t96;
t122 = -t11 + t135;
t118 = t88 * t42;
t117 = g(1) * t179 + g(2) * t178;
t113 = t191 * t98 - t195 * t95;
t10 = -pkin(7) * t87 + t113;
t23 = -t37 + t184;
t114 = t150 * t23 - t10 - t134;
t112 = -pkin(7) * qJDD(4) + (t23 + t37 + t184) * qJD(4);
t111 = t19 * t27 + t190 * t29 - t42 * t7 - t45 * t8;
t35 = qJD(2) * t98 + t141 * qJD(3);
t110 = -qJDD(4) * t47 + (-t150 * t46 - t23 - t35) * qJD(4);
t109 = -qJDD(4) * t95 + 0.2e1 * t98 * t192;
t108 = -t117 + 0.2e1 * t151;
t107 = t120 - t135;
t106 = pkin(7) * t100 - t122 + t174 + t185;
t2 = -t24 * t157 + qJDD(4) * pkin(4) - t10 * t94 + (t146 + t167) * pkin(8);
t105 = t17 * t156 + t18 * t27 + (-t17 * t88 - t2) * t93 - g(3) * t80 - t134 * t81;
t104 = t100 * t47 - t46 * t87 + t122 - t175;
t103 = t113 + t134;
t3 = -t121 * pkin(8) + t10 * t97 - t24 * t158;
t102 = g(3) * t81 + t130 * qJD(5) - t134 * t80 + t18 * t29 + t96 * t2 - t93 * t3;
t101 = qJD(1) ^ 2;
t54 = qJDD(4) * t97 - t100 * t94;
t53 = qJDD(4) * t94 + t100 * t97;
t49 = t97 * t142;
t48 = t94 * t142;
t39 = t46 + t183;
t30 = -0.2e1 * t94 * t146 - t87 * t90;
t22 = t36 - t148;
t21 = t161 * t192 - t94 * t164;
t13 = t97 * t140 - t35 * t94;
t12 = t94 * t140 + t35 * t97;
t9 = -t27 ^ 2 + t29 ^ 2;
t5 = -t29 * t88 - t8;
t4 = t27 * t88 + t7;
t1 = [qJDD(1), t199, t117, 0.2e1 * t160 + t196, t108 + 0.2e1 * t152, t189 * pkin(1) + (t108 + t152) * qJ(2), t87, -t141 * t87 + t107 + t175, t150 * t35 + t162 * t87 + t103, -t30, -0.2e1 * t21, -t53, -t54, 0, -t104 * t97 + t110 * t94, t104 * t94 + t110 * t97, -t131, -t111, t129, t128, 0, t22 * t27 + t39 * t8 + (-t126 * qJD(5) - t12 * t93 + t13 * t96) * t88 + t127 * t86 + t188, -t22 * t29 + t39 * t7 - (t127 * qJD(5) + t12 * t96 + t13 * t93) * t88 - t126 * t86 + t187; 0, 0, 0, -qJDD(1), -t101, -t101 * qJ(2) - t189, 0, -t197 * t95 - t171, -t197 * t98 + t87 * t95, 0, 0, 0, 0, 0, t109 * t94 - t198 * t97, t109 * t97 + t198 * t94, 0, 0, 0, 0, 0, (-t8 + t194) * t98 + ((t94 * t156 + t93 * t158 + t116) * t88 - t172 - t150 * t27) * t95, (-t150 * t118 - t7) * t98 + (-(-t94 * t155 - t96 * t158 + t200 * t93) * t88 + t173 + t150 * t29) * t95; 0, 0, 0, 0, 0, 0, -t87, -t107 - t174, -t150 * t37 - t103, t30, 0.2e1 * t21, t53, t54, 0, -t106 * t97 + t112 * t94, t106 * t94 + t112 * t97, t131, t111, -t129, -t128, 0, t78 * t8 + (-t124 * qJD(5) + t48 * t93 - t49 * t96) * t88 + t125 * t86 + t37 * t190 + t136 * t27 - t188, t78 * t7 - (t125 * qJD(5) - t48 * t96 - t49 * t93) * t88 - t124 * t86 - t37 * t118 - t136 * t29 - t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * t197 * t97, t161 * t197, -t167, -t164, qJDD(4), g(3) * t97 + t114 * t94, -g(3) * t94 + t114 * t97, -t176, t9, t4, t5, t86, -(t16 * t93 - t177) * t88 + (-t88 * t156 + t27 * t170 + t96 * t86) * pkin(4) + t102, (-qJD(5) * t15 - t16 * t88 - t3) * t96 + (-t88 * t155 - t29 * t170 - t93 * t86) * pkin(4) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t9, t4, t5, t86, -t130 * t88 + t102, (-t3 + (-qJD(5) + t88) * t15) * t96 + t105;];
tau_reg = t1;
