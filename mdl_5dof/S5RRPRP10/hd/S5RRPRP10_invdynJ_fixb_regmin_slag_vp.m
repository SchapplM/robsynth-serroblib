% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:08
% EndTime: 2019-12-31 20:11:14
% DurationCPUTime: 2.30s
% Computational Cost: add. (1648->323), mult. (3464->412), div. (0->0), fcn. (1997->6), ass. (0->189)
t112 = sin(qJ(2));
t168 = qJD(1) * qJD(2);
t153 = t112 * t168;
t115 = cos(qJ(2));
t166 = t115 * qJDD(1);
t236 = -t153 + t166;
t224 = pkin(3) + pkin(6);
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t117 = -pkin(2) - pkin(7);
t97 = t112 * qJ(3);
t151 = -pkin(1) - t97;
t127 = t117 * t115 + t151;
t30 = t127 * qJD(1);
t180 = qJD(1) * t112;
t91 = pkin(6) * t180;
t231 = qJD(3) + t91;
t182 = pkin(3) * t180 + t231;
t33 = t117 * qJD(2) + t182;
t13 = t111 * t33 + t114 * t30;
t200 = qJ(3) * t115;
t140 = pkin(7) * t112 - t200;
t174 = qJD(3) * t112;
t124 = qJD(2) * t140 - t174;
t81 = pkin(2) * t153;
t11 = qJD(1) * t124 + qJDD(1) * t127 + t81;
t152 = t115 * t168;
t167 = t112 * qJDD(1);
t128 = t152 + t167;
t80 = pkin(6) * t152;
t88 = pkin(6) * t167;
t157 = qJDD(3) + t80 + t88;
t22 = t128 * pkin(3) + t117 * qJDD(2) + t157;
t150 = -t11 * t111 + t114 * t22;
t122 = -qJD(4) * t13 + t150;
t178 = qJD(2) * t111;
t179 = qJD(1) * t115;
t55 = t114 * t179 + t178;
t19 = qJD(4) * t55 - t114 * qJDD(2) + t236 * t111;
t54 = qJDD(4) + t128;
t154 = t111 * t179;
t176 = qJD(2) * t114;
t57 = -t154 + t176;
t1 = pkin(4) * t54 + qJ(5) * t19 - qJD(5) * t57 + t122;
t104 = g(3) * t115;
t116 = cos(qJ(1));
t190 = t112 * t116;
t113 = sin(qJ(1));
t192 = t112 * t113;
t162 = g(1) * t190 + g(2) * t192 - t104;
t172 = qJD(4) * t114;
t165 = -t114 * t11 - t111 * t22 - t33 * t172;
t173 = qJD(4) * t111;
t20 = -qJD(4) * t154 + qJDD(2) * t111 + (qJD(2) * qJD(4) + t236) * t114;
t2 = -qJ(5) * t20 - qJD(5) * t55 - t173 * t30 - t165;
t12 = -t111 * t30 + t114 * t33;
t7 = -qJ(5) * t57 + t12;
t82 = qJD(4) + t180;
t6 = pkin(4) * t82 + t7;
t8 = -qJ(5) * t55 + t13;
t235 = -(t6 * t82 - t2) * t111 + (t8 * t82 + t1) * t114 - t162;
t213 = t55 * t82;
t234 = t19 - t213;
t212 = t57 * t82;
t233 = -t20 + t212;
t40 = t114 * t54;
t232 = -t82 * t173 + t40;
t142 = g(1) * t116 + g(2) * t113;
t230 = pkin(4) * t20 + qJDD(5);
t186 = t114 * t115;
t185 = t114 * t116;
t43 = -t111 * t113 + t112 * t185;
t188 = t113 * t114;
t45 = t111 * t116 + t112 * t188;
t227 = -g(1) * t43 - g(2) * t45 + g(3) * t186;
t107 = qJD(2) * qJ(3);
t92 = pkin(6) * t179;
t70 = -t92 - t107;
t93 = pkin(3) * t179;
t41 = t93 - t70;
t226 = t117 * t54 + t41 * t82;
t225 = t57 ^ 2;
t223 = -t7 + t6;
t219 = pkin(4) * t111;
t218 = g(1) * t113;
t215 = g(2) * t116;
t214 = g(3) * t112;
t101 = t115 * pkin(2);
t95 = pkin(2) * t180;
t37 = qJD(1) * t140 + t95;
t64 = t92 + t93;
t210 = t111 * t64 + t114 * t37;
t183 = qJ(5) - t117;
t194 = t111 * t112;
t49 = t114 * t64;
t209 = -qJD(5) * t114 + t173 * t183 + t111 * t37 - t49 - (pkin(4) * t115 - qJ(5) * t194) * qJD(1);
t67 = t183 * t114;
t208 = -qJ(5) * t114 * t180 - qJD(4) * t67 - qJD(5) * t111 - t210;
t202 = t101 + t97;
t69 = -pkin(1) - t202;
t51 = -pkin(7) * t115 + t69;
t71 = t224 * t112;
t207 = t111 * t71 + t114 * t51;
t206 = t111 * t54;
t205 = t114 * t19;
t204 = t114 * t57;
t72 = t224 * t115;
t201 = pkin(6) * qJDD(2);
t199 = qJD(2) * t55;
t198 = qJD(2) * t57;
t197 = qJD(4) * t30;
t196 = qJDD(2) * pkin(2);
t110 = -qJ(5) - pkin(7);
t195 = t110 * t115;
t193 = t111 * t115;
t191 = t112 * t114;
t119 = qJD(1) ^ 2;
t189 = t112 * t119;
t187 = t113 * t115;
t184 = t115 * t116;
t108 = t112 ^ 2;
t109 = t115 ^ 2;
t181 = t108 - t109;
t177 = qJD(2) * t112;
t175 = qJD(2) * t115;
t171 = qJD(4) * t115;
t170 = qJD(4) * t117;
t169 = qJD(5) * t115;
t164 = pkin(4) * t193;
t105 = qJDD(2) * qJ(3);
t106 = qJD(2) * qJD(3);
t89 = pkin(6) * t166;
t161 = t105 + t106 + t89;
t159 = t111 * t190;
t158 = t115 * t189;
t87 = pkin(4) * t114 + pkin(3);
t156 = t224 * qJD(2);
t155 = t111 * t171;
t149 = qJ(5) * t115 - t51;
t148 = t116 * pkin(1) + pkin(2) * t184 + t113 * pkin(6) + qJ(3) * t190;
t147 = -t88 + t162;
t146 = -qJD(2) * pkin(2) + qJD(3);
t145 = pkin(3) * t166 + t161;
t144 = qJD(1) * t156;
t118 = qJD(2) ^ 2;
t143 = pkin(6) * t118 + t215;
t68 = t146 + t91;
t139 = t112 * t70 + t115 * t68;
t138 = t111 * t82;
t137 = t151 - t101;
t42 = t137 * qJD(1);
t135 = t42 * t180 + qJDD(3) - t147;
t133 = -t172 * t82 - t206;
t132 = t142 * t115;
t131 = -0.2e1 * pkin(1) * t168 - t201;
t94 = pkin(2) * t177;
t27 = t94 + t124;
t65 = t224 * t175;
t130 = t111 * t65 + t114 * t27 + t71 * t172 - t173 * t51;
t129 = -qJ(3) * t175 - t174;
t126 = 0.2e1 * qJDD(1) * pkin(1) - t143;
t125 = t201 + (-qJD(1) * t69 - t42) * qJD(2);
t23 = -t112 * t144 + t145;
t123 = -t132 + t23 - t214;
t21 = qJD(1) * t129 + qJDD(1) * t137 + t81;
t39 = t129 + t94;
t121 = qJD(1) * t39 + qJDD(1) * t69 + t143 + t21;
t31 = pkin(6) * t153 - t161;
t36 = t157 - t196;
t120 = qJD(2) * t139 + t36 * t112 - t31 * t115;
t102 = t116 * pkin(6);
t85 = g(1) * t187;
t79 = qJ(3) * t184;
t77 = qJ(3) * t187;
t66 = t183 * t111;
t63 = t112 * t156;
t61 = -qJ(3) * t179 + t95;
t60 = t114 * t71;
t53 = t55 ^ 2;
t50 = t114 * t65;
t46 = -t111 * t192 + t185;
t44 = t159 + t188;
t24 = pkin(4) * t55 + qJD(5) + t41;
t16 = -qJ(5) * t186 + t207;
t15 = pkin(4) * t112 + t111 * t149 + t60;
t5 = t23 + t230;
t4 = -t114 * t169 + (t112 * t176 + t155) * qJ(5) + t130;
t3 = pkin(4) * t175 + t50 + t149 * t172 + (-qJ(5) * t177 - qJD(4) * t71 + t169 - t27) * t111;
t9 = [qJDD(1), -t215 + t218, t142, qJDD(1) * t108 + 0.2e1 * t112 * t152, 0.2e1 * t112 * t166 - 0.2e1 * t168 * t181, qJDD(2) * t112 + t115 * t118, qJDD(2) * t115 - t112 * t118, 0, t112 * t131 + t115 * t126 + t85, t131 * t115 + (-t126 - t218) * t112, (t108 + t109) * qJDD(1) * pkin(6) + t120 - t142, t112 * t125 + t115 * t121 - t85, t125 * t115 + (-t121 + t218) * t112, pkin(6) * t120 - g(1) * t102 - g(2) * t148 - t137 * t218 + t21 * t69 + t42 * t39, -t171 * t204 + (t115 * t19 + t177 * t57) * t111, (-t111 * t55 + t204) * t177 + (t111 * t20 + t205 + (t111 * t57 + t114 * t55) * qJD(4)) * t115, (t178 * t82 - t19) * t112 + (t133 + t198) * t115, (t176 * t82 - t20) * t112 + (-t199 - t232) * t115, t112 * t54 + t175 * t82, (-t111 * t27 + t50) * t82 + (-t111 * t51 + t60) * t54 + t150 * t112 - t63 * t55 + t72 * t20 + t23 * t186 - g(1) * t46 - g(2) * t44 + (t12 * t115 - t191 * t41) * qJD(2) + (-t13 * t112 - t41 * t193 - t207 * t82) * qJD(4), -t130 * t82 - t207 * t54 - t63 * t57 - t72 * t19 + g(1) * t45 - g(2) * t43 + ((qJD(2) * t41 + t197) * t111 + t165) * t112 + (-qJD(2) * t13 - t23 * t111 - t172 * t41) * t115, t15 * t19 - t16 * t20 - t3 * t57 - t4 * t55 + t85 + (-t111 * t6 + t114 * t8) * t177 + (-t215 + t1 * t111 - t114 * t2 + (t111 * t8 + t114 * t6) * qJD(4)) * t115, t2 * t16 + t8 * t4 + t1 * t15 + t6 * t3 + t5 * (pkin(4) * t186 + t72) - g(1) * (t116 * t87 + t102) - g(2) * (pkin(4) * t159 - t110 * t184 + t148) + (-g(1) * (-pkin(4) * t194 + t137 + t195) - g(2) * t87) * t113 + (-pkin(4) * t155 + (-pkin(6) - t87) * t177) * t24; 0, 0, 0, -t158, t181 * t119, t167, t166, qJDD(2), pkin(1) * t189 + t147, t214 - t89 + (pkin(1) * t119 + t142) * t115, (-pkin(2) * t112 + t200) * qJDD(1) + ((-t70 - t107) * t112 + (t146 - t68) * t115) * qJD(1), -t179 * t61 + t135 - 0.2e1 * t196, 0.2e1 * t105 + 0.2e1 * t106 + t89 + (qJD(1) * t61 - g(3)) * t112 + (qJD(1) * t42 - t142) * t115, -t31 * qJ(3) - t70 * qJD(3) - t36 * pkin(2) - t42 * t61 - g(1) * (-pkin(2) * t190 + t79) - g(2) * (-pkin(2) * t192 + t77) - g(3) * t202 - t139 * qJD(1) * pkin(6), -t138 * t57 - t205, (-t20 - t212) * t114 + (t19 + t213) * t111, (-t115 * t57 - t194 * t82) * qJD(1) + t232, (t115 * t55 - t191 * t82) * qJD(1) + t133, -t82 * t179, -t12 * t179 + qJ(3) * t20 - t49 * t82 + t182 * t55 + t226 * t114 + ((t37 - t170) * t82 + t123) * t111, -qJ(3) * t19 + t210 * t82 + t13 * t179 + t182 * t57 - t226 * t111 + (-t170 * t82 + t123) * t114, -t19 * t67 + t20 * t66 - t208 * t55 - t209 * t57 - t235, -t2 * t66 - t1 * t67 + t5 * (qJ(3) + t219) - g(1) * (t116 * t164 + t79) - g(2) * (t113 * t164 + t77) - g(3) * (-t195 + t202) + t208 * t8 + t209 * t6 + (pkin(4) * t172 + t231) * t24 + (t24 * t87 * qJD(1) - g(3) * t219 + t142 * (pkin(2) - t110)) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, qJDD(2) + t158, -t108 * t119 - t118, qJD(2) * t70 + t135 - t196 + t80, 0, 0, 0, 0, 0, -t138 * t82 - t199 + t40, -t114 * t82 ^ 2 - t198 - t206, t233 * t111 + t234 * t114, -qJD(2) * t24 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t53 + t225, -t234, t233, t54, t13 * t82 - t41 * t57 + t122 + t227, g(1) * t44 - g(2) * t46 + t12 * t82 + t41 * t55 + (t197 - t104) * t111 + t165, pkin(4) * t19 - t223 * t55, t223 * t8 + (-t24 * t57 + t1 + t227) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 - t225, t55 * t8 + t57 * t6 - t132 + (-g(3) - t144) * t112 + t145 + t230;];
tau_reg = t9;
