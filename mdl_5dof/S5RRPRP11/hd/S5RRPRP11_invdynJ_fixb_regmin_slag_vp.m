% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP11
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
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:07
% EndTime: 2019-12-31 20:14:14
% DurationCPUTime: 2.60s
% Computational Cost: add. (2176->376), mult. (4460->458), div. (0->0), fcn. (2570->6), ass. (0->190)
t112 = sin(qJ(2));
t174 = qJD(1) * qJD(2);
t162 = t112 * t174;
t115 = cos(qJ(2));
t172 = t115 * qJDD(1);
t238 = -t162 + t172;
t226 = pkin(3) + pkin(6);
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t152 = g(1) * t116 + g(2) * t113;
t219 = g(3) * t112;
t126 = -t152 * t115 - t219;
t227 = pkin(2) + pkin(7);
t179 = qJD(4) * t227;
t177 = t112 * qJD(1);
t83 = qJD(4) + t177;
t120 = t83 * t179 + t126;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t178 = t111 * qJD(2);
t185 = qJD(1) * t115;
t57 = t114 * t185 + t178;
t20 = t57 * qJD(4) - t114 * qJDD(2) + t238 * t111;
t163 = t111 * t185;
t21 = -qJD(4) * t163 + t111 * qJDD(2) + (qJD(2) * qJD(4) + t238) * t114;
t106 = qJDD(2) * qJ(3);
t107 = qJD(2) * qJD(3);
t90 = pkin(6) * t172;
t169 = t106 + t107 + t90;
t184 = qJD(2) * t112;
t64 = t226 * t184;
t26 = pkin(3) * t172 - qJD(1) * t64 + t169;
t175 = t114 * qJD(2);
t59 = -t163 + t175;
t3 = t21 * pkin(4) + t20 * qJ(5) - t59 * qJD(5) + t26;
t237 = t3 + t120;
t216 = t57 * t83;
t236 = t20 - t216;
t214 = t59 * t83;
t235 = -t21 + t214;
t143 = t111 * t83;
t201 = qJD(2) * t57;
t161 = t115 * t174;
t173 = t112 * qJDD(1);
t133 = t161 + t173;
t56 = qJDD(4) + t133;
t45 = t114 * t56;
t234 = -t83 * t143 - t201 + t45;
t182 = qJD(4) * t111;
t233 = -t83 * t182 + t45;
t92 = pkin(6) * t177;
t232 = qJD(3) + t92;
t105 = g(3) * t115;
t202 = qJ(3) * t115;
t147 = pkin(7) * t112 - t202;
t176 = t112 * qJD(3);
t124 = t147 * qJD(2) - t176;
t98 = t112 * qJ(3);
t160 = -pkin(1) - t98;
t132 = -t227 * t115 + t160;
t82 = pkin(2) * t162;
t12 = t124 * qJD(1) + t132 * qJDD(1) + t82;
t181 = qJD(4) * t114;
t81 = pkin(6) * t161;
t89 = pkin(6) * t173;
t165 = qJDD(3) + t81 + t89;
t25 = t133 * pkin(3) - t227 * qJDD(2) + t165;
t188 = pkin(3) * t177 + t232;
t39 = -t227 * qJD(2) + t188;
t171 = t111 * t25 + t114 * t12 + t39 * t181;
t36 = t132 * qJD(1);
t199 = qJD(4) * t36;
t190 = t116 * t111;
t193 = t113 * t114;
t49 = t112 * t190 + t193;
t189 = t116 * t114;
t194 = t113 * t111;
t51 = -t112 * t194 + t189;
t231 = -g(1) * t49 + g(2) * t51 + (-t199 + t105) * t111 + t171;
t200 = qJD(2) * t59;
t210 = t111 * t56;
t228 = t83 ^ 2;
t230 = t114 * t228 + t200 + t210;
t229 = t59 ^ 2;
t225 = t56 * pkin(4);
t14 = t111 * t39 + t114 * t36;
t8 = t83 * qJ(5) + t14;
t224 = t8 * t83;
t223 = g(1) * t113;
t220 = g(2) * t116;
t102 = t115 * pkin(2);
t218 = t115 * pkin(7);
t217 = t14 * t83;
t215 = t59 * t57;
t96 = pkin(2) * t177;
t42 = t147 * qJD(1) + t96;
t93 = pkin(6) * t185;
t94 = pkin(3) * t185;
t65 = t93 + t94;
t213 = t111 * t65 + t114 * t42;
t149 = pkin(4) * t114 + qJ(5) * t111;
t141 = -pkin(3) - t149;
t212 = t149 * qJD(4) - t114 * qJD(5) - t141 * t177 + t232;
t204 = t102 + t98;
t69 = -pkin(1) - t204;
t54 = t69 - t218;
t71 = t226 * t112;
t211 = t111 * t71 + t114 * t54;
t209 = t112 * t83;
t208 = t114 * t59;
t207 = t227 * t56;
t206 = t20 * t114;
t205 = t56 * qJ(5);
t108 = qJD(2) * qJ(3);
t70 = -t93 - t108;
t203 = pkin(6) * qJDD(2);
t198 = qJDD(2) * pkin(2);
t197 = t112 * t113;
t196 = t112 * t116;
t119 = qJD(1) ^ 2;
t195 = t112 * t119;
t192 = t113 * t115;
t191 = t115 * t116;
t13 = -t111 * t36 + t114 * t39;
t187 = qJD(5) - t13;
t72 = t226 * t115;
t109 = t112 ^ 2;
t110 = t115 ^ 2;
t186 = t109 - t110;
t183 = qJD(2) * t115;
t180 = qJD(4) * t115;
t170 = g(1) * t196 + g(2) * t197 - t105;
t46 = t94 - t70;
t166 = t115 * t195;
t164 = g(3) * t204;
t159 = -t111 * t12 + t114 * t25 - t36 * t181 - t39 * t182;
t158 = t116 * pkin(1) + pkin(2) * t191 + t113 * pkin(6) + qJ(3) * t196;
t157 = -t89 + t170;
t156 = -qJD(2) * pkin(2) + qJD(3);
t48 = -t112 * t189 + t194;
t50 = t112 * t193 + t190;
t155 = -g(1) * t50 - g(2) * t48;
t154 = -g(1) * t51 - g(2) * t49;
t118 = qJD(2) ^ 2;
t153 = pkin(6) * t118 + t220;
t6 = -t83 * pkin(4) + t187;
t151 = t111 * t6 + t114 * t8;
t148 = t111 * pkin(4) - t114 * qJ(5);
t67 = t156 + t92;
t145 = t112 * t70 + t115 * t67;
t142 = t160 - t102;
t47 = t142 * qJD(1);
t139 = t47 * t177 + qJDD(3) - t157;
t137 = -t83 * t181 - t210;
t136 = -0.2e1 * pkin(1) * t174 - t203;
t95 = pkin(2) * t184;
t34 = t95 + t124;
t66 = t226 * t183;
t135 = t111 * t66 + t114 * t34 + t71 * t181 - t54 * t182;
t134 = -qJ(3) * t183 - t176;
t130 = 0.2e1 * qJDD(1) * pkin(1) - t153;
t129 = t203 + (-qJD(1) * t69 - t47) * qJD(2);
t15 = t57 * pkin(4) - t59 * qJ(5) + t46;
t128 = t83 * t15 - t207;
t127 = g(1) * t48 - g(2) * t50 + t114 * t105 + t159;
t22 = t134 * qJD(1) + t142 * qJDD(1) + t82;
t44 = t134 + t95;
t123 = qJD(1) * t44 + qJDD(1) * t69 + t153 + t22;
t37 = pkin(6) * t162 - t169;
t41 = t165 - t198;
t122 = t145 * qJD(2) + t41 * t112 - t37 * t115;
t121 = t15 * t59 + qJDD(5) - t127;
t103 = t116 * pkin(6);
t87 = g(1) * t192;
t80 = qJ(3) * t191;
t78 = qJ(3) * t192;
t68 = qJ(3) + t148;
t62 = -qJ(3) * t185 + t96;
t33 = t149 * t115 + t72;
t27 = t59 * pkin(4) + t57 * qJ(5);
t24 = -t112 * pkin(4) + t111 * t54 - t114 * t71;
t23 = t112 * qJ(5) + t211;
t17 = -pkin(4) * t185 + t111 * t42 - t114 * t65;
t16 = qJ(5) * t185 + t213;
t9 = (-t148 * qJD(4) + qJD(5) * t111) * t115 + (-pkin(6) + t141) * t184;
t5 = -pkin(4) * t183 + t211 * qJD(4) + t111 * t34 - t114 * t66;
t4 = qJ(5) * t183 + t112 * qJD(5) + t135;
t2 = qJDD(5) - t159 - t225;
t1 = t83 * qJD(5) - t36 * t182 + t171 + t205;
t7 = [qJDD(1), -t220 + t223, t152, t109 * qJDD(1) + 0.2e1 * t112 * t161, 0.2e1 * t112 * t172 - 0.2e1 * t186 * t174, qJDD(2) * t112 + t118 * t115, qJDD(2) * t115 - t118 * t112, 0, t136 * t112 + t130 * t115 + t87, t136 * t115 + (-t130 - t223) * t112, (t109 + t110) * qJDD(1) * pkin(6) + t122 - t152, t129 * t112 + t123 * t115 - t87, t129 * t115 + (-t123 + t223) * t112, pkin(6) * t122 - g(1) * t103 - g(2) * t158 - t142 * t223 + t22 * t69 + t47 * t44, -t180 * t208 + (t115 * t20 + t184 * t59) * t111, (-t111 * t57 + t208) * t184 + (t111 * t21 + t206 + (t111 * t59 + t114 * t57) * qJD(4)) * t115, (t178 * t83 - t20) * t112 + (t137 + t200) * t115, (t175 * t83 - t21) * t112 + (-t201 - t233) * t115, t56 * t112 + t183 * t83, t159 * t112 + t13 * t183 - t64 * t57 + t72 * t21 + ((-qJD(4) * t71 - t34) * t83 - t54 * t56 - t46 * t180) * t111 + ((-qJD(4) * t54 + t66) * t83 + t71 * t56 + t26 * t115 - t46 * t184) * t114 + t154, -t135 * t83 - t211 * t56 - t64 * t59 - t72 * t20 + ((qJD(2) * t46 + t199) * t111 - t171) * t112 + (-t14 * qJD(2) - t26 * t111 - t181 * t46) * t115 - t155, t33 * t21 - t24 * t56 - t5 * t83 + t9 * t57 + (-t15 * t175 - t2) * t112 + (-qJD(2) * t6 + t3 * t114 - t15 * t182) * t115 + t154, -t24 * t20 - t23 * t21 - t4 * t57 + t5 * t59 + t87 + t151 * t184 + (-t220 - t1 * t114 - t111 * t2 + (t111 * t8 - t114 * t6) * qJD(4)) * t115, t33 * t20 + t23 * t56 + t4 * t83 - t9 * t59 + (-t15 * t178 + t1) * t112 + (qJD(2) * t8 + t3 * t111 + t15 * t181) * t115 + t155, t1 * t23 + t8 * t4 + t3 * t33 + t15 * t9 + t2 * t24 + t6 * t5 - g(1) * (t116 * pkin(3) + t51 * pkin(4) + t50 * qJ(5) + t103) - g(2) * (t49 * pkin(4) + pkin(7) * t191 + t48 * qJ(5) + t158) + (-g(1) * (t142 - t218) - g(2) * pkin(3)) * t113; 0, 0, 0, -t166, t186 * t119, t173, t172, qJDD(2), pkin(1) * t195 + t157, t219 - t90 + (pkin(1) * t119 + t152) * t115, (-pkin(2) * t112 + t202) * qJDD(1) + ((-t70 - t108) * t112 + (t156 - t67) * t115) * qJD(1), -t62 * t185 + t139 - 0.2e1 * t198, 0.2e1 * t106 + 0.2e1 * t107 + t90 + (qJD(1) * t62 - g(3)) * t112 + (qJD(1) * t47 - t152) * t115, -t37 * qJ(3) - t70 * qJD(3) - t41 * pkin(2) - t47 * t62 - g(1) * (-pkin(2) * t196 + t80) - g(2) * (-pkin(2) * t197 + t78) - t164 - t145 * qJD(1) * pkin(6), -t143 * t59 - t206, (-t21 - t214) * t114 + (t20 + t216) * t111, (-t111 * t209 - t115 * t59) * qJD(1) + t233, (-t114 * t209 + t115 * t57) * qJD(1) + t137, -t83 * t185, -t13 * t185 + qJ(3) * t21 + t188 * t57 + (-t207 + (t46 - t65) * t83) * t114 + (t26 + (t42 + t179) * t83 + t126) * t111, -qJ(3) * t20 + t213 * t83 + t14 * t185 + t188 * t59 + (-t46 * t83 + t207) * t111 + (t26 + t120) * t114, t111 * t237 + t128 * t114 + t17 * t83 + t6 * t185 + t68 * t21 + t212 * t57, t16 * t57 - t17 * t59 + (-t8 * t177 - t227 * t20 + t2 + (t227 * t57 - t8) * qJD(4)) * t114 + (-t6 * t177 + t227 * t21 - t1 + (-t227 * t59 - t6) * qJD(4)) * t111 + t170, t128 * t111 - t114 * t237 - t16 * t83 - t8 * t185 + t68 * t20 - t212 * t59, t3 * t68 - t8 * t16 - t6 * t17 - g(1) * t80 - g(2) * t78 - t164 + t212 * t15 + (-g(3) * pkin(7) - t148 * t152) * t115 + (-g(3) * t148 + t152 * t227) * t112 - (qJD(4) * t151 + t1 * t111 - t2 * t114) * t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, qJDD(2) + t166, -t109 * t119 - t118, t70 * qJD(2) + t139 - t198 + t81, 0, 0, 0, 0, 0, t234, -t230, t234, t111 * t235 + t114 * t236, t230, -t15 * qJD(2) + (-t2 + t224) * t114 + (t6 * t83 + t1) * t111 - t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t57 ^ 2 + t229, -t236, t235, t56, -t46 * t59 + t127 + t217, t13 * t83 + t46 * t57 - t231, -t27 * t57 - t121 + t217 + 0.2e1 * t225, pkin(4) * t20 - t21 * qJ(5) + (-t14 + t8) * t59 + (t6 - t187) * t57, 0.2e1 * t205 - t15 * t57 + t27 * t59 + (0.2e1 * qJD(5) - t13) * t83 + t231, t1 * qJ(5) - t2 * pkin(4) - t15 * t27 - t6 * t14 - g(1) * (-t48 * pkin(4) + t49 * qJ(5)) - g(2) * (t50 * pkin(4) - t51 * qJ(5)) + t187 * t8 + t149 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 + t215, -t236, -t228 - t229, t121 - t224 - t225;];
tau_reg = t7;
