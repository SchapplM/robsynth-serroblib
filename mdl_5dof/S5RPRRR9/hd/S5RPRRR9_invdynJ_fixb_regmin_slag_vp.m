% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:03
% EndTime: 2019-12-31 19:08:10
% DurationCPUTime: 2.65s
% Computational Cost: add. (3692->294), mult. (9090->384), div. (0->0), fcn. (7375->14), ass. (0->173)
t139 = cos(qJ(5));
t192 = qJD(5) * t139;
t136 = sin(qJ(4));
t140 = cos(qJ(4));
t134 = cos(pkin(9));
t141 = cos(qJ(3));
t202 = t141 * t134;
t133 = sin(pkin(9));
t137 = sin(qJ(3));
t205 = t137 * t133;
t93 = -t202 + t205;
t85 = t93 * qJD(1);
t94 = t133 * t141 + t134 * t137;
t86 = t94 * qJD(1);
t59 = t136 * t86 + t140 * t85;
t249 = t139 * t59;
t257 = t192 + t249;
t131 = pkin(9) + qJ(3);
t127 = qJ(4) + t131;
t119 = sin(t127);
t138 = sin(qJ(1));
t142 = cos(qJ(1));
t165 = g(1) * t142 + g(2) * t138;
t256 = t165 * t119;
t132 = qJD(3) + qJD(4);
t208 = t59 * t132;
t194 = qJD(4) * t140;
t195 = qJD(4) * t136;
t182 = qJD(1) * t205;
t189 = t134 * qJDD(1);
t190 = t133 * qJDD(1);
t196 = qJD(3) * t141;
t184 = qJD(1) * t134 * t196 + t137 * t189 + t141 * t190;
t66 = -qJD(3) * t182 + t184;
t162 = t137 * t190 - t141 * t189;
t88 = t94 * qJD(3);
t67 = qJD(1) * t88 + t162;
t25 = -t136 * t67 + t140 * t66 - t194 * t85 - t195 * t86;
t255 = t25 + t208;
t128 = qJDD(3) + qJDD(4);
t135 = sin(qJ(5));
t161 = -t136 * t85 + t140 * t86;
t193 = qJD(5) * t135;
t13 = t128 * t135 + t132 * t192 + t139 * t25 - t161 * t193;
t52 = t132 * t135 + t139 * t161;
t14 = qJD(5) * t52 - t128 * t139 + t135 * t25;
t50 = -t132 * t139 + t135 * t161;
t254 = t13 * t139 - t135 * t14 - t257 * t50;
t11 = t13 * t135;
t253 = t257 * t52 + t11;
t26 = qJD(4) * t161 + t136 * t66 + t140 * t67;
t23 = qJDD(5) + t26;
t17 = t135 * t23;
t221 = t52 * t161;
t247 = qJD(5) + t59;
t54 = t247 * t192;
t252 = t247 * t249 + t17 - t221 + t54;
t219 = pkin(6) + qJ(2);
t106 = t219 * t133;
t95 = qJD(1) * t106;
t107 = t219 * t134;
t96 = qJD(1) * t107;
t159 = t137 * t95 - t141 * t96;
t49 = -pkin(7) * t85 - t159;
t214 = t136 * t49;
t236 = -t137 * t96 - t141 * t95;
t48 = -pkin(7) * t86 + t236;
t47 = qJD(3) * pkin(3) + t48;
t29 = t140 * t47 - t214;
t27 = -pkin(4) * t132 - t29;
t251 = t27 * t59;
t120 = cos(t127);
t224 = g(3) * t120;
t191 = qJD(1) * qJD(2);
t229 = qJDD(1) * t219 + t191;
t74 = t229 * t133;
t75 = t229 * t134;
t174 = -t137 * t75 - t141 * t74;
t21 = qJDD(3) * pkin(3) - t66 * pkin(7) + qJD(3) * t159 + t174;
t160 = -t137 * t74 + t141 * t75;
t24 = -t67 * pkin(7) + qJD(3) * t236 + t160;
t210 = t140 * t49;
t30 = t136 * t47 + t210;
t231 = qJD(4) * t30 + t136 * t24 - t140 * t21;
t3 = -t128 * pkin(4) + t231;
t250 = t3 + t224;
t248 = t161 * t59;
t209 = t161 * t132;
t245 = -t26 + t209;
t243 = t161 ^ 2 - t59 ^ 2;
t114 = g(3) * t119;
t230 = (qJD(4) * t47 + t24) * t140 + t136 * t21 - t49 * t195;
t121 = -pkin(2) * t134 - pkin(1);
t100 = qJD(1) * t121 + qJD(2);
t70 = t85 * pkin(3) + t100;
t242 = t120 * t165 + t70 * t59 + t114 - t230;
t240 = pkin(4) * t161;
t220 = t161 * t50;
t122 = pkin(3) * t136 + pkin(8);
t239 = (pkin(3) * t86 + pkin(8) * t59 + qJD(5) * t122 + t240) * t247;
t238 = (pkin(8) * t247 + t240) * t247;
t237 = t247 * t161;
t235 = qJ(2) * qJDD(1);
t28 = pkin(8) * t132 + t30;
t31 = t59 * pkin(4) - pkin(8) * t161 + t70;
t6 = -t135 * t28 + t139 * t31;
t234 = t139 * t256 - t161 * t6 + t27 * t193;
t7 = t135 * t31 + t139 * t28;
t233 = t135 * t250 + t161 * t7 + t27 * t192;
t232 = -t161 * t70 - t224 - t231 + t256;
t180 = t128 * pkin(8) + qJD(5) * t31 + t230;
t172 = -t106 * t141 - t107 * t137;
t55 = -pkin(7) * t94 + t172;
t216 = -t106 * t137 + t107 * t141;
t56 = -pkin(7) * t93 + t216;
t36 = t136 * t55 + t140 * t56;
t68 = t136 * t94 + t140 * t93;
t69 = -t136 * t93 + t140 * t94;
t72 = pkin(3) * t93 + t121;
t37 = pkin(4) * t68 - pkin(8) * t69 + t72;
t87 = t93 * qJD(3);
t40 = -qJD(4) * t68 - t136 * t88 - t140 * t87;
t35 = t136 * t56 - t140 * t55;
t151 = qJD(2) * t202 - t106 * t196 + (-qJD(2) * t133 - qJD(3) * t107) * t137;
t42 = -t88 * pkin(7) + t151;
t145 = -qJD(2) * t94 - qJD(3) * t216;
t43 = t87 * pkin(7) + t145;
t8 = -qJD(4) * t35 + t136 * t43 + t140 * t42;
t228 = -(qJD(5) * t37 + t8) * t247 - t180 * t68 - t36 * t23 + t27 * t40 + t3 * t69;
t227 = pkin(3) * t88;
t223 = t27 * t69;
t222 = t37 * t23;
t215 = t135 * t52;
t206 = qJDD(1) * pkin(1);
t204 = t138 * t135;
t203 = t138 * t139;
t201 = t142 * t135;
t200 = t142 * t139;
t197 = t133 ^ 2 + t134 ^ 2;
t185 = t69 * t193;
t99 = qJDD(1) * t121 + qJDD(2);
t53 = t67 * pkin(3) + t99;
t5 = t26 * pkin(4) - t25 * pkin(8) + t53;
t179 = qJD(5) * t28 - t5;
t170 = t135 * t247;
t169 = t197 * qJD(1) ^ 2;
t167 = 0.2e1 * t197;
t32 = t136 * t48 + t210;
t166 = pkin(3) * t195 - t32;
t164 = g(1) * t138 - g(2) * t142;
t163 = t23 * t69 + t247 * t40;
t18 = t139 * t23;
t158 = t18 - (t135 * t59 + t193) * t247;
t157 = -t180 + t114;
t155 = -pkin(8) * t23 + t247 * t29 + t251;
t153 = -t164 - t206;
t124 = qJDD(2) - t206;
t152 = -t124 - t153;
t33 = t140 * t48 - t214;
t149 = -t122 * t23 + t251 + (-pkin(3) * t194 + t33) * t247;
t148 = t167 * t191 - t165;
t126 = cos(t131);
t125 = sin(t131);
t123 = -pkin(3) * t140 - pkin(4);
t82 = t120 * t200 + t204;
t81 = -t120 * t201 + t203;
t80 = -t120 * t203 + t201;
t79 = t120 * t204 + t200;
t41 = qJD(4) * t69 - t136 * t87 + t140 * t88;
t15 = pkin(4) * t41 - pkin(8) * t40 + t227;
t9 = qJD(4) * t36 + t136 * t42 - t140 * t43;
t4 = t139 * t5;
t1 = [qJDD(1), t164, t165, t152 * t134, -t152 * t133, t167 * t235 + t148, (-t124 + t164) * pkin(1) + (t197 * t235 + t148) * qJ(2), t66 * t94 - t86 * t87, -t66 * t93 - t67 * t94 + t85 * t87 - t86 * t88, -qJD(3) * t87 + qJDD(3) * t94, -qJD(3) * t88 - qJDD(3) * t93, 0, qJD(3) * t145 + qJDD(3) * t172 + t100 * t88 + t121 * t67 + t126 * t164 + t99 * t93, -qJD(3) * t151 - qJDD(3) * t216 - t100 * t87 + t121 * t66 - t125 * t164 + t99 * t94, t161 * t40 + t25 * t69, -t161 * t41 - t25 * t68 - t26 * t69 - t40 * t59, t128 * t69 + t132 * t40, -t128 * t68 - t132 * t41, 0, t120 * t164 - t35 * t128 - t9 * t132 + t227 * t59 + t72 * t26 + t70 * t41 + t53 * t68, -t119 * t164 - t36 * t128 - t8 * t132 + t161 * t227 + t72 * t25 + t70 * t40 + t53 * t69, -t52 * t185 + (t13 * t69 + t40 * t52) * t139, (-t139 * t50 - t215) * t40 + (-t11 - t139 * t14 + (t135 * t50 - t139 * t52) * qJD(5)) * t69, t13 * t68 + t139 * t163 - t185 * t247 + t52 * t41, -t135 * t163 - t14 * t68 - t50 * t41 - t54 * t69, t23 * t68 + t247 * t41, -g(1) * t80 - g(2) * t82 + t35 * t14 + t4 * t68 + t6 * t41 + t9 * t50 + (t15 * t247 + t222 + (-t247 * t36 - t28 * t68 + t223) * qJD(5)) * t139 + t228 * t135, -g(1) * t79 - g(2) * t81 + t35 * t13 - t7 * t41 + t9 * t52 + (-(-qJD(5) * t36 + t15) * t247 - t222 + t179 * t68 - qJD(5) * t223) * t135 + t228 * t139; 0, 0, 0, -t189, t190, -t169, -qJ(2) * t169 + qJDD(2) + t153, 0, 0, 0, 0, 0, 0.2e1 * t86 * qJD(3) + t162, (-t85 - t182) * qJD(3) + t184, 0, 0, 0, 0, 0, t26 + t209, t25 - t208, 0, 0, 0, 0, 0, t158 - t220, -t139 * t247 ^ 2 - t17 - t221; 0, 0, 0, 0, 0, 0, 0, t86 * t85, -t85 ^ 2 + t86 ^ 2, (t85 - t182) * qJD(3) + t184, -t162, qJDD(3), -g(3) * t126 - t100 * t86 + t125 * t165 + t174, g(3) * t125 + t100 * t85 + t126 * t165 - t160, t248, t243, t255, t245, t128, t32 * t132 + (t128 * t140 - t132 * t195 - t59 * t86) * pkin(3) + t232, t33 * t132 + (-t128 * t136 - t132 * t194 - t161 * t86) * pkin(3) + t242, t253, -t215 * t247 + t254, t252, t158 + t220, -t237, t123 * t14 + t166 * t50 + (-t250 - t239) * t139 + t149 * t135 + t234, t123 * t13 + t166 * t52 + t149 * t139 + (-t256 + t239) * t135 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t243, t255, t245, t128, t30 * t132 + t232, t29 * t132 + t242, t253, -t170 * t52 + t254, t252, -t170 * t247 + t18 + t220, -t237, -pkin(4) * t14 - t30 * t50 + t155 * t135 + (-t250 - t238) * t139 + t234, -pkin(4) * t13 - t30 * t52 + t155 * t139 + (-t256 + t238) * t135 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t50, -t50 ^ 2 + t52 ^ 2, t247 * t50 + t13, t247 * t52 - t14, t23, -g(1) * t81 + g(2) * t79 + t135 * t157 - t192 * t28 + t247 * t7 - t27 * t52 + t4, g(1) * t82 - g(2) * t80 + t135 * t179 + t139 * t157 + t247 * t6 + t27 * t50;];
tau_reg = t1;
