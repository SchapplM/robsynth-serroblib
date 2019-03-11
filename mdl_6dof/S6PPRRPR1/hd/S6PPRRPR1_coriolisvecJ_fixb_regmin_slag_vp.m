% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:28
% EndTime: 2019-03-08 18:47:35
% DurationCPUTime: 2.53s
% Computational Cost: add. (2753->313), mult. (7618->496), div. (0->0), fcn. (6662->14), ass. (0->174)
t134 = cos(pkin(6));
t119 = qJD(1) * t134 + qJD(2);
t128 = sin(pkin(12));
t130 = sin(pkin(6));
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t132 = cos(pkin(12));
t133 = cos(pkin(7));
t203 = t132 * t133;
t146 = (t128 * t140 + t137 * t203) * t130;
t129 = sin(pkin(7));
t208 = t129 * t137;
t61 = qJD(1) * t146 + t119 * t208;
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t167 = pkin(4) * t136 - qJ(5) * t139;
t89 = t167 * qJD(4) - qJD(5) * t136;
t233 = t61 - t89;
t197 = qJD(3) * t139;
t121 = -qJD(6) + t197;
t232 = qJD(6) + t121;
t127 = sin(pkin(13));
t131 = cos(pkin(13));
t190 = t131 * qJD(4);
t199 = qJD(3) * t136;
t103 = t127 * t199 - t190;
t195 = qJD(4) * t127;
t105 = t131 * t199 + t195;
t135 = sin(qJ(6));
t138 = cos(qJ(6));
t159 = t103 * t135 - t105 * t138;
t231 = t121 * t159;
t207 = t129 * t140;
t225 = (-t128 * t137 + t140 * t203) * t130;
t230 = t134 * t207 + t225;
t229 = qJD(1) * t225 + t119 * t207;
t194 = qJD(4) * t136;
t188 = pkin(9) * t194;
t118 = t127 * t188;
t209 = t127 * t139;
t200 = qJD(1) * t130;
t181 = t132 * t200;
t173 = t133 * t181;
t182 = t128 * t200;
t59 = -t137 * t182 + (t119 * t129 + t173) * t140;
t221 = t233 * t131 - t59 * t209 - t118;
t204 = t131 * t139;
t228 = -t233 * t127 - t59 * t204;
t58 = qJD(3) * pkin(9) + t61;
t79 = t119 * t133 - t129 * t181;
t227 = -t136 * t58 + t139 * t79;
t198 = qJD(3) * t137;
t179 = t129 * t198;
t196 = qJD(3) * t140;
t54 = t119 * t179 + t173 * t198 + t182 * t196;
t226 = qJD(3) * t61 - t54;
t224 = t159 * qJD(6);
t223 = pkin(10) + qJ(5);
t53 = t229 * qJD(3);
t16 = t139 * t53 + (qJD(5) + t227) * qJD(4);
t38 = t89 * qJD(3) + t54;
t7 = t127 * t38 + t131 * t16;
t158 = pkin(5) * t136 - pkin(10) * t204;
t150 = t158 * qJD(4);
t222 = -t150 + t221;
t205 = t131 * t136;
t220 = (-pkin(9) * t205 - pkin(10) * t209) * qJD(4) + t228;
t174 = t131 * t188;
t219 = t174 - t228;
t35 = t136 * t79 + t139 * t58;
t33 = qJD(4) * qJ(5) + t35;
t114 = -pkin(4) * t139 - qJ(5) * t136 - pkin(3);
t45 = t114 * qJD(3) - t59;
t13 = t127 * t45 + t131 * t33;
t109 = t167 * qJD(3);
t25 = t127 * t109 + t131 * t227;
t202 = t138 * t131;
t210 = t127 * t135;
t107 = -t202 + t210;
t151 = t107 * t139;
t218 = qJD(3) * t151 - t107 * qJD(6);
t108 = t127 * t138 + t131 * t135;
t152 = t108 * t139;
t217 = -qJD(3) * t152 + t108 * qJD(6);
t216 = qJD(3) * pkin(3);
t193 = qJD(4) * t139;
t17 = t136 * t53 + t58 * t193 + t79 * t194;
t215 = t127 * t17;
t214 = t131 * t17;
t189 = qJD(3) * qJD(4);
t177 = t139 * t189;
t191 = qJD(6) * t138;
t212 = -t103 * t191 + t177 * t202;
t142 = qJD(3) ^ 2;
t206 = t129 * t142;
t81 = pkin(9) * t204 + t127 * t114;
t201 = t136 ^ 2 - t139 ^ 2;
t192 = qJD(6) * t136;
t185 = t137 * t206;
t184 = pkin(5) * t127 + pkin(9);
t6 = -t127 * t16 + t131 * t38;
t4 = qJD(3) * t150 + t6;
t170 = t127 * t177;
t5 = -pkin(10) * t170 + t7;
t183 = -t135 * t5 + t138 * t4;
t180 = t127 * t197;
t178 = t129 * t196;
t176 = t136 * t189;
t12 = -t127 * t33 + t131 * t45;
t24 = t131 * t109 - t127 * t227;
t175 = -qJD(4) * pkin(4) + qJD(5);
t172 = t136 * t178;
t171 = t139 * t178;
t169 = t135 * t4 + t138 * t5;
t8 = -pkin(5) * t197 - pkin(10) * t105 + t12;
t9 = -pkin(10) * t103 + t13;
t168 = t135 * t9 - t138 * t8;
t2 = t135 * t8 + t138 * t9;
t69 = t134 * t208 + t146;
t90 = -t129 * t130 * t132 + t133 * t134;
t44 = t136 * t90 + t139 * t69;
t22 = -t127 * t44 - t131 * t230;
t23 = -t127 * t230 + t131 * t44;
t166 = -t135 * t23 + t138 * t22;
t165 = t135 * t22 + t138 * t23;
t102 = t131 * t114;
t67 = -pkin(10) * t205 + t102 + (-pkin(9) * t127 - pkin(5)) * t139;
t76 = -pkin(10) * t127 * t136 + t81;
t164 = -t135 * t76 + t138 * t67;
t163 = t135 * t67 + t138 * t76;
t95 = t133 * t136 + t139 * t208;
t72 = -t127 * t95 - t131 * t207;
t73 = -t127 * t207 + t131 * t95;
t162 = -t135 * t73 + t138 * t72;
t161 = t135 * t72 + t138 * t73;
t43 = t136 * t69 - t139 * t90;
t141 = qJD(4) ^ 2;
t157 = pkin(9) * t141 - t226;
t57 = -t59 - t216;
t156 = qJD(4) * (t57 + t59 - t216);
t94 = -t133 * t139 + t136 * t208;
t116 = t223 * t127;
t154 = pkin(10) * t180 + qJD(5) * t131 - qJD(6) * t116 - t25;
t117 = t223 * t131;
t153 = t158 * qJD(3) + qJD(5) * t127 + qJD(6) * t117 + t24;
t31 = t175 - t227;
t145 = qJD(4) * t152;
t144 = -qJ(5) * t194 + (t175 - t31) * t139;
t39 = (-qJD(6) * t105 - t170) * t135 + t212;
t40 = qJD(3) * t145 - t224;
t123 = -pkin(5) * t131 - pkin(4);
t110 = t184 * t136;
t99 = t184 * t193;
t91 = t138 * t103;
t87 = t107 * t136;
t86 = t108 * t136;
t80 = -pkin(9) * t209 + t102;
t75 = t95 * qJD(4) + t172;
t74 = -t94 * qJD(4) + t171;
t64 = t105 * t135 + t91;
t63 = t69 * qJD(3);
t62 = t230 * qJD(3);
t56 = t191 * t205 - t192 * t210 + t145;
t55 = -qJD(4) * t151 - t108 * t192;
t52 = t127 * t179 + t131 * t74;
t51 = -t127 * t74 + t131 * t179;
t29 = pkin(5) * t180 + t35;
t26 = pkin(5) * t103 + t31;
t21 = -t43 * qJD(4) + t139 * t62;
t20 = t44 * qJD(4) + t136 * t62;
t14 = pkin(5) * t170 + t17;
t11 = t127 * t63 + t131 * t21;
t10 = -t127 * t21 + t131 * t63;
t1 = [0, 0, 0, -t63 * qJD(3), -t62 * qJD(3), 0, 0, 0, 0, 0, -qJD(4) * t20 + (-t139 * t63 - t194 * t230) * qJD(3), -qJD(4) * t21 + (t136 * t63 - t193 * t230) * qJD(3), t103 * t20 + (-t10 * t139 + (t136 * t22 + t43 * t209) * qJD(4)) * qJD(3), t105 * t20 + (t11 * t139 + (-t136 * t23 + t43 * t204) * qJD(4)) * qJD(3), -t10 * t105 - t103 * t11 + (-t127 * t23 - t131 * t22) * t177, t10 * t12 + t11 * t13 + t17 * t43 + t20 * t31 + t22 * t6 + t23 * t7, 0, 0, 0, 0, 0 -(-qJD(6) * t165 + t10 * t138 - t11 * t135) * t121 + t166 * t176 + t20 * t64 + t43 * t40 (qJD(6) * t166 + t10 * t135 + t11 * t138) * t121 - t165 * t176 - t20 * t159 + t43 * t39; 0, 0, 0, -t185, -t140 * t206, 0, 0, 0, 0, 0, -t139 * t185 + (-t75 - t172) * qJD(4), t136 * t185 + (-t74 - t171) * qJD(4), t103 * t75 + (-t139 * t51 + (t136 * t72 + t94 * t209) * qJD(4)) * qJD(3), t105 * t75 + (t139 * t52 + (-t136 * t73 + t94 * t204) * qJD(4)) * qJD(3), -t103 * t52 - t105 * t51 + (-t127 * t73 - t131 * t72) * t177, t12 * t51 + t13 * t52 + t17 * t94 + t31 * t75 + t6 * t72 + t7 * t73, 0, 0, 0, 0, 0 -(-qJD(6) * t161 - t135 * t52 + t138 * t51) * t121 + t162 * t176 + t75 * t64 + t94 * t40 (qJD(6) * t162 + t135 * t51 + t138 * t52) * t121 - t161 * t176 - t75 * t159 + t94 * t39; 0, 0, 0, t226 (-t229 + t59) * qJD(3), 0.2e1 * t139 * t176, -0.2e1 * t201 * t189, t141 * t139, -t141 * t136, 0, t136 * t156 - t157 * t139, t157 * t136 + t139 * t156 (-t103 * t59 + t215 + (qJD(3) * t80 + t12) * qJD(4)) * t136 + (-t6 + (pkin(9) * t103 + t127 * t31) * qJD(4) + (t118 + t221) * qJD(3)) * t139 (-t105 * t59 + t214 + (-qJD(3) * t81 - t13) * qJD(4)) * t136 + (t7 + (pkin(9) * t105 + t131 * t31) * qJD(4) + (t174 - t219) * qJD(3)) * t139 (-t127 * t7 - t131 * t6) * t136 + t221 * t105 + t219 * t103 + (-t12 * t131 - t127 * t13 + (-t127 * t81 - t131 * t80) * qJD(3)) * t193, -t136 * t31 * t59 + t6 * t80 + t7 * t81 - t219 * t13 - t221 * t12 + (t136 * t17 + t31 * t193) * pkin(9), -t159 * t55 - t39 * t87, t159 * t56 - t39 * t86 + t40 * t87 - t55 * t64, -t121 * t55 - t139 * t39 + (-qJD(3) * t87 - t159) * t194, t121 * t56 + t139 * t40 + (-qJD(3) * t86 - t64) * t194 (-t121 - t197) * t194, -t183 * t139 + t99 * t64 + t110 * t40 + t14 * t86 + t26 * t56 + (t220 * t135 + t222 * t138) * t121 + (t121 * t163 + t139 * t2) * qJD(6) + (-t59 * t64 + (qJD(3) * t164 - t168) * qJD(4)) * t136, t169 * t139 - t99 * t159 + t110 * t39 - t14 * t87 + t26 * t55 + (-t222 * t135 + t220 * t138) * t121 + (t121 * t164 - t139 * t168) * qJD(6) + (t59 * t159 + (-qJD(3) * t163 - t2) * qJD(4)) * t136; 0, 0, 0, 0, 0, -t139 * t142 * t136, t201 * t142, 0, 0, 0, qJD(4) * t35 - t57 * t199 - t17 (-qJD(3) * t57 - t53) * t139, -t103 * t35 - t214 + (-t12 * t136 + t127 * t144 + t139 * t24) * qJD(3), -t105 * t35 + t215 + (t13 * t136 + t131 * t144 - t139 * t25) * qJD(3), t103 * t25 + t105 * t24 + (-qJD(5) * t103 + t12 * t197 + t7) * t131 + (qJD(5) * t105 + t13 * t197 - t6) * t127, -pkin(4) * t17 - t12 * t24 - t13 * t25 - t31 * t35 + (-t12 * t127 + t13 * t131) * qJD(5) + (-t127 * t6 + t131 * t7) * qJ(5), t108 * t39 - t159 * t218, -t107 * t39 - t108 * t40 + t159 * t217 - t218 * t64, -t218 * t121 + (qJD(4) * t108 + t159) * t199, t217 * t121 + (-qJD(4) * t107 + t64) * t199, t121 * t199, t14 * t107 + t123 * t40 - t29 * t64 + t217 * t26 + (t135 * t154 + t138 * t153) * t121 + ((-t116 * t138 - t117 * t135) * qJD(4) + t168) * t199, t14 * t108 + t123 * t39 + t29 * t159 + t218 * t26 + (-t135 * t153 + t138 * t154) * t121 + (-(-t116 * t135 + t117 * t138) * qJD(4) + t2) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t105 + t195) * t197 (t103 + t190) * t197, -t103 ^ 2 - t105 ^ 2, t103 * t13 + t105 * t12 + t17, 0, 0, 0, 0, 0, t40 + t231, t91 * t121 + (-t170 + (-qJD(6) + t121) * t105) * t135 + t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159 * t64, t159 ^ 2 - t64 ^ 2, -t121 * t64 + t39, -t108 * t177 + t224 + t231, t176, t26 * t159 - t232 * t2 + t183, t232 * t168 + t26 * t64 - t169;];
tauc_reg  = t1;
