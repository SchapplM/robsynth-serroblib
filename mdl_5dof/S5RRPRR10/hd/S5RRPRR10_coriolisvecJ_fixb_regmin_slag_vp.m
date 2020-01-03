% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:39
% EndTime: 2019-12-31 20:26:47
% DurationCPUTime: 2.44s
% Computational Cost: add. (4056->336), mult. (12549->488), div. (0->0), fcn. (10113->10), ass. (0->177)
t132 = sin(pkin(10));
t134 = cos(pkin(10));
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t156 = t132 * t141 + t134 * t138;
t133 = sin(pkin(5));
t198 = qJD(1) * t133;
t110 = t156 * t198;
t140 = cos(qJ(4));
t135 = cos(pkin(5));
t190 = t135 * qJD(1);
t171 = qJD(2) + t190;
t120 = t140 * t171;
t137 = sin(qJ(4));
t77 = t110 * t137 - t120;
t76 = qJD(5) + t77;
t202 = t134 * t141;
t183 = t133 * t202;
t119 = qJD(1) * t183;
t182 = t138 * t198;
t107 = -t132 * t182 + t119;
t104 = qJD(4) - t107;
t197 = qJD(2) * t133;
t181 = t138 * t197;
t228 = pkin(2) * t181;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t145 = -t140 * t110 - t137 * t171;
t43 = -t104 * t139 - t136 * t145;
t227 = t104 * t43;
t113 = t156 * t133;
t108 = qJD(2) * t113;
t102 = qJD(1) * t108;
t225 = pkin(1) * t138;
t188 = t135 * t225;
t203 = t133 * t141;
t219 = pkin(7) + qJ(3);
t105 = t203 * t219 + t188;
t204 = t133 * t138;
t82 = -qJD(2) * t105 - qJD(3) * t204;
t143 = qJD(1) * t82;
t224 = pkin(1) * t141;
t187 = t135 * t224;
t124 = qJD(1) * t187;
t122 = t124 * qJD(2);
t179 = t219 * t138;
t144 = (-qJD(2) * t179 + qJD(3) * t141) * t133;
t73 = qJD(1) * t144 + t122;
t29 = t132 * t143 + t134 * t73;
t167 = qJD(1) * t181;
t103 = qJD(2) * t119 - t132 * t167;
t121 = pkin(2) * t167;
t51 = pkin(3) * t102 - pkin(8) * t103 + t121;
t175 = t137 * t29 - t140 * t51;
t96 = t105 * qJD(1);
t215 = t134 * t96;
t169 = t133 * t179;
t80 = qJD(2) * pkin(2) + t124 + (t135 * pkin(2) - t169) * qJD(1);
t39 = t132 * t80 + t215;
t32 = pkin(8) * t171 + t39;
t166 = (-pkin(2) * t141 - pkin(1)) * t133;
t154 = qJD(1) * t166;
t114 = qJD(3) + t154;
t50 = -t107 * pkin(3) - t110 * pkin(8) + t114;
t18 = t137 * t50 + t140 * t32;
t6 = -t102 * pkin(4) + qJD(4) * t18 + t175;
t226 = (-t145 * pkin(4) + pkin(9) * t76) * t76 + t6;
t195 = qJD(4) * t137;
t41 = qJD(4) * t120 + t103 * t140 - t110 * t195;
t45 = t104 * t136 - t139 * t145;
t15 = qJD(5) * t45 - t102 * t139 + t136 * t41;
t42 = -qJD(4) * t145 + t137 * t103;
t28 = t132 * t73 - t134 * t143;
t10 = pkin(4) * t42 - pkin(9) * t41 + t28;
t13 = pkin(9) * t104 + t18;
t84 = t132 * t96;
t38 = t134 * t80 - t84;
t31 = -pkin(3) * t171 - t38;
t16 = pkin(4) * t77 + pkin(9) * t145 + t31;
t163 = t13 * t136 - t139 * t16;
t193 = qJD(4) * t140;
t151 = t137 * t51 + t140 * t29 + t193 * t50 - t195 * t32;
t5 = pkin(9) * t102 + t151;
t1 = -qJD(5) * t163 + t136 * t10 + t139 * t5;
t142 = qJD(1) ^ 2;
t223 = t43 * t76;
t222 = t45 * t76;
t206 = t107 * t140;
t69 = -t110 * t139 + t136 * t206;
t221 = t69 * t76;
t70 = t110 * t136 + t139 * t206;
t220 = t70 * t76;
t94 = (pkin(2) + t224) * t135 - t169;
t62 = t105 * t134 + t132 * t94;
t53 = pkin(8) * t135 + t62;
t112 = t132 * t204 - t183;
t68 = t112 * pkin(3) - t113 * pkin(8) + t166;
t159 = t137 * t68 + t140 * t53;
t95 = -qJD(1) * t169 + t124;
t55 = t134 * t95 - t84;
t64 = pkin(2) * t182 + pkin(3) * t110 - pkin(8) * t107;
t218 = t137 * t64 + t140 * t55;
t217 = t107 * t45;
t216 = t110 * t77;
t214 = t136 * t42;
t213 = t139 * t42;
t174 = t139 * t76;
t191 = qJD(5) * t139;
t192 = qJD(5) * t136;
t14 = t102 * t136 + t104 * t191 + t139 * t41 + t145 * t192;
t212 = t14 * t136;
t211 = t77 * t104;
t210 = t145 * t104;
t209 = t145 * t110;
t54 = t132 * t95 + t215;
t208 = -t54 + t104 * (pkin(4) * t137 - pkin(9) * t140);
t207 = t104 * t137;
t129 = t133 ^ 2;
t205 = t129 * t142;
t201 = t137 * t102;
t199 = t138 ^ 2 - t141 ^ 2;
t196 = qJD(4) * t136;
t194 = qJD(4) * t139;
t189 = qJD(1) * qJD(2);
t186 = t76 * t196;
t185 = t76 * t194;
t184 = t138 * t205;
t128 = -pkin(2) * t134 - pkin(3);
t180 = t133 * t135 * t142;
t178 = t129 * t189;
t125 = qJD(2) * t187;
t81 = t125 + t144;
t33 = t132 * t81 - t134 * t82;
t61 = -t105 * t132 + t134 * t94;
t115 = -pkin(4) * t140 - pkin(9) * t137 + t128;
t173 = pkin(9) * t110 - qJD(5) * t115 + t218;
t172 = t104 * t140;
t170 = qJD(2) + 0.2e1 * t190;
t168 = t141 * t178;
t4 = t13 * t139 + t136 * t16;
t24 = pkin(9) * t112 + t159;
t52 = -pkin(3) * t135 - t61;
t87 = t113 * t137 - t135 * t140;
t88 = t113 * t140 + t135 * t137;
t25 = pkin(4) * t87 - pkin(9) * t88 + t52;
t162 = t136 * t25 + t139 * t24;
t161 = -t136 * t24 + t139 * t25;
t17 = -t137 * t32 + t140 * t50;
t34 = t132 * t82 + t134 * t81;
t109 = (-t132 * t138 + t202) * t197;
t65 = pkin(3) * t108 - pkin(8) * t109 + t228;
t160 = -t137 * t34 + t140 * t65;
t158 = -t137 * t53 + t140 * t68;
t157 = t112 * t139 - t136 * t88;
t67 = t112 * t136 + t139 * t88;
t155 = t102 * t140 - t104 * t195 + t107 * t207;
t153 = -t191 * t76 - t214;
t152 = t192 * t76 - t213;
t150 = t137 * t65 + t140 * t34 + t193 * t68 - t195 * t53;
t149 = -pkin(7) * t203 - t188;
t127 = pkin(2) * t132 + pkin(8);
t148 = -t127 * t102 + t104 * t31;
t147 = t149 * t135;
t12 = -pkin(4) * t104 - t17;
t146 = -pkin(9) * t42 + (t12 + t17) * t76;
t2 = -qJD(5) * t4 + t10 * t139 - t136 * t5;
t60 = -qJD(4) * t87 + t109 * t140;
t59 = qJD(4) * t88 + t109 * t137;
t36 = t45 * t195;
t23 = -pkin(4) * t112 - t158;
t22 = qJD(5) * t157 + t108 * t136 + t60 * t139;
t21 = qJD(5) * t67 - t108 * t139 + t60 * t136;
t19 = -pkin(4) * t110 + t137 * t55 - t140 * t64;
t11 = pkin(4) * t59 - pkin(9) * t60 + t33;
t8 = -t108 * pkin(4) + qJD(4) * t159 - t160;
t7 = pkin(9) * t108 + t150;
t3 = [0, 0, 0, 0.2e1 * t138 * t168, -0.2e1 * t199 * t178, t170 * t141 * t197, -t170 * t181, 0, t149 * qJD(2) ^ 2 + 0.2e1 * (-t129 * t225 + t147) * t189, -0.2e1 * pkin(1) * t168 - (-pkin(7) * t181 + t125) * t171 - (-pkin(7) * t167 + t122) * t135, -t102 * t62 - t103 * t61 + t107 * t34 - t108 * t39 - t109 * t38 + t110 * t33 - t112 * t29 + t113 * t28, -t28 * t61 + t29 * t62 - t38 * t33 + t39 * t34 + (t114 + t154) * t228, -t145 * t60 + t41 * t88, t145 * t59 - t41 * t87 - t42 * t88 - t60 * t77, t102 * t88 + t104 * t60 - t108 * t145 + t112 * t41, -t102 * t87 - t104 * t59 - t108 * t77 - t112 * t42, t102 * t112 + t104 * t108, t160 * t104 + t158 * t102 - t175 * t112 + t17 * t108 + t33 * t77 + t52 * t42 + t28 * t87 + t31 * t59 + (-t104 * t159 - t112 * t18) * qJD(4), -t102 * t159 - t104 * t150 - t18 * t108 - t112 * t151 - t145 * t33 + t28 * t88 + t31 * t60 + t52 * t41, t14 * t67 + t22 * t45, t14 * t157 - t15 * t67 - t21 * t45 - t22 * t43, t14 * t87 + t22 * t76 + t42 * t67 + t45 * t59, -t15 * t87 + t157 * t42 - t21 * t76 - t43 * t59, t42 * t87 + t59 * t76, (-qJD(5) * t162 + t139 * t11 - t136 * t7) * t76 + t161 * t42 + t2 * t87 - t163 * t59 + t8 * t43 + t23 * t15 - t6 * t157 + t12 * t21, -(qJD(5) * t161 + t136 * t11 + t139 * t7) * t76 - t162 * t42 - t1 * t87 - t4 * t59 + t8 * t45 + t23 * t14 + t6 * t67 + t12 * t22; 0, 0, 0, -t141 * t184, t199 * t205, -t141 * t180, t138 * t180, 0, pkin(1) * t184 - t142 * t147, t205 * t224 + (-pkin(7) * t182 + t124) * t190, (t39 - t54) * t110 + (t38 - t55) * t107 + (-t102 * t132 - t103 * t134) * pkin(2), t38 * t54 - t39 * t55 + (-t114 * t182 + t132 * t29 - t134 * t28) * pkin(2), t41 * t137 - t145 * t172, (t41 - t211) * t140 + (-t42 + t210) * t137, t104 * t172 + t201 + t209, t155 + t216, -t104 * t110, -t17 * t110 + t128 * t42 - t54 * t77 + (-t28 + (-qJD(4) * t127 - t64) * t104) * t140 + (t55 * t104 + t148) * t137, t18 * t110 + t128 * t41 + t28 * t137 + t54 * t145 + (t127 * t195 + t218) * t104 + t148 * t140, t14 * t139 * t137 + (-t137 * t192 + t139 * t193 - t70) * t45, t70 * t43 + t45 * t69 + (-t136 * t45 - t139 * t43) * t193 + (-t212 - t139 * t15 + (t136 * t43 - t139 * t45) * qJD(5)) * t137, -t220 + t36 + (-t14 + t185) * t140 + (-t152 - t217) * t137, t221 + (t15 - t186) * t140 + (t153 - t227) * t137, -t42 * t140 + t207 * t76, t115 * t213 - t12 * t69 - t19 * t43 + (t136 * t173 + t139 * t208) * t76 + (t12 * t196 - t2 + (qJD(4) * t43 + t153) * t127) * t140 + (t12 * t191 + t163 * t107 + t127 * t15 + t6 * t136 + (t127 * t136 * t76 - t163) * qJD(4)) * t137, -t115 * t214 - t12 * t70 - t19 * t45 + (-t136 * t208 + t139 * t173) * t76 + (t12 * t194 + t1 + (qJD(4) * t45 + t152) * t127) * t140 + (-t12 * t192 + t4 * t107 + t127 * t14 + t6 * t139 + (t127 * t174 - t4) * qJD(4)) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 ^ 2 - t110 ^ 2, -t107 * t39 + t110 * t38 + t121, 0, 0, 0, 0, 0, t155 - t216, -t104 ^ 2 * t140 - t201 + t209, 0, 0, 0, 0, 0, t221 + (-t15 - t186) * t140 + (t153 + t227) * t137, t220 + t36 + (-t14 - t185) * t140 + (t152 - t217) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145 * t77, t145 ^ 2 - t77 ^ 2, t41 + t211, -t42 - t210, t102, t31 * t145 - t175 + (-qJD(4) + t104) * t18, t104 * t17 + t31 * t77 - t151, t174 * t45 + t212, (t14 - t223) * t139 + (-t15 - t222) * t136, t145 * t45 + t174 * t76 + t214, -t136 * t76 ^ 2 - t145 * t43 + t213, t76 * t145, -pkin(4) * t15 + t146 * t136 - t139 * t226 - t145 * t163 - t18 * t43, -pkin(4) * t14 + t136 * t226 + t146 * t139 - t145 * t4 - t18 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t14 + t223, -t15 + t222, t42, -t12 * t45 + t4 * t76 + t2, t12 * t43 - t163 * t76 - t1;];
tauc_reg = t3;
