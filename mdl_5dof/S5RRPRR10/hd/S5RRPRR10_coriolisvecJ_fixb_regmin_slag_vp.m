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
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 22:02:43
% EndTime: 2021-01-15 22:02:59
% DurationCPUTime: 2.76s
% Computational Cost: add. (4158->359), mult. (12937->519), div. (0->0), fcn. (10383->10), ass. (0->182)
t133 = sin(pkin(10));
t135 = cos(pkin(5));
t138 = sin(qJ(2));
t232 = pkin(1) * t138;
t193 = t135 * t232;
t134 = sin(pkin(5));
t141 = cos(qJ(2));
t210 = t134 * t141;
t226 = pkin(7) + qJ(3);
t104 = t226 * t210 + t193;
t202 = qJD(3) * t138;
t143 = -t104 * qJD(2) - t134 * t202;
t235 = t143 * t133;
t215 = cos(pkin(10));
t151 = t133 * t141 + t215 * t138;
t205 = qJD(1) * t134;
t109 = t151 * t205;
t140 = cos(qJ(4));
t195 = t135 * qJD(1);
t171 = qJD(2) + t195;
t121 = t140 * t171;
t137 = sin(qJ(4));
t77 = t137 * t109 - t121;
t76 = qJD(5) + t77;
t177 = t215 * t141;
t166 = t134 * t177;
t120 = qJD(1) * t166;
t188 = t138 * t205;
t106 = t133 * t188 - t120;
t103 = qJD(4) + t106;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t146 = -t140 * t109 - t137 * t171;
t43 = -t139 * t103 - t136 * t146;
t234 = t103 * t43;
t113 = t151 * t134;
t108 = qJD(2) * t113;
t101 = qJD(1) * t108;
t231 = pkin(1) * t141;
t192 = t135 * t231;
t125 = qJD(1) * t192;
t123 = t125 * qJD(2);
t184 = t138 * t226;
t145 = (-qJD(2) * t184 + qJD(3) * t141) * t134;
t73 = qJD(1) * t145 + t123;
t183 = t215 * t73;
t29 = qJD(1) * t235 + t183;
t204 = qJD(2) * t138;
t187 = t134 * t204;
t167 = qJD(1) * t187;
t102 = qJD(2) * t120 - t133 * t167;
t122 = pkin(2) * t167;
t51 = t101 * pkin(3) - t102 * pkin(8) + t122;
t178 = t137 * t29 - t140 * t51;
t95 = t104 * qJD(1);
t182 = t215 * t95;
t169 = t134 * t184;
t80 = qJD(2) * pkin(2) + t125 + (t135 * pkin(2) - t169) * qJD(1);
t39 = t133 * t80 + t182;
t32 = t171 * pkin(8) + t39;
t116 = (-pkin(2) * t141 - pkin(1)) * t134;
t206 = qJD(1) * t116;
t114 = qJD(3) + t206;
t50 = t106 * pkin(3) - t109 * pkin(8) + t114;
t18 = t137 * t50 + t140 * t32;
t6 = -t101 * pkin(4) + t18 * qJD(4) + t178;
t233 = (-t146 * pkin(4) + t76 * pkin(9)) * t76 + t6;
t200 = qJD(4) * t137;
t41 = qJD(4) * t121 + t140 * t102 - t109 * t200;
t45 = t136 * t103 - t139 * t146;
t15 = t45 * qJD(5) - t139 * t101 + t136 * t41;
t42 = -t146 * qJD(4) + t137 * t102;
t75 = t215 * t143;
t28 = -qJD(1) * t75 + t133 * t73;
t10 = t42 * pkin(4) - t41 * pkin(9) + t28;
t13 = t103 * pkin(9) + t18;
t83 = t133 * t95;
t38 = t215 * t80 - t83;
t31 = -t171 * pkin(3) - t38;
t16 = t77 * pkin(4) + pkin(9) * t146 + t31;
t163 = t136 * t13 - t139 * t16;
t198 = qJD(4) * t140;
t153 = t137 * t51 + t140 * t29 + t50 * t198 - t32 * t200;
t5 = t101 * pkin(9) + t153;
t1 = -t163 * qJD(5) + t136 * t10 + t139 * t5;
t142 = qJD(1) ^ 2;
t230 = t43 * t76;
t229 = t45 * t76;
t213 = t106 * t140;
t69 = -t139 * t109 - t136 * t213;
t228 = t69 * t76;
t70 = t136 * t109 - t139 * t213;
t227 = t70 * t76;
t93 = (pkin(2) + t231) * t135 - t169;
t62 = t215 * t104 + t133 * t93;
t53 = t135 * pkin(8) + t62;
t211 = t133 * t138;
t112 = t134 * t211 - t166;
t68 = t112 * pkin(3) - t113 * pkin(8) + t116;
t159 = t137 * t68 + t140 * t53;
t94 = -qJD(1) * t169 + t125;
t55 = t215 * t94 - t83;
t172 = pkin(2) * t188;
t64 = t109 * pkin(3) + t106 * pkin(8) + t172;
t225 = t137 * t64 + t140 * t55;
t224 = t106 * t45;
t223 = t109 * t77;
t222 = t136 * t42;
t221 = t139 * t42;
t176 = t139 * t76;
t196 = qJD(5) * t139;
t197 = qJD(5) * t136;
t14 = t136 * t101 + t103 * t196 + t139 * t41 + t146 * t197;
t220 = t14 * t136;
t219 = t77 * t103;
t218 = t146 * t103;
t217 = t146 * t109;
t54 = t133 * t94 + t182;
t216 = -t54 + t103 * (pkin(4) * t137 - pkin(9) * t140);
t214 = t103 * t137;
t130 = t134 ^ 2;
t212 = t130 * t142;
t209 = t137 * t101;
t207 = t138 ^ 2 - t141 ^ 2;
t203 = qJD(2) * t141;
t201 = qJD(4) * t136;
t199 = qJD(4) * t139;
t194 = qJD(1) * qJD(2);
t191 = t76 * t201;
t190 = t76 * t199;
t189 = t138 * t212;
t186 = t134 * t135 * t142;
t181 = t130 * t194;
t126 = qJD(2) * t192;
t81 = t126 + t145;
t33 = t133 * t81 - t75;
t129 = -t215 * pkin(2) - pkin(3);
t115 = -t140 * pkin(4) - t137 * pkin(9) + t129;
t175 = t109 * pkin(9) - qJD(5) * t115 + t225;
t174 = t103 * t140;
t173 = pkin(2) * t187;
t170 = qJD(2) + 0.2e1 * t195;
t168 = t141 * t181;
t61 = -t133 * t104 + t215 * t93;
t4 = t139 * t13 + t136 * t16;
t24 = t112 * pkin(9) + t159;
t52 = -t135 * pkin(3) - t61;
t86 = t137 * t113 - t135 * t140;
t87 = t140 * t113 + t135 * t137;
t25 = t86 * pkin(4) - t87 * pkin(9) + t52;
t162 = t136 * t25 + t139 * t24;
t161 = -t136 * t24 + t139 * t25;
t17 = -t137 * t32 + t140 * t50;
t34 = t215 * t81 + t235;
t111 = (t177 - t211) * t134 * qJD(2);
t65 = t108 * pkin(3) - t111 * pkin(8) + t173;
t160 = -t137 * t34 + t140 * t65;
t158 = -t137 * t53 + t140 * t68;
t157 = t139 * t112 - t136 * t87;
t67 = t136 * t112 + t139 * t87;
t156 = t140 * t101 - t103 * t200 - t106 * t214;
t155 = -t76 * t196 - t222;
t154 = t76 * t197 - t221;
t152 = t137 * t65 + t140 * t34 + t68 * t198 - t53 * t200;
t150 = -pkin(7) * t210 - t193;
t128 = t133 * pkin(2) + pkin(8);
t149 = -t128 * t101 + t103 * t31;
t148 = t150 * t135;
t12 = -t103 * pkin(4) - t17;
t147 = -pkin(9) * t42 + (t12 + t17) * t76;
t2 = -t4 * qJD(5) + t139 * t10 - t136 * t5;
t60 = t87 * qJD(4) + t137 * t111;
t59 = -t86 * qJD(4) + t140 * t111;
t36 = t45 * t200;
t23 = -t112 * pkin(4) - t158;
t22 = t67 * qJD(5) - t139 * t108 + t136 * t59;
t21 = t157 * qJD(5) + t136 * t108 + t139 * t59;
t19 = -t109 * pkin(4) + t137 * t55 - t140 * t64;
t11 = t60 * pkin(4) - t59 * pkin(9) + t33;
t8 = -t108 * pkin(4) + t159 * qJD(4) - t160;
t7 = t108 * pkin(9) + t152;
t3 = [0, 0, 0, 0.2e1 * t138 * t168, -0.2e1 * t207 * t181, t170 * t134 * t203, -t170 * t187, 0, t150 * qJD(2) ^ 2 + 0.2e1 * (-t130 * t232 + t148) * t194, -0.2e1 * pkin(1) * t168 - (-pkin(7) * t187 + t126) * t171 - (-pkin(7) * t167 + t123) * t135, -t33 * t171 - t28 * t135 + t116 * t101 + t114 * t108 + (qJD(1) * t112 + t106) * t173, -t34 * t171 - t29 * t135 + t116 * t102 + t114 * t111 + (qJD(1) * t113 + t109) * t173, -t62 * t101 - t61 * t102 - t34 * t106 - t39 * t108 + t33 * t109 - t38 * t111 - t29 * t112 + t28 * t113, -t28 * t61 + t29 * t62 - t38 * t33 + t39 * t34 + (t114 + t206) * t173, -t146 * t59 + t41 * t87, t146 * t60 - t41 * t86 - t87 * t42 - t59 * t77, t87 * t101 + t59 * t103 - t108 * t146 + t41 * t112, -t86 * t101 - t60 * t103 - t77 * t108 - t42 * t112, t101 * t112 + t103 * t108, t160 * t103 + t158 * t101 - t178 * t112 + t17 * t108 + t33 * t77 + t52 * t42 + t28 * t86 + t31 * t60 + (-t103 * t159 - t112 * t18) * qJD(4), -t159 * t101 - t152 * t103 - t18 * t108 - t153 * t112 - t146 * t33 + t28 * t87 + t31 * t59 + t52 * t41, t14 * t67 + t45 * t21, t14 * t157 - t67 * t15 - t21 * t43 - t45 * t22, t14 * t86 + t21 * t76 + t67 * t42 + t45 * t60, -t15 * t86 + t157 * t42 - t22 * t76 - t43 * t60, t42 * t86 + t76 * t60, (-qJD(5) * t162 + t139 * t11 - t136 * t7) * t76 + t161 * t42 + t2 * t86 - t163 * t60 + t8 * t43 + t23 * t15 - t6 * t157 + t12 * t22, -(qJD(5) * t161 + t136 * t11 + t139 * t7) * t76 - t162 * t42 - t1 * t86 - t4 * t60 + t8 * t45 + t23 * t14 + t6 * t67 + t12 * t21; 0, 0, 0, -t141 * t189, t207 * t212, -t141 * t186, t138 * t186, 0, pkin(1) * t189 - t142 * t148, t212 * t231 + (-pkin(7) * t188 + t125) * t195, -t106 * t172 - t114 * t109 + t54 * t171 - t28, -t183 + t55 * qJD(2) + t114 * t106 + ((t133 * pkin(1) * t204 + t55) * t135 + (-t133 * (-t226 * t203 - t202) - t138 * pkin(2) * t109) * t134) * qJD(1), (t39 - t54) * t109 + (-t38 + t55) * t106 + (-t101 * t133 - t215 * t102) * pkin(2), t38 * t54 - t39 * t55 + (-t114 * t188 + t133 * t29 - t215 * t28) * pkin(2), t41 * t137 - t146 * t174, (t41 - t219) * t140 + (-t42 + t218) * t137, t103 * t174 + t209 + t217, t156 + t223, -t103 * t109, -t17 * t109 + t129 * t42 - t54 * t77 + (-t28 + (-qJD(4) * t128 - t64) * t103) * t140 + (t55 * t103 + t149) * t137, t18 * t109 + t129 * t41 + t28 * t137 + t54 * t146 + (t128 * t200 + t225) * t103 + t149 * t140, t14 * t139 * t137 + (-t137 * t197 + t139 * t198 - t70) * t45, t70 * t43 + t45 * t69 + (-t136 * t45 - t139 * t43) * t198 + (-t220 - t139 * t15 + (t136 * t43 - t139 * t45) * qJD(5)) * t137, -t227 + t36 + (-t14 + t190) * t140 + (-t154 + t224) * t137, t228 + (t15 - t191) * t140 + (t155 - t234) * t137, -t42 * t140 + t76 * t214, t115 * t221 - t12 * t69 - t19 * t43 + (t136 * t175 + t139 * t216) * t76 + (t12 * t201 - t2 + (qJD(4) * t43 + t155) * t128) * t140 + (t12 * t196 - t163 * t106 + t128 * t15 + t6 * t136 + (t128 * t136 * t76 - t163) * qJD(4)) * t137, -t115 * t222 - t12 * t70 - t19 * t45 + (-t136 * t216 + t139 * t175) * t76 + (t12 * t199 + t1 + (qJD(4) * t45 + t154) * t128) * t140 + (-t12 * t197 - t4 * t106 + t128 * t14 + t6 * t139 + (t128 * t176 - t4) * qJD(4)) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * qJD(2) + (t109 * t135 + t108) * qJD(1), -t106 * t171 + t102, -t106 ^ 2 - t109 ^ 2, t39 * t106 + t38 * t109 + t122, 0, 0, 0, 0, 0, t156 - t223, -t103 ^ 2 * t140 - t209 + t217, 0, 0, 0, 0, 0, t228 + (-t15 - t191) * t140 + (t155 + t234) * t137, t227 + t36 + (-t14 - t190) * t140 + (t154 + t224) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146 * t77, t146 ^ 2 - t77 ^ 2, t41 + t219, -t42 - t218, t101, t31 * t146 - t178 + (-qJD(4) + t103) * t18, t17 * t103 + t31 * t77 - t153, t45 * t176 + t220, (t14 - t230) * t139 + (-t15 - t229) * t136, t146 * t45 + t76 * t176 + t222, -t76 ^ 2 * t136 - t146 * t43 + t221, t76 * t146, -pkin(4) * t15 + t147 * t136 - t233 * t139 - t146 * t163 - t18 * t43, -pkin(4) * t14 + t233 * t136 + t147 * t139 - t146 * t4 - t18 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t14 + t230, -t15 + t229, t42, -t12 * t45 + t4 * t76 + t2, t12 * t43 - t163 * t76 - t1;];
tauc_reg = t3;
