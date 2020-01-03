% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:46
% DurationCPUTime: 2.77s
% Computational Cost: add. (4018->342), mult. (9109->460), div. (0->0), fcn. (5686->6), ass. (0->189)
t122 = cos(qJ(4));
t123 = cos(qJ(2));
t192 = t123 * t122;
t119 = sin(qJ(4));
t120 = sin(qJ(2));
t193 = t120 * t119;
t67 = t192 + t193;
t57 = t67 * qJD(1);
t234 = qJD(5) + t57;
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t175 = qJD(2) - qJD(4);
t155 = t121 * t175;
t185 = qJD(1) * t123;
t186 = qJD(1) * t120;
t60 = -t119 * t185 + t122 * t186;
t37 = t118 * t60 + t155;
t239 = t234 * t37;
t39 = -t118 * t175 + t121 * t60;
t238 = t234 * t39;
t177 = qJD(1) * qJD(2);
t168 = t123 * t177;
t96 = pkin(6) * t168;
t140 = -pkin(7) * t168 + t96;
t181 = qJD(4) * t122;
t182 = qJD(4) * t119;
t124 = -pkin(2) - pkin(3);
t170 = t124 * qJD(2);
t103 = pkin(6) * t186;
t73 = pkin(7) * t186 - t103;
t235 = qJD(3) - t73;
t47 = t170 + t235;
t113 = qJD(2) * qJD(3);
t184 = qJD(2) * t120;
t221 = pkin(6) - pkin(7);
t74 = t221 * t184;
t50 = -qJD(1) * t74 + t113;
t114 = qJD(2) * qJ(3);
t104 = pkin(6) * t185;
t149 = -pkin(7) * t185 + t104;
t61 = t114 + t149;
t136 = -t119 * t140 - t122 * t50 - t47 * t181 + t61 * t182;
t62 = -qJD(1) * pkin(1) - pkin(2) * t185 - qJ(3) * t186;
t46 = pkin(3) * t185 - t62;
t22 = t57 * pkin(4) - t60 * pkin(8) + t46;
t28 = t119 * t47 + t122 * t61;
t25 = -t175 * pkin(8) + t28;
t6 = t118 * t22 + t121 * t25;
t107 = t120 * qJD(3);
t207 = qJ(3) * t168 + qJD(1) * t107;
t176 = qJD(1) * qJD(4);
t167 = t123 * t176;
t208 = t122 * t167 + t176 * t193;
t183 = qJD(2) * t123;
t134 = t119 * t183 + t120 * t181;
t169 = t120 * t177;
t30 = t134 * qJD(1) - t119 * t167 - t122 * t169;
t8 = t208 * pkin(8) + t30 * pkin(4) + (-pkin(8) * t192 + (-t119 * pkin(8) + t124) * t120) * t177 + t207;
t2 = -qJD(5) * t6 + t118 * t136 + t121 * t8;
t237 = -t234 * t6 - t2;
t141 = t118 * t25 - t121 * t22;
t1 = -t141 * qJD(5) + t118 * t8 - t121 * t136;
t236 = t234 * t141 + t1;
t201 = t121 * t30;
t226 = t118 * t234;
t233 = t226 * t234 - t37 * t60 - t201;
t205 = t118 * t30;
t225 = t121 * t234;
t232 = t225 * t234 - t39 * t60 + t205;
t135 = t67 * qJD(2);
t130 = -qJD(1) * t135 + t208;
t195 = qJD(5) * t39;
t19 = -t118 * t130 + t195;
t197 = t19 * t121;
t231 = t226 * t37 - t197;
t180 = qJD(5) * t118;
t18 = qJD(5) * t155 + t121 * t130 + t60 * t180;
t198 = t18 * t118;
t230 = t225 * t39 - t198;
t229 = (t18 + t239) * t121 + (t19 + t238) * t118;
t33 = t119 * t149 + t122 * t73;
t76 = -t119 * qJ(3) + t122 * t124;
t48 = t122 * qJD(3) + t76 * qJD(4);
t211 = t48 - t33;
t27 = -t119 * t61 + t122 * t47;
t24 = t175 * pkin(4) - t27;
t228 = t234 * t24;
t227 = t60 * t175 + t30;
t77 = t122 * qJ(3) + t119 * t124;
t222 = t175 ^ 2;
t220 = t141 * t60;
t219 = t6 * t60;
t166 = t119 * t50 - t122 * t140;
t13 = t28 * qJD(4) + t166;
t84 = t221 * t120;
t85 = t221 * t123;
t40 = t119 * t85 - t122 * t84;
t218 = t13 * t40;
t68 = -t123 * t119 + t120 * t122;
t217 = t13 * t68;
t216 = t30 * t67;
t215 = t39 * t37;
t214 = t234 * t60;
t213 = t60 * t57;
t212 = t68 * t30;
t210 = t77 * qJD(4) + t235 * t119 + t122 * t149;
t206 = qJD(2) * pkin(2);
t204 = t118 * t37;
t203 = t118 * t39;
t202 = t119 * t60;
t200 = t121 * t37;
t199 = t121 * t39;
t196 = qJD(4) * t234;
t194 = qJD(5) * t234;
t126 = qJD(1) ^ 2;
t191 = t123 * t126;
t125 = qJD(2) ^ 2;
t190 = t125 * t120;
t108 = t125 * t123;
t188 = qJ(3) * t183 + t107;
t115 = t120 ^ 2;
t187 = t123 ^ 2 - t115;
t179 = qJD(5) * t121;
t178 = t122 * qJD(2);
t174 = -t28 + t210;
t79 = -t123 * pkin(2) - t120 * qJ(3) - pkin(1);
t173 = t57 ^ 2 - t60 ^ 2;
t172 = t68 * t180;
t171 = t68 * t179;
t165 = qJD(1) * t79 + t62;
t160 = -0.2e1 * pkin(1) * t177;
t159 = qJD(3) - t206;
t158 = t27 * t175;
t157 = t57 * t175;
t65 = t123 * pkin(3) - t79;
t154 = qJD(2) * t85;
t153 = t120 * t170;
t152 = pkin(8) * t194 + t13;
t151 = t120 * t168;
t31 = t60 * pkin(4) + t57 * pkin(8);
t72 = -pkin(8) + t77;
t150 = t72 * t194 - t13;
t148 = -t118 * t6 + t121 * t141;
t147 = -t118 * t141 - t121 * t6;
t35 = -t67 * qJD(4) + t135;
t146 = t24 * t35 + t217;
t145 = t234 * t35 + t212;
t26 = t67 * pkin(4) - t68 * pkin(8) + t65;
t41 = t119 * t84 + t122 * t85;
t16 = -t118 * t41 + t121 * t26;
t17 = t118 * t26 + t121 * t41;
t139 = -t46 * t60 - t166;
t100 = qJ(3) * t185;
t53 = t124 * t186 + t100;
t45 = pkin(2) * t169 - t207;
t54 = pkin(2) * t184 - t188;
t138 = -pkin(6) * t125 - qJD(1) * t54 - t45;
t137 = -pkin(8) * t30 + t228;
t43 = t153 + t188;
t133 = t46 * t57 + t136;
t132 = -t234 * t48 - t30 * t72 - t228;
t131 = -qJD(2) * t57 + t208;
t129 = t148 * qJD(5) + t1 * t121 - t2 * t118;
t75 = -pkin(6) * t169 + t113;
t78 = t103 + t159;
t81 = t104 + t114;
t127 = t75 * t123 + (t123 * t78 + (-t81 + t104) * t120) * qJD(2);
t93 = t120 * t191;
t83 = -0.2e1 * t151;
t82 = 0.2e1 * t151;
t80 = t187 * t126;
t71 = pkin(4) - t76;
t70 = pkin(2) * t186 - t100;
t66 = t187 * t177;
t59 = t118 * t186 + t121 * t178;
t56 = -t118 * t178 + t121 * t186;
t36 = qJD(1) * t153 + t207;
t34 = -t120 * t178 - t123 * t182 + t134;
t23 = -t31 + t53;
t21 = t41 * qJD(4) - t119 * t74 - t122 * t154;
t20 = -t40 * qJD(4) + t119 * t154 - t122 * t74;
t15 = t118 * t31 + t121 * t27;
t14 = -t118 * t27 + t121 * t31;
t11 = t34 * pkin(4) - t35 * pkin(8) + t43;
t10 = t118 * t23 + t121 * t33;
t9 = -t118 * t33 + t121 * t23;
t4 = -t17 * qJD(5) + t121 * t11 - t118 * t20;
t3 = t16 * qJD(5) + t118 * t11 + t121 * t20;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0.2e1 * t66, t108, t83, -t190, 0, -pkin(6) * t108 + t120 * t160, pkin(6) * t190 + t123 * t160, 0, 0, t82, t108, -0.2e1 * t66, 0, t190, t83, t138 * t123 + t165 * t184, t127, t138 * t120 - t165 * t183, pkin(6) * t127 + t45 * t79 + t62 * t54, -t130 * t68 + t60 * t35, t130 * t67 - t60 * t34 - t35 * t57 - t212, -t35 * t175, t57 * t34 + t216, t34 * t175, 0, t21 * t175 + t65 * t30 + t46 * t34 + t36 * t67 + t43 * t57, -t65 * t130 + t20 * t175 + t46 * t35 + t36 * t68 + t43 * t60, -t131 * t40 + t136 * t67 - t20 * t57 + t21 * t60 - t27 * t35 - t28 * t34 - t41 * t30 + t217, -t136 * t41 + t28 * t20 - t27 * t21 + t36 * t65 + t46 * t43 + t218, -t39 * t172 + (-t18 * t68 + t35 * t39) * t121, (-t200 - t203) * t35 + (t198 - t197 + (-t199 + t204) * qJD(5)) * t68, t121 * t145 - t172 * t234 - t18 * t67 + t39 * t34, t37 * t171 + (t19 * t68 + t35 * t37) * t118, -t118 * t145 - t171 * t234 - t19 * t67 - t37 * t34, t234 * t34 + t216, t118 * t146 - t141 * t34 + t16 * t30 + t171 * t24 + t40 * t19 + t2 * t67 + t21 * t37 + t234 * t4, -t1 * t67 + t121 * t146 - t17 * t30 - t172 * t24 - t40 * t18 + t21 * t39 - t234 * t3 - t6 * t34, t16 * t18 - t17 * t19 - t3 * t37 - t4 * t39 + t148 * t35 + (qJD(5) * t147 - t1 * t118 - t121 * t2) * t68, t1 * t17 - t141 * t4 + t2 * t16 + t24 * t21 + t6 * t3 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t80, 0, t93, 0, 0, t126 * pkin(1) * t120, pkin(1) * t191, 0, 0, -t93, 0, t80, 0, 0, t93, (-t120 * t62 + t123 * t70) * qJD(1), ((t81 - t114) * t120 + (t159 - t78) * t123) * qJD(1), 0.2e1 * t113 + (t120 * t70 + t123 * t62) * qJD(1), t75 * qJ(3) + t81 * qJD(3) - t62 * t70 + (t120 * t81 + (-t78 - t206) * t123) * qJD(1) * pkin(6), -t213, t173, t157 + t131, t213, t227, 0, t210 * qJD(2) - t174 * qJD(4) - t53 * t57 - t139, t211 * t175 - t53 * t60 - t133, -t77 * t30 + t76 * t131 + t174 * t60 + (t27 - t211) * t57, -t13 * t76 - t136 * t77 - t210 * t27 + t211 * t28 - t46 * t53, -t230, t229, -t232, -t231, t233, t214, t118 * t132 - t121 * t150 + t71 * t19 + t210 * t37 - t234 * t9 - t220, t10 * t234 + t118 * t150 + t121 * t132 - t71 * t18 + t210 * t39 - t219, t10 * t37 + t9 * t39 + (-t19 * t72 - t37 * t48 - t141 * t57 - t1 + (t39 * t72 - t141) * qJD(5)) * t121 + (-t18 * t72 + t39 * t48 + t57 * t6 + t2 + (t37 * t72 + t6) * qJD(5)) * t118, -t6 * t10 + t129 * t72 + t13 * t71 + t141 * t9 - t147 * t48 + t210 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, 0, -t115 * t126 - t125, -t81 * qJD(2) + t62 * t186 + t96, 0, 0, 0, 0, 0, 0, -t119 * t222 - t57 * t186, -t122 * t222 - t60 * t186, -t119 * t30 + t122 * t208 + (-t122 * t57 + t202) * qJD(4) - qJD(2) * t202, -t46 * t186 + (-t175 * t28 - t13) * t122 + (-t136 + t158) * t119, 0, 0, 0, 0, 0, 0, -t56 * t234 + (-t118 * t196 - t19) * t122 + (-t175 * t37 - t179 * t234 - t205) * t119, t59 * t234 + (-t121 * t196 + t18) * t122 + (-t175 * t39 + t180 * t234 - t201) * t119, t59 * t37 + t56 * t39 + (-t200 + t203) * t181 + (-t198 - t197 + (t199 + t204) * qJD(5)) * t119, t141 * t56 - t6 * t59 + (-qJD(4) * t147 - t13) * t122 + (-t175 * t24 + t129) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t173, -t157 - t130, -t213, -t227, 0, -t28 * qJD(2) + t139, -t158 + t133, 0, 0, t230, -t229, t232, t231, -t233, -t214, -pkin(4) * t19 + t118 * t137 - t121 * t152 - t14 * t234 - t28 * t37 + t220, pkin(4) * t18 + t118 * t152 + t121 * t137 + t15 * t234 - t28 * t39 + t219, t14 * t39 + t15 * t37 + ((-t19 + t195) * pkin(8) + t236) * t121 + ((qJD(5) * t37 - t18) * pkin(8) + t237) * t118, -t13 * pkin(4) + pkin(8) * t129 + t14 * t141 - t6 * t15 - t24 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t37 ^ 2 + t39 ^ 2, -t18 + t239, -t215, -t19 + t238, t30, -t24 * t39 - t237, t24 * t37 - t236, 0, 0;];
tauc_reg = t5;
