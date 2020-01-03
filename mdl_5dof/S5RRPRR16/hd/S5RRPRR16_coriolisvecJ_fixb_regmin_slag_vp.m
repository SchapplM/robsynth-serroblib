% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR16_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:34
% DurationCPUTime: 2.63s
% Computational Cost: add. (2629->348), mult. (7046->507), div. (0->0), fcn. (5186->8), ass. (0->178)
t139 = cos(qJ(2));
t132 = sin(pkin(5));
t136 = sin(qJ(2));
t206 = qJD(1) * t136;
t189 = t132 * t206;
t133 = cos(pkin(5));
t207 = qJD(1) * t133;
t195 = pkin(1) * t207;
t83 = pkin(7) * t189 - t139 * t195;
t211 = qJD(3) + t83;
t121 = qJD(2) + t207;
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t208 = qJD(1) * t132;
t188 = t139 * t208;
t75 = t121 * t135 + t138 * t188;
t74 = qJD(5) + t75;
t108 = qJD(4) + t189;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t172 = t135 * t188;
t77 = t121 * t138 - t172;
t34 = -t108 * t137 + t134 * t77;
t235 = t108 * t34;
t230 = pkin(1) * t136;
t126 = t133 * t230;
t214 = t132 * t139;
t234 = pkin(7) * t214 + t126;
t212 = pkin(3) * t189 + t211;
t196 = qJD(1) * qJD(2);
t183 = t132 * t196;
t110 = t139 * t183;
t140 = -pkin(2) - pkin(8);
t29 = t121 * t140 + t212;
t182 = -qJ(3) * t136 - pkin(1);
t63 = (t139 * t140 + t182) * t132;
t46 = qJD(1) * t63;
t17 = t135 * t29 + t138 * t46;
t171 = t136 * t183;
t103 = pkin(2) * t171;
t166 = pkin(8) * t136 - qJ(3) * t139;
t203 = qJD(3) * t136;
t142 = (qJD(2) * t166 - t203) * t132;
t31 = qJD(1) * t142 + t103;
t194 = pkin(1) * qJD(2) * t133;
t174 = qJD(1) * t194;
t78 = pkin(7) * t110 + t136 * t174;
t50 = pkin(3) * t110 + t78;
t143 = -qJD(4) * t17 - t135 * t31 + t138 * t50;
t6 = -pkin(4) * t110 - t143;
t233 = (pkin(4) * t77 + pkin(9) * t74) * t74 + t6;
t36 = t108 * t134 + t137 * t77;
t38 = -qJD(4) * t75 + t135 * t171;
t15 = qJD(5) * t36 - t110 * t137 + t134 * t38;
t215 = t132 * t136;
t123 = pkin(7) * t215;
t190 = -pkin(1) * t139 - pkin(2);
t49 = pkin(3) * t215 + t123 + (-pkin(8) + t190) * t133;
t159 = t135 * t49 + t138 * t63;
t205 = qJD(2) * t136;
t187 = t132 * t205;
t115 = pkin(2) * t187;
t43 = t115 + t142;
t231 = pkin(3) + pkin(7);
t68 = (t214 * t231 + t126) * qJD(2);
t232 = -qJD(4) * t159 - t135 * t43 + t138 * t68;
t210 = pkin(7) * t171 - t139 * t174;
t54 = -t121 * qJD(3) + t210;
t32 = -pkin(3) * t171 - t54;
t199 = qJD(4) * t138;
t39 = -qJD(4) * t172 + t121 * t199 - t138 * t171;
t10 = pkin(4) * t39 - pkin(9) * t38 + t32;
t12 = pkin(9) * t108 + t17;
t111 = t121 * qJ(3);
t84 = pkin(7) * t188 + t136 * t195;
t67 = pkin(3) * t188 + t84;
t41 = t111 + t67;
t18 = pkin(4) * t75 - pkin(9) * t77 + t41;
t165 = t12 * t134 - t137 * t18;
t201 = qJD(4) * t135;
t151 = t135 * t50 + t138 * t31 + t199 * t29 - t201 * t46;
t5 = pkin(9) * t110 + t151;
t1 = -qJD(5) * t165 + t10 * t134 + t137 * t5;
t229 = t34 * t74;
t228 = t36 * t74;
t168 = pkin(4) * t138 + pkin(9) * t135;
t227 = (-pkin(3) - t168) * t189 - qJD(4) * t168 - t211;
t120 = pkin(2) * t189;
t65 = t166 * t208 + t120;
t226 = t135 * t67 + t138 * t65;
t225 = t108 * t75;
t224 = t108 * t77;
t197 = qJD(5) * t137;
t198 = qJD(5) * t134;
t14 = t108 * t197 + t110 * t134 + t137 * t38 - t198 * t77;
t223 = t134 * t14;
t222 = t134 * t39;
t221 = t137 * t39;
t177 = t137 * t74;
t220 = t140 * t74;
t219 = t108 * t135;
t218 = t108 * t138;
t217 = t108 * t140;
t129 = t132 ^ 2;
t216 = t129 * qJD(1) ^ 2;
t213 = t136 * t137;
t130 = t136 ^ 2;
t209 = -t139 ^ 2 + t130;
t204 = qJD(2) * t139;
t202 = qJD(4) * t134;
t200 = qJD(4) * t137;
t193 = t136 * t218;
t192 = t139 * t216;
t191 = t135 * t214;
t79 = -qJ(3) * t133 - t234;
t186 = t132 * t204;
t185 = t108 * t199;
t184 = t129 * t196;
t97 = pkin(4) * t135 - pkin(9) * t138 + qJ(3);
t178 = pkin(9) * t188 - qJD(5) * t97 + t226;
t176 = t121 + t207;
t175 = 0.2e1 * t184;
t173 = t136 * t192;
t62 = pkin(3) * t214 - t79;
t170 = t121 * t84 - t78;
t169 = -0.2e1 * pkin(1) * t184;
t4 = t12 * t137 + t134 * t18;
t85 = t234 * qJD(2);
t164 = t121 * t85 + t133 * t78;
t22 = pkin(9) * t215 + t159;
t86 = t133 * t135 + t138 * t214;
t87 = t133 * t138 - t191;
t25 = pkin(4) * t86 - pkin(9) * t87 + t62;
t163 = t134 * t25 + t137 * t22;
t162 = -t134 * t22 + t137 * t25;
t16 = -t135 * t46 + t138 * t29;
t160 = -t135 * t63 + t138 * t49;
t158 = -t135 * t65 + t138 * t67;
t119 = t139 * t194;
t157 = -pkin(7) * t187 + t119;
t156 = -t121 * t188 + t110;
t155 = -t197 * t74 - t222;
t154 = t198 * t74 - t221;
t153 = t132 * t213 - t134 * t87;
t60 = t134 * t215 + t137 * t87;
t152 = t136 * t41 + t140 * t204;
t150 = t135 * t68 + t138 * t43 + t199 * t49 - t201 * t63;
t149 = t108 * t74;
t148 = t108 * t36;
t80 = (-pkin(2) * t139 + t182) * t132;
t127 = t133 * qJD(3);
t48 = -t187 * t231 + t119 + t127;
t145 = (-qJ(3) * t204 - t203) * t132;
t11 = -pkin(4) * t108 - t16;
t144 = -pkin(9) * t39 + (t11 + t16) * t74;
t2 = -qJD(5) * t4 + t10 * t137 - t134 * t5;
t95 = t138 * t110;
t82 = -qJ(3) * t188 + t120;
t81 = t133 * t190 + t123;
t73 = -t127 - t157;
t72 = (t134 * t139 + t135 * t213) * t208;
t71 = t134 * t135 * t189 - t137 * t188;
t70 = qJD(1) * t80;
t69 = t115 + t145;
t64 = -t111 - t84;
t61 = -pkin(2) * t121 + t211;
t58 = -qJD(4) * t191 + t133 * t199 - t138 * t187;
t57 = -qJD(4) * t86 + t135 * t187;
t52 = qJD(1) * t145 + t103;
t47 = t70 * t189;
t23 = -pkin(4) * t188 - t158;
t21 = -pkin(4) * t215 - t160;
t20 = qJD(5) * t153 + t134 * t186 + t137 * t57;
t19 = qJD(5) * t60 + t134 * t57 - t137 * t186;
t13 = pkin(4) * t58 - pkin(9) * t57 + t48;
t8 = -pkin(4) * t186 - t232;
t7 = pkin(9) * t186 + t150;
t3 = [0, 0, 0, t136 * t139 * t175, -t209 * t175, t176 * t186, -t176 * t187, 0, t136 * t169 - t164, -t121 * t157 + t133 * t210 + t139 * t169, (t136 * t78 - t139 * t54 + (t136 * t64 + t139 * t61) * qJD(2) + (t136 * t85 - t139 * t73 + (t136 * t79 + t139 * t81) * qJD(2)) * qJD(1)) * t132, (-t70 * t205 + t139 * t52 + (t139 * t69 - t205 * t80) * qJD(1)) * t132 + t164, -t121 * t73 - t133 * t54 + (-t70 * t204 - t136 * t52 + (-t136 * t69 - t204 * t80) * qJD(1)) * t132, t52 * t80 + t54 * t79 + t61 * t85 + t64 * t73 + t69 * t70 + t78 * t81, t38 * t87 + t57 * t77, -t38 * t86 - t39 * t87 - t57 * t75 - t58 * t77, t108 * t57 + (t136 * t38 + (qJD(1) * t87 + t77) * t204) * t132, -t108 * t58 + (-t136 * t39 + (-qJD(1) * t86 - t75) * t204) * t132, (t108 * t132 + t129 * t206) * t204, t232 * t108 + t48 * t75 + t62 * t39 + t32 * t86 + t41 * t58 + (t143 * t136 + (qJD(1) * t160 + t16) * t204) * t132, -t150 * t108 + t48 * t77 + t62 * t38 + t32 * t87 + t41 * t57 + (-t151 * t136 + (-qJD(1) * t159 - t17) * t204) * t132, t14 * t60 + t20 * t36, t14 * t153 - t15 * t60 - t19 * t36 - t20 * t34, t14 * t86 + t20 * t74 + t36 * t58 + t39 * t60, -t15 * t86 + t153 * t39 - t19 * t74 - t34 * t58, t39 * t86 + t58 * t74, (-qJD(5) * t163 + t13 * t137 - t134 * t7) * t74 + t162 * t39 + t2 * t86 - t165 * t58 + t8 * t34 + t21 * t15 - t6 * t153 + t11 * t19, -(qJD(5) * t162 + t13 * t134 + t137 * t7) * t74 - t163 * t39 - t1 * t86 - t4 * t58 + t8 * t36 + t21 * t14 + t6 * t60 + t11 * t20; 0, 0, 0, -t173, t209 * t216, t156, (-qJD(2) + t121) * t189, 0, t216 * t230 + t170, pkin(1) * t192 - t121 * t83 + t210, ((-qJ(3) * qJD(2) - t64 - t84) * t136 + (-pkin(2) * qJD(2) + t211 - t61) * t139) * t208, -t188 * t82 - t170 + t47, t211 * t121 + (t136 * t82 + t139 * t70) * t208 - t54, -pkin(2) * t78 - qJ(3) * t54 - t211 * t64 - t61 * t84 - t70 * t82, t138 * t38 - t219 * t77, (-t39 - t224) * t138 + (-t38 + t225) * t135, -t108 * t201 + t95 + (-t136 * t219 - t139 * t77) * t208, -t185 + (-t193 + (-qJD(2) * t135 + t75) * t139) * t208, -t108 * t188, qJ(3) * t39 + t32 * t135 - t158 * t108 + t212 * t75 + (-t135 * t217 + t138 * t41) * qJD(4) + (t138 * t152 - t16 * t139) * t208, qJ(3) * t38 + t32 * t138 + t226 * t108 + t212 * t77 + (-t135 * t41 - t138 * t217) * qJD(4) + (-t135 * t152 + t17 * t139) * t208, t137 * t138 * t14 + (-t135 * t200 - t138 * t198 - t72) * t36, t34 * t72 + t36 * t71 + (t134 * t36 + t137 * t34) * t201 + (-t223 - t137 * t15 + (t134 * t34 - t137 * t36) * qJD(5)) * t138, -t72 * t74 + (-t200 * t74 + t14) * t135 + (t148 - t154) * t138, t71 * t74 + (t202 * t74 - t15) * t135 + (t155 - t235) * t138, t135 * t39 + t218 * t74, t97 * t221 - t11 * t71 - t23 * t34 + (t134 * t178 - t137 * t227) * t74 + (-t11 * t202 + t2 + (qJD(4) * t34 + t155) * t140) * t135 + (-t165 * t189 + t11 * t197 + t6 * t134 - t140 * t15 + (-t134 * t220 - t165) * qJD(4)) * t138, -t97 * t222 - t11 * t72 - t23 * t36 + (t134 * t227 + t137 * t178) * t74 + (-t11 * t200 - t1 + (qJD(4) * t36 + t154) * t140) * t135 + (-t4 * t189 - t11 * t198 + t6 * t137 - t140 * t14 + (-t137 * t220 - t4) * qJD(4)) * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t173, -t121 ^ 2 - t130 * t216, t121 * t64 + t47 + t78, 0, 0, 0, 0, 0, -t108 * t219 - t121 * t75 + t95, -t185 - t121 * t77 + (-t135 * t204 - t193) * t208, 0, 0, 0, 0, 0, -t121 * t177 + (-t134 * t149 - t15) * t138 + (t155 + t235) * t135, t134 * t121 * t74 + (-t137 * t149 - t14) * t138 + (t148 + t154) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t75 ^ 2 + t77 ^ 2, t38 + t225, -t39 + t224, t110, t108 * t17 - t41 * t77 + t143, t108 * t16 + t41 * t75 - t151, t177 * t36 + t223, (t14 - t229) * t137 + (-t15 - t228) * t134, t177 * t74 - t36 * t77 + t222, -t134 * t74 ^ 2 + t34 * t77 + t221, -t74 * t77, -pkin(4) * t15 + t144 * t134 - t137 * t233 + t165 * t77 - t17 * t34, -pkin(4) * t14 + t134 * t233 + t144 * t137 - t17 * t36 + t4 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t34, -t34 ^ 2 + t36 ^ 2, t14 + t229, -t15 + t228, t39, -t11 * t36 + t4 * t74 + t2, t11 * t34 - t165 * t74 - t1;];
tauc_reg = t3;
