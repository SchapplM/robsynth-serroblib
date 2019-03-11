% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:09
% EndTime: 2019-03-09 03:26:14
% DurationCPUTime: 1.86s
% Computational Cost: add. (3878->318), mult. (8275->423), div. (0->0), fcn. (5516->6), ass. (0->155)
t119 = sin(qJ(3));
t121 = cos(qJ(3));
t188 = sin(pkin(9));
t189 = cos(pkin(9));
t128 = t189 * t119 + t188 * t121;
t81 = t128 * qJD(1);
t213 = qJD(5) + t81;
t118 = sin(qJ(5));
t120 = cos(qJ(5));
t175 = qJD(1) * t119;
t122 = -pkin(1) - pkin(7);
t97 = t122 * qJD(1) + qJD(2);
t77 = -qJ(4) * t175 + t119 * t97;
t158 = t189 * t77;
t174 = qJD(1) * t121;
t78 = -qJ(4) * t174 + t121 * t97;
t73 = qJD(3) * pkin(3) + t78;
t39 = t188 * t73 + t158;
t33 = qJD(3) * pkin(8) + t39;
t88 = -t188 * t119 + t189 * t121;
t84 = t88 * qJD(1);
t91 = pkin(3) * t175 + qJD(1) * qJ(2) + qJD(4);
t40 = t81 * pkin(4) - t84 * pkin(8) + t91;
t12 = t118 * t40 + t120 * t33;
t8 = qJ(6) * t213 + t12;
t215 = t213 * t8;
t167 = t120 * qJD(3);
t170 = qJD(5) * t118;
t151 = qJD(3) * t188;
t137 = qJD(1) * t151;
t152 = qJD(3) * t189;
t138 = qJD(1) * t152;
t195 = t119 * t138 + t121 * t137;
t34 = -qJD(5) * t167 + t120 * t195 + t84 * t170;
t62 = t118 * qJD(3) + t120 * t84;
t83 = -t119 * t152 - t121 * t151;
t141 = -t88 * t34 + t83 * t62;
t149 = -qJD(5) * t128 - qJD(1);
t74 = t119 * t137 - t121 * t138;
t68 = t120 * t74;
t82 = t119 * t151 - t121 * t152;
t214 = -t128 * t68 + (t118 * t149 - t120 * t82) * t213 + t141;
t70 = t188 * t77;
t38 = t189 * t73 - t70;
t32 = -qJD(3) * pkin(4) - t38;
t60 = t118 * t84 - t167;
t13 = t60 * pkin(5) - t62 * qJ(6) + t32;
t106 = t188 * pkin(3) + pkin(8);
t194 = t106 * t74;
t212 = t13 * t213 + t194;
t210 = t62 ^ 2;
t114 = qJD(1) * qJD(2);
t209 = 0.2e1 * t114;
t208 = pkin(5) * t74;
t207 = t12 * t213;
t166 = t121 * qJD(4);
t173 = qJD(3) * t119;
t126 = -t97 * t173 + (qJ(4) * t173 - t166) * qJD(1);
t168 = t119 * qJD(4);
t172 = qJD(3) * t121;
t59 = t97 * t172 + (-qJ(4) * t172 - t168) * qJD(1);
t22 = -t189 * t126 + t188 * t59;
t206 = t22 * t88;
t181 = qJ(4) - t122;
t92 = t181 * t119;
t93 = t181 * t121;
t57 = -t188 * t93 - t189 * t92;
t205 = t57 * t74;
t204 = t60 * t81;
t203 = t62 * t60;
t202 = t62 * t84;
t201 = t84 * t60;
t200 = t128 * t74;
t136 = pkin(5) * t118 - qJ(6) * t120;
t45 = t188 * t78 + t158;
t199 = t118 * qJD(6) - t213 * t136 + t45;
t169 = qJD(5) * t120;
t159 = t118 * t195;
t35 = t62 * qJD(5) - t159;
t198 = -t118 * t35 - t60 * t169;
t46 = t189 * t78 - t70;
t48 = pkin(3) * t174 + t84 * pkin(4) + t81 * pkin(8);
t197 = t118 * t48 + t120 * t46;
t182 = t119 * pkin(3) + qJ(2);
t53 = pkin(4) * t128 - t88 * pkin(8) + t182;
t196 = t118 * t53 + t120 * t57;
t66 = t118 * t74;
t193 = t120 * t213;
t191 = t34 * t118;
t190 = t74 * qJ(6);
t187 = qJD(5) * t88;
t123 = qJD(3) ^ 2;
t186 = t123 * t119;
t185 = t123 * t121;
t124 = qJD(1) ^ 2;
t184 = t124 * qJ(2);
t183 = t124 * t121;
t11 = -t118 * t33 + t120 * t40;
t179 = qJD(6) - t11;
t164 = qJD(1) * qJD(3);
t157 = t121 * t164;
t178 = pkin(3) * t157 + t114;
t177 = t119 ^ 2 - t121 ^ 2;
t176 = -t123 - t124;
t171 = qJD(5) * t106;
t165 = pkin(3) * t172 + qJD(2);
t163 = 0.2e1 * qJD(1);
t162 = t213 * t171;
t161 = t88 * t170;
t160 = t88 * t169;
t156 = t213 * t62;
t23 = t188 * t126 + t189 * t59;
t29 = -t74 * pkin(4) + t195 * pkin(8) + t178;
t153 = t118 * t23 - t120 * t29 + t33 * t169 + t40 * t170;
t150 = t118 * t213;
t4 = t35 * pkin(5) + t34 * qJ(6) - t62 * qJD(6) + t22;
t148 = -t4 - t162;
t132 = t118 * t29 + t120 * t23 + t40 * t169 - t33 * t170;
t1 = qJD(6) * t213 + t132 - t190;
t147 = t1 * t128 - t8 * t82;
t2 = t153 + t208;
t7 = -pkin(5) * t213 + t179;
t146 = t128 * t2 - t7 * t82;
t107 = -t189 * pkin(3) - pkin(4);
t145 = t13 * t83 + t4 * t88;
t144 = -t118 * t8 + t120 * t7;
t143 = t32 * t83 + t206;
t142 = -t128 * t34 - t62 * t82;
t140 = -t128 * t35 + t60 * t82;
t139 = -t213 * t83 + t74 * t88;
t75 = t181 * t173 - t166;
t76 = -qJD(3) * t93 - t168;
t42 = t188 * t76 - t189 * t75;
t56 = -t188 * t92 + t189 * t93;
t135 = t169 * t213 + t81 * t193 - t66;
t134 = -t68 + (-t118 * t81 - t170) * t213;
t133 = t13 * t62 + t153;
t43 = t188 * t75 + t189 * t76;
t47 = -t82 * pkin(4) - t83 * pkin(8) + t165;
t131 = t118 * t47 + t120 * t43 + t53 * t169 - t57 * t170;
t130 = t213 * t32 + t194;
t127 = -t128 * t23 - t38 * t83 + t39 * t82 + t206;
t125 = t128 * t66 - t88 * t35 - t83 * t60 + (t118 * t82 + t149 * t120) * t213;
t86 = -t120 * pkin(5) - t118 * qJ(6) + t107;
t27 = t62 * pkin(5) + t60 * qJ(6);
t18 = t136 * t88 + t56;
t16 = -pkin(5) * t128 + t118 * t57 - t120 * t53;
t15 = qJ(6) * t128 + t196;
t14 = t213 * t60 - t34;
t10 = -t84 * pkin(5) + t118 * t46 - t120 * t48;
t9 = t84 * qJ(6) + t197;
t6 = (pkin(5) * t83 + qJ(6) * t187) * t118 + (-qJ(6) * t83 + (pkin(5) * qJD(5) - qJD(6)) * t88) * t120 + t42;
t5 = t82 * pkin(5) + t196 * qJD(5) + t118 * t43 - t120 * t47;
t3 = -t82 * qJ(6) + qJD(6) * t128 + t131;
t17 = [0, 0, 0, 0, t209, qJ(2) * t209, -0.2e1 * t119 * t157, 0.2e1 * t177 * t164, -t186, -t185, 0, -t122 * t186 + (qJ(2) * t172 + qJD(2) * t119) * t163, -t122 * t185 + (-qJ(2) * t173 + qJD(2) * t121) * t163, -t56 * t195 + t42 * t84 - t43 * t81 + t127 + t205, t165 * t91 + t178 * t182 + t22 * t56 + t23 * t57 - t38 * t42 + t39 * t43, t120 * t141 - t161 * t62 (-t118 * t62 - t120 * t60) * t83 + (t191 - t120 * t35 + (t118 * t60 - t120 * t62) * qJD(5)) * t88, -t120 * t139 - t161 * t213 + t142, t118 * t139 - t160 * t213 + t140, -t213 * t82 - t200, -t153 * t128 - t11 * t82 + t42 * t60 + t56 * t35 + ((-qJD(5) * t57 + t47) * t213 - t53 * t74 + t32 * t187) * t120 + ((-qJD(5) * t53 - t43) * t213 + t205 + t143) * t118, t12 * t82 + t143 * t120 - t128 * t132 - t131 * t213 - t32 * t161 + t196 * t74 - t56 * t34 + t42 * t62, t118 * t145 + t13 * t160 + t16 * t74 + t18 * t35 - t213 * t5 + t6 * t60 - t146, -t15 * t35 - t16 * t34 - t3 * t60 + t5 * t62 + t144 * t83 + (-t1 * t118 + t2 * t120 + (-t118 * t7 - t120 * t8) * qJD(5)) * t88, -t120 * t145 + t13 * t161 - t15 * t74 + t18 * t34 + t213 * t3 - t6 * t62 + t147, t1 * t15 + t13 * t6 + t2 * t16 + t4 * t18 + t8 * t3 + t7 * t5; 0, 0, 0, 0, -t124, -t184, 0, 0, 0, 0, 0, t176 * t119, t176 * t121, t88 * t195 + t82 * t81 - t83 * t84 + t200, -t91 * qJD(1) - t127, 0, 0, 0, 0, 0, t125, -t214, t125 (-t149 * t62 + t140) * t120 + (-t149 * t60 + t142) * t118, t214 (-t149 * t7 + t147) * t120 + (t149 * t8 + t146) * t118 - t145; 0, 0, 0, 0, 0, 0, t119 * t183, -t177 * t124, 0, 0, 0, -qJ(2) * t183, t119 * t184 (t39 - t45) * t84 - (-t46 + t38) * t81 + (t188 * t74 + t189 * t195) * pkin(3), t38 * t45 - t39 * t46 + (-t174 * t91 + t188 * t23 - t189 * t22) * pkin(3), t120 * t156 - t191 (-t34 - t204) * t120 - t62 * t150 + t198, t135 - t202, t134 + t201, -t213 * t84, t107 * t35 - t11 * t84 - t45 * t60 + (-t22 + (-t48 - t171) * t213) * t120 + (t213 * t46 + t130) * t118, -t107 * t34 + t197 * t213 + t12 * t84 - t45 * t62 + (t22 + t162) * t118 + t130 * t120, t10 * t213 + t212 * t118 + t148 * t120 - t199 * t60 + t86 * t35 + t7 * t84, -t10 * t62 + t9 * t60 + (-t106 * t35 + t7 * t81 + t1 + (t106 * t62 + t7) * qJD(5)) * t120 + (-t106 * t34 - t8 * t81 + t2 + (t106 * t60 - t8) * qJD(5)) * t118, t148 * t118 - t212 * t120 + t199 * t62 - t213 * t9 + t86 * t34 - t8 * t84, -t7 * t10 + t4 * t86 - t8 * t9 - t199 * t13 + (qJD(5) * t144 + t1 * t120 + t2 * t118) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81 ^ 2 - t84 ^ 2, t38 * t84 + t39 * t81 + t178, 0, 0, 0, 0, 0, t134 - t201, -t193 * t213 - t202 + t66, -t150 * t213 - t201 - t68 (t34 - t204) * t120 + t118 * t156 + t198, t135 + t202, -t13 * t84 + (-t2 + t215) * t120 + (t213 * t7 + t1) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t60 ^ 2 + t210, t14, t159 + (-qJD(5) + t213) * t62, -t74, -t32 * t62 - t153 + t207, t11 * t213 + t32 * t60 - t132, -t27 * t60 - t133 + t207 - 0.2e1 * t208, pkin(5) * t34 - t35 * qJ(6) + (-t12 + t8) * t62 + (t7 - t179) * t60, -0.2e1 * t190 - t13 * t60 + t27 * t62 + (0.2e1 * qJD(6) - t11) * t213 + t132, -t2 * pkin(5) + t1 * qJ(6) - t7 * t12 - t13 * t27 + t179 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 + t203, t14, -t213 ^ 2 - t210, t133 + t208 - t215;];
tauc_reg  = t17;
