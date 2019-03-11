% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:40
% EndTime: 2019-03-08 19:01:44
% DurationCPUTime: 2.01s
% Computational Cost: add. (2990->269), mult. (7870->414), div. (0->0), fcn. (6856->14), ass. (0->175)
t123 = cos(qJ(5));
t124 = cos(qJ(4));
t183 = qJD(3) * t124;
t101 = t123 * t183;
t119 = sin(qJ(5));
t120 = sin(qJ(4));
t185 = qJD(3) * t120;
t163 = t119 * t185;
t83 = -t101 + t163;
t80 = qJD(6) + t83;
t221 = qJD(6) - t80;
t213 = pkin(9) + pkin(10);
t117 = cos(pkin(6));
t113 = sin(pkin(7));
t125 = cos(qJ(3));
t190 = t113 * t125;
t112 = sin(pkin(13));
t114 = sin(pkin(6));
t121 = sin(qJ(3));
t115 = cos(pkin(13));
t116 = cos(pkin(7));
t188 = t115 * t116;
t216 = (-t112 * t121 + t125 * t188) * t114;
t220 = t117 * t190 + t216;
t102 = qJD(1) * t117 + qJD(2);
t219 = qJD(1) * t216 + t102 * t190;
t106 = pkin(4) * t119 + pkin(11);
t176 = pkin(4) * t185;
t85 = -t119 * t183 - t123 * t185;
t64 = -pkin(5) * t85 + pkin(11) * t83;
t218 = (qJD(6) * t106 + t176 + t64) * t80;
t186 = qJD(1) * t114;
t166 = t115 * t186;
t151 = t116 * t166;
t184 = qJD(3) * t121;
t165 = t113 * t184;
t167 = t112 * t186;
t182 = qJD(3) * t125;
t51 = t102 * t165 + t151 * t184 + t167 * t182;
t131 = (t112 * t125 + t121 * t188) * t114;
t191 = t113 * t121;
t59 = qJD(1) * t131 + t102 * t191;
t217 = qJD(3) * t59 - t51;
t178 = t120 * qJD(4);
t175 = pkin(4) * t178;
t148 = -t59 + t175;
t160 = t213 * qJD(3) + t59;
t78 = t102 * t116 - t113 * t166;
t36 = -t160 * t120 + t124 * t78;
t50 = t219 * qJD(3);
t22 = t36 * qJD(4) + t124 * t50;
t35 = qJD(4) * pkin(4) + t36;
t215 = (qJD(5) * t35 + t22) * t123;
t109 = qJD(4) + qJD(5);
t58 = -t121 * t167 + (t102 * t113 + t151) * t125;
t37 = t120 * t78 + t160 * t124;
t196 = t119 * t37;
t13 = t123 * t35 - t196;
t11 = -pkin(5) * t109 - t13;
t181 = qJD(5) * t119;
t23 = -t37 * qJD(4) - t120 * t50;
t157 = t119 * t23 - t37 * t181;
t2 = t157 + t215;
t94 = t213 * t120;
t95 = t213 * t124;
t76 = t119 * t95 + t123 * t94;
t88 = t119 * t120 - t123 * t124;
t168 = qJD(4) * t213;
t91 = t120 * t168;
t92 = t124 * t168;
t203 = t76 * qJD(5) + t119 * t92 + t123 * t91 - t88 * t58;
t193 = t123 * t37;
t14 = t119 * t35 + t193;
t158 = t119 * t22 - t123 * t23;
t3 = t14 * qJD(5) + t158;
t108 = -pkin(4) * t124 - pkin(3);
t49 = t108 * qJD(3) - t58;
t34 = pkin(5) * t83 + pkin(11) * t85 + t49;
t89 = t119 * t124 + t120 * t123;
t69 = t109 * t89;
t61 = t69 * qJD(3);
t65 = pkin(5) * t88 - pkin(11) * t89 + t108;
t68 = t109 * t88;
t77 = -t119 * t94 + t123 * t95;
t214 = -(qJD(6) * t34 + t2) * t88 - t11 * t68 + t3 * t89 + (-qJD(6) * t65 + t203) * t80 - t77 * t61;
t118 = sin(qJ(6));
t177 = qJD(3) * qJD(4);
t161 = t124 * t177;
t60 = qJD(5) * t101 - t109 * t163 + t123 * t161;
t122 = cos(qJ(6));
t74 = t109 * t118 - t122 * t85;
t40 = t74 * qJD(6) + t118 * t60;
t212 = pkin(4) * t123;
t211 = t11 * t83;
t210 = t11 * t89;
t209 = t61 * t89;
t208 = t65 * t61;
t72 = -t122 * t109 - t118 * t85;
t207 = t72 * t80;
t206 = t74 * t80;
t205 = t80 * t85;
t204 = t85 * t83;
t202 = pkin(5) * t69 + pkin(11) * t68 + t148;
t201 = t77 * qJD(5) - t119 * t91 + t123 * t92 - t89 * t58;
t200 = qJD(3) * pkin(3);
t179 = qJD(6) * t122;
t180 = qJD(6) * t118;
t39 = t109 * t179 + t122 * t60 + t85 * t180;
t199 = t118 * t39;
t197 = t118 * t61;
t194 = t122 * t61;
t127 = qJD(3) ^ 2;
t189 = t113 * t127;
t187 = t120 ^ 2 - t124 ^ 2;
t12 = pkin(11) * t109 + t14;
t144 = t118 * t12 - t122 * t34;
t174 = t11 * t180 - t144 * t85;
t173 = t89 * t180;
t172 = t80 * t179;
t169 = t121 * t189;
t164 = t113 * t182;
t162 = -pkin(4) * t109 - t35;
t52 = -t58 - t200;
t156 = -t52 * qJD(3) - t50;
t154 = t122 * t80;
t47 = qJD(3) * t175 + t51;
t8 = t118 * t34 + t12 * t122;
t152 = t11 * t179 + t3 * t118 - t8 * t85;
t150 = t120 * t164;
t149 = t124 * t164;
t15 = t119 * t36 + t193;
t147 = pkin(4) * t181 - t15;
t146 = -t68 * t80 + t209;
t67 = t117 * t191 + t131;
t79 = -t113 * t114 * t115 + t116 * t117;
t43 = -t120 * t67 + t124 * t79;
t44 = t120 * t79 + t124 * t67;
t26 = t119 * t43 + t123 * t44;
t143 = -t118 * t220 + t122 * t26;
t142 = -t118 * t26 - t122 * t220;
t25 = t119 * t44 - t123 * t43;
t81 = t116 * t124 - t120 * t191;
t82 = t116 * t120 + t124 * t191;
t55 = t119 * t82 - t123 * t81;
t56 = t119 * t81 + t123 * t82;
t141 = t49 * t85 - t158;
t140 = t49 * t83 - t157;
t126 = qJD(4) ^ 2;
t139 = pkin(9) * t126 - t217;
t138 = qJD(4) * (t52 + t58 - t200);
t137 = -t118 * t56 - t122 * t190;
t136 = t118 * t190 - t122 * t56;
t16 = t123 * t36 - t196;
t129 = -t106 * t61 + t211 + (-qJD(5) * t212 + t16) * t80;
t107 = -pkin(5) - t212;
t71 = -t82 * qJD(4) - t150;
t70 = t81 * qJD(4) + t149;
t63 = t67 * qJD(3);
t62 = t220 * qJD(3);
t48 = -t83 ^ 2 + t85 ^ 2;
t46 = (-qJD(3) * t89 - t85) * t109;
t45 = t109 * t83 + t60;
t30 = t43 * qJD(4) + t124 * t62;
t29 = -t44 * qJD(4) - t120 * t62;
t28 = t56 * qJD(5) + t119 * t70 - t123 * t71;
t27 = -t55 * qJD(5) + t119 * t71 + t123 * t70;
t24 = pkin(5) * t61 - pkin(11) * t60 + t47;
t21 = t122 * t24;
t19 = t80 * t154 + t74 * t85 + t197;
t18 = -t80 ^ 2 * t118 - t72 * t85 + t194;
t17 = t74 * t154 + t199;
t6 = (t39 - t207) * t122 + (-t40 - t206) * t118;
t5 = t26 * qJD(5) + t119 * t30 - t123 * t29;
t4 = -t25 * qJD(5) + t119 * t29 + t123 * t30;
t1 = [0, 0, 0, -t63 * qJD(3), -t62 * qJD(3), 0, 0, 0, 0, 0, qJD(4) * t29 + (-t124 * t63 - t178 * t220) * qJD(3), -qJD(4) * t30 + (-qJD(4) * t124 * t220 + t120 * t63) * qJD(3), 0, 0, 0, 0, 0, -t109 * t5 - t220 * t61 + t63 * t83, -t109 * t4 - t220 * t60 - t63 * t85, 0, 0, 0, 0, 0 (-qJD(6) * t143 - t118 * t4 + t122 * t63) * t80 + t142 * t61 + t5 * t72 + t25 * t40 -(qJD(6) * t142 + t118 * t63 + t122 * t4) * t80 - t143 * t61 + t5 * t74 + t25 * t39; 0, 0, 0, -t169, -t125 * t189, 0, 0, 0, 0, 0, -t124 * t169 + (t71 - t150) * qJD(4), t120 * t169 + (-t70 - t149) * qJD(4), 0, 0, 0, 0, 0, -t109 * t28 + (-t125 * t61 + t83 * t184) * t113, -t109 * t27 + (-t125 * t60 - t85 * t184) * t113, 0, 0, 0, 0, 0 (qJD(6) * t136 - t118 * t27 + t122 * t165) * t80 + t137 * t61 + t28 * t72 + t55 * t40 -(qJD(6) * t137 + t118 * t165 + t122 * t27) * t80 + t136 * t61 + t28 * t74 + t55 * t39; 0, 0, 0, t217 (-t219 + t58) * qJD(3), 0.2e1 * t120 * t161, -0.2e1 * t187 * t177, t126 * t124, -t126 * t120, 0, t120 * t138 - t139 * t124, t139 * t120 + t124 * t138, t60 * t89 + t68 * t85, -t60 * t88 + t68 * t83 + t69 * t85 - t209, -t68 * t109, -t69 * t109, 0, t108 * t61 - t201 * t109 + t148 * t83 + t47 * t88 + t49 * t69, t108 * t60 + t203 * t109 - t148 * t85 + t47 * t89 - t49 * t68, -t74 * t173 + (t39 * t89 - t68 * t74) * t122 -(-t118 * t74 - t122 * t72) * t68 + (-t199 - t122 * t40 + (t118 * t72 - t122 * t74) * qJD(6)) * t89, t122 * t146 - t80 * t173 + t39 * t88 + t69 * t74, -t118 * t146 - t89 * t172 - t40 * t88 - t69 * t72, t61 * t88 + t69 * t80, t21 * t88 + t76 * t40 - t144 * t69 + t201 * t72 + (t208 + t202 * t80 + (-t12 * t88 - t77 * t80 + t210) * qJD(6)) * t122 + t214 * t118, t76 * t39 - t8 * t69 + t201 * t74 + (-t208 - (-qJD(6) * t12 + t24) * t88 - qJD(6) * t210 + (qJD(6) * t77 - t202) * t80) * t118 + t214 * t122; 0, 0, 0, 0, 0, -t120 * t127 * t124, t187 * t127, 0, 0, 0, t156 * t120, t156 * t124, -t204, t48, t45, t46, 0, -t83 * t176 + t109 * t15 + (t162 * t119 - t193) * qJD(5) + t141, t85 * t176 + t109 * t16 + (t162 * qJD(5) - t22) * t123 + t140, t17, t6, t19, t18, t205, t107 * t40 + t147 * t72 + (-t3 - t218) * t122 + t129 * t118 + t174, t107 * t39 + t118 * t218 + t129 * t122 + t147 * t74 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t48, t45, t46, 0, t141 + (-qJD(5) + t109) * t14, t109 * t13 + t140 - t215, t17, t6, t19, t18, t205, -pkin(5) * t40 - t3 * t122 - (-t118 * t13 + t122 * t64) * t80 - t14 * t72 + t118 * t211 + (-t172 - t197) * pkin(11) + t174, -pkin(5) * t39 + (t118 * t64 + t122 * t13) * t80 - t14 * t74 + t122 * t211 + (t80 * t180 - t194) * pkin(11) + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t72 ^ 2 + t74 ^ 2, t39 + t207, t206 - t40, t61, -t11 * t74 - t118 * t2 - t221 * t8 + t21, t11 * t72 - t118 * t24 - t122 * t2 + t144 * t221;];
tauc_reg  = t1;
