% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:30
% DurationCPUTime: 2.09s
% Computational Cost: add. (2314->295), mult. (4864->375), div. (0->0), fcn. (2552->6), ass. (0->163)
t71 = -cos(pkin(9)) * pkin(1) - pkin(2);
t60 = qJD(1) * t71;
t108 = cos(qJ(3));
t106 = sin(qJ(3));
t70 = sin(pkin(9)) * pkin(1) + pkin(7);
t59 = t70 * qJD(1);
t187 = t106 * t59;
t143 = -t108 * qJD(2) + t187;
t208 = qJD(4) + t143;
t31 = -qJD(3) * pkin(3) + t208;
t169 = qJD(1) * t108;
t105 = sin(qJ(6));
t107 = cos(qJ(6));
t159 = t107 * qJD(3);
t55 = t105 * t169 - t159;
t170 = qJD(1) * t106;
t69 = qJD(6) + t170;
t121 = t55 * t69;
t162 = qJD(6) * t108;
t148 = t105 * t162;
t209 = t106 * t159 + t148;
t27 = t209 * qJD(1) - qJD(6) * t159;
t204 = t27 - t121;
t168 = qJD(3) * t105;
t56 = t107 * t169 + t168;
t196 = t56 * t69;
t158 = qJD(1) * qJD(3);
t146 = t106 * t158;
t176 = qJD(6) * t56;
t28 = -t105 * t146 + t176;
t116 = t28 - t196;
t109 = -pkin(3) - pkin(4);
t97 = -pkin(8) + t109;
t117 = pkin(5) * t108 + t97 * t106;
t115 = t117 * qJD(3);
t145 = t108 * t158;
t85 = t106 * qJD(4);
t192 = qJ(4) * t145 + qJD(1) * t85;
t12 = qJD(1) * t115 + t192;
t131 = pkin(5) * t106 + pkin(8) * t108;
t36 = -pkin(3) * t169 - qJ(4) * t170 + t60;
t21 = pkin(4) * t169 + qJD(5) - t36;
t13 = t131 * qJD(1) + t21;
t76 = qJ(5) * t170;
t120 = t208 - t76;
t15 = t97 * qJD(3) + t120;
t126 = t105 * t15 - t107 * t13;
t166 = qJD(3) * t108;
t149 = qJ(5) * t166;
t160 = t106 * qJD(5);
t157 = qJD(2) * qJD(3);
t39 = t106 * t157 + t59 * t166;
t16 = (-t149 - t160) * qJD(1) + t39;
t1 = -t126 * qJD(6) + t105 * t12 + t107 * t16;
t207 = t126 * t69 + t1;
t6 = t105 * t13 + t107 * t15;
t2 = -qJD(6) * t6 - t105 * t16 + t107 * t12;
t206 = -t6 * t69 - t2;
t161 = t106 * qJD(2);
t30 = -(qJ(5) * qJD(1) - t59) * t108 + t161;
t98 = qJD(3) * qJ(4);
t205 = t30 + t98;
t44 = t108 * t59 + t161;
t35 = t44 + t98;
t177 = -qJ(5) + t70;
t151 = qJD(3) * t109;
t203 = (-t143 + t187) * qJD(3);
t100 = t108 ^ 2;
t186 = t106 * t69;
t202 = -(qJD(1) * t100 - t186) * t159 + t69 * t148;
t129 = t105 * t6 - t107 * t126;
t201 = -qJD(6) * t129 + t1 * t107 - t2 * t105;
t165 = qJD(5) * t108;
t167 = qJD(3) * t106;
t75 = t108 * t157;
t38 = -t59 * t167 + t75;
t96 = qJD(3) * qJD(4);
t26 = t38 + t96;
t66 = qJ(5) * t146;
t14 = qJD(1) * t165 - t26 - t66;
t51 = t177 * t108;
t198 = t14 * t51;
t197 = t55 * t56;
t195 = t27 * t106 - t56 * t166;
t194 = t209 * t55;
t147 = t107 * t162;
t193 = t105 * t100 * t158 + t69 * t147;
t191 = t75 + 0.2e1 * t96;
t190 = qJ(4) * t166 + t85;
t188 = t105 * t14;
t185 = t107 * t14;
t184 = t107 * t69;
t183 = t108 * t205;
t182 = t108 * t55;
t181 = t27 * t105;
t180 = t28 * t106;
t179 = t39 * t106;
t178 = t39 * t108;
t111 = qJD(1) ^ 2;
t175 = t100 * t111;
t110 = qJD(3) ^ 2;
t88 = t110 * t106;
t89 = t110 * t108;
t19 = qJD(3) * pkin(5) + t205;
t174 = t19 * qJD(3);
t29 = t143 - t76;
t173 = qJD(4) + t29;
t171 = -qJD(5) - t21;
t164 = qJD(6) * t105;
t163 = qJD(6) * t107;
t156 = t105 * t186;
t155 = t106 * t184;
t154 = t56 * t167;
t153 = t69 * t164;
t152 = t69 * t163;
t49 = -t108 * pkin(3) - t106 * qJ(4) + t71;
t42 = t108 * pkin(4) - t49;
t142 = qJD(1) * t42 + t21;
t135 = t106 * t151;
t22 = qJD(1) * t135 + t192;
t32 = t135 + t190;
t141 = qJD(1) * t32 + t22;
t140 = qJD(1) * t49 + t36;
t139 = t171 * t106;
t136 = t105 * t154;
t134 = t106 * t145;
t133 = t44 * qJD(3) - t39;
t128 = t105 * t126 + t107 * t6;
t25 = t131 + t42;
t50 = t177 * t106;
t9 = -t105 * t50 + t107 * t25;
t10 = t105 * t25 + t107 * t50;
t125 = 0.2e1 * qJD(3) * t60;
t37 = pkin(3) * t146 - t192;
t47 = pkin(3) * t167 - t190;
t119 = -qJD(1) * t47 - t110 * t70 - t37;
t118 = -t106 * t19 - t97 * t166;
t114 = t128 * qJD(6) + t1 * t105 + t107 * t2;
t113 = t179 + t26 * t108 + (-t106 * t35 + t108 * t31) * qJD(3);
t112 = t179 + t38 * t108 + (-t106 * t44 + t108 * t143) * qJD(3);
t104 = qJ(4) + pkin(5);
t99 = t106 ^ 2;
t84 = t99 * t111;
t79 = qJ(4) * t169;
t68 = t106 * t111 * t108;
t65 = -t84 - t110;
t62 = -0.2e1 * t134;
t61 = 0.2e1 * t134;
t58 = -t84 + t175;
t57 = pkin(3) * t170 - t79;
t54 = (t100 - t99) * t158;
t48 = 0.2e1 * t54;
t45 = t109 * t170 + t79;
t34 = t177 * t166 - t160;
t33 = t167 * t177 + t165;
t20 = qJD(1) * t117 + t79;
t18 = t151 + t120;
t17 = t115 + t190;
t8 = t105 * t20 + t107 * t30;
t7 = -t105 * t30 + t107 * t20;
t4 = -t10 * qJD(6) - t105 * t34 + t107 * t17;
t3 = t9 * qJD(6) + t105 * t17 + t107 * t34;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t48, t89, t62, -t88, 0, t106 * t125 - t70 * t89, t108 * t125 + t70 * t88, t112, t112 * t70, t61, t89, -0.2e1 * t54, 0, t88, t62, t119 * t108 + t140 * t167, t113, t119 * t106 - t140 * t166, t113 * t70 + t36 * t47 + t37 * t49, t62, t48, -t88, t61, t89, 0, t141 * t106 + (t142 * t108 - t33) * qJD(3), -t141 * t108 + (t142 * t106 + t34) * qJD(3), -t106 * t16 + t108 * t14 + (t106 * t205 - t108 * t18) * qJD(3) + (-t106 * t34 + t108 * t33 + (t106 * t51 - t108 * t50) * qJD(3)) * qJD(1), t16 * t50 + t18 * t34 - t205 * t33 + t21 * t32 + t22 * t42 - t198, -t56 * t148 + (-t108 * t27 - t154) * t107, t136 + (t181 + (-t28 - t176) * t107) * t108 + t194, t195 + t202, t55 * t147 + (t108 * t28 - t55 * t167) * t105, t180 + (-t156 + t182) * qJD(3) + t193 (t69 + t170) * t166, -t28 * t51 + t33 * t55 + t4 * t69 + (t19 * t168 + t2) * t106 + (-t19 * t163 + t188 + (qJD(1) * t9 - t126) * qJD(3)) * t108, t27 * t51 - t3 * t69 + t33 * t56 + (t159 * t19 - t1) * t106 + (t19 * t164 + t185 + (-qJD(1) * t10 - t6) * qJD(3)) * t108, t10 * t28 + t108 * t114 - t129 * t167 - t27 * t9 + t3 * t55 + t4 * t56, t1 * t10 - t126 * t4 - t19 * t33 + t2 * t9 + t3 * t6 - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t89, 0, t38 * t106 - t178 + (t106 * t143 + t108 * t44) * qJD(3), 0, 0, 0, 0, 0, 0, -t88, 0, t89, t26 * t106 - t178 + (t106 * t31 + t108 * t35) * qJD(3), 0, 0, 0, 0, 0, 0, t89, t88, 0, -t14 * t106 - t16 * t108 + (t106 * t18 + t183) * qJD(3), 0, 0, 0, 0, 0, 0, -t180 + (-t156 - t182) * qJD(3) + t193, t195 - t202, -t136 + (-t181 + (-t28 + t176) * t107) * t108 + t194 (qJD(3) * t128 - t14) * t106 + (t174 - t201) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t58, 0, t68, 0, 0, -t60 * t170 + t133, -t60 * t169 + t203 - t75, 0, 0, -t68, 0, t58, 0, 0, t68 (-t106 * t36 + t108 * t57) * qJD(1) + t133, 0, -t203 + (t106 * t57 + t108 * t36) * qJD(1) + t191, -t39 * pkin(3) + t26 * qJ(4) + t208 * t35 - t31 * t44 - t36 * t57, t68, -t58, 0, -t68, 0, 0, t66 + (t29 - t187) * qJD(3) + (-t106 * t45 + t171 * t108) * qJD(1) + t191, -t30 * qJD(3) + ((-qJ(5) * qJD(3) + t45) * t108 + t139) * qJD(1) + t39 (-t173 + t18 - t151) * t169, -t14 * qJ(4) + t16 * t109 + t173 * t205 - t18 * t30 - t21 * t45, t56 * t184 - t181 (-t27 - t121) * t107 + (-t28 - t196) * t105, -t152 + (-t155 + (t56 - t168) * t108) * qJD(1), t105 * t121 - t28 * t107, t153 + (t156 + (-t55 - t159) * t108) * qJD(1), -t69 * t169, -t104 * t28 - t185 - t69 * t7 - t173 * t55 + (-t105 * t19 - t97 * t184) * qJD(6) + (t105 * t118 + t108 * t126) * qJD(1), t104 * t27 + t188 + t69 * t8 - t173 * t56 + (t105 * t69 * t97 - t107 * t19) * qJD(6) + (t107 * t118 + t108 * t6) * qJD(1), -t55 * t8 - t56 * t7 + (-t126 * t170 + t28 * t97 - t1 + (-t56 * t97 - t126) * qJD(6)) * t107 + (t6 * t170 + t27 * t97 + t2 + (-t55 * t97 + t6) * qJD(6)) * t105, -t104 * t14 + t126 * t7 + t173 * t19 + t201 * t97 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, t65, -qJD(3) * t35 + t36 * t170 + t39, 0, 0, 0, 0, 0, 0, t65, t68, 0, -t205 * qJD(3) + (t139 - t149) * qJD(1) + t39, 0, 0, 0, 0, 0, 0, -t152 + qJD(3) * t55 + (-t105 * t166 - t155) * qJD(1), t153 + qJD(3) * t56 + (-t108 * t159 + t156) * qJD(1), t204 * t105 + t116 * t107, t206 * t105 + t207 * t107 - t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145, 0.2e1 * t146, -t84 - t175 (t183 + (t18 + t151) * t106) * qJD(1) + t192, 0, 0, 0, 0, 0, 0, -t153 + (-t156 + (-t55 + t159) * t108) * qJD(1), -t152 + (-t155 + (-t56 - t168) * t108) * qJD(1), t116 * t105 - t107 * t204 (t106 * t128 + t108 * t19) * qJD(1) + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t55 ^ 2 + t56 ^ 2, t204, -t197, t116, t145, t19 * t56 - t206, -t19 * t55 - t207, 0, 0;];
tauc_reg  = t5;
