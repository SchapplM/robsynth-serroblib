% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:02
% EndTime: 2019-03-09 01:10:12
% DurationCPUTime: 3.78s
% Computational Cost: add. (4471->390), mult. (13079->709), div. (0->0), fcn. (13262->14), ass. (0->198)
t151 = sin(qJ(4));
t252 = -0.4e1 * t151;
t149 = sin(qJ(6));
t150 = sin(qJ(5));
t154 = cos(qJ(6));
t155 = cos(qJ(5));
t121 = t149 * t155 + t150 * t154;
t106 = t121 * t151;
t156 = cos(qJ(4));
t176 = -pkin(4) * t156 - pkin(11) * t151;
t125 = -pkin(3) + t176;
t226 = t155 * t156;
t135 = pkin(10) * t226;
t99 = t150 * t125 + t135;
t215 = qJD(4) * t156;
t194 = t155 * t215;
t213 = qJD(5) * t150;
t251 = -t151 * t213 + t194;
t143 = t155 ^ 2;
t222 = t150 ^ 2 - t143;
t183 = t222 * qJD(5);
t250 = qJD(5) + qJD(6);
t211 = qJD(5) * t156;
t192 = t150 * t211;
t216 = qJD(4) * t155;
t161 = t151 * t216 + t192;
t175 = pkin(4) * t151 - pkin(11) * t156;
t164 = qJD(4) * t175;
t212 = qJD(5) * t155;
t71 = pkin(10) * t161 - t125 * t212 - t150 * t164;
t145 = sin(pkin(7));
t249 = 0.2e1 * t145;
t248 = pkin(11) + pkin(12);
t152 = sin(qJ(3));
t247 = pkin(2) * t152;
t147 = cos(pkin(7));
t236 = t145 * t152;
t109 = -t156 * t147 + t151 * t236;
t246 = pkin(5) * t109;
t245 = pkin(5) * t154;
t244 = pkin(10) * t145;
t243 = pkin(10) * t150;
t157 = cos(qJ(3));
t100 = pkin(9) * t236 + (-pkin(2) * t157 - pkin(3)) * t147;
t110 = t147 * t151 + t156 * t236;
t67 = t109 * pkin(4) - t110 * pkin(11) + t100;
t235 = t145 * t157;
t205 = pkin(9) * t235;
t101 = t205 + (pkin(10) + t247) * t147;
t102 = (-pkin(3) * t157 - pkin(10) * t152 - pkin(2)) * t145;
t242 = t156 * t101 + t151 * t102;
t69 = -pkin(11) * t235 + t242;
t31 = t150 * t67 + t155 * t69;
t219 = qJD(3) * t152;
t190 = t145 * t219;
t218 = qJD(3) * t157;
t189 = t145 * t218;
t81 = -qJD(4) * t109 + t156 * t189;
t83 = t110 * t150 + t155 * t235;
t44 = -qJD(5) * t83 + t150 * t190 + t155 * t81;
t241 = t150 * t44;
t27 = -pkin(12) * t83 + t31;
t240 = t154 * t27;
t233 = t150 * t151;
t85 = -pkin(12) * t233 + t99;
t239 = t154 * t85;
t238 = t44 * t155;
t128 = t248 * t150;
t237 = t128 * t149;
t231 = t151 * t155;
t153 = sin(qJ(2));
t230 = t152 * t153;
t158 = cos(qJ(2));
t229 = t152 * t158;
t228 = t153 * t157;
t227 = t154 * t155;
t225 = t157 * t158;
t139 = qJD(4) * t151;
t196 = t150 * t139;
t223 = pkin(10) * t196 + t155 * t164;
t142 = t151 ^ 2;
t221 = -t156 ^ 2 + t142;
t146 = sin(pkin(6));
t220 = qJD(2) * t146;
t217 = qJD(4) * t150;
t214 = qJD(4) * t157;
t210 = qJD(6) * t149;
t209 = qJD(6) * t154;
t208 = -0.2e1 * pkin(3) * qJD(4);
t207 = -0.2e1 * pkin(4) * qJD(5);
t206 = t156 * t243;
t204 = pkin(5) * t213;
t203 = pkin(10) * t215;
t202 = pkin(5) * t210;
t201 = pkin(5) * t209;
t200 = t150 * t235;
t105 = (t147 * t247 + t205) * qJD(3);
t82 = qJD(4) * t110 + t151 * t189;
t159 = t82 * pkin(4) - t81 * pkin(11) + t105;
t103 = (pkin(3) * t152 - pkin(10) * t157) * t145 * qJD(3);
t104 = -t147 * pkin(2) * t218 + pkin(9) * t190;
t38 = t101 * t139 - t102 * t215 - t151 * t103 + t156 * t104;
t33 = pkin(11) * t190 - t38;
t13 = -t31 * qJD(5) - t150 * t33 + t155 * t159;
t7 = pkin(5) * t82 - pkin(12) * t44 + t13;
t12 = -t150 * t159 - t155 * t33 - t67 * t212 + t69 * t213;
t45 = -qJD(5) * t200 + t110 * t212 + t150 * t81 - t155 * t190;
t9 = -pkin(12) * t45 - t12;
t199 = -t149 * t9 + t154 * t7;
t140 = t145 ^ 2;
t197 = t140 * t218;
t195 = t150 * t215;
t191 = t155 * t211;
t188 = t153 * t220;
t187 = t150 * t212;
t186 = t151 * t215;
t48 = (pkin(5) * t151 - pkin(12) * t226) * qJD(4) + (-t135 + (pkin(12) * t151 - t125) * t150) * qJD(5) + t223;
t160 = t151 * t212 + t195;
t54 = -pkin(12) * t160 - t71;
t185 = -t149 * t54 + t154 * t48;
t30 = -t150 * t69 + t155 * t67;
t184 = -t151 * t101 + t102 * t156;
t182 = t221 * qJD(4);
t181 = 0.2e1 * t186;
t180 = t140 * t188;
t179 = t145 * t188;
t178 = t152 * t197;
t177 = t150 * t194;
t68 = pkin(4) * t235 - t184;
t84 = t110 * t155 - t200;
t26 = -pkin(12) * t84 + t246 + t30;
t11 = t149 * t26 + t240;
t148 = cos(pkin(6));
t108 = -t145 * t146 * t158 + t148 * t147;
t165 = t147 * t229 + t228;
t78 = t146 * t165 + t148 * t236;
t66 = t108 * t151 + t156 * t78;
t166 = t147 * t225 - t230;
t77 = -t146 * t166 - t148 * t235;
t36 = -t150 * t66 + t155 * t77;
t37 = t150 * t77 + t155 * t66;
t23 = t149 * t36 + t154 * t37;
t119 = t155 * t125;
t76 = -pkin(12) * t231 + t119 + (-pkin(5) - t243) * t156;
t50 = t149 * t76 + t239;
t52 = t149 * t84 + t154 * t83;
t53 = -t149 * t83 + t154 * t84;
t174 = -t150 * t84 - t155 * t83;
t173 = t108 * t156 - t151 * t78;
t129 = t248 * t155;
t87 = t129 * t154 - t237;
t120 = t149 * t150 - t227;
t39 = -t101 * t215 - t102 * t139 + t103 * t156 + t151 * t104;
t1 = -t149 * t7 - t154 * t9 - t26 * t209 + t27 * t210;
t61 = t148 * t189 + (t166 * qJD(3) + (-t147 * t230 + t225) * qJD(2)) * t146;
t28 = qJD(4) * t66 + t61 * t151 - t156 * t179;
t172 = t150 * t28 - t173 * t212;
t171 = -t155 * t28 - t173 * t213;
t34 = -pkin(4) * t190 - t39;
t170 = t34 * t150 + t68 * t212;
t169 = -t34 * t155 + t68 * t213;
t168 = t109 * t212 + t150 * t82;
t167 = t109 * t213 - t155 * t82;
t20 = -t149 * t48 - t154 * t54 - t76 * t209 + t85 * t210;
t163 = t151 * t214 + t156 * t219;
t162 = t151 * t219 - t156 * t214;
t138 = -pkin(5) * t155 - pkin(4);
t132 = -0.2e1 * t186;
t122 = (pkin(5) * t150 + pkin(10)) * t151;
t107 = t120 * t151;
t98 = t119 - t206;
t94 = pkin(5) * t160 + t203;
t86 = -t128 * t154 - t129 * t149;
t80 = t250 * t121;
t79 = t250 * t120;
t73 = 0.2e1 * t109 * t82;
t72 = -t99 * qJD(5) + t223;
t70 = t109 * t139 - t156 * t82;
t60 = t148 * t190 + (t165 * qJD(3) + (t147 * t228 + t229) * qJD(2)) * t146;
t59 = -t87 * qJD(6) + (-t248 * t227 + t237) * qJD(5);
t58 = t121 * qJD(5) * t248 + t128 * t209 + t129 * t210;
t56 = -t210 * t233 + (t250 * t231 + t195) * t154 + t251 * t149;
t55 = -t250 * t106 - t120 * t215;
t49 = -t149 * t85 + t154 * t76;
t40 = pkin(5) * t83 + t68;
t29 = qJD(4) * t173 + t151 * t179 + t61 * t156;
t24 = pkin(5) * t45 + t34;
t22 = -t149 * t37 + t154 * t36;
t21 = -qJD(6) * t50 + t185;
t19 = qJD(6) * t53 + t149 * t44 + t154 * t45;
t18 = -qJD(6) * t52 - t149 * t45 + t154 * t44;
t17 = qJD(5) * t36 + t60 * t150 + t29 * t155;
t16 = -qJD(5) * t37 - t29 * t150 + t60 * t155;
t10 = -t149 * t27 + t154 * t26;
t4 = -qJD(6) * t23 - t149 * t17 + t154 * t16;
t3 = -t149 * t16 - t154 * t17 - t36 * t209 + t37 * t210;
t2 = -qJD(6) * t11 + t199;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t188, -t158 * t220, 0, 0, 0, 0, 0, t108 * t190 - t147 * t60 - t157 * t180, t108 * t189 - t147 * t61 + t152 * t180, 0, 0, 0, 0, 0, t60 * t109 + t77 * t82 + (t157 * t28 + t173 * t219) * t145, t60 * t110 + t77 * t81 + (t157 * t29 - t219 * t66) * t145, 0, 0, 0, 0, 0, t109 * t16 - t173 * t45 + t28 * t83 + t36 * t82, -t109 * t17 - t173 * t44 + t28 * t84 - t37 * t82, 0, 0, 0, 0, 0, t109 * t4 - t173 * t19 + t22 * t82 + t28 * t52, t109 * t3 - t173 * t18 - t23 * t82 + t28 * t53; 0, 0, 0, 0, 0.2e1 * t178, 0.2e1 * (-t152 ^ 2 + t157 ^ 2) * t140 * qJD(3), 0.2e1 * t147 * t189, -0.2e1 * t147 * t190, 0, -0.2e1 * pkin(2) * t140 * t219 - 0.2e1 * t105 * t147, -0.2e1 * pkin(2) * t197 + 0.2e1 * t104 * t147, 0.2e1 * t110 * t81, -0.2e1 * t109 * t81 - 0.2e1 * t110 * t82 (t110 * t219 - t157 * t81) * t249 (-t109 * t219 + t157 * t82) * t249, -0.2e1 * t178, 0.2e1 * t100 * t82 + 0.2e1 * t105 * t109 + 0.2e1 * (-t39 * t157 + t184 * t219) * t145, 0.2e1 * t100 * t81 + 0.2e1 * t105 * t110 + 0.2e1 * (-t38 * t157 - t242 * t219) * t145, 0.2e1 * t84 * t44, -0.2e1 * t44 * t83 - 0.2e1 * t45 * t84, 0.2e1 * t109 * t44 + 0.2e1 * t82 * t84, -0.2e1 * t109 * t45 - 0.2e1 * t82 * t83, t73, 0.2e1 * t109 * t13 + 0.2e1 * t30 * t82 + 0.2e1 * t34 * t83 + 0.2e1 * t45 * t68, 0.2e1 * t109 * t12 - 0.2e1 * t31 * t82 + 0.2e1 * t34 * t84 + 0.2e1 * t44 * t68, 0.2e1 * t53 * t18, -0.2e1 * t18 * t52 - 0.2e1 * t19 * t53, 0.2e1 * t109 * t18 + 0.2e1 * t53 * t82, -0.2e1 * t109 * t19 - 0.2e1 * t52 * t82, t73, 0.2e1 * t10 * t82 + 0.2e1 * t109 * t2 + 0.2e1 * t19 * t40 + 0.2e1 * t24 * t52, 0.2e1 * t1 * t109 - 0.2e1 * t11 * t82 + 0.2e1 * t18 * t40 + 0.2e1 * t24 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t61, 0, 0, 0, 0, 0, t139 * t77 - t156 * t60, t151 * t60 + t215 * t77, 0, 0, 0, 0, 0 (-t173 * t217 - t16) * t156 + (qJD(4) * t36 + t172) * t151 (-t173 * t216 + t17) * t156 + (-qJD(4) * t37 - t171) * t151, 0, 0, 0, 0, 0, t106 * t28 + t139 * t22 - t156 * t4 - t173 * t56, -t107 * t28 - t139 * t23 - t156 * t3 - t173 * t55; 0, 0, 0, 0, 0, 0, t189, -t190, 0, -t105, t104, t110 * t215 + t151 * t81, -t151 * t82 + t81 * t156 + (-t109 * t156 - t110 * t151) * qJD(4), t162 * t145, t163 * t145, 0, -pkin(3) * t82 + t100 * t139 - t105 * t156 - t162 * t244, -pkin(3) * t81 + t100 * t215 + t105 * t151 - t163 * t244, t84 * t194 + (-t213 * t84 + t238) * t151, t174 * t215 + (-t241 - t155 * t45 + (t150 * t83 - t155 * t84) * qJD(5)) * t151 (t109 * t216 - t44) * t156 + (qJD(4) * t84 - t167) * t151 (-t109 * t217 + t45) * t156 + (-qJD(4) * t83 - t168) * t151, t70, t109 * t72 + t82 * t98 + (-t13 + (pkin(10) * t83 + t150 * t68) * qJD(4)) * t156 + (pkin(10) * t45 + qJD(4) * t30 + t170) * t151, t109 * t71 - t82 * t99 + (-t12 + (pkin(10) * t84 + t155 * t68) * qJD(4)) * t156 + (pkin(10) * t44 - qJD(4) * t31 - t169) * t151, -t107 * t18 + t53 * t55, -t106 * t18 + t107 * t19 - t52 * t55 - t53 * t56, -t107 * t82 + t109 * t55 + t139 * t53 - t156 * t18, -t106 * t82 - t109 * t56 - t139 * t52 + t156 * t19, t70, t10 * t139 + t106 * t24 + t109 * t21 + t122 * t19 - t156 * t2 + t40 * t56 + t49 * t82 + t52 * t94, -t1 * t156 - t107 * t24 + t109 * t20 - t11 * t139 + t122 * t18 + t40 * t55 - t50 * t82 + t53 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -0.2e1 * t182, 0, 0, 0, t151 * t208, t156 * t208, -0.2e1 * t142 * t187 + 0.2e1 * t143 * t186, 0.2e1 * t142 * t183 + t177 * t252, 0.2e1 * t151 * t192 + 0.2e1 * t216 * t221, -0.2e1 * t150 * t182 + 0.2e1 * t151 * t191, t132, 0.2e1 * t98 * t139 - 0.2e1 * t72 * t156 + 0.2e1 * (t142 * t212 + t150 * t181) * pkin(10), -0.2e1 * t99 * t139 - 0.2e1 * t71 * t156 + 0.2e1 * (-t142 * t213 + t155 * t181) * pkin(10), -0.2e1 * t107 * t55, -0.2e1 * t106 * t55 + 0.2e1 * t107 * t56, -0.2e1 * t107 * t139 - 0.2e1 * t156 * t55, -0.2e1 * t106 * t139 + 0.2e1 * t156 * t56, t132, 0.2e1 * t106 * t94 + 0.2e1 * t122 * t56 + 0.2e1 * t139 * t49 - 0.2e1 * t156 * t21, -0.2e1 * t107 * t94 + 0.2e1 * t122 * t55 - 0.2e1 * t139 * t50 - 0.2e1 * t156 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, 0, 0, 0, 0, t171, t172, 0, 0, 0, 0, 0, t120 * t28 - t173 * t80, t121 * t28 + t173 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t82, t190, t39, t38, t212 * t84 + t241, qJD(5) * t174 - t150 * t45 + t238, t168, -t167, 0, -pkin(4) * t45 - pkin(11) * t168 + t169, -pkin(4) * t44 + pkin(11) * t167 + t170, t121 * t18 - t53 * t79, -t120 * t18 - t121 * t19 + t52 * t79 - t53 * t80, -t109 * t79 + t121 * t82, -t109 * t80 - t120 * t82, 0, t109 * t59 + t120 * t24 + t138 * t19 + t204 * t52 + t40 * t80 + t82 * t86, t109 * t58 + t121 * t24 + t138 * t18 + t204 * t53 - t40 * t79 - t82 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t139, 0, -t203, pkin(10) * t139, -t151 * t183 + t177, t187 * t252 - t215 * t222, -t191 + t196, t161, 0 (pkin(11) * t226 + (-pkin(4) * t155 + t243) * t151) * qJD(5) + (t150 * t176 - t135) * qJD(4) (pkin(10) * t231 + t150 * t175) * qJD(5) + (t155 * t176 + t206) * qJD(4), t107 * t79 + t121 * t55, t106 * t79 + t107 * t80 - t120 * t55 - t121 * t56, t121 * t139 + t156 * t79, -t120 * t139 + t156 * t80, 0, t106 * t204 + t120 * t94 + t122 * t80 + t138 * t56 + t139 * t86 - t156 * t59, -t107 * t204 + t121 * t94 - t122 * t79 + t138 * t55 - t139 * t87 - t156 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t187, -0.2e1 * t183, 0, 0, 0, t150 * t207, t155 * t207, -0.2e1 * t121 * t79, 0.2e1 * t120 * t79 - 0.2e1 * t121 * t80, 0, 0, 0, 0.2e1 * t120 * t204 + 0.2e1 * t138 * t80, 0.2e1 * t121 * t204 - 0.2e1 * t138 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t45, t82, t13, t12, 0, 0, t18, -t19, t82, t82 * t245 + (-t240 + (-t26 - t246) * t149) * qJD(6) + t199 (-t109 * t209 - t149 * t82) * pkin(5) + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, -t160, t139, t72, t71, 0, 0, t55, -t56, t139, t139 * t245 + (-t239 + (pkin(5) * t156 - t76) * t149) * qJD(6) + t185 (-t139 * t149 + t156 * t209) * pkin(5) + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, -t213, 0, -pkin(11) * t212, pkin(11) * t213, 0, 0, -t79, -t80, 0, t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t202, -0.2e1 * t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, t82, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t56, t139, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t80, 0, t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, -t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
