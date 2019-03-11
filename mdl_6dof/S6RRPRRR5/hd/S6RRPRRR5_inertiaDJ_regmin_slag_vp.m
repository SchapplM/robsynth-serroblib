% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:47
% EndTime: 2019-03-09 13:46:59
% DurationCPUTime: 3.80s
% Computational Cost: add. (6572->356), mult. (18507->642), div. (0->0), fcn. (18948->12), ass. (0->198)
t160 = sin(qJ(4));
t260 = -0.4e1 * t160;
t155 = sin(pkin(6));
t149 = t155 ^ 2;
t161 = sin(qJ(2));
t224 = qJD(2) * t161;
t259 = t149 * t224;
t163 = cos(qJ(5));
t159 = sin(qJ(5));
t219 = qJD(5) * t159;
t154 = sin(pkin(12));
t156 = cos(pkin(12));
t165 = cos(qJ(2));
t223 = qJD(2) * t165;
t195 = t155 * t223;
t196 = t155 * t224;
t108 = -t154 * t196 + t156 * t195;
t110 = (t154 * t165 + t156 * t161) * t155;
t157 = cos(pkin(6));
t164 = cos(qJ(4));
t96 = t110 * t164 + t157 * t160;
t74 = qJD(4) * t96 + t108 * t160;
t95 = t110 * t160 - t157 * t164;
t172 = -t163 * t74 + t219 * t95;
t258 = t172 * t160;
t145 = -pkin(2) * t156 - pkin(3);
t181 = -pkin(4) * t164 - pkin(10) * t160;
t127 = t145 + t181;
t144 = pkin(2) * t154 + pkin(9);
t230 = t163 * t164;
t132 = t144 * t230;
t92 = t127 * t159 + t132;
t252 = pkin(1) * t157;
t213 = t161 * t252;
t238 = t155 * t165;
t250 = pkin(8) + qJ(3);
t104 = t238 * t250 + t213;
t239 = t155 * t161;
t257 = -qJD(2) * t104 - qJD(3) * t239;
t221 = qJD(4) * t164;
t189 = t163 * t221;
t256 = -t160 * t219 + t189;
t158 = sin(qJ(6));
t162 = cos(qJ(6));
t131 = t158 * t163 + t159 * t162;
t255 = qJD(5) + qJD(6);
t217 = qJD(5) * t164;
t198 = t159 * t217;
t222 = qJD(4) * t163;
t120 = t160 * t222 + t198;
t180 = pkin(4) * t160 - pkin(10) * t164;
t170 = t180 * qJD(4);
t218 = qJD(5) * t163;
t63 = t120 * t144 - t127 * t218 - t159 * t170;
t254 = pkin(10) + pkin(11);
t253 = t95 * pkin(5);
t251 = pkin(5) * t162;
t109 = t154 * t239 - t156 * t238;
t188 = t250 * t161;
t98 = (pkin(1) * t165 + pkin(2)) * t157 - t155 * t188;
t77 = t104 * t156 + t154 * t98;
t62 = pkin(9) * t157 + t77;
t205 = -pkin(2) * t165 - pkin(1);
t82 = t109 * pkin(3) - t110 * pkin(9) + t155 * t205;
t248 = t160 * t82 + t164 * t62;
t35 = pkin(10) * t109 + t248;
t76 = -t104 * t154 + t156 * t98;
t61 = -pkin(3) * t157 - t76;
t39 = pkin(4) * t95 - pkin(10) * t96 + t61;
t22 = t159 * t39 + t163 * t35;
t114 = t131 * t160;
t192 = t159 * t221;
t216 = qJD(6) * t158;
t233 = t160 * t163;
t235 = t159 * t160;
t67 = -t216 * t235 + (t233 * t255 + t192) * t162 + t256 * t158;
t249 = -t114 * t74 - t67 * t95;
t177 = t109 * t163 - t159 * t96;
t17 = pkin(11) * t177 + t22;
t247 = t162 * t17;
t88 = -pkin(11) * t235 + t92;
t246 = t162 * t88;
t107 = qJD(2) * t110;
t242 = qJD(4) * t95;
t75 = t108 * t164 - t242;
t33 = qJD(5) * t177 + t107 * t159 + t75 * t163;
t245 = t33 * t159;
t244 = t33 * t163;
t243 = qJD(4) * t177;
t241 = t144 * t159;
t240 = t144 * t164;
t138 = t254 * t159;
t237 = t158 * t138;
t234 = t160 * t107;
t231 = t162 * t163;
t229 = t164 * t107;
t148 = qJD(4) * t160;
t193 = t159 * t148;
t227 = t144 * t193 + t163 * t170;
t152 = t163 ^ 2;
t226 = t159 ^ 2 - t152;
t151 = t160 ^ 2;
t225 = -t164 ^ 2 + t151;
t220 = qJD(5) * t151;
t215 = qJD(6) * t162;
t214 = -0.2e1 * pkin(4) * qJD(5);
t212 = 0.2e1 * qJD(4) * t145;
t211 = pkin(5) * t219;
t210 = pkin(5) * t216;
t209 = pkin(5) * t215;
t208 = t159 * t242;
t207 = t95 * t222;
t206 = t159 * t240;
t141 = t223 * t252;
t89 = t141 + (-qJD(2) * t188 + qJD(3) * t165) * t155;
t55 = t154 * t89 - t156 * t257;
t167 = pkin(4) * t74 - pkin(10) * t75 + t55;
t56 = t154 * t257 + t156 * t89;
t140 = pkin(2) * t196;
t79 = pkin(3) * t107 - pkin(9) * t108 + t140;
t23 = t148 * t62 - t160 * t79 - t164 * t56 - t221 * t82;
t19 = pkin(10) * t107 - t23;
t9 = -qJD(5) * t22 - t159 * t19 + t163 * t167;
t5 = t74 * pkin(5) - t33 * pkin(11) + t9;
t81 = t109 * t159 + t163 * t96;
t32 = qJD(5) * t81 - t107 * t163 + t75 * t159;
t8 = -t159 * t167 - t163 * t19 - t218 * t39 + t219 * t35;
t7 = -pkin(11) * t32 - t8;
t204 = -t158 * t7 + t162 * t5;
t202 = t149 * t223;
t201 = t144 * t221;
t200 = t144 * t220;
t197 = t163 * t217;
t191 = t159 * t218;
t190 = t160 * t221;
t41 = t158 * t81 - t162 * t177;
t13 = -qJD(6) * t41 - t158 * t32 + t162 * t33;
t42 = t158 * t177 + t162 * t81;
t187 = -t13 * t164 + t148 * t42;
t50 = (pkin(5) * t160 - pkin(11) * t230) * qJD(4) + (-t132 + (pkin(11) * t160 - t127) * t159) * qJD(5) + t227;
t168 = t160 * t218 + t192;
t54 = -pkin(11) * t168 - t63;
t186 = -t158 * t54 + t162 * t50;
t21 = -t159 * t35 + t163 * t39;
t185 = -t160 * t62 + t164 * t82;
t184 = t226 * qJD(5);
t183 = t225 * qJD(4);
t182 = t159 * t189;
t130 = t158 * t159 - t231;
t115 = t130 * t160;
t94 = t255 * t131;
t66 = t158 * t192 + t160 * t94 - t162 * t189;
t179 = t115 * t74 + t66 * t95;
t16 = -pkin(11) * t81 + t21 + t253;
t11 = t158 * t16 + t247;
t117 = t163 * t127;
t86 = -pkin(11) * t233 + t117 + (-pkin(5) - t241) * t164;
t49 = t158 * t86 + t246;
t178 = -t159 * t81 + t163 * t177;
t24 = -t148 * t82 - t160 * t56 + t164 * t79 - t221 * t62;
t139 = t254 * t163;
t102 = t139 * t162 - t237;
t34 = -pkin(4) * t109 - t185;
t1 = -t158 * t5 - t16 * t215 - t162 * t7 + t17 * t216;
t14 = qJD(6) * t42 + t158 * t33 + t162 * t32;
t176 = t14 * t164 - t148 * t41;
t20 = -pkin(4) * t107 - t24;
t175 = t20 * t159 + t218 * t34;
t174 = -t20 * t163 + t219 * t34;
t173 = t159 * t74 + t218 * t95;
t171 = t130 * t148 - t164 * t94;
t25 = -t158 * t50 - t162 * t54 - t215 * t86 + t216 * t88;
t169 = t109 * t221 + t234;
t147 = -pkin(5) * t163 - pkin(4);
t142 = -0.2e1 * t190;
t122 = t193 - t197;
t118 = (pkin(5) * t159 + t144) * t160;
t113 = (-pkin(8) * t238 - t213) * qJD(2);
t112 = pkin(8) * t196 - t141;
t101 = -t138 * t162 - t139 * t158;
t97 = pkin(5) * t168 + t201;
t93 = t255 * t130;
t91 = t117 - t206;
t84 = -t109 * t148 + t229;
t83 = t131 * t148 + t164 * t93;
t71 = -t102 * qJD(6) + (-t231 * t254 + t237) * qJD(5);
t70 = qJD(5) * t131 * t254 + t138 * t215 + t139 * t216;
t65 = t81 * t148;
t64 = -qJD(5) * t92 + t227;
t48 = -t158 * t88 + t162 * t86;
t45 = 0.2e1 * t95 * t74;
t43 = t148 * t95 - t164 * t74;
t27 = -pkin(5) * t177 + t34;
t26 = -qJD(6) * t49 + t186;
t12 = pkin(5) * t32 + t20;
t10 = -t158 * t17 + t16 * t162;
t2 = -qJD(6) * t11 + t204;
t3 = [0, 0, 0, 0.2e1 * t161 * t202, 0.2e1 * (-t161 ^ 2 + t165 ^ 2) * t149 * qJD(2), 0.2e1 * t157 * t195, -0.2e1 * t157 * t196, 0, -0.2e1 * pkin(1) * t259 + 0.2e1 * t113 * t157, -0.2e1 * pkin(1) * t202 + 0.2e1 * t112 * t157, -0.2e1 * t107 * t77 - 0.2e1 * t108 * t76 - 0.2e1 * t109 * t56 + 0.2e1 * t110 * t55, 0.2e1 * pkin(2) * t205 * t259 - 0.2e1 * t76 * t55 + 0.2e1 * t77 * t56, 0.2e1 * t96 * t75, -0.2e1 * t74 * t96 - 0.2e1 * t75 * t95, 0.2e1 * t107 * t96 + 0.2e1 * t109 * t75, -0.2e1 * t107 * t95 - 0.2e1 * t109 * t74, 0.2e1 * t109 * t107, 0.2e1 * t107 * t185 + 0.2e1 * t109 * t24 + 0.2e1 * t55 * t95 + 0.2e1 * t61 * t74, -0.2e1 * t107 * t248 + 0.2e1 * t109 * t23 + 0.2e1 * t55 * t96 + 0.2e1 * t61 * t75, 0.2e1 * t81 * t33, 0.2e1 * t177 * t33 - 0.2e1 * t32 * t81, 0.2e1 * t33 * t95 + 0.2e1 * t74 * t81, 0.2e1 * t177 * t74 - 0.2e1 * t32 * t95, t45, -0.2e1 * t177 * t20 + 0.2e1 * t21 * t74 + 0.2e1 * t32 * t34 + 0.2e1 * t9 * t95, 0.2e1 * t20 * t81 - 0.2e1 * t22 * t74 + 0.2e1 * t33 * t34 + 0.2e1 * t8 * t95, 0.2e1 * t42 * t13, -0.2e1 * t13 * t41 - 0.2e1 * t14 * t42, 0.2e1 * t13 * t95 + 0.2e1 * t42 * t74, -0.2e1 * t14 * t95 - 0.2e1 * t41 * t74, t45, 0.2e1 * t10 * t74 + 0.2e1 * t12 * t41 + 0.2e1 * t14 * t27 + 0.2e1 * t2 * t95, 0.2e1 * t1 * t95 - 0.2e1 * t11 * t74 + 0.2e1 * t12 * t42 + 0.2e1 * t13 * t27; 0, 0, 0, 0, 0, t195, -t196, 0, t113, t112 (-t107 * t154 - t108 * t156) * pkin(2) (t154 * t56 - t156 * t55) * pkin(2), t160 * t75 + t221 * t96, -t160 * t74 + t75 * t164 + (-t160 * t96 - t164 * t95) * qJD(4), t169, t84, 0, -t144 * t234 + t145 * t74 - t55 * t164 + (-t109 * t240 + t160 * t61) * qJD(4), -t144 * t229 + t145 * t75 + t55 * t160 + (t109 * t144 * t160 + t164 * t61) * qJD(4), t81 * t189 + (-t219 * t81 + t244) * t160, t178 * t221 + (-t245 - t163 * t32 + (-t159 * t177 - t163 * t81) * qJD(5)) * t160, t65 + (-t33 + t207) * t164 - t258 (t32 - t208) * t164 + (-t173 + t243) * t160, t43, t64 * t95 + t91 * t74 + (-t9 + (-t144 * t177 + t159 * t34) * qJD(4)) * t164 + (qJD(4) * t21 + t144 * t32 + t175) * t160, t63 * t95 - t92 * t74 + (-t8 + (t144 * t81 + t163 * t34) * qJD(4)) * t164 + (-qJD(4) * t22 + t144 * t33 - t174) * t160, -t115 * t13 - t42 * t66, -t114 * t13 + t115 * t14 + t41 * t66 - t42 * t67, -t179 + t187, t176 + t249, t43, t10 * t148 + t114 * t12 + t118 * t14 - t164 * t2 + t26 * t95 + t27 * t67 + t41 * t97 + t48 * t74, -t1 * t164 - t11 * t148 - t115 * t12 + t118 * t13 + t25 * t95 - t27 * t66 + t42 * t97 - t49 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t190, -0.2e1 * t183, 0, 0, 0, t160 * t212, t164 * t212, -0.2e1 * t151 * t191 + 0.2e1 * t152 * t190, t182 * t260 + 0.2e1 * t220 * t226, 0.2e1 * t160 * t198 + 0.2e1 * t222 * t225, -0.2e1 * t159 * t183 + 0.2e1 * t160 * t197, t142, 0.2e1 * t163 * t200 - 0.2e1 * t64 * t164 + 0.2e1 * (t91 + 0.2e1 * t206) * t148, -0.2e1 * t159 * t200 - 0.2e1 * t63 * t164 + 0.2e1 * (-t92 + 0.2e1 * t132) * t148, 0.2e1 * t115 * t66, 0.2e1 * t114 * t66 + 0.2e1 * t115 * t67, -0.2e1 * t115 * t148 + 0.2e1 * t164 * t66, -0.2e1 * t114 * t148 + 0.2e1 * t164 * t67, t142, 0.2e1 * t114 * t97 + 0.2e1 * t118 * t67 + 0.2e1 * t148 * t48 - 0.2e1 * t164 * t26, -0.2e1 * t115 * t97 - 0.2e1 * t118 * t66 - 0.2e1 * t148 * t49 - 0.2e1 * t164 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, 0, 0, 0, 0, t84, -t169, 0, 0, 0, 0, 0 (-t32 - t208) * t164 + (-t173 - t243) * t160, t65 + (-t33 - t207) * t164 + t258, 0, 0, 0, 0, 0, -t176 + t249, t179 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t74, t107, t24, t23, t218 * t81 + t245, qJD(5) * t178 - t159 * t32 + t244, t173, -t172, 0, -pkin(4) * t32 - pkin(10) * t173 + t174, -pkin(4) * t33 + pkin(10) * t172 + t175, t13 * t131 - t42 * t93, -t13 * t130 - t131 * t14 + t41 * t93 - t42 * t94, t131 * t74 - t93 * t95, -t130 * t74 - t94 * t95, 0, t101 * t74 + t12 * t130 + t14 * t147 + t211 * t41 + t27 * t94 + t71 * t95, -t102 * t74 + t12 * t131 + t13 * t147 + t211 * t42 - t27 * t93 + t70 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -t148, 0, -t201, t144 * t148, -t160 * t184 + t182, t191 * t260 - t221 * t226, t122, t120, 0 (pkin(10) * t230 + (-pkin(4) * t163 + t241) * t160) * qJD(5) + (t159 * t181 - t132) * qJD(4) (t144 * t233 + t159 * t180) * qJD(5) + (t163 * t181 + t206) * qJD(4), t115 * t93 - t131 * t66, t114 * t93 + t115 * t94 + t130 * t66 - t131 * t67, t83, -t171, 0, t101 * t148 + t114 * t211 + t118 * t94 + t130 * t97 + t147 * t67 - t164 * t71, -t102 * t148 - t115 * t211 - t118 * t93 + t131 * t97 - t147 * t66 - t164 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t221, 0, 0, 0, 0, 0, -t120, t122, 0, 0, 0, 0, 0, t171, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t191, -0.2e1 * t184, 0, 0, 0, t159 * t214, t163 * t214, -0.2e1 * t131 * t93, 0.2e1 * t130 * t93 - 0.2e1 * t131 * t94, 0, 0, 0, 0.2e1 * t130 * t211 + 0.2e1 * t147 * t94, 0.2e1 * t131 * t211 - 0.2e1 * t147 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t74, t9, t8, 0, 0, t13, -t14, t74, t74 * t251 + (-t247 + (-t16 - t253) * t158) * qJD(6) + t204 (-t158 * t74 - t215 * t95) * pkin(5) + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, -t168, t148, t64, t63, 0, 0, -t66, -t67, t148, t148 * t251 + (-t246 + (pkin(5) * t164 - t86) * t158) * qJD(6) + t186 (-t148 * t158 + t164 * t215) * pkin(5) + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t256, 0, 0, 0, 0, 0, -t67, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t219, 0, -pkin(10) * t218, pkin(10) * t219, 0, 0, -t93, -t94, 0, t71, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t210, -0.2e1 * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, t74, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t67, t148, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t94, 0, t71, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, -t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
