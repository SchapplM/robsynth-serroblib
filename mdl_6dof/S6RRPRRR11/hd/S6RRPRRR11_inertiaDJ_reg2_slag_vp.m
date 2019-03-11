% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:50
% EndTime: 2019-03-09 14:33:09
% DurationCPUTime: 7.86s
% Computational Cost: add. (8386->377), mult. (17310->673), div. (0->0), fcn. (15340->8), ass. (0->208)
t118 = sin(qJ(6));
t244 = cos(qJ(6));
t121 = cos(qJ(4));
t119 = sin(qJ(4));
t245 = cos(qJ(5));
t209 = t245 * t119;
t243 = sin(qJ(5));
t86 = t121 * t243 + t209;
t207 = t243 * t119;
t208 = t245 * t121;
t87 = t208 - t207;
t165 = t118 * t86 - t244 * t87;
t166 = t118 * t87 + t244 * t86;
t194 = t243 * qJD(5);
t154 = qJD(4) * t243 + t194;
t195 = t245 * qJD(5);
t155 = qJD(4) * t245 + t195;
t236 = t155 * t121;
t269 = t119 * t154 - t236;
t150 = t86 * qJD(5);
t62 = qJD(4) * t86 + t150;
t17 = qJD(6) * t165 + t118 * t62 + t244 * t269;
t197 = qJD(6) * t244;
t227 = qJD(6) * t118;
t19 = t118 * t269 - t197 * t86 - t227 * t87 - t244 * t62;
t291 = (t118 * t17 - t244 * t19 - (t118 * t165 + t166 * t244) * qJD(6)) * pkin(5);
t259 = t166 * t17;
t271 = t165 * t19;
t266 = -t259 - t271;
t122 = cos(qJ(2));
t120 = sin(qJ(2));
t199 = qJD(2) * t245;
t181 = t120 * t199;
t198 = qJD(2) * t243;
t182 = t120 * t198;
t92 = t119 * t182;
t160 = t121 * t181 - t92;
t129 = t122 * t62 + t160;
t140 = t122 * t154 + t181;
t37 = t140 * t119 + t121 * t182 - t122 * t236;
t96 = t122 * t208;
t171 = t122 * t207 - t96;
t70 = t86 * t122;
t44 = t118 * t171 - t244 * t70;
t14 = qJD(6) * t44 + t118 * t37 - t129 * t244;
t153 = t244 * t171;
t43 = -t118 * t70 - t153;
t290 = t14 * t166 - t17 * t43;
t13 = -qJD(6) * t153 - t118 * t129 - t227 * t70 - t244 * t37;
t289 = t13 * t165 + t19 * t44;
t109 = qJD(2) * t122;
t288 = -t109 * t166 + t120 * t17;
t287 = -t109 * t165 + t120 * t19;
t220 = pkin(5) * t109;
t247 = pkin(3) + pkin(7);
t163 = t121 * qJ(3) - t119 * t247;
t156 = t163 * t120;
t164 = t119 * qJ(3) + t121 * t247;
t158 = pkin(4) + t164;
t248 = pkin(2) + pkin(8);
t225 = pkin(9) + t248;
t173 = t122 * t225 + pkin(1);
t226 = t120 * qJD(3);
t204 = t119 * t226;
t255 = t120 * t225;
t126 = t204 + (t121 * t173 + t156) * qJD(4) + (-t119 * t255 + t122 * t158) * qJD(2);
t200 = t121 * t226;
t212 = t247 * t120;
t238 = t120 * qJ(3);
t254 = t163 * t122;
t127 = -t200 + (t121 * t212 + (t173 + t238) * t119) * qJD(4) + (t121 * t255 - t254) * qJD(2);
t136 = t119 * t173 + t120 * t158;
t41 = t243 * t136;
t237 = t121 * t122;
t214 = t248 * t122;
t162 = -t214 - t238;
t61 = t121 * (-pkin(1) + t162) + t119 * t212;
t48 = -pkin(9) * t237 + t61;
t9 = -qJD(5) * t41 + t126 * t245 - t127 * t243 - t195 * t48;
t123 = -t37 * pkin(10) + t220 + t9;
t42 = t245 * t136;
t8 = -qJD(5) * t42 - t126 * t243 - t127 * t245 + t194 * t48;
t125 = pkin(10) * t129 - t8;
t28 = -t243 * t48 + t42;
t151 = pkin(5) * t120 + pkin(10) * t70 + t28;
t143 = t244 * t151;
t29 = t245 * t48 + t41;
t16 = pkin(10) * t171 + t29;
t1 = -qJD(6) * t143 - t118 * t123 - t125 * t244 + t16 * t227;
t11 = -t118 * t16 + t143;
t144 = t118 * t151;
t12 = t16 * t244 + t144;
t273 = -t118 * t125 + t123 * t244;
t2 = -qJD(6) * t12 + t273;
t286 = -t1 * t166 + t11 * t19 - t12 * t17 - t165 * t2;
t178 = t225 * t245;
t169 = t121 * t178;
t177 = t225 * t243;
t63 = t119 * t177 - t169;
t139 = -t87 * pkin(10) + t63;
t134 = t244 * t139;
t83 = t121 * t177;
t64 = -t119 * t178 - t83;
t45 = -pkin(10) * t86 + t64;
t25 = -t118 * t45 + t134;
t135 = t118 * t139;
t26 = t244 * t45 + t135;
t167 = qJD(5) * t177;
t168 = qJD(5) * t178;
t230 = qJD(4) * t119;
t176 = t225 * t230;
t33 = qJD(4) * t83 + t119 * t168 + t121 * t167 + t176 * t245;
t128 = t62 * pkin(10) + t33;
t32 = qJD(4) * t169 - t119 * t167 + t121 * t168 - t176 * t243;
t130 = -pkin(10) * t269 + t32;
t6 = -qJD(6) * t134 - t118 * t128 + t130 * t244 + t227 * t45;
t7 = -qJD(6) * t135 + t118 * t130 + t128 * t244 - t197 * t45;
t285 = t165 * t7 + t166 * t6 + t17 * t26 - t19 * t25;
t222 = t245 * pkin(4);
t184 = t222 + pkin(5);
t159 = t244 * t184;
t190 = pkin(4) * t195;
t49 = -qJD(6) * t159 - t244 * t190 + (qJD(6) * t243 + t194) * t118 * pkin(4);
t183 = t244 * t243;
t131 = (qJD(5) + qJD(6)) * (-t118 * t245 - t183) * pkin(4);
t219 = pkin(5) * t227;
t50 = t131 - t219;
t221 = t243 * pkin(4);
t76 = -t118 * t221 + t159;
t77 = pkin(4) * t183 + t118 * t184;
t284 = t165 * t50 + t166 * t49 + t17 * t77 - t19 * t76;
t242 = t87 * t62;
t281 = t269 * t86;
t282 = 0.2e1 * t242 + 0.2e1 * t281;
t280 = t171 * t269;
t279 = t269 * t64 + t32 * t86 - t33 * t87 + t63 * t62;
t278 = -t269 * t29 - t28 * t62 - t8 * t86 + t9 * t87;
t276 = -t109 * t86 + t120 * t269;
t189 = pkin(4) * t194;
t275 = t189 * t87 - t190 * t86 + t221 * t269 + t222 * t62;
t148 = t121 * t154;
t213 = t248 * t120;
t117 = t122 ^ 2;
t192 = qJD(2) * (t120 ^ 2 - t117);
t114 = t119 ^ 2;
t116 = t121 ^ 2;
t235 = t114 - t116;
t191 = t235 * qJD(4);
t180 = pkin(1) + t214;
t60 = t119 * t180 + t120 * t164;
t231 = qJD(3) * t122;
t95 = t247 * t122;
t252 = qJD(2) * t162 + t95 * qJD(4) + t231;
t250 = t122 * (-qJD(4) - qJD(5));
t249 = 0.2e1 * qJD(3);
t246 = pkin(5) * t86;
t110 = t119 * pkin(4);
t240 = t62 * t120;
t239 = qJ(3) * t122;
t104 = qJ(3) + t110;
t233 = qJD(2) * t120;
t232 = qJD(2) * t121;
t229 = qJD(4) * t121;
t228 = qJD(4) * t122;
t100 = pkin(4) * t229 + qJD(3);
t224 = -0.2e1 * pkin(1) * qJD(2);
t223 = t8 * t243;
t218 = pkin(7) * t233;
t217 = pkin(7) * t109;
t216 = t119 * t248;
t215 = t121 * t248;
t78 = pkin(4) * t237 + t95;
t210 = qJD(4) * t248;
t206 = t119 * t228;
t205 = t121 * t228;
t203 = t120 * t109;
t202 = t120 * t232;
t201 = t119 * t229;
t193 = -0.2e1 * t219;
t188 = pkin(5) * t197;
t186 = t119 * t202;
t185 = t117 * t201;
t179 = t37 * t87 + t62 * t70;
t175 = -pkin(2) * t122 - t238;
t174 = -t119 * t60 + t121 * t61;
t89 = t247 * t233;
t146 = -t89 + (-t239 + t213) * qJD(4);
t145 = qJD(2) * t175 + t231;
t124 = -qJD(6) * t144 - t16 * t197 + t273;
t113 = qJ(3) * t249;
t99 = -0.2e1 * t203;
t98 = 0.2e1 * t203;
t91 = -pkin(1) + t175;
t85 = -0.2e1 * t192;
t81 = -t109 * t119 - t120 * t229;
t80 = t109 * t121 - t120 * t230;
t73 = -t226 + (pkin(2) * t120 - t239) * qJD(2);
t69 = t104 + t246;
t66 = -t122 * t191 - t186;
t65 = -pkin(4) * t206 + (-pkin(4) * t121 - t247) * t233;
t51 = -pkin(5) * t171 + t78;
t47 = -pkin(5) * t269 + t100;
t38 = t109 * t87 - t240;
t31 = t204 + (t121 * t180 + t156) * qJD(4) + (-t119 * t213 + t122 * t164) * qJD(2);
t30 = t200 - t60 * qJD(4) + (-t121 * t213 + t254) * qJD(2);
t27 = t92 * pkin(5) + ((-pkin(5) * t245 - pkin(4)) * t121 - t247) * t233 + (-pkin(5) * t150 + (-t110 - t246) * qJD(4)) * t122;
t10 = qJD(4) * t174 - t30 * t119 + t31 * t121;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t85, 0, t99, 0, 0, t120 * t224, t122 * t224, 0, 0, 0, 0, 0, t98, t85, t99, 0, 0.2e1 * t122 * t73 - 0.2e1 * t233 * t91, -0.2e1 * t109 * t91 - 0.2e1 * t120 * t73, 0.2e1 * t91 * t73, -0.2e1 * t114 * t203 + 0.2e1 * t185, -0.2e1 * t117 * t191 - 0.4e1 * t122 * t186, 0.2e1 * t119 * t192 - 0.2e1 * t120 * t205, -0.2e1 * t116 * t203 - 0.2e1 * t185, 0.2e1 * t120 * t206 + 0.2e1 * t121 * t192, t98, 0.2e1 * (-t232 * t95 + t31) * t120 + 0.2e1 * (qJD(2) * t60 - t89 * t121 - t230 * t95) * t122, 0.2e1 * (qJD(2) * t119 * t95 + t30) * t120 + 0.2e1 * (-qJD(2) * t61 + t89 * t119 - t229 * t95) * t122, 0.2e1 * t174 * t233 + 0.2e1 * (t119 * t31 + t121 * t30 + (t119 * t61 + t121 * t60) * qJD(4)) * t122, -0.2e1 * t30 * t61 + 0.2e1 * t31 * t60 - 0.2e1 * t89 * t95, -0.2e1 * t70 * t37, -0.2e1 * t37 * t96 - 0.2e1 * t70 * t160 + 0.2e1 * (-t70 * t148 + (-t155 * t70 + t243 * t37) * t119) * t122, -0.2e1 * t109 * t70 + 0.2e1 * t120 * t37, 0.2e1 * t171 * t129, 0.2e1 * t160 * t120 + 0.2e1 * (qJD(2) * t171 + t240) * t122, t98, 0.2e1 * t9 * t120 + 0.2e1 * t65 * t96 - 0.2e1 * t78 * t160 + 0.2e1 * (t28 * qJD(2) - t78 * t148 + (-t155 * t78 - t243 * t65) * t119) * t122, -0.2e1 * t109 * t29 + 0.2e1 * t120 * t8 + 0.2e1 * t37 * t78 - 0.2e1 * t65 * t70, 0.2e1 * t8 * t96 + 0.2e1 * t29 * t160 + 0.2e1 * t9 * t70 - 0.2e1 * t28 * t37 + 0.2e1 * (t29 * t148 + (t155 * t29 - t223) * t119) * t122, 0.2e1 * t28 * t9 - 0.2e1 * t29 * t8 + 0.2e1 * t65 * t78, -0.2e1 * t44 * t13, 0.2e1 * t13 * t43 - 0.2e1 * t14 * t44, 0.2e1 * t109 * t44 - 0.2e1 * t120 * t13, 0.2e1 * t43 * t14, -0.2e1 * t109 * t43 - 0.2e1 * t120 * t14, t98, 0.2e1 * t109 * t11 + 0.2e1 * t120 * t2 + 0.2e1 * t14 * t51 + 0.2e1 * t27 * t43, 0.2e1 * t1 * t120 - 0.2e1 * t109 * t12 - 0.2e1 * t13 * t51 + 0.2e1 * t27 * t44, 0.2e1 * t1 * t43 + 0.2e1 * t11 * t13 - 0.2e1 * t12 * t14 - 0.2e1 * t2 * t44, -0.2e1 * t1 * t12 + 0.2e1 * t11 * t2 + 0.2e1 * t27 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, -t233, 0, -t217, t218, 0, 0, 0, -t109, t233, 0, 0, 0, t145, t217, -t218, t145 * pkin(7), -t66, 0.4e1 * t122 * t201 - t233 * t235, t80, t66, t81, 0, t119 * t146 + t121 * t252, -t119 * t252 + t121 * t146, -t10, t30 * t216 - t31 * t215 - t89 * qJ(3) + t95 * qJD(3) + (-t215 * t61 + t216 * t60) * qJD(4), t179, t129 * t87 - t171 * t62 - t269 * t70 - t37 * t86, t38, t280 + (-t140 * t121 + t209 * t250 + t92) * t86, t276, 0, t33 * t120 + t63 * t109 - t100 * t171 + t104 * (t250 * t86 - t160) + t65 * t86 - t78 * t269, -t100 * t70 + t104 * t37 - t109 * t64 + t120 * t32 - t62 * t78 + t65 * t87, t129 * t64 - t171 * t32 + t33 * t70 - t63 * t37 - t278, t100 * t78 + t104 * t65 + t28 * t33 - t29 * t32 + t63 * t9 - t64 * t8, t289, t13 * t166 + t14 * t165 + t17 * t44 - t19 * t43, t287, t290, t288, 0, t109 * t25 + t120 * t7 + t14 * t69 + t166 * t27 - t17 * t51 + t43 * t47, -t109 * t26 + t120 * t6 - t13 * t69 - t165 * t27 + t19 * t51 + t44 * t47, t13 * t25 - t14 * t26 + t43 * t6 - t44 * t7 - t286, -t1 * t26 + t11 * t7 - t12 * t6 + t2 * t25 + t27 * t69 + t47 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t113, -0.2e1 * t201, 0.2e1 * t191, 0, 0.2e1 * t201, 0, 0, 0.2e1 * qJ(3) * t229 + 0.2e1 * qJD(3) * t119, -0.2e1 * qJ(3) * t230 + 0.2e1 * qJD(3) * t121, 0, t113, -0.2e1 * t242, 0.2e1 * t269 * t87 + 0.2e1 * t62 * t86, 0, -0.2e1 * t281, 0, 0, 0.2e1 * t100 * t86 - 0.2e1 * t104 * t269, 0.2e1 * t100 * t87 - 0.2e1 * t104 * t62, 0.2e1 * t279, 0.2e1 * t100 * t104 - 0.2e1 * t32 * t64 + 0.2e1 * t33 * t63, -0.2e1 * t271, -0.2e1 * t165 * t17 - 0.2e1 * t166 * t19, 0, -0.2e1 * t259, 0, 0, 0.2e1 * t166 * t47 - 0.2e1 * t17 * t69, -0.2e1 * t165 * t47 + 0.2e1 * t19 * t69, 0.2e1 * t285, 0.2e1 * t25 * t7 - 0.2e1 * t26 * t6 + 0.2e1 * t47 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, t217, 0, 0, 0, 0, 0, 0, t80, t81, 0, t10, 0, 0, 0, 0, 0, 0, t38, t276, t129 * t86 - t179 - t280, t278, 0, 0, 0, 0, 0, 0, t287, t288, -t289 - t290, t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, -t279, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t266, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t233 - t205, 0, t202 + t206, t109, t31, t30, 0, 0, 0, 0, t37, 0, t129, t109, pkin(4) * t122 * t199 - t120 * t189 + t9 (-t120 * t195 - t122 * t198) * pkin(4) + t8 (t129 * t243 + t171 * t195 - t194 * t70 - t245 * t37) * pkin(4) (t245 * t9 - t223 + (-t243 * t28 + t245 * t29) * qJD(5)) * pkin(4), 0, 0, -t13, 0, -t14, t109, t109 * t76 + t50 * t120 + t124, -t109 * t77 + t120 * t49 + t1, t13 * t76 - t14 * t77 + t43 * t49 - t44 * t50, -t1 * t77 + t11 * t50 - t12 * t49 + t2 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, 0, -t229, 0, t119 * t210, t121 * t210, 0, 0, 0, 0, -t62, 0, t269, 0, t33, t32, t275 (t245 * t33 - t243 * t32 + (-t243 * t63 + t245 * t64) * qJD(5)) * pkin(4), 0, 0, t19, 0, t17, 0, t7, t6, t284, t25 * t50 - t26 * t49 - t6 * t77 + t7 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, -t229, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t269, 0, -t275, 0, 0, 0, 0, 0, 0, t19, t17, 0, -t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t189, -0.2e1 * t190, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t50, 0.2e1 * t49, 0, -0.2e1 * t49 * t77 + 0.2e1 * t50 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t129, t109, t9, t8, 0, 0, 0, 0, -t13, 0, -t14, t109, -t120 * t219 + t220 * t244 + t124 (-t109 * t118 - t120 * t197) * pkin(5) + t1 (t244 * t13 - t118 * t14 + (t118 * t44 - t244 * t43) * qJD(6)) * pkin(5) (t244 * t2 - t1 * t118 + (-t11 * t118 + t12 * t244) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, t269, 0, t33, t32, 0, 0, 0, 0, t19, 0, t17, 0, t7, t6, t291 (t244 * t7 - t118 * t6 + (-t118 * t25 + t244 * t26) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t269, 0, 0, 0, 0, 0, 0, 0, 0, t19, t17, 0, -t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t190, 0, 0, 0, 0, 0, 0, 0, 0, t131 + t193, -t188 + t49, 0 (t244 * t50 - t118 * t49 + (-t118 * t76 + t244 * t77) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -0.2e1 * t188, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, -t14, t109, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t17, 0, t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, -t188, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
