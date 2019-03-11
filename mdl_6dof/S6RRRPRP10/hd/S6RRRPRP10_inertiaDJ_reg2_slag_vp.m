% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:36:56
% EndTime: 2019-03-09 17:37:19
% DurationCPUTime: 7.99s
% Computational Cost: add. (10620->543), mult. (28622->948), div. (0->0), fcn. (27708->10), ass. (0->243)
t150 = sin(pkin(11));
t152 = cos(pkin(11));
t151 = sin(pkin(6));
t154 = sin(qJ(2));
t273 = qJD(2) * t154;
t241 = t151 * t273;
t153 = sin(qJ(3));
t281 = t151 * t154;
t247 = qJD(3) * t281;
t122 = t153 * t247;
t155 = cos(qJ(3));
t285 = cos(pkin(6));
t230 = t285 * qJD(3);
t156 = cos(qJ(2));
t272 = qJD(2) * t156;
t240 = t151 * t272;
t302 = -t155 * (t230 + t240) + t122;
t164 = t150 * t241 - t152 * t302;
t307 = t164 * t150;
t162 = t164 * t153;
t271 = qJD(3) * t155;
t182 = t285 * t153 + t155 * t281;
t280 = t151 * t156;
t88 = -t150 * t280 + t152 * t182;
t306 = t88 * t271 + t162;
t226 = t152 * t241;
t83 = t302 * t150;
t197 = t83 + t226;
t276 = t182 * t150 + t152 * t280;
t305 = t197 * t153 - t276 * t271;
t294 = sin(qJ(5));
t233 = qJD(5) * t294;
t295 = cos(qJ(5));
t234 = qJD(5) * t295;
t108 = t150 * t233 - t152 * t234;
t252 = pkin(1) * t285;
t304 = pkin(8) * t281 - t156 * t252;
t110 = t153 * t281 - t285 * t155;
t143 = qJD(3) * t153;
t179 = t182 * qJD(3);
t89 = t153 * t240 + t179;
t60 = t110 * t143 - t155 * t89;
t117 = t294 * t150 - t295 * t152;
t149 = t155 ^ 2;
t303 = (t153 ^ 2 - t149) * qJD(3);
t106 = t304 * qJD(2);
t278 = t155 * t106;
t102 = (-pkin(2) * t156 - pkin(9) * t154 - pkin(1)) * t151;
t138 = t154 * t252;
t112 = pkin(8) * t280 + t138;
t196 = t285 * pkin(9) + t112;
t68 = t155 * t102 - t153 * t196;
t301 = qJD(3) * t68 - t278;
t147 = t152 ^ 2;
t300 = -0.2e1 * t108;
t251 = t295 * t150;
t118 = t294 * t152 + t251;
t109 = t118 * qJD(5);
t299 = 0.2e1 * t109;
t298 = 0.2e1 * t151;
t297 = 2 * qJD(6);
t296 = t89 * pkin(5);
t293 = pkin(3) * t153;
t292 = pkin(9) * t151;
t291 = pkin(9) * t155;
t144 = t153 * pkin(9);
t290 = qJ(4) + pkin(10);
t101 = -t285 * pkin(2) + t304;
t57 = t110 * pkin(3) - qJ(4) * t182 + t101;
t69 = t153 * t102 + t155 * t196;
t58 = -qJ(4) * t280 + t69;
t34 = -t150 * t58 + t152 * t57;
t177 = pkin(4) * t110 - pkin(10) * t88 + t34;
t24 = t294 * t177;
t35 = t150 * t57 + t152 * t58;
t31 = -t276 * pkin(10) + t35;
t10 = t295 * t31 + t24;
t209 = -pkin(3) * t155 - qJ(4) * t153;
t121 = -pkin(2) + t209;
t253 = pkin(9) * t150 + pkin(4);
t172 = -t253 * t155 + (-pkin(10) * t153 + t121) * t152;
t79 = t294 * t172;
t283 = t150 * t153;
t136 = t152 * t291;
t99 = t150 * t121 + t136;
t90 = -pkin(10) * t283 + t99;
t49 = t295 * t90 + t79;
t289 = t152 * t89;
t218 = pkin(2) * t154 - pkin(9) * t156;
t184 = qJD(3) * t196;
t256 = t102 * t143 - t153 * t106 + t155 * t184;
t274 = qJD(2) * t151;
t37 = (-t154 * pkin(3) - t155 * t218) * t274 + t256;
t287 = t37 * t150;
t286 = t89 * qJ(6);
t284 = qJ(4) * t155;
t282 = t150 * t155;
t279 = t152 * t153;
t277 = t155 * t156;
t142 = pkin(9) * t271;
t243 = t150 * t271;
t113 = pkin(4) * t243 + t142;
t119 = pkin(4) * t283 + t144;
t270 = qJD(3) * t156;
t269 = qJD(4) * t152;
t268 = qJD(4) * t153;
t267 = qJD(4) * t155;
t266 = t101 * qJD(3);
t265 = t110 * qJD(6);
t264 = t155 * qJD(6);
t263 = t156 * qJD(4);
t215 = t294 * t276;
t27 = -qJD(5) * t215 + t294 * t164 - t295 * t197 + t88 * t234;
t216 = t295 * t276;
t50 = t294 * t88 + t216;
t262 = 0.2e1 * t50 * t27;
t261 = -0.2e1 * pkin(2) * qJD(3);
t104 = t118 * t153;
t75 = -t108 * t153 + t118 * t271;
t260 = 0.2e1 * t104 * t75;
t70 = 0.2e1 * t110 * t89;
t259 = t117 * t299;
t258 = pkin(9) * t282;
t257 = pkin(9) * t279;
t255 = pkin(5) * t143;
t141 = -t152 * pkin(4) - pkin(3);
t146 = t151 ^ 2;
t248 = t146 * t272;
t246 = qJ(6) * t143;
t244 = t153 * t270;
t242 = t150 * t268;
t239 = t152 * t271;
t238 = t152 * t268;
t237 = t153 * t271;
t125 = t290 * t152;
t198 = t290 * t251;
t232 = t295 * qJD(4);
t235 = qJD(4) * t294;
t71 = qJD(5) * t198 + t125 * t233 + t150 * t235 - t152 * t232;
t220 = t290 * t294;
t72 = t125 * t234 + t152 * t235 + (-qJD(5) * t220 + t232) * t150;
t91 = t294 * t125 + t198;
t92 = t295 * t125 - t150 * t220;
t236 = -t92 * t71 + t72 * t91;
t231 = qJD(2) * t285;
t229 = 0.2e1 * t237;
t227 = t154 * t248;
t225 = t150 * t239;
t145 = t150 ^ 2;
t222 = 0.2e1 * (t145 + t147) * qJD(4);
t221 = -0.2e1 * t303;
t59 = pkin(3) * t280 - t68;
t219 = t154 * t231;
t26 = qJD(5) * t216 - t295 * t164 - t294 * t197 + t88 * t233;
t51 = t295 * t88 - t215;
t217 = t26 * t50 - t27 * t51;
t213 = t104 * t27 + t50 * t75;
t212 = t110 * t27 + t50 * t89;
t211 = t110 * t71 - t89 * t92;
t210 = -t110 * t72 - t89 * t91;
t208 = -t284 + t293;
t105 = t117 * t153;
t74 = t153 * t109 - t295 * t239 + t294 * t243;
t207 = t104 * t74 + t105 * t75;
t206 = t109 * t50 + t117 * t27;
t167 = t89 * pkin(3) - qJD(4) * t182;
t160 = -(t155 * t230 - t122) * qJ(4) + t167;
t163 = t102 * t271 - t151 * t263 - t153 * t184 - t278;
t178 = t154 * qJ(4) + t153 * t218;
t195 = pkin(8) * t156 - qJ(4) * t277;
t19 = -t150 * t163 + t152 * t160 + (t152 * t138 + (-t150 * t178 + t152 * t195) * t151) * qJD(2);
t159 = t150 * t160 + t152 * t163;
t200 = t150 * t138;
t20 = (t200 + (t150 * t195 + t152 * t178) * t151) * qJD(2) + t159;
t205 = -t19 * t150 + t20 * t152;
t86 = -t238 + (pkin(9) * t283 + t152 * t208) * qJD(3);
t87 = -t242 + (t150 * t208 - t257) * qJD(3);
t204 = -t150 * t86 + t152 * t87;
t185 = t218 * t274;
t39 = -t153 * t185 - t301;
t40 = t155 * t185 - t256;
t203 = -t40 * t153 - t39 * t155;
t202 = t104 * t109 + t117 * t75;
t201 = t109 * t110 + t117 * t89;
t199 = t108 * t117 - t109 * t118;
t194 = -t92 * t143 - t155 * t71;
t193 = -t91 * t143 + t155 * t72;
t192 = -t290 * t155 + t293;
t191 = t104 * t143 - t155 * t75;
t190 = -t109 * t155 + t117 * t143;
t189 = -t26 * t91 - t92 * t27 + t71 * t50 + t51 * t72;
t188 = t71 * t104 - t105 * t72 - t74 * t91 - t92 * t75;
t187 = t197 * t152;
t157 = -t150 * ((qJD(2) * t178 - t263) * t151 + t301) + t152 * (pkin(1) * t219 + pkin(8) * t240 + qJ(4) * t302 + t167) - t164 * pkin(10) + t89 * pkin(4);
t158 = t83 * pkin(10) + (t200 + (t150 * (pkin(8) - t284) * t156 + (-t156 * t144 + (t153 * pkin(2) + t290) * t154) * t152) * t151) * qJD(2) + t159;
t168 = t295 * t177;
t3 = -qJD(5) * t168 - t294 * t157 - t295 * t158 + t233 * t31;
t161 = -t238 + (t192 * t152 + t253 * t153) * qJD(3);
t165 = -t242 + (t150 * t192 - t257) * qJD(3);
t166 = t295 * t172;
t28 = -qJD(5) * t166 - t294 * t161 - t295 * t165 + t233 * t90;
t183 = t153 * t273 - t155 * t270;
t41 = t276 * pkin(4) + t59;
t180 = -0.2e1 * t108 * t91 - 0.2e1 * t92 * t109 + 0.2e1 * t71 * t117 + 0.2e1 * t118 * t72;
t176 = t104 * t26 + t105 * t27 + t50 * t74 - t51 * t75;
t175 = t108 * t50 - t109 * t51 + t117 * t26 - t118 * t27;
t174 = t104 * t108 + t105 * t109 + t117 * t74 - t118 * t75;
t171 = t104 * t89 + t110 * t75 + t50 * t143 - t155 * t27;
t9 = -t294 * t31 + t168;
t48 = -t294 * t90 + t166;
t32 = -t83 * pkin(4) + (pkin(9) * t277 + (-t155 * pkin(2) + t141) * t154) * t274 + t256;
t29 = -qJD(5) * t79 + t295 * t161 - t294 * t165 - t90 * t234;
t4 = -qJD(5) * t24 + t295 * t157 - t294 * t158 - t31 * t234;
t130 = -0.2e1 * t237;
t120 = -0.2e1 * t227;
t107 = t112 * qJD(2);
t98 = t121 * t152 - t258;
t85 = t118 * t300;
t80 = pkin(5) * t117 - qJ(6) * t118 + t141;
t77 = t108 * t155 + t118 * t143;
t62 = pkin(5) * t109 + qJ(6) * t108 - qJD(6) * t118;
t61 = pkin(5) * t104 + qJ(6) * t105 + t119;
t53 = 0.2e1 * t105 * t74;
t46 = -0.2e1 * t105 * t143 + 0.2e1 * t155 * t74;
t45 = t155 * pkin(5) - t48;
t43 = -qJ(6) * t155 + t49;
t42 = -t108 * t110 + t118 * t89;
t38 = t105 * t108 - t118 * t74;
t33 = pkin(5) * t75 + qJ(6) * t74 + qJD(6) * t105 + t113;
t25 = -t29 - t255;
t22 = t246 - t28 - t264;
t18 = t50 * pkin(5) - t51 * qJ(6) + t41;
t17 = -0.2e1 * t51 * t26;
t16 = -t108 * t51 - t118 * t26;
t14 = -0.2e1 * t110 * t26 + 0.2e1 * t51 * t89;
t13 = t105 * t26 - t51 * t74;
t8 = -t110 * pkin(5) - t9;
t7 = qJ(6) * t110 + t10;
t6 = -t105 * t89 - t110 * t74 + t51 * t143 + t155 * t26;
t5 = t27 * pkin(5) + t26 * qJ(6) - t51 * qJD(6) + t32;
t2 = -t4 - t296;
t1 = t265 - t3 + t286;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t227, 0.2e1 * (-t154 ^ 2 + t156 ^ 2) * t146 * qJD(2), 0.2e1 * t231 * t280, t120, -0.2e1 * t151 * t219, 0, -0.2e1 * t146 * pkin(1) * t273 - 0.2e1 * t107 * t285, -0.2e1 * pkin(1) * t248 + 0.2e1 * t106 * t285 (-t106 * t156 + t107 * t154 + (-t112 * t154 + t156 * t304) * qJD(2)) * t298, -0.2e1 * t106 * t112 + 0.2e1 * t107 * t304, -0.2e1 * t182 * t302, 0.2e1 * t110 * t302 - 0.2e1 * t182 * t89, 0.2e1 * t182 * t241 + 0.2e1 * t280 * t302, t70 (-t110 * t273 + t156 * t89) * t298, t120, 0.2e1 * t101 * t89 + 0.2e1 * t107 * t110 + 0.2e1 * (-t156 * t40 + t273 * t68) * t151, -0.2e1 * t101 * t302 + 0.2e1 * t107 * t182 - 0.2e1 * t241 * t69 - 0.2e1 * t280 * t39, 0.2e1 * t39 * t110 - 0.2e1 * t182 * t40 + 0.2e1 * t302 * t68 - 0.2e1 * t69 * t89, 0.2e1 * t101 * t107 - 0.2e1 * t39 * t69 + 0.2e1 * t40 * t68, 0.2e1 * t88 * t164, -0.2e1 * t164 * t276 + 0.2e1 * t197 * t88, 0.2e1 * t110 * t164 + 0.2e1 * t88 * t89, -0.2e1 * t276 * t197, 0.2e1 * t110 * t197 - 0.2e1 * t276 * t89, t70, 0.2e1 * t19 * t110 - 0.2e1 * t197 * t59 + 0.2e1 * t276 * t37 + 0.2e1 * t34 * t89, -0.2e1 * t20 * t110 + 0.2e1 * t164 * t59 - 0.2e1 * t35 * t89 + 0.2e1 * t37 * t88, -0.2e1 * t164 * t34 - 0.2e1 * t19 * t88 + 0.2e1 * t197 * t35 - 0.2e1 * t20 * t276, 0.2e1 * t19 * t34 + 0.2e1 * t20 * t35 + 0.2e1 * t37 * t59, t17, 0.2e1 * t217, t14, t262, -0.2e1 * t212, t70, 0.2e1 * t110 * t4 + 0.2e1 * t27 * t41 + 0.2e1 * t32 * t50 + 0.2e1 * t89 * t9, -0.2e1 * t10 * t89 + 0.2e1 * t110 * t3 - 0.2e1 * t26 * t41 + 0.2e1 * t32 * t51, -0.2e1 * t10 * t27 + 0.2e1 * t26 * t9 + 0.2e1 * t3 * t50 - 0.2e1 * t4 * t51, -0.2e1 * t10 * t3 + 0.2e1 * t32 * t41 + 0.2e1 * t4 * t9, t17, t14, -0.2e1 * t217, t70, 0.2e1 * t212, t262, -0.2e1 * t110 * t2 + 0.2e1 * t18 * t27 + 0.2e1 * t5 * t50 - 0.2e1 * t8 * t89, -0.2e1 * t1 * t50 + 0.2e1 * t2 * t51 - 0.2e1 * t26 * t8 - 0.2e1 * t27 * t7, 0.2e1 * t1 * t110 + 0.2e1 * t18 * t26 - 0.2e1 * t5 * t51 + 0.2e1 * t7 * t89, 0.2e1 * t1 * t7 + 0.2e1 * t18 * t5 + 0.2e1 * t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, -t241, 0, -t107, t106, 0, 0, t149 * t247 + (-t122 + (0.2e1 * t230 + t240) * t155) * t153, -t110 * t271 - t155 * t302 + (-t179 - t89) * t153, t183 * t151, t60 (t155 * t273 + t244) * t151, 0, -pkin(2) * t89 - t107 * t155 + t153 * t266 - t183 * t292, pkin(2) * t302 + t107 * t153 + t155 * t266 - t241 * t291 - t244 * t292, t60 * pkin(9) + t142 * t182 - t143 * t69 - t144 * t302 - t271 * t68 + t203, -pkin(2) * t107 + ((-t69 * t153 - t68 * t155) * qJD(3) + t203) * pkin(9), t306 * t152, -t150 * t162 + t305 * t152 - t243 * t88, t110 * t239 + t143 * t88 - t155 * t164 + t279 * t89, -t305 * t150, -t89 * t283 - t197 * t155 + (-t110 * t282 - t153 * t276) * qJD(3), t60, t86 * t110 - t19 * t155 + t98 * t89 + (-pkin(9) * t197 + t287) * t153 + (t34 * t153 + (pkin(9) * t276 + t59 * t150) * t155) * qJD(3), t306 * pkin(9) - t87 * t110 - t143 * t35 + t20 * t155 + t239 * t59 + t279 * t37 - t99 * t89, -t164 * t98 - t19 * t279 + t197 * t99 - t20 * t283 - t239 * t34 - t243 * t35 - t276 * t87 - t86 * t88, t19 * t98 + t20 * t99 + t34 * t86 + t35 * t87 + (t37 * t153 + t271 * t59) * pkin(9), t13, t176, t6, t213, -t171, t60, t104 * t32 + t110 * t29 + t113 * t50 + t119 * t27 + t143 * t9 - t155 * t4 + t41 * t75 + t48 * t89, -t10 * t143 - t105 * t32 + t110 * t28 + t113 * t51 - t119 * t26 - t155 * t3 - t41 * t74 - t49 * t89, -t10 * t75 + t104 * t3 + t105 * t4 + t26 * t48 - t27 * t49 + t28 * t50 - t29 * t51 + t74 * t9, -t10 * t28 + t113 * t41 + t119 * t32 + t29 * t9 - t3 * t49 + t4 * t48, t13, t6, -t176, t60, t171, t213, t104 * t5 - t110 * t25 - t143 * t8 + t155 * t2 + t18 * t75 + t27 * t61 + t33 * t50 - t45 * t89, -t1 * t104 - t105 * t2 - t22 * t50 + t25 * t51 - t26 * t45 - t27 * t43 - t7 * t75 - t74 * t8, -t1 * t155 + t105 * t5 + t110 * t22 + t143 * t7 + t18 * t74 + t26 * t61 - t33 * t51 + t43 * t89, t1 * t43 + t18 * t33 + t2 * t45 + t22 * t7 + t25 * t8 + t5 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t221, 0, t130, 0, 0, t153 * t261, t155 * t261, 0, 0, t147 * t229, -0.4e1 * t153 * t225, 0.2e1 * t152 * t303, t145 * t229, t150 * t221, t130, -0.2e1 * t155 * t86 + 0.2e1 * (t98 + 0.2e1 * t258) * t143, 0.2e1 * t155 * t87 + 0.2e1 * (-t99 + 0.2e1 * t136) * t143, 0.2e1 * (-t150 * t87 - t152 * t86) * t153 + 0.2e1 * (-t150 * t99 - t152 * t98) * t271, 0.2e1 * pkin(9) ^ 2 * t237 + 0.2e1 * t86 * t98 + 0.2e1 * t87 * t99, t53, 0.2e1 * t207, t46, t260, -0.2e1 * t191, t130, 0.2e1 * t104 * t113 + 0.2e1 * t119 * t75 + 0.2e1 * t143 * t48 - 0.2e1 * t155 * t29, -0.2e1 * t105 * t113 - 0.2e1 * t119 * t74 - 0.2e1 * t143 * t49 - 0.2e1 * t155 * t28, 0.2e1 * t104 * t28 + 0.2e1 * t105 * t29 + 0.2e1 * t48 * t74 - 0.2e1 * t49 * t75, 0.2e1 * t113 * t119 - 0.2e1 * t28 * t49 + 0.2e1 * t29 * t48, t53, t46, -0.2e1 * t207, t130, 0.2e1 * t191, t260, 0.2e1 * t104 * t33 - 0.2e1 * t143 * t45 + 0.2e1 * t155 * t25 + 0.2e1 * t61 * t75, -0.2e1 * t104 * t22 - 0.2e1 * t105 * t25 - 0.2e1 * t43 * t75 - 0.2e1 * t45 * t74, 0.2e1 * t105 * t33 + 0.2e1 * t143 * t43 - 0.2e1 * t155 * t22 + 0.2e1 * t61 * t74, 0.2e1 * t22 * t43 + 0.2e1 * t25 * t45 + 0.2e1 * t33 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t302, 0, -t89, t241, t40, t39, 0, 0, t307, -t302 * t147 + (t83 + 0.2e1 * t226) * t150, t150 * t89, t187, t289, 0, pkin(3) * t197 - t37 * t152 + (-qJ(4) * t89 - qJD(4) * t110) * t150, -pkin(3) * t164 - qJ(4) * t289 - t110 * t269 + t287, t150 * qJD(4) * t88 - t269 * t276 + t205 + (t187 + t307) * qJ(4), -pkin(3) * t37 + (-t150 * t34 + t152 * t35) * qJD(4) + t205 * qJ(4), t16, t175, t42, t206, -t201, 0, t109 * t41 + t117 * t32 + t141 * t27 + t210, -t108 * t41 + t118 * t32 - t141 * t26 + t211, -t10 * t109 + t108 * t9 + t117 * t3 - t118 * t4 + t189, -t10 * t71 + t141 * t32 - t3 * t92 - t4 * t91 - t72 * t9, t16, t42, -t175, 0, t201, t206, t109 * t18 + t117 * t5 + t27 * t80 + t50 * t62 + t210, -t1 * t117 - t108 * t8 - t109 * t7 + t118 * t2 + t189, t108 * t18 - t118 * t5 + t26 * t80 - t51 * t62 - t211, t1 * t92 + t18 * t62 + t2 * t91 + t5 * t80 - t7 * t71 + t72 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, 0, -t143, 0, -t142, pkin(9) * t143, 0, 0, t225 (-t145 + t147) * t271, t150 * t143, -t225, t152 * t143, 0, t150 * t267 + (t150 * t209 - t136) * qJD(3), t152 * t267 + (t152 * t209 + t258) * qJD(3), t204, -pkin(3) * t142 + (-t150 * t98 + t152 * t99) * qJD(4) + t204 * qJ(4), t38, t174, t77, t202, -t190, 0, t109 * t119 + t113 * t117 + t141 * t75 + t193, -t108 * t119 + t113 * t118 - t141 * t74 + t194, t108 * t48 - t109 * t49 + t117 * t28 - t118 * t29 + t188, t113 * t141 - t28 * t92 - t29 * t91 - t48 * t72 - t49 * t71, t38, t77, -t174, 0, t190, t202, t104 * t62 + t109 * t61 + t117 * t33 + t75 * t80 + t193, -t108 * t45 - t109 * t43 - t117 * t22 + t118 * t25 + t188, t105 * t62 + t108 * t61 - t118 * t33 + t74 * t80 - t194, t22 * t92 + t25 * t91 + t33 * t80 - t43 * t71 + t45 * t72 + t61 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, qJ(4) * t222, t85, 0.2e1 * t199, 0, t259, 0, 0, t141 * t299, t141 * t300, t180, 0.2e1 * t236, t85, 0, -0.2e1 * t199, 0, 0, t259, 0.2e1 * t109 * t80 + 0.2e1 * t117 * t62, t180, 0.2e1 * t108 * t80 - 0.2e1 * t118 * t62, 0.2e1 * t62 * t80 + 0.2e1 * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t164, 0, t37, 0, 0, 0, 0, 0, 0, t27, -t26, 0, t32, 0, 0, 0, 0, 0, 0, t27, 0, t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t239, 0, t142, 0, 0, 0, 0, 0, 0, t75, -t74, 0, t113, 0, 0, 0, 0, 0, 0, t75, 0, t74, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t108, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, t108, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, t89, t4, t3, 0, 0, 0, -t26, 0, t89, t27, 0, t4 + 0.2e1 * t296, pkin(5) * t26 - qJ(6) * t27 - qJD(6) * t50, 0.2e1 * t265 - t3 + 0.2e1 * t286, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, -t75, t143, t29, t28, 0, 0, 0, -t74, 0, t143, t75, 0, t29 + 0.2e1 * t255, pkin(5) * t74 - qJ(6) * t75 - qJD(6) * t104, 0.2e1 * t246 - t28 - 0.2e1 * t264, -pkin(5) * t25 + qJ(6) * t22 + qJD(6) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t109, 0, -t72, t71, 0, 0, 0, -t108, 0, 0, t109, 0, -t72, pkin(5) * t108 - qJ(6) * t109 - qJD(6) * t117, -t71, -pkin(5) * t72 - qJ(6) * t71 + qJD(6) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, qJ(6) * t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t74, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
