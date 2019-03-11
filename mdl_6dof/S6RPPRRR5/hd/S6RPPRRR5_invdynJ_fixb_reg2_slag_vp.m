% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:08
% EndTime: 2019-03-09 02:29:15
% DurationCPUTime: 4.22s
% Computational Cost: add. (5843->451), mult. (10698->540), div. (0->0), fcn. (6748->10), ass. (0->247)
t160 = sin(qJ(6));
t164 = cos(qJ(6));
t161 = sin(qJ(5));
t165 = cos(qJ(5));
t252 = qJD(1) * qJD(5);
t166 = cos(qJ(4));
t271 = t161 * t166;
t162 = sin(qJ(4));
t253 = qJD(1) * qJD(4);
t235 = t162 * t253;
t247 = t166 * qJDD(1);
t315 = t235 - t247;
t313 = -t162 * t252 - t315;
t234 = t166 * t253;
t248 = t162 * qJDD(1);
t319 = t234 + t248;
t179 = -t161 * t319 + t165 * t313 - t252 * t271;
t246 = qJD(4) + qJD(5);
t221 = t164 * t246;
t244 = qJDD(4) + qJDD(5);
t256 = qJD(6) * t160;
t261 = qJD(1) * t166;
t262 = qJD(1) * t162;
t83 = -t161 * t262 + t165 * t261;
t24 = -qJD(6) * t221 - t160 * t244 - t164 * t179 + t256 * t83;
t86 = t162 * t165 + t271;
t81 = t86 * qJD(1);
t317 = qJD(6) + t81;
t62 = t160 * t83 - t221;
t310 = t317 * t62;
t322 = -t24 + t310;
t133 = qJ(2) * qJD(1) + qJD(3);
t108 = -pkin(7) * qJD(1) + t133;
t72 = (-pkin(8) * qJD(1) + t108) * t162;
t281 = t165 * t72;
t73 = -pkin(8) * t261 + t166 * t108;
t70 = qJD(4) * pkin(4) + t73;
t48 = t161 * t70 + t281;
t40 = pkin(9) * t246 + t48;
t159 = pkin(1) + qJ(3);
t305 = qJD(1) * t159;
t109 = -qJD(2) + t305;
t85 = pkin(4) * t262 + t109;
t46 = t81 * pkin(5) - t83 * pkin(9) + t85;
t20 = -t160 * t40 + t164 * t46;
t321 = t20 * t317;
t21 = t160 * t46 + t164 * t40;
t320 = t21 * t317;
t163 = sin(qJ(1));
t167 = cos(qJ(1));
t265 = g(1) * t167 + g(2) * t163;
t318 = -t109 * qJD(1) - t265;
t225 = t160 * t317;
t190 = t313 * t161;
t220 = t246 * t166;
t42 = (qJD(1) * t220 + t248) * t165 + t190;
t38 = qJDD(6) + t42;
t284 = t164 * t38;
t314 = t225 * t317 - t284;
t311 = t246 * t81;
t154 = t162 ^ 2;
t155 = t166 ^ 2;
t264 = t154 + t155;
t309 = t108 * t264;
t156 = qJ(4) + qJ(5);
t137 = cos(t156);
t147 = g(1) * t163;
t306 = -g(2) * t167 + t147;
t308 = t306 * t137;
t307 = t319 * pkin(4);
t136 = sin(t156);
t230 = -pkin(5) * t136 + t137 * pkin(9);
t143 = t162 * pkin(4);
t218 = -t159 - t143;
t260 = qJD(4) * t162;
t151 = qJD(1) * qJD(2);
t152 = qJ(2) * qJDD(1);
t231 = qJDD(3) + t151 + t152;
t88 = -pkin(7) * qJDD(1) + t231;
t80 = t166 * t88;
t50 = qJDD(4) * pkin(4) + pkin(8) * t315 - t108 * t260 + t80;
t259 = qJD(4) * t166;
t53 = -pkin(8) * t319 + t108 * t259 + t162 * t88;
t229 = t161 * t53 - t165 * t50;
t15 = -qJD(5) * t48 - t229;
t158 = -pkin(7) + qJ(2);
t304 = qJD(4) * (qJD(2) + t109 + t305) + qJDD(4) * t158;
t258 = qJD(5) * t161;
t59 = -t161 * t260 - t162 * t258 + t165 * t220;
t303 = -t244 * t86 - t246 * t59;
t302 = t81 ^ 2;
t301 = t83 ^ 2;
t141 = 0.2e1 * t151;
t299 = pkin(5) * t137;
t123 = g(3) * t136;
t124 = g(3) * t137;
t298 = g(3) * t162;
t257 = qJD(5) * t165;
t14 = t161 * t50 + t165 * t53 + t70 * t257 - t72 * t258;
t12 = pkin(9) * t244 + t14;
t149 = qJDD(1) * qJ(3);
t150 = qJD(3) * qJD(1);
t157 = qJDD(1) * pkin(1);
t245 = t157 - qJDD(2);
t214 = -t149 - t150 - t245;
t66 = -t214 + t307;
t17 = pkin(5) * t42 - pkin(9) * t179 + t66;
t2 = qJD(6) * t20 + t164 * t12 + t160 * t17;
t1 = t2 * t164;
t285 = t161 * t72;
t47 = t165 * t70 - t285;
t39 = -pkin(5) * t246 - t47;
t297 = t39 * t81;
t296 = t62 * t83;
t64 = t160 * t246 + t164 * t83;
t295 = t64 * t62;
t294 = t64 * t83;
t293 = t317 * t83;
t292 = t83 * t81;
t291 = pkin(8) - t158;
t51 = t161 * t73 + t281;
t290 = t48 - t51;
t58 = t246 * t86;
t87 = -t161 * t162 + t165 * t166;
t289 = t87 * t244 - t58 * t246;
t288 = t160 * t21;
t287 = t160 * t38;
t286 = t160 * t62;
t283 = t164 * t64;
t282 = t164 * t81;
t280 = t24 * t160;
t226 = t160 * t179 - t164 * t244;
t25 = qJD(6) * t64 + t226;
t279 = t25 * t164;
t278 = qJD(1) * t85;
t277 = qJD(6) * t317;
t276 = t136 * t163;
t275 = t136 * t167;
t170 = qJD(1) ^ 2;
t274 = t154 * t170;
t273 = t160 * t163;
t272 = t160 * t167;
t270 = t163 * t164;
t269 = t164 * t167;
t267 = t167 * pkin(1) + t163 * qJ(2);
t169 = qJD(4) ^ 2;
t263 = -t169 - t170;
t255 = qJD(6) * t164;
t111 = pkin(4) * t259 + qJD(3);
t251 = qJDD(1) * t159;
t249 = qJDD(4) * t162;
t243 = pkin(4) * t257;
t126 = pkin(4) * t161 + pkin(9);
t242 = t126 * t277;
t241 = t87 * t256;
t240 = t87 * t255;
t239 = t166 * t170 * t162;
t35 = t39 * t256;
t36 = t39 * t255;
t238 = g(1) * t275 + g(2) * t276 + t124;
t237 = t167 * qJ(3) + t267;
t140 = t167 * qJ(2);
t168 = -pkin(8) - pkin(7);
t236 = g(1) * (t167 * t168 + t140);
t94 = t291 * t166;
t232 = qJDD(2) - t306;
t228 = t317 * t64;
t227 = t264 * t88;
t224 = t164 * t317;
t223 = t83 * t246;
t222 = qJD(6) * t86 + qJD(1);
t219 = 0.2e1 * t152 + t141 - t265;
t217 = t164 * t123 - t20 * t83 + t35;
t216 = t162 * t234;
t215 = -t157 + t232;
t213 = pkin(4) * t258 - t51;
t56 = pkin(5) * t83 + pkin(9) * t81;
t212 = t2 * t86 + t21 * t59;
t16 = t164 * t17;
t3 = -qJD(6) * t21 - t160 * t12 + t16;
t211 = -t20 * t59 - t3 * t86;
t13 = -pkin(5) * t244 - t15;
t207 = t13 * t87 - t39 * t58;
t206 = -t87 * t24 - t58 * t64;
t205 = -t24 * t86 + t64 * t59;
t204 = t87 * t25 - t58 * t62;
t203 = -t25 * t86 - t62 * t59;
t202 = -t317 * t58 + t38 * t87;
t201 = t179 * t87 - t58 * t83;
t200 = t42 * t86 + t59 * t81;
t198 = -t159 * t163 + t140;
t197 = t164 * t20 + t288;
t196 = t160 * t20 - t164 * t21;
t55 = pkin(5) * t86 - pkin(9) * t87 - t218;
t93 = t291 * t162;
t61 = -t161 * t94 - t165 * t93;
t30 = -t160 * t61 + t164 * t55;
t31 = t160 * t55 + t164 * t61;
t60 = -t161 * t93 + t165 * t94;
t194 = t167 * t143 + t163 * t168 + t237;
t193 = t109 * qJD(3) - t159 * t214;
t192 = t149 - t215;
t191 = -t20 * t282 - t81 * t288 + t1 - t238;
t189 = -qJD(6) * t46 - t12 + t124;
t188 = t13 * t160 + t21 * t83 + t36 + (g(1) * t272 + g(2) * t273) * t137;
t187 = t265 * t137;
t186 = -t150 - t192;
t183 = t224 * t317 + t287;
t182 = -t187 - t13;
t181 = -g(1) * (pkin(9) * t275 + t167 * t299) - g(2) * (pkin(9) * t276 + t163 * t299);
t180 = t166 * qJD(2) + t260 * t291;
t178 = t85 * t81 - t14 + t238;
t177 = -t126 * t38 - t243 * t317 + t297;
t176 = -qJD(6) * t197 - t3 * t160;
t175 = t14 * t86 + t15 * t87 - t47 * t58 + t48 * t59 - t265;
t174 = t176 + t1;
t173 = -t85 * t83 + t123 - t187 - t229;
t172 = -t158 * t169 + t150 - t214 + t251 + t306;
t171 = (-t246 * t261 - t248) * t165 - t190;
t135 = t155 * t170;
t134 = qJDD(4) * t166;
t127 = -pkin(4) * t165 - pkin(5);
t77 = t136 * t269 - t273;
t76 = -t136 * t272 - t270;
t75 = -t136 * t270 - t272;
t74 = t136 * t273 - t269;
t71 = t162 * qJD(2) - qJD(4) * t94;
t54 = pkin(4) * t261 + t56;
t52 = t165 * t73 - t285;
t45 = t301 - t302;
t34 = t223 + t171;
t33 = t179 + t311;
t32 = pkin(5) * t59 + pkin(9) * t58 + t111;
t29 = qJD(5) * t61 + t161 * t71 - t165 * t180;
t28 = -qJD(5) * t60 + t161 * t180 + t165 * t71;
t27 = t160 * t56 + t164 * t47;
t26 = -t160 * t47 + t164 * t56;
t23 = t160 * t54 + t164 * t52;
t22 = -t160 * t52 + t164 * t54;
t10 = t183 - t294;
t9 = t296 - t314;
t8 = t225 * t62 - t279;
t7 = t224 * t64 - t280;
t6 = -qJD(6) * t31 - t160 * t28 + t164 * t32;
t5 = qJD(6) * t30 + t160 * t32 + t164 * t28;
t4 = (-t24 - t310) * t164 + (-t25 - t228) * t160;
t11 = [0, 0, 0, 0, 0, qJDD(1), t306, t265, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t157 + t232, t219, t245 * pkin(1) - g(1) * (-pkin(1) * t163 + t140) - g(2) * t267 + (t141 + t152) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t219, 0.2e1 * t150 + t192 + t251, -g(1) * t198 - g(2) * t237 + qJ(2) * t231 + t133 * qJD(2) + t193, qJDD(1) * t155 - 0.2e1 * t216, -0.2e1 * t162 * t247 + 0.2e1 * (t154 - t155) * t253, -t162 * t169 + t134, qJDD(1) * t154 + 0.2e1 * t216, -t166 * t169 - t249, 0, t172 * t162 + t304 * t166, -t304 * t162 + t172 * t166, t265 + t264 * (-qJDD(1) * t158 - t151 - t88) -g(1) * (-t167 * pkin(7) + t198) - g(2) * (-pkin(7) * t163 + t237) + t158 * t227 + qJD(2) * t309 + t193, t201, -t179 * t86 - t42 * t87 + t58 * t81 - t59 * t83, t289, t200, t303, 0, t111 * t81 + t136 * t306 - t218 * t42 - t244 * t60 - t246 * t29 + t85 * t59 + t66 * t86, t111 * t83 - t179 * t218 - t244 * t61 - t246 * t28 - t85 * t58 + t66 * t87 + t308, t179 * t60 - t28 * t81 + t29 * t83 - t42 * t61 - t175, -g(2) * t194 + t85 * t111 + t14 * t61 - t15 * t60 + t48 * t28 - t47 * t29 - t236 + (-t147 - t66) * t218, t164 * t206 - t241 * t64 (t160 * t64 + t164 * t62) * t58 + (t280 - t279 + (-t283 + t286) * qJD(6)) * t87, t164 * t202 - t241 * t317 + t205, t160 * t204 + t240 * t62, -t160 * t202 - t240 * t317 + t203, t317 * t59 + t38 * t86, -g(1) * t75 - g(2) * t77 + t160 * t207 + t60 * t25 + t29 * t62 + t30 * t38 + t317 * t6 + t36 * t87 - t211, -g(1) * t74 - g(2) * t76 + t164 * t207 - t60 * t24 + t29 * t64 - t31 * t38 - t317 * t5 - t35 * t87 - t212, t30 * t24 - t31 * t25 - t5 * t62 - t6 * t64 + t197 * t58 - t308 + (qJD(6) * t196 - t2 * t160 - t3 * t164) * t87, t2 * t31 + t21 * t5 + t3 * t30 + t20 * t6 + t13 * t60 + t39 * t29 - t236 - g(2) * (-t167 * t230 + t194) - (t230 + t218) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t170, -qJ(2) * t170 + t215, 0, 0, 0, 0, 0, 0, 0, -t170, -qJDD(1), -qJD(1) * t133 + t186, 0, 0, 0, 0, 0, 0, -0.2e1 * t234 - t248, 0.2e1 * t235 - t247, t135 + t274, -qJD(1) * t309 + t186, 0, 0, 0, 0, 0, 0, -t223 + t171, -t179 + t311, t301 + t302, -t47 * t83 - t48 * t81 + t186 - t307, 0, 0, 0, 0, 0, 0, t296 + t314, t183 + t294, t322 * t164 + (t25 - t228) * t160, t39 * t83 + (-t3 - t320) * t164 + (-t2 + t321) * t160 - t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t170, t318 + t231, 0, 0, 0, 0, 0, 0, t162 * t263 + t134, t166 * t263 - t249, -t264 * qJDD(1), t227 + t318, 0, 0, 0, 0, 0, 0, -qJD(1) * t81 + t289, -qJD(1) * t83 + t303, -t200 - t201, t175 - t278, 0, 0, 0, 0, 0, 0, -t86 * t287 + (-t160 * t59 - t164 * t222) * t317 - t204, -t86 * t284 + (t160 * t222 - t164 * t59) * t317 - t206 (t222 * t64 + t203) * t164 + (t222 * t62 + t205) * t160 (-t20 * t222 + t212) * t164 + (-t21 * t222 + t211) * t160 - t207 - t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, t135 - t274, t247, -t239, -t248, qJDD(4), t166 * t318 + t298 + t80, g(3) * t166 + (-t318 - t88) * t162, 0, 0, t292, t45, t33, -t292, t34, t244, t51 * qJD(4) - t290 * qJD(5) + (t165 * t244 - t246 * t258 - t261 * t81) * pkin(4) + t173, t52 * t246 + (-t161 * t244 - t246 * t257 - t261 * t83) * pkin(4) + t178, t290 * t83 + (-t47 + t52) * t81 + (-t161 * t42 - t165 * t179 + (t161 * t83 - t165 * t81) * qJD(5)) * pkin(4), t47 * t51 - t48 * t52 + (t298 + t14 * t161 + t15 * t165 + (-t161 * t47 + t165 * t48) * qJD(5) + (-t265 - t278) * t166) * pkin(4), t7, t4, t10, t8, t9, -t293, t127 * t25 - t22 * t317 + t213 * t62 + t177 * t160 + (t182 - t242) * t164 + t217, -t127 * t24 + t23 * t317 + t213 * t64 + (t242 - t123) * t160 + t177 * t164 + t188, t22 * t64 + t23 * t62 + (-t62 * t243 - t126 * t25 + (t126 * t64 - t20) * qJD(6)) * t164 + (t64 * t243 - t126 * t24 - t3 + (t126 * t62 - t21) * qJD(6)) * t160 + t191, t13 * t127 - t21 * t23 - t20 * t22 - t39 * t51 - g(3) * (-t143 + t230) + (-t265 * t166 + (t161 * t39 - t165 * t196) * qJD(5)) * pkin(4) + t174 * t126 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, t45, t33, -t292, t34, t244, t48 * qJD(4) + t173, t246 * t47 + t178, 0, 0, t7, t4, t10, t8, t9, -t293, -pkin(5) * t25 - t26 * t317 - t48 * t62 + (-pkin(9) * t38 + t297) * t160 + (-pkin(9) * t277 + t182) * t164 + t217, -t160 * t123 + t39 * t282 + pkin(5) * t24 + t27 * t317 - t48 * t64 + (t256 * t317 - t284) * pkin(9) + t188, t26 * t64 + t27 * t62 + (-t280 - t279 + (t283 + t286) * qJD(6)) * pkin(9) + t176 + t191, -t13 * pkin(5) + pkin(9) * t174 - g(3) * t230 - t20 * t26 - t21 * t27 - t39 * t48 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t62 ^ 2 + t64 ^ 2, t322, -t295, -t226 + (-qJD(6) + t317) * t64, t38, -g(1) * t76 + g(2) * t74 + t160 * t189 - t255 * t40 - t39 * t64 + t16 + t320, g(1) * t77 - g(2) * t75 + t321 + t39 * t62 + (qJD(6) * t40 - t17) * t160 + t189 * t164, 0, 0;];
tau_reg  = t11;
