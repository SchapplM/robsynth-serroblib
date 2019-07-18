% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:28:57
% EndTime: 2019-05-06 04:29:14
% DurationCPUTime: 6.77s
% Computational Cost: add. (47284->461), mult. (94415->642), div. (0->0), fcn. (65679->10), ass. (0->310)
t264 = sin(qJ(6));
t266 = sin(qJ(4));
t271 = cos(qJ(4));
t272 = cos(qJ(3));
t313 = qJD(1) * t272;
t238 = -qJD(3) * t271 + t266 * t313;
t240 = qJD(3) * t266 + t271 * t313;
t265 = sin(qJ(5));
t270 = cos(qJ(5));
t214 = t238 * t270 + t240 * t265;
t216 = -t238 * t265 + t240 * t270;
t269 = cos(qJ(6));
t179 = t214 * t269 + t216 * t264;
t181 = -t214 * t264 + t216 * t269;
t142 = t181 * t179;
t308 = qJD(1) * qJD(3);
t254 = t272 * t308;
t267 = sin(qJ(3));
t255 = t267 * qJDD(1);
t242 = -t255 - t254;
t237 = qJDD(4) - t242;
t234 = qJDD(5) + t237;
t232 = qJDD(6) + t234;
t349 = -t142 + t232;
t357 = t264 * t349;
t183 = t216 * t214;
t347 = -t183 + t234;
t356 = t265 * t347;
t221 = t240 * t238;
t346 = -t221 + t237;
t355 = t266 * t346;
t354 = t269 * t349;
t353 = t270 * t347;
t352 = t271 * t346;
t275 = qJD(1) ^ 2;
t345 = pkin(7) + pkin(1);
t351 = t345 * t275;
t257 = t272 * qJDD(1);
t300 = t267 * t308;
t243 = t257 - t300;
t296 = -qJDD(3) * t271 + t266 * t243;
t209 = -qJD(4) * t240 - t296;
t285 = -t266 * qJDD(3) - t271 * t243;
t210 = -qJD(4) * t238 - t285;
t297 = -t209 * t270 + t265 * t210;
t158 = -qJD(5) * t216 - t297;
t159 = -qJD(5) * t214 + t209 * t265 + t210 * t270;
t102 = -qJD(6) * t179 + t158 * t264 + t159 * t269;
t253 = qJD(1) * t267 + qJD(4);
t251 = qJD(5) + t253;
t247 = qJD(6) + t251;
t166 = t247 * t179;
t350 = t102 - t166;
t200 = t251 * t214;
t348 = t159 - t200;
t139 = t159 + t200;
t229 = t253 * t238;
t192 = t210 + t229;
t298 = -t158 * t269 + t264 * t159;
t83 = (qJD(6) - t247) * t181 + t298;
t136 = (qJD(5) - t251) * t216 + t297;
t188 = (qJD(4) - t253) * t240 + t296;
t177 = t179 ^ 2;
t178 = t181 ^ 2;
t212 = t214 ^ 2;
t213 = t216 ^ 2;
t235 = t238 ^ 2;
t236 = t240 ^ 2;
t246 = t247 ^ 2;
t250 = t251 ^ 2;
t252 = t253 ^ 2;
t307 = qJD(2) * qJD(1);
t259 = 0.2e1 * t307;
t261 = qJDD(1) * qJ(2);
t268 = sin(qJ(1));
t273 = cos(qJ(1));
t290 = t273 * g(1) + t268 * g(2);
t283 = -t261 + t290;
t278 = t259 - t283;
t287 = -t243 + t300;
t288 = -t242 + t254;
t187 = pkin(3) * t288 + pkin(8) * t287 + t278 - t351;
t299 = t268 * g(1) - g(2) * t273;
t289 = qJDD(2) - t299;
t315 = t275 * qJ(2);
t277 = t289 - t315;
t228 = -qJDD(1) * t345 + t277;
t218 = t272 * g(3) - t228 * t267;
t274 = qJD(3) ^ 2;
t293 = pkin(3) * t267 - pkin(8) * t272;
t280 = t275 * t293;
t197 = -t274 * pkin(3) + qJDD(3) * pkin(8) - t267 * t280 - t218;
t153 = -t187 * t271 + t266 * t197;
t118 = pkin(4) * t346 - pkin(9) * t192 - t153;
t154 = t187 * t266 + t197 * t271;
t292 = pkin(4) * t253 - pkin(9) * t240;
t121 = -t235 * pkin(4) + t209 * pkin(9) - t253 * t292 + t154;
t76 = -t118 * t270 + t265 * t121;
t77 = t118 * t265 + t121 * t270;
t44 = t265 * t77 - t270 * t76;
t344 = pkin(4) * t44;
t97 = -t136 * t265 - t139 * t270;
t343 = pkin(4) * t97;
t217 = g(3) * t267 + t228 * t272;
t196 = qJDD(3) * pkin(3) + pkin(8) * t274 - t272 * t280 + t217;
t146 = pkin(4) * t209 + pkin(9) * t235 - t240 * t292 + t196;
t291 = pkin(5) * t251 - pkin(10) * t216;
t88 = pkin(5) * t158 + pkin(10) * t212 - t216 * t291 + t146;
t342 = t264 * t88;
t58 = pkin(5) * t347 - pkin(10) * t139 - t76;
t64 = -t212 * pkin(5) + t158 * pkin(10) - t251 * t291 + t77;
t35 = t264 * t64 - t269 * t58;
t36 = t264 * t58 + t269 * t64;
t18 = t264 * t36 - t269 * t35;
t341 = t265 * t18;
t340 = t266 * t44;
t339 = t269 * t88;
t338 = t270 * t18;
t337 = t271 * t44;
t336 = qJDD(1) * pkin(1);
t335 = t247 * t264;
t334 = t247 * t269;
t333 = t251 * t265;
t332 = t251 * t270;
t331 = t253 * t266;
t330 = t253 * t271;
t262 = t267 ^ 2;
t329 = t262 * t275;
t263 = t272 ^ 2;
t328 = t263 * t275;
t128 = t142 + t232;
t327 = t264 * t128;
t326 = t265 * t146;
t169 = t183 + t234;
t325 = t265 * t169;
t324 = t266 * t196;
t204 = t221 + t237;
t323 = t266 * t204;
t302 = t267 * t275 * t272;
t322 = t267 * (qJDD(3) + t302);
t321 = t269 * t128;
t320 = t270 * t146;
t319 = t270 * t169;
t318 = t271 * t196;
t317 = t271 * t204;
t316 = t272 * (qJDD(3) - t302);
t314 = t262 + t263;
t311 = qJD(4) + t253;
t17 = pkin(5) * t18;
t19 = t264 * t35 + t269 * t36;
t8 = t19 * t265 + t338;
t306 = pkin(4) * t8 + t17;
t305 = t267 * t142;
t304 = t267 * t183;
t303 = t267 * t221;
t86 = t102 + t166;
t52 = -t264 * t83 - t269 * t86;
t54 = t264 * t86 - t269 * t83;
t28 = t265 * t54 + t270 * t52;
t50 = pkin(5) * t52;
t301 = pkin(4) * t28 + t50;
t45 = t265 * t76 + t270 * t77;
t115 = t153 * t266 + t154 * t271;
t134 = -t246 - t177;
t94 = t134 * t264 + t354;
t295 = pkin(5) * t94 - t35;
t195 = -t213 - t250;
t144 = t195 * t270 - t325;
t294 = pkin(4) * t144 - t77;
t286 = t153 * t271 - t154 * t266;
t186 = t272 * t217 - t267 * t218;
t284 = qJ(2) + t293;
t174 = -t250 - t212;
t125 = t174 * t265 + t353;
t282 = pkin(4) * t125 - t76;
t161 = -t178 - t246;
t108 = t161 * t269 - t327;
t281 = pkin(5) * t108 - t36;
t95 = t134 * t269 - t357;
t59 = t265 * t95 + t270 * t94;
t279 = pkin(4) * t59 + t295;
t109 = -t161 * t264 - t321;
t65 = t108 * t270 + t109 * t265;
t276 = pkin(4) * t65 + t281;
t245 = t314 * qJDD(1);
t244 = t257 - 0.2e1 * t300;
t241 = t255 + 0.2e1 * t254;
t230 = -t277 + t336;
t227 = -t236 + t252;
t226 = t235 - t252;
t225 = t283 - 0.2e1 * t307 + t351;
t223 = -t322 + t272 * (-t274 - t328);
t222 = t267 * (-t274 - t329) + t316;
t220 = t236 - t235;
t219 = -t236 - t252;
t211 = -t252 - t235;
t202 = t235 + t236;
t199 = -t213 + t250;
t198 = t212 - t250;
t193 = t238 * t311 + t285;
t191 = t210 - t229;
t189 = -t240 * t311 - t296;
t182 = t213 - t212;
t176 = -t219 * t266 - t317;
t175 = t219 * t271 - t323;
t172 = t211 * t271 - t355;
t171 = t211 * t266 + t352;
t165 = -t178 + t246;
t164 = t177 - t246;
t163 = (-t214 * t270 + t216 * t265) * t251;
t162 = (-t214 * t265 - t216 * t270) * t251;
t160 = -t212 - t213;
t156 = -t188 * t271 + t192 * t266;
t151 = t198 * t270 - t325;
t150 = -t199 * t265 + t353;
t149 = t198 * t265 + t319;
t148 = t199 * t270 + t356;
t147 = t176 * t267 + t193 * t272;
t145 = -t195 * t265 - t319;
t143 = t172 * t267 + t189 * t272;
t141 = t178 - t177;
t135 = (qJD(5) + t251) * t216 + t297;
t133 = t159 * t270 - t216 * t333;
t132 = t159 * t265 + t216 * t332;
t131 = -t158 * t265 + t214 * t332;
t130 = t158 * t270 + t214 * t333;
t126 = t174 * t270 - t356;
t124 = t156 * t267 + t202 * t272;
t123 = (-t179 * t269 + t181 * t264) * t247;
t122 = (-t179 * t264 - t181 * t269) * t247;
t119 = -t177 - t178;
t113 = t164 * t269 - t327;
t112 = -t165 * t264 + t354;
t111 = t164 * t264 + t321;
t110 = t165 * t269 + t357;
t106 = -pkin(9) * t144 - t320;
t105 = -t144 * t266 + t145 * t271;
t104 = t144 * t271 + t145 * t266;
t103 = t115 * t267 + t196 * t272;
t101 = -qJD(6) * t181 - t298;
t100 = -pkin(9) * t125 - t326;
t99 = -t136 * t270 + t139 * t265;
t98 = -t135 * t270 - t265 * t348;
t96 = -t135 * t265 + t270 * t348;
t92 = -t125 * t266 + t126 * t271;
t91 = t125 * t271 + t126 * t266;
t90 = -t122 * t265 + t123 * t270;
t89 = t122 * t270 + t123 * t265;
t82 = (qJD(6) + t247) * t181 + t298;
t81 = t102 * t269 - t181 * t335;
t80 = t102 * t264 + t181 * t334;
t79 = -t101 * t264 + t179 * t334;
t78 = t101 * t269 + t179 * t335;
t74 = -pkin(4) * t348 + pkin(9) * t145 - t326;
t73 = t105 * t267 - t272 * t348;
t72 = -pkin(4) * t135 + pkin(9) * t126 + t320;
t71 = -t135 * t272 + t267 * t92;
t70 = -t111 * t265 + t113 * t270;
t69 = -t110 * t265 + t112 * t270;
t68 = t111 * t270 + t113 * t265;
t67 = t110 * t270 + t112 * t265;
t66 = -t108 * t265 + t109 * t270;
t63 = -pkin(10) * t108 - t339;
t62 = -t266 * t97 + t271 * t99;
t61 = t266 * t99 + t271 * t97;
t60 = -t265 * t94 + t270 * t95;
t56 = -pkin(10) * t94 - t342;
t55 = -t160 * t272 + t267 * t62;
t53 = -t264 * t350 - t269 * t82;
t51 = -t264 * t82 + t269 * t350;
t49 = -t265 * t80 + t270 * t81;
t48 = -t265 * t78 + t270 * t79;
t47 = t265 * t81 + t270 * t80;
t46 = t265 * t79 + t270 * t78;
t43 = pkin(4) * t146 + pkin(9) * t45;
t42 = -pkin(5) * t350 + pkin(10) * t109 - t342;
t41 = -t266 * t65 + t271 * t66;
t40 = t266 * t66 + t271 * t65;
t39 = -pkin(5) * t82 + pkin(10) * t95 + t339;
t38 = -pkin(9) * t97 - t44;
t37 = -pkin(4) * t160 + pkin(9) * t99 + t45;
t33 = -t266 * t59 + t271 * t60;
t32 = t266 * t60 + t271 * t59;
t31 = t267 * t41 - t272 * t350;
t30 = -t265 * t52 + t270 * t54;
t29 = -t265 * t51 + t270 * t53;
t27 = t265 * t53 + t270 * t51;
t26 = t267 * t33 - t272 * t82;
t25 = t271 * t45 - t340;
t24 = t266 * t45 + t337;
t23 = t146 * t272 + t25 * t267;
t22 = -pkin(9) * t65 - t265 * t42 + t270 * t63;
t21 = -pkin(9) * t59 - t265 * t39 + t270 * t56;
t20 = -pkin(4) * t350 + pkin(9) * t66 + t265 * t63 + t270 * t42;
t16 = -pkin(4) * t82 + pkin(9) * t60 + t265 * t56 + t270 * t39;
t15 = pkin(5) * t88 + pkin(10) * t19;
t14 = -t266 * t28 + t271 * t30;
t13 = t266 * t30 + t271 * t28;
t12 = -t119 * t272 + t14 * t267;
t11 = -pkin(10) * t52 - t18;
t10 = -pkin(5) * t119 + pkin(10) * t54 + t19;
t9 = t19 * t270 - t341;
t7 = -pkin(9) * t28 - t10 * t265 + t11 * t270;
t6 = -pkin(4) * t119 + pkin(9) * t30 + t10 * t270 + t11 * t265;
t5 = -t266 * t8 + t271 * t9;
t4 = t266 * t9 + t271 * t8;
t3 = -pkin(9) * t8 - pkin(10) * t338 - t15 * t265;
t2 = t267 * t5 + t272 * t88;
t1 = pkin(4) * t88 + pkin(9) * t9 - pkin(10) * t341 + t15 * t270;
t34 = [0, 0, 0, 0, 0, qJDD(1), t299, t290, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t289 - 0.2e1 * t336, t259 + 0.2e1 * t261 - t290, pkin(1) * t230 + qJ(2) * (-t275 * pkin(1) + t278), -t287 * t272, -t241 * t272 - t244 * t267, t316 - t267 * (t274 - t328), t288 * t267, t272 * (-t274 + t329) - t322, 0, qJ(2) * t241 - t222 * t345 - t267 * t225, qJ(2) * t244 - t223 * t345 - t272 * t225, t245 * t345 - t314 * t315 - t186, -qJ(2) * t225 - t186 * t345, t272 * (t210 * t271 - t240 * t331) + t303, t272 * (t189 * t271 - t191 * t266) + t267 * t220, t272 * (-t227 * t266 + t352) + t267 * t192, t272 * (-t209 * t266 + t238 * t330) - t303, t272 * (t226 * t271 - t323) - t267 * t188, t267 * t237 + t272 * (-t238 * t271 + t240 * t266) * t253, t272 * (-pkin(8) * t171 - t324) - t267 * (-pkin(3) * t171 + t153) + qJ(2) * t171 - t345 * t143, t272 * (-pkin(8) * t175 - t318) - t267 * (-pkin(3) * t175 + t154) + qJ(2) * t175 - t345 * t147, t272 * t286 + t284 * (-t188 * t266 - t192 * t271) - t345 * t124, -t103 * t345 - t284 * t286, t272 * (-t132 * t266 + t133 * t271) + t304, t272 * (-t266 * t96 + t271 * t98) + t267 * t182, t272 * (-t148 * t266 + t150 * t271) + t267 * t139, t272 * (-t130 * t266 + t131 * t271) - t304, t272 * (-t149 * t266 + t151 * t271) - t267 * t136, t272 * (-t162 * t266 + t163 * t271) + t267 * t234, t272 * (-pkin(8) * t91 + t100 * t271 - t266 * t72) - t267 * (-pkin(3) * t91 - t282) + qJ(2) * t91 - t345 * t71, t272 * (-pkin(8) * t104 + t106 * t271 - t266 * t74) - t267 * (-pkin(3) * t104 - t294) + qJ(2) * t104 - t345 * t73, t272 * (-pkin(8) * t61 - t266 * t37 + t271 * t38) - t267 * (-pkin(3) * t61 - t343) + qJ(2) * t61 - t345 * t55, t272 * (-pkin(8) * t24 - pkin(9) * t337 - t266 * t43) - t267 * (-pkin(3) * t24 - t344) + qJ(2) * t24 - t345 * t23, t272 * (-t266 * t47 + t271 * t49) + t305, t272 * (-t266 * t27 + t271 * t29) + t267 * t141, t272 * (-t266 * t67 + t271 * t69) + t267 * t86, t272 * (-t266 * t46 + t271 * t48) - t305, t272 * (-t266 * t68 + t271 * t70) - t267 * t83, t272 * (-t266 * t89 + t271 * t90) + t267 * t232, t272 * (-pkin(8) * t32 - t16 * t266 + t21 * t271) - t267 * (-pkin(3) * t32 - t279) + qJ(2) * t32 - t345 * t26, t272 * (-pkin(8) * t40 - t20 * t266 + t22 * t271) - t267 * (-pkin(3) * t40 - t276) + qJ(2) * t40 - t345 * t31, t272 * (-pkin(8) * t13 - t266 * t6 + t271 * t7) - t267 * (-pkin(3) * t13 - t301) + qJ(2) * t13 - t345 * t12, t272 * (-pkin(8) * t4 - t1 * t266 + t271 * t3) - t267 * (-pkin(3) * t4 - t306) + qJ(2) * t4 - t345 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t275, -t230, 0, 0, 0, 0, 0, 0, t222, t223, -t245, t186, 0, 0, 0, 0, 0, 0, t143, t147, t124, t103, 0, 0, 0, 0, 0, 0, t71, t73, t55, t23, 0, 0, 0, 0, 0, 0, t26, t31, t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, (-t262 + t263) * t275, t257, -t302, -t255, qJDD(3), t217, t218, 0, 0, t210 * t266 + t240 * t330, t189 * t266 + t191 * t271, t227 * t271 + t355, t209 * t271 + t238 * t331, t226 * t266 + t317, (-t238 * t266 - t240 * t271) * t253, pkin(3) * t189 + pkin(8) * t172 + t318, pkin(3) * t193 + pkin(8) * t176 - t324, pkin(3) * t202 + pkin(8) * t156 + t115, pkin(3) * t196 + pkin(8) * t115, t132 * t271 + t133 * t266, t266 * t98 + t271 * t96, t148 * t271 + t150 * t266, t130 * t271 + t131 * t266, t149 * t271 + t151 * t266, t162 * t271 + t163 * t266, -pkin(3) * t135 + pkin(8) * t92 + t100 * t266 + t271 * t72, -pkin(3) * t348 + pkin(8) * t105 + t106 * t266 + t271 * t74, -pkin(3) * t160 + pkin(8) * t62 + t266 * t38 + t271 * t37, pkin(3) * t146 + pkin(8) * t25 - pkin(9) * t340 + t271 * t43, t266 * t49 + t271 * t47, t266 * t29 + t27 * t271, t266 * t69 + t271 * t67, t266 * t48 + t271 * t46, t266 * t70 + t271 * t68, t266 * t90 + t271 * t89, -pkin(3) * t82 + pkin(8) * t33 + t16 * t271 + t21 * t266, -pkin(3) * t350 + pkin(8) * t41 + t20 * t271 + t22 * t266, -pkin(3) * t119 + pkin(8) * t14 + t266 * t7 + t271 * t6, pkin(3) * t88 + pkin(8) * t5 + t1 * t271 + t266 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t220, t192, -t221, -t188, t237, -t153, -t154, 0, 0, t183, t182, t139, -t183, -t136, t234, t282, t294, t343, t344, t142, t141, t86, -t142, -t83, t232, t279, t276, t301, t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t182, t139, -t183, -t136, t234, -t76, -t77, 0, 0, t142, t141, t86, -t142, -t83, t232, t295, t281, t50, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t141, t86, -t142, -t83, t232, -t35, -t36, 0, 0;];
tauJ_reg  = t34;