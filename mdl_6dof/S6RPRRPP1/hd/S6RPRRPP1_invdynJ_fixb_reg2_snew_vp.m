% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:17:27
% EndTime: 2019-05-05 21:17:47
% DurationCPUTime: 8.31s
% Computational Cost: add. (19388->414), mult. (39003->530), div. (0->0), fcn. (26375->10), ass. (0->271)
t241 = sin(pkin(9));
t243 = cos(pkin(9));
t247 = sin(qJ(3));
t250 = cos(qJ(3));
t246 = sin(qJ(4));
t249 = cos(qJ(4));
t303 = qJD(1) * t247;
t208 = -t249 * qJD(3) + t246 * t303;
t297 = qJD(1) * qJD(3);
t289 = t250 * t297;
t296 = t247 * qJDD(1);
t214 = t289 + t296;
t268 = -t246 * qJDD(3) - t249 * t214;
t181 = -qJD(4) * t208 - t268;
t240 = sin(pkin(10));
t242 = cos(pkin(10));
t209 = qJD(3) * t246 + t249 * t303;
t269 = t249 * qJDD(3) - t246 * t214;
t258 = -qJD(4) * t209 + t269;
t255 = t242 * t181 + t240 * t258;
t186 = t242 * t208 + t209 * t240;
t227 = qJD(1) * t250 - qJD(4);
t316 = t186 * t227;
t351 = t255 + t316;
t231 = t247 * t297;
t295 = t250 * qJDD(1);
t215 = -t231 + t295;
t207 = -qJDD(4) + t215;
t188 = -t208 * t240 + t209 * t242;
t315 = t188 * t186;
t260 = t207 - t315;
t324 = t260 * t240;
t185 = t188 ^ 2;
t340 = t227 ^ 2;
t353 = -t185 - t340;
t82 = t242 * t353 + t324;
t323 = t260 * t242;
t84 = -t240 * t353 + t323;
t47 = t246 * t82 - t249 * t84;
t37 = -t247 * t351 + t250 * t47;
t53 = t246 * t84 + t249 * t82;
t425 = pkin(1) * (t241 * t37 + t243 * t53) + pkin(2) * t53 + pkin(7) * t37;
t421 = pkin(3) * t53;
t420 = pkin(8) * t53;
t418 = pkin(3) * t351 + pkin(8) * t47;
t417 = t247 * t47 + t250 * t351;
t168 = t188 * t227;
t284 = t181 * t240 - t242 * t258;
t105 = -t284 - t168;
t342 = t186 ^ 2;
t160 = t342 - t340;
t89 = t160 * t240 - t323;
t93 = t160 * t242 + t324;
t416 = t250 * t105 + t247 * (t246 * t89 - t249 * t93);
t352 = t185 - t342;
t355 = -t168 + t284;
t64 = -t355 * t240 + t242 * t351;
t326 = t351 * t240;
t66 = t355 * t242 + t326;
t415 = t247 * (t246 * t64 + t249 * t66) + t250 * t352;
t350 = -t316 + t255;
t372 = t105 * t242 + t350 * t240;
t373 = t105 * t240 - t242 * t350;
t386 = t246 * t372 + t249 * t373;
t115 = -t342 - t185;
t387 = -t246 * t373 + t249 * t372;
t402 = t115 * t247 + t250 * t387;
t414 = pkin(1) * (t241 * t402 - t243 * t386) + pkin(7) * t402 - pkin(2) * t386;
t413 = pkin(4) * t82;
t411 = pkin(3) * t386;
t410 = pkin(8) * t386;
t409 = qJ(5) * t82;
t408 = qJ(5) * t84;
t406 = -pkin(3) * t115 + pkin(8) * t387;
t405 = t246 * t66 - t249 * t64;
t403 = -t115 * t250 + t247 * t387;
t401 = t246 * t93 + t249 * t89;
t128 = t207 + t315;
t322 = t128 * t240;
t349 = -t340 - t342;
t358 = t242 * t349 + t322;
t121 = t242 * t128;
t359 = t240 * t349 - t121;
t370 = t246 * t358 + t249 * t359;
t371 = -t246 * t359 + t249 * t358;
t388 = t247 * t355 + t250 * t371;
t399 = pkin(1) * (t241 * t388 - t243 * t370) + pkin(7) * t388 - pkin(2) * t370;
t397 = pkin(3) * t370;
t337 = pkin(4) * t373;
t396 = pkin(8) * t370;
t393 = qJ(5) * t373;
t391 = -pkin(3) * t355 + pkin(8) * t371;
t390 = -pkin(4) * t115 + qJ(5) * t372;
t389 = t247 * t371 - t250 * t355;
t161 = -t185 + t340;
t374 = t242 * t161 - t322;
t375 = -t161 * t240 - t121;
t385 = t246 * t375 + t249 * t374;
t384 = t247 * (-t246 * t374 + t249 * t375) - t250 * t350;
t336 = pkin(4) * t359;
t378 = qJ(6) * t351;
t306 = -g(3) + qJDD(2);
t230 = t250 * t306;
t252 = qJD(1) ^ 2;
t248 = sin(qJ(1));
t251 = cos(qJ(1));
t287 = t248 * g(1) - g(2) * t251;
t210 = qJDD(1) * pkin(1) + t287;
t276 = g(1) * t251 + g(2) * t248;
t211 = -pkin(1) * t252 - t276;
t304 = t241 * t210 + t243 * t211;
t175 = -pkin(2) * t252 + qJDD(1) * pkin(7) + t304;
t277 = -pkin(3) * t250 - pkin(8) * t247;
t280 = t252 * t277 + t175;
t339 = qJD(3) ^ 2;
t141 = -qJDD(3) * pkin(3) - t339 * pkin(8) + t280 * t247 - t230;
t193 = -pkin(4) * t227 - qJ(5) * t209;
t341 = t208 ^ 2;
t77 = -t258 * pkin(4) - t341 * qJ(5) + t209 * t193 + qJDD(5) + t141;
t381 = pkin(5) * t284 - t378 + t77;
t380 = qJ(5) * t358;
t379 = qJ(5) * t359;
t314 = t209 * t208;
t259 = -t207 - t314;
t357 = t246 * t259;
t356 = t249 * t259;
t302 = qJD(5) * t186;
t178 = -0.2e1 * t302;
t300 = qJD(6) * t227;
t354 = t178 - 0.2e1 * t300;
t196 = t208 * t227;
t151 = t181 - t196;
t143 = pkin(5) * t186 - qJ(6) * t188;
t283 = t243 * t210 - t241 * t211;
t174 = -qJDD(1) * pkin(2) - t252 * pkin(7) - t283;
t272 = -t215 + t231;
t273 = t214 + t289;
t134 = pkin(3) * t272 - pkin(8) * t273 + t174;
t285 = t247 * t306;
t142 = -t339 * pkin(3) + qJDD(3) * pkin(8) + t280 * t250 + t285;
t79 = -t249 * t134 + t246 * t142;
t61 = t259 * pkin(4) - t151 * qJ(5) - t79;
t80 = t246 * t134 + t249 * t142;
t71 = -t341 * pkin(4) + qJ(5) * t258 + t227 * t193 + t80;
t334 = t240 * t61 + t242 * t71;
t281 = -t207 * qJ(6) - t186 * t143 + t334;
t347 = -t413 - pkin(5) * (t353 + t340) - qJ(6) * t260 + t281;
t147 = (qJD(4) + t227) * t209 - t269;
t264 = (t186 * t240 + t188 * t242) * t227;
t313 = t227 * t240;
t159 = t188 * t313;
t312 = t227 * t242;
t294 = t186 * t312;
t274 = -t159 + t294;
t346 = t246 * t274 + t249 * t264;
t266 = t240 * t284 - t294;
t275 = -t186 * t313 - t242 * t284;
t345 = t246 * t266 + t249 * t275;
t197 = t250 * t207;
t344 = t247 * (-t246 * t264 + t249 * t274) + t197;
t293 = t250 * t315;
t343 = t247 * (-t246 * t275 + t249 * t266) + t293;
t206 = t209 ^ 2;
t286 = t240 * t71 - t242 * t61;
t301 = qJD(5) * t188;
t31 = t286 + 0.2e1 * t301;
t32 = t178 + t334;
t16 = t240 * t32 - t242 * t31;
t338 = pkin(4) * t16;
t335 = pkin(5) * t242;
t333 = t16 * t246;
t332 = t16 * t249;
t331 = t240 * t77;
t330 = t242 * t77;
t329 = qJ(6) * t242;
t321 = t141 * t246;
t320 = t141 * t249;
t170 = t207 - t314;
t318 = t170 * t246;
t317 = t170 * t249;
t311 = t227 * t246;
t310 = t227 * t249;
t226 = t250 * t252 * t247;
t220 = qJDD(3) + t226;
t309 = t247 * t220;
t219 = -t226 + qJDD(3);
t307 = t250 * t219;
t299 = qJD(4) - t227;
t292 = t250 * t314;
t290 = pkin(1) * t241 + pkin(7);
t288 = -qJ(6) * t240 - pkin(4);
t17 = t240 * t31 + t242 * t32;
t51 = t246 * t79 + t249 * t80;
t156 = t247 * t175 - t230;
t157 = t250 * t175 + t285;
t114 = t247 * t156 + t250 * t157;
t282 = (0.2e1 * qJD(5) + t143) * t188;
t279 = -t334 + t413;
t100 = t242 * t255 + t159;
t99 = -t188 * t312 + t240 * t255;
t278 = t247 * (t100 * t249 - t246 * t99) - t293;
t271 = t246 * t80 - t249 * t79;
t267 = t281 + t354;
t265 = t207 * pkin(5) - qJ(6) * t340 + qJDD(6) + t286;
t24 = -pkin(5) * t340 + t267;
t25 = t282 + t265;
t12 = t24 * t240 - t242 * t25;
t263 = pkin(4) * t12 - pkin(5) * t25 + qJ(6) * t24;
t262 = -pkin(5) * t350 + qJ(6) * t105 + t337;
t261 = -pkin(1) * t243 - pkin(2) + t277;
t254 = pkin(5) * t128 - qJ(6) * t349 + t265 - t336;
t253 = 0.2e1 * qJD(6) * t188 - t381;
t237 = t250 ^ 2;
t236 = t247 ^ 2;
t234 = t237 * t252;
t232 = t236 * t252;
t224 = -t234 - t339;
t223 = -t232 - t339;
t218 = t232 + t234;
t217 = (t236 + t237) * qJDD(1);
t216 = -0.2e1 * t231 + t295;
t213 = 0.2e1 * t289 + t296;
t195 = -t206 + t340;
t194 = -t340 + t341;
t192 = -t223 * t247 - t307;
t191 = t224 * t250 - t309;
t190 = t206 - t341;
t189 = -t206 - t340;
t182 = -t340 - t341;
t180 = -0.2e1 * t301;
t179 = 0.2e1 * t302;
t169 = t206 + t341;
t152 = t299 * t208 + t268;
t150 = t181 + t196;
t148 = -t299 * t209 + t269;
t139 = -t189 * t246 + t317;
t138 = t189 * t249 + t318;
t126 = t182 * t249 - t357;
t125 = t182 * t246 + t356;
t112 = -t147 * t249 + t151 * t246;
t86 = t139 * t250 - t152 * t247;
t81 = t126 * t250 - t148 * t247;
t60 = t100 * t246 + t249 * t99;
t52 = t330 - t409;
t49 = t331 - t379;
t41 = -pkin(4) * t351 + t331 + t408;
t40 = (-pkin(5) * t227 - 0.2e1 * qJD(6)) * t188 + t381;
t39 = -pkin(4) * t355 - t330 + t380;
t29 = (-t355 + t168) * pkin(5) + t253;
t28 = pkin(5) * t168 + t253 + t378;
t23 = -qJ(6) * t115 + t25;
t22 = (-t115 - t340) * pkin(5) + t267;
t21 = -t240 * t29 - t329 * t355 - t379;
t20 = -pkin(5) * t326 + t242 * t28 + t409;
t19 = t242 * t29 + t288 * t355 + t380;
t18 = -t408 + t240 * t28 + (pkin(4) + t335) * t351;
t15 = -pkin(4) * t77 + qJ(5) * t17;
t14 = -t16 - t393;
t13 = t24 * t242 + t240 * t25;
t11 = t17 + t390;
t10 = -t22 * t240 + t23 * t242 - t393;
t9 = t22 * t242 + t23 * t240 + t390;
t8 = t17 * t249 - t333;
t7 = t17 * t246 + t332;
t6 = -qJ(5) * t12 + (pkin(5) * t240 - t329) * t40;
t5 = t247 * t77 + t250 * t8;
t4 = qJ(5) * t13 + (t288 - t335) * t40;
t3 = -t12 * t246 + t13 * t249;
t2 = t12 * t249 + t13 * t246;
t1 = t247 * t40 + t250 * t3;
t26 = [0, 0, 0, 0, 0, qJDD(1), t287, t276, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t243 - t241 * t252) + t283, pkin(1) * (-qJDD(1) * t241 - t243 * t252) - t304, 0, pkin(1) * (t241 * t304 + t243 * t283), t273 * t247, t213 * t250 + t216 * t247, t309 + t250 * (-t232 + t339), -t272 * t250, t247 * (t234 - t339) + t307, 0, -t250 * t174 + pkin(2) * t216 + pkin(7) * t191 + pkin(1) * (t191 * t241 + t216 * t243), t247 * t174 - pkin(2) * t213 + pkin(7) * t192 + pkin(1) * (t192 * t241 - t213 * t243), pkin(2) * t218 + pkin(7) * t217 + pkin(1) * (t217 * t241 + t218 * t243) + t114, -pkin(2) * t174 + pkin(7) * t114 + pkin(1) * (t114 * t241 - t174 * t243), t247 * (t181 * t249 + t209 * t311) - t292, t247 * (t148 * t249 - t150 * t246) - t250 * t190, t247 * (-t195 * t246 + t356) - t250 * t151, t247 * (-t208 * t310 - t246 * t258) + t292, t247 * (t194 * t249 + t318) + t250 * t147, t197 + t247 * (t208 * t249 - t209 * t246) * t227, t247 * (-pkin(8) * t125 + t321) + t250 * (-pkin(3) * t125 + t79) - pkin(2) * t125 + pkin(7) * t81 + pkin(1) * (-t125 * t243 + t241 * t81), t247 * (-pkin(8) * t138 + t320) + t250 * (-pkin(3) * t138 + t80) - pkin(2) * t138 + pkin(7) * t86 + pkin(1) * (-t138 * t243 + t241 * t86), -t247 * t271 + t290 * (t112 * t250 - t169 * t247) + t261 * (-t147 * t246 - t151 * t249), t290 * (t141 * t247 + t250 * t51) + t261 * t271, t278, -t415, t384, t343, -t416, t344, t247 * (-t246 * t39 + t249 * t49 - t396) + t250 * (t31 - t336 - t397) + t399, t247 * (-t246 * t41 + t249 * t52 - t420) + t250 * (t178 - t279 - t421) - t425, t247 * (-t11 * t246 + t14 * t249 - t410) + t250 * (-t337 - t411) + t414, t247 * (-pkin(8) * t7 - qJ(5) * t332 - t15 * t246) + t250 * (-pkin(3) * t7 - t338) - pkin(2) * t7 + pkin(7) * t5 + pkin(1) * (t241 * t5 - t243 * t7), t278, t384, t415, t344, t416, t343, t247 * (-t19 * t246 + t21 * t249 - t396) + t250 * (t254 + t282 - t397) + t399, t247 * (t10 * t249 - t246 * t9 - t410) + t250 * (-t262 - t411) + t414, t247 * (-t18 * t246 + t20 * t249 + t420) + t250 * (t179 + 0.2e1 * t300 - t347 + t421) + t425, t247 * (-pkin(8) * t2 - t246 * t4 + t249 * t6) + t250 * (-pkin(3) * t2 - t263) - pkin(2) * t2 + pkin(7) * t1 + pkin(1) * (t1 * t241 - t2 * t243); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, 0, 0, 0, 0, 0, 0, t220 * t250 + t224 * t247, -t219 * t247 + t223 * t250, 0, -t156 * t250 + t157 * t247, 0, 0, 0, 0, 0, 0, t126 * t247 + t148 * t250, t139 * t247 + t152 * t250, t112 * t247 + t169 * t250, -t141 * t250 + t247 * t51, 0, 0, 0, 0, 0, 0, t389, -t417, t403, t247 * t8 - t250 * t77, 0, 0, 0, 0, 0, 0, t389, t403, t417, t247 * t3 - t250 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, t232 - t234, t296, t226, t295, qJDD(3), -t156, -t157, 0, 0, t181 * t246 - t209 * t310, t148 * t246 + t150 * t249, t195 * t249 + t357, -t208 * t311 + t249 * t258, t194 * t246 - t317, (t208 * t246 + t209 * t249) * t227, pkin(3) * t148 + pkin(8) * t126 - t320, pkin(3) * t152 + pkin(8) * t139 + t321, pkin(3) * t169 + pkin(8) * t112 + t51, -pkin(3) * t141 + pkin(8) * t51, t60, -t405, t385, t345, t401, t346, t246 * t49 + t249 * t39 + t391, t246 * t52 + t249 * t41 - t418, t11 * t249 + t14 * t246 + t406, -pkin(3) * t77 + pkin(8) * t8 - qJ(5) * t333 + t15 * t249, t60, t385, t405, t346, -t401, t345, t19 * t249 + t21 * t246 + t391, t10 * t246 + t249 * t9 + t406, t18 * t249 + t20 * t246 + t418, -pkin(3) * t40 + pkin(8) * t3 + t246 * t6 + t249 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t190, t151, -t314, -t147, -t207, -t79, -t80, 0, 0, t315, t352, t350, -t315, t105, -t207, t180 - t286 + t336, t179 + t279, t337, t338, t315, t350, -t352, -t207, -t105, -t315, -t143 * t188 + t180 - t254, t262, t347 + t354, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, t351, t115, t77, 0, 0, 0, 0, 0, 0, t355, t115, -t351, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t350, t353, t25;];
tauJ_reg  = t26;
