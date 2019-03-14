% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:09
% EndTime: 2019-03-10 01:23:27
% DurationCPUTime: 7.92s
% Computational Cost: add. (10914->487), mult. (26561->672), div. (0->0), fcn. (19277->8), ass. (0->254)
t245 = cos(qJ(3));
t243 = sin(qJ(2));
t246 = cos(qJ(2));
t344 = t245 * t246;
t277 = pkin(3) * t243 - pkin(9) * t344;
t374 = pkin(8) + pkin(9);
t311 = qJD(3) * t374;
t283 = pkin(2) * t243 - pkin(8) * t246;
t203 = t283 * qJD(1);
t242 = sin(qJ(3));
t328 = qJD(1) * t243;
t309 = t242 * t328;
t331 = pkin(7) * t309 + t245 * t203;
t412 = qJD(1) * t277 + t245 * t311 + t331;
t183 = t242 * t203;
t345 = t243 * t245;
t346 = t242 * t246;
t411 = -t183 - (-pkin(7) * t345 - pkin(9) * t346) * qJD(1) - t242 * t311;
t241 = sin(qJ(4));
t244 = cos(qJ(4));
t279 = t241 * t242 - t244 * t245;
t327 = qJD(1) * t246;
t384 = qJD(4) + qJD(3);
t399 = t384 * t279;
t336 = t279 * t327 - t399;
t410 = t411 * t241 + t244 * t412;
t240 = sin(qJ(5));
t316 = t245 * qJD(2);
t196 = t309 - t316;
t325 = qJD(2) * t242;
t198 = t245 * t328 + t325;
t281 = -t196 * t241 + t244 * t198;
t282 = -t196 * t244 - t241 * t198;
t371 = cos(qJ(5));
t376 = t240 * t282 + t281 * t371;
t375 = t376 ^ 2;
t86 = t240 * t281 - t282 * t371;
t83 = t86 ^ 2;
t405 = -t83 + t375;
t199 = t241 * t245 + t242 * t244;
t165 = t199 * t327;
t268 = t199 * qJD(3);
t253 = -qJD(4) * t199 - t268;
t404 = t165 + t253;
t213 = t374 * t242;
t318 = qJD(4) * t244;
t408 = -t213 * t318 - t241 * t412 + t411 * t244;
t214 = t374 * t245;
t332 = -t241 * t213 + t244 * t214;
t407 = pkin(4) * t328 + t336 * pkin(10) + qJD(4) * t332 + t410;
t349 = t241 * t214;
t368 = pkin(10) * t199;
t406 = (-t349 - t368) * qJD(4) + t408 + (t165 - t268) * pkin(10);
t401 = t376 * t86;
t400 = t86 * qJ(6);
t224 = -qJD(3) + t327;
t217 = -qJD(4) + t224;
t395 = pkin(10) * t281;
t208 = -pkin(2) * t246 - pkin(8) * t243 - pkin(1);
t189 = t208 * qJD(1);
t233 = pkin(7) * t327;
t212 = qJD(2) * pkin(8) + t233;
t149 = t245 * t189 - t212 * t242;
t112 = -pkin(9) * t198 + t149;
t104 = -pkin(3) * t224 + t112;
t348 = t242 * t189;
t150 = t212 * t245 + t348;
t113 = -pkin(9) * t196 + t150;
t107 = t241 * t113;
t55 = t244 * t104 - t107;
t42 = t55 - t395;
t37 = -pkin(4) * t217 + t42;
t394 = pkin(10) * t282;
t109 = t244 * t113;
t56 = t104 * t241 + t109;
t43 = t56 + t394;
t39 = t240 * t43;
t20 = t371 * t37 - t39;
t392 = qJ(6) * t376;
t14 = t20 - t392;
t209 = -qJD(5) + t217;
t13 = -pkin(5) * t209 + t14;
t41 = t371 * t43;
t21 = t240 * t37 + t41;
t15 = t21 - t400;
t403 = -t13 * t86 + t15 * t376;
t315 = qJD(1) * qJD(2);
t302 = t246 * t315;
t322 = qJD(3) * t242;
t307 = t243 * t322;
t314 = qJD(2) * qJD(3);
t158 = -qJD(1) * t307 + (t302 + t314) * t245;
t323 = qJD(2) * t246;
t308 = t242 * t323;
t321 = qJD(3) * t245;
t265 = t243 * t321 + t308;
t159 = qJD(1) * t265 + t242 * t314;
t320 = qJD(4) * t241;
t276 = t241 * t158 + t244 * t159 - t196 * t320 + t198 * t318;
t303 = qJD(5) * t371;
t317 = qJD(5) * t240;
t75 = t244 * t158 - t241 * t159 - t196 * t318 - t198 * t320;
t28 = t240 * t276 + t281 * t317 - t282 * t303 - t371 * t75;
t398 = -t209 * t86 - t28;
t211 = -qJD(2) * pkin(2) + pkin(7) * t328;
t161 = pkin(3) * t196 + t211;
t102 = -pkin(4) * t282 + t161;
t228 = t243 * t315;
t206 = t283 * qJD(2);
t190 = qJD(1) * t206;
t285 = pkin(7) * t228;
t335 = -t245 * t190 - t242 * t285;
t261 = -qJD(3) * t150 - t335;
t65 = pkin(3) * t228 - pkin(9) * t158 + t261;
t275 = t189 * t321 + t242 * t190 - t212 * t322;
t258 = -t245 * t285 + t275;
t80 = -pkin(9) * t159 + t258;
t288 = -t104 * t318 + t113 * t320 - t241 * t65 - t244 * t80;
t12 = -pkin(10) * t276 - t288;
t260 = -qJD(4) * t56 - t241 * t80 + t244 * t65;
t9 = pkin(4) * t228 - pkin(10) * t75 + t260;
t300 = -t371 * t12 - t240 * t9 - t37 * t303 + t43 * t317;
t397 = t102 * t86 + t300;
t264 = t371 * t279;
t359 = -qJD(5) * t264 - t199 * t317 + t240 * t404 + t336 * t371;
t146 = t371 * t199 - t240 * t279;
t358 = t146 * qJD(5) + t240 * t336 - t371 * t404;
t370 = pkin(3) * t242;
t192 = t327 * t370 + t233;
t387 = pkin(3) * t322 - t192;
t257 = -qJD(5) * t21 - t240 * t12 + t371 * t9;
t381 = -t102 * t376 + t257;
t29 = qJD(5) * t376 + t240 * t75 + t371 * t276;
t378 = -t209 * t376 - t29;
t396 = -0.2e1 * t315;
t47 = pkin(5) * t86 + qJD(6) + t102;
t393 = t376 * t47;
t391 = t281 * t282;
t289 = -t244 * t213 - t349;
t122 = t289 - t368;
t123 = -pkin(10) * t279 + t332;
t339 = t240 * t122 + t371 * t123;
t390 = -t404 * pkin(4) + t387;
t230 = pkin(3) * t244 + pkin(4);
t350 = t240 * t241;
t292 = -t112 * t241 - t109;
t48 = t292 - t394;
t340 = t244 * t112 - t107;
t49 = t340 - t395;
t389 = t230 * t303 + (-t241 * t317 + (t371 * t244 - t350) * qJD(4)) * pkin(3) - t240 * t48 - t371 * t49;
t310 = t371 * t241;
t388 = t230 * t317 - (-t241 * t303 + (-t240 * t244 - t310) * qJD(4)) * pkin(3) - t240 * t49 + t371 * t48;
t195 = t245 * t208;
t369 = pkin(7) * t242;
t148 = -pkin(9) * t345 + t195 + (-pkin(3) - t369) * t246;
t226 = pkin(7) * t344;
t330 = t242 * t208 + t226;
t347 = t242 * t243;
t154 = -pkin(9) * t347 + t330;
t337 = t241 * t148 + t244 * t154;
t386 = qJD(5) * t339 + t406 * t240 + t407 * t371;
t326 = qJD(2) * t199;
t263 = t246 * t326;
t249 = t243 * t399 - t263;
t385 = t122 * t303 - t123 * t317 - t407 * t240 + t406 * t371;
t383 = t281 ^ 2 - t282 ^ 2;
t382 = t217 * t282 + t75;
t380 = -t161 * t281 + t260;
t379 = -t161 * t282 + t288;
t377 = -t217 * t281 - t276;
t145 = t199 * t240 + t264;
t373 = -qJ(6) * t358 - qJD(6) * t145 + t385;
t372 = -pkin(5) * t328 - qJ(6) * t359 - t146 * qJD(6) - t386;
t367 = t13 - t14;
t366 = t371 * t42 - t39;
t172 = t279 * t243;
t290 = t244 * t148 - t241 * t154;
t78 = -pkin(4) * t246 + pkin(10) * t172 + t290;
t270 = t199 * t243;
t81 = -pkin(10) * t270 + t337;
t363 = t240 * t78 + t371 * t81;
t362 = t75 * t199;
t361 = t389 + t392;
t360 = -t388 - t400;
t137 = t159 * pkin(3) + pkin(7) * t302;
t357 = t137 * t199;
t356 = t158 * t242;
t355 = t196 * t224;
t354 = t198 * t224;
t353 = t211 * t242;
t352 = t211 * t245;
t351 = t224 * t245;
t248 = qJD(1) ^ 2;
t343 = t246 * t248;
t247 = qJD(2) ^ 2;
t342 = t247 * t243;
t341 = t247 * t246;
t334 = t242 * t206 + t208 * t321;
t324 = qJD(2) * t243;
t333 = t245 * t206 + t324 * t369;
t207 = pkin(3) * t347 + t243 * pkin(7);
t238 = t243 ^ 2;
t329 = -t246 ^ 2 + t238;
t319 = qJD(4) * t242;
t94 = t277 * qJD(2) + (-t226 + (pkin(9) * t243 - t208) * t242) * qJD(3) + t333;
t306 = t246 * t322;
t99 = -t265 * pkin(9) + (-t243 * t316 - t306) * pkin(7) + t334;
t313 = t148 * t318 + t241 * t94 + t244 * t99;
t162 = pkin(3) * t265 + pkin(7) * t323;
t231 = -t245 * pkin(3) - pkin(2);
t304 = t224 * t321;
t299 = -t240 * t42 - t41;
t295 = -t240 * t81 + t371 * t78;
t293 = pkin(1) * t396;
t291 = t371 * t122 - t123 * t240;
t287 = t196 + t316;
t286 = -t198 + t325;
t284 = -pkin(3) * t350 + t371 * t230;
t106 = pkin(3) * t198 + pkin(4) * t281;
t278 = qJD(1) * t238 - t224 * t246;
t267 = t279 * qJD(2);
t100 = t243 * t253 - t246 * t267;
t259 = -qJD(4) * t337 - t241 * t99 + t244 * t94;
t24 = pkin(4) * t324 - pkin(10) * t100 + t259;
t266 = t246 * t316 - t307;
t26 = (-t345 * t384 - t308) * pkin(10) * t244 + (-qJD(4) * t154 + (t243 * t319 - t266) * pkin(10)) * t241 + t313;
t274 = t240 * t24 + t371 * t26 + t78 * t303 - t317 * t81;
t271 = t161 * t199;
t262 = t371 * t270;
t167 = pkin(4) * t279 + t231;
t151 = pkin(4) * t270 + t207;
t50 = pkin(4) * t276 + t137;
t256 = -t363 * qJD(5) + t371 * t24 - t240 * t26;
t116 = -t371 * t172 - t240 * t270;
t11 = t29 * pkin(5) + t50;
t251 = t199 * t384 - t165;
t82 = -pkin(4) * t249 + t162;
t229 = t371 * pkin(4) + pkin(5);
t178 = pkin(3) * t310 + t240 * t230;
t173 = pkin(5) + t284;
t115 = -t172 * t240 + t262;
t52 = -qJ(6) * t145 + t339;
t51 = -qJ(6) * t146 + t291;
t34 = t116 * qJD(5) + t240 * t100 - t371 * t249;
t33 = qJD(5) * t262 - t371 * t100 - t172 * t317 - t240 * t249;
t32 = -qJ(6) * t115 + t363;
t31 = -pkin(5) * t246 - qJ(6) * t116 + t295;
t17 = t366 - t392;
t16 = t299 + t400;
t4 = -qJ(6) * t34 - qJD(6) * t115 + t274;
t3 = pkin(5) * t324 + t33 * qJ(6) - t116 * qJD(6) + t256;
t2 = -qJ(6) * t29 - qJD(6) * t86 - t300;
t1 = pkin(5) * t228 + t28 * qJ(6) - qJD(6) * t376 + t257;
t5 = [0, 0, 0, 0.2e1 * t246 * t228, t329 * t396, t341, -t342, 0, -pkin(7) * t341 + t243 * t293, pkin(7) * t342 + t246 * t293, t158 * t345 + t198 * t266 (-t196 * t245 - t198 * t242) * t323 + (-t356 - t159 * t245 + (t196 * t242 - t198 * t245) * qJD(3)) * t243, t224 * t307 - t158 * t246 + (t198 * t243 + t245 * t278) * qJD(2), t243 * t304 + t159 * t246 + (-t196 * t243 - t242 * t278) * qJD(2) (-t224 - t327) * t324 -(-t208 * t322 + t333) * t224 + (t211 * t321 + pkin(7) * t159 + (qJD(1) * t195 + t149) * qJD(2)) * t243 + ((pkin(7) * t196 + t353) * qJD(2) + (t348 + (pkin(7) * t224 + t212) * t245) * qJD(3) + t335) * t246 (-pkin(7) * t306 + t334) * t224 + t275 * t246 + (pkin(7) * t158 - t211 * t322) * t243 + ((pkin(7) * t198 + t352) * t246 + (-pkin(7) * t351 - qJD(1) * t330 - t150) * t243) * qJD(2), t100 * t281 - t172 * t75, t100 * t282 + t172 * t276 - t281 * t263 + (-t362 + t281 * (-t244 * t321 - t245 * t318 + (t319 + t322) * t241)) * t243, -t100 * t217 - t75 * t246 + (-qJD(1) * t172 + t281) * t324 (t217 * t326 + t276) * t246 + (-t399 * t217 + (-qJD(1) * t270 + t282) * qJD(2)) * t243 (-t217 - t327) * t324, -t259 * t217 - t162 * t282 + t207 * t276 + (qJD(2) * t271 - t260) * t246 + (t357 - t161 * t399 + (qJD(1) * t290 + t55) * qJD(2)) * t243 (-t154 * t320 + t313) * t217 - t288 * t246 + t162 * t281 + t207 * t75 - t137 * t172 + t161 * t100 + (-qJD(1) * t337 - t56) * t324, -t116 * t28 - t33 * t376, t115 * t28 - t116 * t29 + t33 * t86 - t34 * t376, t209 * t33 + t246 * t28 + (qJD(1) * t116 + t376) * t324, t209 * t34 + t246 * t29 + (-qJD(1) * t115 - t86) * t324 (-t209 - t327) * t324, t102 * t34 + t50 * t115 + t151 * t29 + t20 * t324 - t209 * t256 + t228 * t295 - t246 * t257 + t82 * t86, t274 * t209 - t300 * t246 + t82 * t376 - t151 * t28 + t50 * t116 - t102 * t33 + (-qJD(1) * t363 - t21) * t324, -t1 * t116 - t115 * t2 + t13 * t33 - t15 * t34 + t28 * t31 - t29 * t32 - t3 * t376 - t4 * t86, t2 * t32 + t15 * t4 + t1 * t31 + t13 * t3 + t11 * (t115 * pkin(5) + t151) + t47 * (t34 * pkin(5) + t82); 0, 0, 0, -t243 * t343, t329 * t248, 0, 0, 0, t248 * pkin(1) * t243, pkin(1) * t343, -t198 * t351 + t356 (t158 + t355) * t245 + (-t159 + t354) * t242, -t304 + (t224 * t344 + t243 * t286) * qJD(1), t224 * t322 + (-t224 * t346 + t243 * t287) * qJD(1), t224 * t328, -pkin(2) * t159 + t331 * t224 + (pkin(8) * t351 + t353) * qJD(3) + ((-pkin(8) * t325 - t149) * t243 + (-pkin(7) * t287 - t353) * t246) * qJD(1), -pkin(2) * t158 - t183 * t224 + (-pkin(8) * t224 * t242 + t352) * qJD(3) + (-t211 * t344 + (-pkin(8) * t316 + t150) * t243 + (t224 * t345 + t246 * t286) * pkin(7)) * qJD(1), t281 * t336 + t362, -t199 * t276 - t75 * t279 + t404 * t281 + t336 * t282, -t336 * t217 + (-t281 + t326) * t328, t251 * t217 + (-t267 - t282) * t328, t217 * t328, t231 * t276 + t137 * t279 + t192 * t282 - t161 * t165 + t410 * t217 + (t217 * t332 + t271) * qJD(4) + (-t282 * t370 + t271) * qJD(3) + (qJD(2) * t289 - t55) * t328, t357 + t231 * t75 + (-t214 * t320 + t408) * t217 + t336 * t161 + t387 * t281 + (-qJD(2) * t332 + t56) * t328, -t28 * t146 + t359 * t376, t145 * t28 - t146 * t29 - t358 * t376 - t359 * t86, -t359 * t209 + (qJD(2) * t146 - t376) * t328, t358 * t209 + (-qJD(2) * t145 + t86) * t328, t209 * t328, t102 * t358 + t50 * t145 + t167 * t29 - t20 * t328 + t209 * t386 + t228 * t291 + t390 * t86, t50 * t146 - t167 * t28 + t390 * t376 + t385 * t209 + t359 * t102 + (-qJD(2) * t339 + t21) * t328, -t1 * t146 - t359 * t13 - t145 * t2 - t358 * t15 + t28 * t51 - t29 * t52 - t372 * t376 - t373 * t86, t2 * t52 + t1 * t51 + t11 * (t145 * pkin(5) + t167) + (pkin(4) * t251 + pkin(5) * t358 + t387) * t47 + t373 * t15 + t372 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198 * t196, -t196 ^ 2 + t198 ^ 2, t158 - t355, -t159 - t354, t228, -t150 * t224 - t198 * t211 + t261, -t149 * t224 + t196 * t211 - t258, -t391, t383, t382, t377, t228, t292 * t217 + (t198 * t282 + t217 * t320 + t228 * t244) * pkin(3) + t380, -t340 * t217 + (-t198 * t281 + t217 * t318 - t228 * t241) * pkin(3) + t379, t401, t405, t398, t378, t228, -t106 * t86 + t209 * t388 + t284 * t228 + t381, -t106 * t376 - t178 * t228 + t209 * t389 + t397, t173 * t28 - t178 * t29 - t360 * t376 - t361 * t86 + t403, t2 * t178 + t1 * t173 - t47 * (pkin(5) * t376 + t106) + t361 * t15 + t360 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, t383, t382, t377, t228, -t217 * t56 + t380, -t217 * t55 + t379, t401, t405, t398, t378, t228, t299 * t209 + (t209 * t317 + t371 * t228 - t281 * t86) * pkin(4) + t381, -t366 * t209 + (t209 * t303 - t228 * t240 - t281 * t376) * pkin(4) + t397, t16 * t376 + t17 * t86 + t229 * t28 + (-t240 * t29 + (t240 * t376 - t371 * t86) * qJD(5)) * pkin(4) + t403, -pkin(5) * t393 + t1 * t229 - t13 * t16 - t15 * t17 + (-t47 * t281 + t2 * t240 + (-t13 * t240 + t371 * t15) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t401, t405, t398, t378, t228, -t21 * t209 + t381, -t20 * t209 + t397, pkin(5) * t28 - t367 * t86, t367 * t15 + (t1 - t393) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 - t375, t13 * t376 + t15 * t86 + t11;];
tauc_reg  = t5;