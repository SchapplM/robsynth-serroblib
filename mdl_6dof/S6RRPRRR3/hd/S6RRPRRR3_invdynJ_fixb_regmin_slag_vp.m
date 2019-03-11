% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:23
% EndTime: 2019-03-09 13:25:43
% DurationCPUTime: 8.41s
% Computational Cost: add. (11511->541), mult. (27083->726), div. (0->0), fcn. (21487->18), ass. (0->282)
t262 = sin(pkin(11));
t268 = sin(qJ(2));
t343 = qJD(1) * t268;
t263 = cos(pkin(11));
t273 = cos(qJ(2));
t358 = t263 * t273;
t210 = qJD(1) * t358 - t262 * t343;
t398 = qJD(4) + qJD(5);
t423 = t210 - t398;
t223 = t262 * t273 + t263 * t268;
t212 = t223 * qJD(1);
t267 = sin(qJ(4));
t272 = cos(qJ(4));
t336 = t272 * qJD(2);
t172 = t212 * t267 - t336;
t174 = qJD(2) * t267 + t212 * t272;
t266 = sin(qJ(5));
t271 = cos(qJ(5));
t108 = t271 * t172 + t174 * t266;
t265 = sin(qJ(6));
t270 = cos(qJ(6));
t299 = t172 * t266 - t271 * t174;
t300 = t108 * t265 + t270 * t299;
t58 = t270 * t108 - t265 * t299;
t415 = t300 * t58;
t225 = t266 * t267 - t271 * t272;
t347 = t423 * t225;
t356 = t266 * t272;
t226 = t267 * t271 + t356;
t346 = t423 * t226;
t411 = t300 ^ 2 - t58 ^ 2;
t200 = qJD(4) - t210;
t194 = qJD(5) + t200;
t190 = qJD(6) + t194;
t335 = qJD(1) * qJD(2);
t325 = t273 * t335;
t326 = t268 * t335;
t157 = qJDD(1) * t223 - t262 * t326 + t263 * t325;
t342 = qJD(4) * t267;
t100 = qJD(4) * t336 + t267 * qJDD(2) + t272 * t157 - t212 * t342;
t101 = qJD(4) * t174 - t272 * qJDD(2) + t267 * t157;
t279 = qJD(5) * t299 - t266 * t100 - t271 * t101;
t337 = qJD(6) * t270;
t338 = qJD(6) * t265;
t339 = qJD(5) * t271;
t340 = qJD(5) * t266;
t37 = t271 * t100 - t266 * t101 - t172 * t339 - t174 * t340;
t8 = -t108 * t337 + t265 * t279 + t270 * t37 + t299 * t338;
t409 = t190 * t58 + t8;
t261 = qJ(4) + qJ(5);
t257 = qJ(6) + t261;
t246 = sin(t257);
t247 = cos(t257);
t274 = cos(qJ(1));
t258 = qJ(2) + pkin(11);
t251 = cos(t258);
t269 = sin(qJ(1));
t365 = t251 * t269;
t177 = t246 * t274 - t247 * t365;
t364 = t251 * t274;
t179 = t246 * t269 + t247 * t364;
t249 = pkin(2) * t273 + pkin(1);
t233 = -qJD(1) * t249 + qJD(3);
t120 = -t210 * pkin(3) - t212 * pkin(8) + t233;
t391 = qJ(3) + pkin(7);
t235 = t391 * t268;
t227 = qJD(1) * t235;
t385 = qJD(2) * pkin(2);
t219 = -t227 + t385;
t236 = t391 * t273;
t228 = qJD(1) * t236;
t359 = t263 * t228;
t155 = t262 * t219 + t359;
t143 = qJD(2) * pkin(8) + t155;
t81 = t272 * t120 - t143 * t267;
t64 = -pkin(9) * t174 + t81;
t44 = pkin(4) * t200 + t64;
t82 = t120 * t267 + t143 * t272;
t65 = -pkin(9) * t172 + t82;
t55 = t271 * t65;
t27 = t266 * t44 + t55;
t417 = pkin(10) * t108;
t21 = t27 - t417;
t19 = t21 * t338;
t250 = sin(t258);
t395 = g(3) * t250;
t215 = t262 * t228;
t154 = t219 * t263 - t215;
t142 = -qJD(2) * pkin(3) - t154;
t102 = pkin(4) * t172 + t142;
t48 = pkin(5) * t108 + t102;
t408 = g(1) * t179 - g(2) * t177 + t247 * t395 + t48 * t58 + t19;
t176 = t246 * t365 + t247 * t274;
t178 = -t246 * t364 + t247 * t269;
t211 = t223 * qJD(2);
t333 = t273 * qJDD(1);
t334 = t268 * qJDD(1);
t302 = -t262 * t334 + t263 * t333;
t151 = qJD(1) * t211 + qJDD(4) - t302;
t146 = qJDD(5) + t151;
t156 = -qJD(2) * t212 + t302;
t288 = pkin(2) * t326 - qJDD(1) * t249 + qJDD(3);
t83 = -t156 * pkin(3) - t157 * pkin(8) + t288;
t80 = t272 * t83;
t318 = qJD(2) * t391;
t209 = -t268 * qJD(3) - t273 * t318;
t150 = qJDD(2) * pkin(2) + qJD(1) * t209 - qJDD(1) * t235;
t208 = t273 * qJD(3) - t268 * t318;
t159 = qJD(1) * t208 + qJDD(1) * t236;
t95 = t262 * t150 + t263 * t159;
t93 = qJDD(2) * pkin(8) + t95;
t12 = t151 * pkin(4) - t100 * pkin(9) - qJD(4) * t82 - t267 * t93 + t80;
t341 = qJD(4) * t272;
t292 = t120 * t341 - t143 * t342 + t267 * t83 + t272 * t93;
t17 = -pkin(9) * t101 + t292;
t323 = t271 * t12 - t266 * t17;
t281 = -qJD(5) * t27 + t323;
t2 = t146 * pkin(5) - t37 * pkin(10) + t281;
t316 = -t266 * t12 - t271 * t17 - t44 * t339 + t65 * t340;
t3 = pkin(10) * t279 - t316;
t330 = t270 * t2 - t265 * t3;
t422 = -g(1) * t178 + g(2) * t176 + t246 * t395 + t48 * t300 + t330;
t280 = qJD(6) * t300 - t265 * t37 + t270 * t279;
t403 = -t190 * t300 + t280;
t242 = pkin(2) * t262 + pkin(8);
t392 = pkin(9) + t242;
t317 = qJD(4) * t392;
t132 = pkin(2) * t343 + pkin(3) * t212 - pkin(8) * t210;
t162 = -t227 * t263 - t215;
t349 = t267 * t132 + t272 * t162;
t371 = t210 * t267;
t421 = -pkin(9) * t371 + t267 * t317 + t349;
t124 = t272 * t132;
t420 = pkin(4) * t212 - t162 * t267 + t124 + (-pkin(9) * t210 + t317) * t272;
t306 = g(1) * t274 + g(2) * t269;
t286 = -g(3) * t251 + t250 * t306;
t94 = t150 * t263 - t262 * t159;
t92 = -qJDD(2) * pkin(3) - t94;
t419 = -qJD(4) * t242 * t200 + t286 - t92;
t418 = t342 - t371;
t416 = pkin(10) * t299;
t163 = -t225 * t265 + t226 * t270;
t387 = qJD(6) * t163 + t265 * t347 - t270 * t346;
t413 = t299 * t108;
t328 = t223 * t341;
t222 = t262 * t268 - t358;
t214 = t222 * qJD(2);
t354 = t267 * t214;
t412 = t328 - t354;
t410 = -t108 ^ 2 + t299 ^ 2;
t407 = t108 * t194 + t37;
t53 = t266 * t65;
t26 = t271 * t44 - t53;
t20 = t26 + t416;
t18 = pkin(5) * t194 + t20;
t383 = t270 * t21;
t7 = t265 * t18 + t383;
t406 = -qJD(6) * t7 + t422;
t256 = cos(t261);
t361 = t256 * t269;
t255 = sin(t261);
t362 = t255 * t274;
t181 = -t251 * t361 + t362;
t360 = t256 * t274;
t363 = t255 * t269;
t183 = t251 * t360 + t363;
t405 = g(1) * t183 - g(2) * t181 + t102 * t108 + t256 * t395 + t316;
t180 = t251 * t363 + t360;
t182 = -t251 * t362 + t361;
t404 = -g(1) * t182 + g(2) * t180 + t102 * t299 + t255 * t395 + t281;
t402 = -t194 * t299 + t279;
t401 = t420 * t271;
t138 = t226 * t223;
t160 = -t227 * t262 + t359;
t309 = pkin(4) * t418 - t160;
t220 = t392 * t267;
t221 = t392 * t272;
t345 = -t266 * t220 + t271 * t221;
t399 = t220 * t339 + t221 * t340 + t420 * t266 + t271 * t421;
t137 = qJDD(6) + t146;
t161 = t270 * t225 + t226 * t265;
t388 = -qJD(6) * t161 + t265 * t346 + t270 * t347;
t397 = -t163 * t137 - t190 * t388;
t396 = -t226 * t146 - t194 * t347;
t393 = g(3) * t273;
t390 = t271 * t64 - t53;
t153 = pkin(3) * t222 - pkin(8) * t223 - t249;
t141 = t272 * t153;
t171 = -t235 * t262 + t236 * t263;
t368 = t223 * t272;
t74 = pkin(4) * t222 - pkin(9) * t368 - t171 * t267 + t141;
t164 = t272 * t171;
t348 = t267 * t153 + t164;
t369 = t223 * t267;
t85 = -pkin(9) * t369 + t348;
t386 = t266 * t74 + t271 * t85;
t384 = t212 * t58;
t382 = t300 * t212;
t381 = -pkin(5) * t346 + t309;
t380 = t100 * t267;
t379 = t108 * t212;
t378 = t299 * t212;
t377 = t137 * t266;
t375 = t172 * t200;
t374 = t172 * t212;
t373 = t174 * t200;
t372 = t174 * t212;
t370 = t214 * t272;
t248 = pkin(4) * t271 + pkin(5);
t366 = t248 * t137;
t357 = t266 * t270;
t355 = t267 * t151;
t353 = t267 * t269;
t352 = t267 * t274;
t351 = t269 * t272;
t136 = t272 * t151;
t350 = t272 * t274;
t259 = t268 ^ 2;
t344 = -t273 ^ 2 + t259;
t332 = t268 * t385;
t243 = -pkin(2) * t263 - pkin(3);
t329 = t223 * t342;
t324 = qJD(6) * t18 + t3;
t133 = pkin(3) * t211 + pkin(8) * t214 + t332;
t125 = t272 * t133;
t131 = t208 * t263 + t209 * t262;
t32 = pkin(9) * t370 + t211 * pkin(4) - t267 * t131 + t125 + (-t164 + (pkin(9) * t223 - t153) * t267) * qJD(4);
t291 = t272 * t131 + t267 * t133 + t153 * t341 - t171 * t342;
t40 = -pkin(9) * t412 + t291;
t321 = -t266 * t40 + t271 * t32;
t320 = -t266 * t64 - t55;
t319 = -t266 * t85 + t271 * t74;
t314 = -qJD(4) * t120 - t93;
t130 = t208 * t262 - t263 * t209;
t312 = -t271 * t220 - t221 * t266;
t170 = t263 * t235 + t236 * t262;
t311 = t200 * t272;
t310 = -t161 * t137 - t190 * t387;
t115 = -pkin(10) * t226 + t312;
t308 = -pkin(10) * t346 - qJD(6) * t115 + t399;
t116 = -pkin(10) * t225 + t345;
t307 = pkin(5) * t212 + pkin(10) * t347 + t345 * qJD(5) + qJD(6) * t116 - t266 * t421 + t401;
t305 = g(1) * t269 - g(2) * t274;
t304 = -t143 * t341 + t80;
t303 = -t225 * t146 + t194 * t346;
t128 = pkin(4) * t369 + t170;
t139 = t225 * t223;
t86 = t270 * t138 - t139 * t265;
t87 = -t138 * t265 - t139 * t270;
t234 = -pkin(4) * t272 + t243;
t91 = pkin(4) * t412 + t130;
t298 = -t200 * t418 + t136;
t297 = -0.2e1 * pkin(1) * t335 - pkin(7) * qJDD(2);
t296 = t266 * t32 + t271 * t40 + t74 * t339 - t340 * t85;
t294 = -t329 - t370;
t289 = t142 * t200 - t242 * t151;
t41 = pkin(4) * t101 + t92;
t275 = qJD(2) ^ 2;
t283 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t275 + t305;
t276 = qJD(1) ^ 2;
t282 = pkin(1) * t276 - pkin(7) * qJDD(1) + t306;
t204 = t251 * t350 + t353;
t203 = -t251 * t352 + t351;
t202 = -t251 * t351 + t352;
t201 = t251 * t353 + t350;
t175 = pkin(5) * t225 + t234;
t89 = pkin(5) * t138 + t128;
t88 = pkin(4) * t174 - pkin(5) * t299;
t47 = -t214 * t356 - t266 * t329 - t340 * t369 + (t368 * t398 - t354) * t271;
t46 = -t138 * t398 + t225 * t214;
t36 = pkin(5) * t47 + t91;
t29 = -pkin(10) * t138 + t386;
t28 = pkin(5) * t222 + pkin(10) * t139 + t319;
t23 = t390 + t416;
t22 = t320 + t417;
t15 = qJD(6) * t87 + t265 * t46 + t270 * t47;
t14 = -qJD(6) * t86 - t265 * t47 + t270 * t46;
t13 = -pkin(5) * t279 + t41;
t6 = t18 * t270 - t21 * t265;
t5 = -pkin(10) * t47 + t296;
t4 = t211 * pkin(5) - t46 * pkin(10) - qJD(5) * t386 + t321;
t1 = [qJDD(1), t305, t306, qJDD(1) * t259 + 0.2e1 * t268 * t325, 0.2e1 * t268 * t333 - 0.2e1 * t335 * t344, qJDD(2) * t268 + t273 * t275, qJDD(2) * t273 - t268 * t275, 0, t268 * t297 + t273 * t283, -t268 * t283 + t273 * t297, t130 * t212 + t131 * t210 + t154 * t214 - t155 * t211 + t156 * t171 + t157 * t170 - t222 * t95 - t223 * t94 - t306, t95 * t171 + t155 * t131 - t94 * t170 - t154 * t130 - t288 * t249 + t233 * t332 - g(1) * (-t249 * t269 + t274 * t391) - g(2) * (t249 * t274 + t269 * t391) t100 * t368 + t174 * t294 -(-t172 * t272 - t174 * t267) * t214 + (-t380 - t101 * t272 + (t172 * t267 - t174 * t272) * qJD(4)) * t223, t100 * t222 + t136 * t223 + t174 * t211 + t200 * t294, -t101 * t222 - t172 * t211 - t200 * t412 - t223 * t355, t151 * t222 + t200 * t211 (-t171 * t341 + t125) * t200 + t141 * t151 + t304 * t222 + t81 * t211 + t130 * t172 + t170 * t101 + t142 * t328 - g(1) * t202 - g(2) * t204 + ((-qJD(4) * t153 - t131) * t200 - t171 * t151 + t314 * t222 + t92 * t223 - t142 * t214) * t267, -g(1) * t201 - g(2) * t203 + t170 * t100 + t130 * t174 + t142 * t294 - t151 * t348 - t200 * t291 - t82 * t211 - t222 * t292 + t368 * t92, -t139 * t37 - t299 * t46, -t108 * t46 - t138 * t37 - t139 * t279 + t299 * t47, -t139 * t146 + t194 * t46 - t211 * t299 + t222 * t37, -t108 * t211 - t138 * t146 - t194 * t47 + t222 * t279, t146 * t222 + t194 * t211, t91 * t108 - t128 * t279 + t41 * t138 + t102 * t47 + t321 * t194 + t319 * t146 + t323 * t222 + t26 * t211 - g(1) * t181 - g(2) * t183 + (-t194 * t386 - t222 * t27) * qJD(5), -g(1) * t180 - g(2) * t182 + t102 * t46 + t128 * t37 - t41 * t139 - t146 * t386 - t194 * t296 - t27 * t211 + t222 * t316 - t299 * t91, -t14 * t300 + t8 * t87, -t14 * t58 + t15 * t300 + t280 * t87 - t8 * t86, t137 * t87 + t14 * t190 - t211 * t300 + t222 * t8, -t137 * t86 - t15 * t190 - t211 * t58 + t222 * t280, t137 * t222 + t190 * t211 (-t265 * t5 + t270 * t4) * t190 + (-t265 * t29 + t270 * t28) * t137 + t330 * t222 + t6 * t211 + t36 * t58 - t89 * t280 + t13 * t86 + t48 * t15 - g(1) * t177 - g(2) * t179 + ((-t265 * t28 - t270 * t29) * t190 - t7 * t222) * qJD(6), -g(1) * t176 - g(2) * t178 + t13 * t87 + t48 * t14 + t19 * t222 - t7 * t211 - t36 * t300 + t89 * t8 + (-(-qJD(6) * t29 + t4) * t190 - t28 * t137 - t2 * t222) * t265 + (-(qJD(6) * t28 + t5) * t190 - t29 * t137 - t324 * t222) * t270; 0, 0, 0, -t268 * t276 * t273, t344 * t276, t334, t333, qJDD(2), t268 * t282 - t393, g(3) * t268 + t273 * t282 (t155 - t160) * t212 + (t154 - t162) * t210 + (t156 * t262 - t157 * t263) * pkin(2), t154 * t160 - t155 * t162 + (-t393 + t262 * t95 + t263 * t94 + (-qJD(1) * t233 + t306) * t268) * pkin(2), t174 * t311 + t380 (t100 - t375) * t272 + (-t101 - t373) * t267, t200 * t311 + t355 - t372, t298 + t374, -t200 * t212, t243 * t101 - t124 * t200 - t160 * t172 - t81 * t212 + (t162 * t200 + t289) * t267 + t419 * t272, t243 * t100 - t160 * t174 + t349 * t200 + t82 * t212 - t267 * t419 + t289 * t272, t37 * t226 - t299 * t347, -t108 * t347 - t37 * t225 + t226 * t279 - t299 * t346, t378 - t396, t303 + t379, -t194 * t212, -t234 * t279 + t41 * t225 + t312 * t146 - t26 * t212 + (-t221 * t339 + (qJD(5) * t220 + t421) * t266 - t401) * t194 + t309 * t108 - t346 * t102 + t286 * t256, t347 * t102 - t345 * t146 + t194 * t399 + t27 * t212 + t41 * t226 + t234 * t37 - t255 * t286 - t299 * t309, t8 * t163 - t300 * t388, -t8 * t161 + t163 * t280 + t300 * t387 - t388 * t58, t382 - t397, t310 + t384, -t190 * t212 (t115 * t270 - t116 * t265) * t137 - t175 * t280 + t13 * t161 - t6 * t212 + t381 * t58 + t387 * t48 + (t265 * t308 - t270 * t307) * t190 + t286 * t247 -(t115 * t265 + t116 * t270) * t137 + t175 * t8 + t13 * t163 + t7 * t212 - t381 * t300 + t388 * t48 + (t265 * t307 + t270 * t308) * t190 - t286 * t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210 ^ 2 - t212 ^ 2, t154 * t212 - t155 * t210 + t288 - t305, 0, 0, 0, 0, 0, t298 - t374, -t200 ^ 2 * t272 - t355 - t372, 0, 0, 0, 0, 0, t303 - t379, t378 + t396, 0, 0, 0, 0, 0, t310 - t384, t382 + t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174 * t172, -t172 ^ 2 + t174 ^ 2, t100 + t375, -t101 + t373, t151, -g(1) * t203 + g(2) * t201 - t142 * t174 + t82 * t200 + (t314 + t395) * t267 + t304, g(1) * t204 - g(2) * t202 + t142 * t172 + t200 * t81 + t272 * t395 - t292, -t413, t410, t407, t402, t146, -t320 * t194 + (-t108 * t174 + t146 * t271 - t194 * t340) * pkin(4) + t404, t390 * t194 + (-t146 * t266 + t174 * t299 - t194 * t339) * pkin(4) + t405, -t415, t411, t409, t403, t137, t270 * t366 - (t22 * t270 - t23 * t265) * t190 - t88 * t58 + (-t265 * t377 + (-t265 * t271 - t357) * t190 * qJD(5)) * pkin(4) + ((-pkin(4) * t357 - t248 * t265) * t190 - t7) * qJD(6) + t422, t88 * t300 + (-t366 - t2 + (t22 - (-qJD(5) - qJD(6)) * t266 * pkin(4)) * t190) * t265 + (-pkin(4) * t377 + (-pkin(4) * t339 - qJD(6) * t248 + t23) * t190 - t324) * t270 + t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, t410, t407, t402, t146, t27 * t194 + t404, t194 * t26 + t405, -t415, t411, t409, t403, t137 -(-t20 * t265 - t383) * t190 + (t137 * t270 - t190 * t338 + t299 * t58) * pkin(5) + t406 (-t190 * t21 - t2) * t265 + (t190 * t20 - t324) * t270 + (-t137 * t265 - t190 * t337 - t299 * t300) * pkin(5) + t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t415, t411, t409, t403, t137, t7 * t190 + t406, t6 * t190 - t265 * t2 - t270 * t324 + t408;];
tau_reg  = t1;
