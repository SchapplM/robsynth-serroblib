% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:49
% EndTime: 2019-03-10 01:56:19
% DurationCPUTime: 11.58s
% Computational Cost: add. (19465->610), mult. (49634->821), div. (0->0), fcn. (39390->10), ass. (0->275)
t262 = sin(qJ(2));
t265 = cos(qJ(3));
t266 = cos(qJ(2));
t257 = sin(pkin(6));
t351 = qJD(1) * t262;
t327 = t257 * t351;
t258 = cos(pkin(6));
t352 = qJD(1) * t258;
t335 = pkin(1) * t352;
t204 = -pkin(8) * t327 + t266 * t335;
t283 = (pkin(2) * t262 - pkin(9) * t266) * t257;
t205 = qJD(1) * t283;
t261 = sin(qJ(3));
t310 = -t204 * t261 + t265 * t205;
t400 = -pkin(10) - pkin(9);
t328 = qJD(3) * t400;
t353 = qJD(1) * t257;
t363 = t265 * t266;
t441 = (pkin(3) * t262 - pkin(10) * t363) * t353 + t310 - t265 * t328;
t350 = qJD(1) * t266;
t326 = t257 * t350;
t303 = t261 * t326;
t355 = t265 * t204 + t261 * t205;
t440 = pkin(10) * t303 + t261 * t328 - t355;
t260 = sin(qJ(4));
t264 = cos(qJ(4));
t216 = t260 * t261 - t264 * t265;
t338 = qJD(3) + qJD(4);
t164 = t338 * t216;
t175 = t216 * t326;
t359 = -t164 + t175;
t217 = t260 * t265 + t261 * t264;
t358 = (-t326 + t338) * t217;
t263 = cos(qJ(5));
t342 = qJD(5) * t263;
t245 = qJD(2) + t352;
t189 = t245 * t265 - t261 * t327;
t190 = t245 * t261 + t265 * t327;
t289 = t264 * t189 - t190 * t260;
t428 = t289 * t263;
t439 = t342 - t428;
t237 = t400 * t261;
t238 = t400 * t265;
t287 = t237 * t264 + t238 * t260;
t438 = qJD(4) * t287 - t260 * t441 + t440 * t264;
t242 = t262 * t335;
t207 = pkin(8) * t326 + t242;
t347 = qJD(3) * t261;
t412 = -t207 + (-t303 + t347) * pkin(3);
t259 = sin(qJ(5));
t236 = -qJD(3) + t326;
t284 = -qJD(4) + t236;
t288 = t189 * t260 + t264 * t190;
t122 = t259 * t288 + t263 * t284;
t124 = -t259 * t284 + t263 * t288;
t423 = -qJD(5) + t289;
t434 = t423 * t259;
t339 = qJD(1) * qJD(2);
t319 = t257 * t339;
t239 = t262 * t319;
t343 = qJD(5) * t259;
t301 = qJD(3) * t327;
t302 = t266 * t319;
t346 = qJD(3) * t265;
t161 = t245 * t346 - t261 * t301 + t265 * t302;
t329 = t245 * t347 + t261 * t302 + t265 * t301;
t272 = t264 * t161 - t260 * t329;
t268 = qJD(4) * t289 + t272;
t424 = -qJD(5) * t284 + t268;
t48 = -t259 * t239 - t263 * t424 + t288 * t343;
t49 = -t263 * t239 + t259 * t424 + t288 * t342;
t437 = -t122 * t439 + t124 * t434 - t259 * t49 - t48 * t263;
t312 = t260 * t161 + t264 * t329;
t85 = qJD(4) * t288 + t312;
t81 = t263 * t85;
t436 = t122 * t288 - t423 * t434 + t81;
t173 = pkin(9) * t245 + t207;
t202 = (-pkin(2) * t266 - pkin(9) * t262 - pkin(1)) * t257;
t184 = qJD(1) * t202;
t137 = -t173 * t261 + t265 * t184;
t112 = -pkin(10) * t190 + t137;
t103 = -pkin(3) * t236 + t112;
t138 = t173 * t265 + t184 * t261;
t113 = pkin(10) * t189 + t138;
t108 = t260 * t113;
t58 = t264 * t103 - t108;
t54 = pkin(4) * t284 - t58;
t30 = t122 * pkin(5) - t124 * qJ(6) + t54;
t435 = t30 * t423;
t429 = pkin(11) * t327 - t438;
t433 = -t358 * pkin(4) + pkin(11) * t359 - t412;
t46 = t48 * t259;
t431 = t124 * t439 - t46;
t383 = t259 * t85 - t342 * t423;
t430 = -t124 * t288 + t423 * t428 + t383;
t427 = t288 * t289;
t177 = t237 * t260 - t238 * t264;
t425 = qJD(4) * t177 + t440 * t260 + t264 * t441;
t422 = t288 ^ 2 - t289 ^ 2;
t99 = pkin(4) * t288 - pkin(11) * t289;
t172 = -pkin(2) * t245 - t204;
t151 = -pkin(3) * t189 + t172;
t344 = qJD(4) * t264;
t345 = qJD(4) * t260;
t206 = qJD(2) * t283;
t197 = qJD(1) * t206;
t365 = t257 * t262;
t246 = pkin(8) * t365;
t398 = pkin(1) * t266;
t208 = (t258 * t398 - t246) * qJD(2);
t198 = qJD(1) * t208;
t270 = -qJD(3) * t138 + t265 * t197 - t261 * t198;
t69 = pkin(3) * t239 - pkin(10) * t161 + t270;
t280 = -t173 * t347 + t184 * t346 + t261 * t197 + t265 * t198;
t76 = -pkin(10) * t329 + t280;
t307 = -t103 * t344 + t113 * t345 - t260 * t69 - t264 * t76;
t421 = -t151 * t289 + t307;
t297 = pkin(5) * t259 - qJ(6) * t263;
t420 = pkin(5) * t343 - qJ(6) * t342 - qJD(6) * t259 - t289 * t297;
t419 = t236 * t289 + t272;
t389 = qJ(6) * t85;
t14 = pkin(11) * t239 - t307;
t199 = pkin(8) * t302 + qJD(2) * t242;
t136 = pkin(3) * t329 + t199;
t33 = t85 * pkin(4) - pkin(11) * t268 + t136;
t109 = t264 * t113;
t59 = t260 * t103 + t109;
t55 = -pkin(11) * t284 + t59;
t70 = -pkin(4) * t289 - pkin(11) * t288 + t151;
t5 = t263 * t14 + t259 * t33 + t70 * t342 - t343 * t55;
t2 = -qJD(6) * t423 + t389 + t5;
t317 = t259 * t14 - t263 * t33 + t55 * t342 + t70 * t343;
t399 = pkin(5) * t85;
t4 = t317 - t399;
t417 = t2 * t263 + t4 * t259;
t415 = pkin(5) * t288;
t384 = pkin(4) * t327 + t425;
t63 = t112 * t260 + t109;
t300 = pkin(3) * t345 - t63;
t414 = qJ(6) * t288;
t413 = t423 * t288;
t252 = -pkin(3) * t265 - pkin(2);
t162 = pkin(4) * t216 - pkin(11) * t217 + t252;
t357 = t259 * t162 + t263 * t177;
t364 = t257 * t266;
t201 = pkin(8) * t364 + (pkin(1) * t262 + pkin(9)) * t258;
t356 = t265 * t201 + t261 * t202;
t411 = -t162 * t342 + t177 * t343 + t259 * t433 + t263 * t429;
t24 = -t259 * t55 + t263 * t70;
t362 = qJD(6) - t24;
t20 = pkin(5) * t423 + t362;
t306 = -t103 * t345 - t113 * t344 - t260 * t76 + t264 * t69;
t15 = -pkin(4) * t239 - t306;
t9 = pkin(5) * t49 + qJ(6) * t48 - qJD(6) * t124 + t15;
t409 = t288 * t20 - t263 * t9 + t30 * t343;
t25 = t259 * t70 + t263 * t55;
t21 = -qJ(6) * t423 + t25;
t408 = -t288 * t21 - t259 * t9;
t407 = -t15 * t263 - t24 * t288 + t54 * t343;
t406 = t15 * t259 + t25 * t288 + t54 * t342;
t404 = -t151 * t288 + t306;
t403 = -t236 * t288 - t312;
t212 = t258 * t261 + t265 * t365;
t311 = -t201 * t261 + t265 * t202;
t120 = -pkin(3) * t364 - pkin(10) * t212 + t311;
t211 = -t258 * t265 + t261 * t365;
t130 = -pkin(10) * t211 + t356;
t324 = qJD(2) * t364;
t167 = -qJD(3) * t211 + t265 * t324;
t271 = -qJD(3) * t356 + t265 * t206 - t208 * t261;
t349 = qJD(2) * t262;
t325 = t257 * t349;
t87 = pkin(3) * t325 - pkin(10) * t167 + t271;
t166 = qJD(3) * t212 + t261 * t324;
t279 = -t201 * t347 + t202 * t346 + t261 * t206 + t265 * t208;
t92 = -pkin(10) * t166 + t279;
t281 = t120 * t344 - t130 * t345 + t260 * t87 + t264 * t92;
t17 = pkin(11) * t325 + t281;
t361 = t260 * t120 + t264 * t130;
t73 = -pkin(11) * t364 + t361;
t156 = t264 * t211 + t212 * t260;
t157 = -t211 * t260 + t212 * t264;
t200 = t246 + (-pkin(2) - t398) * t258;
t160 = pkin(3) * t211 + t200;
t94 = pkin(4) * t156 - pkin(11) * t157 + t160;
t291 = t259 * t94 + t263 * t73;
t209 = t258 * pkin(1) * t349 + pkin(8) * t324;
t152 = pkin(3) * t166 + t209;
t97 = -qJD(4) * t156 - t166 * t260 + t167 * t264;
t98 = qJD(4) * t157 + t264 * t166 + t167 * t260;
t39 = pkin(4) * t98 - pkin(11) * t97 + t152;
t402 = -qJD(5) * t291 - t17 * t259 + t263 * t39;
t401 = t124 ^ 2;
t397 = pkin(3) * t264;
t394 = -qJ(6) * t358 - qJD(6) * t216 + t411;
t393 = t358 * pkin(5) - qJD(5) * t357 + t259 * t429 - t263 * t433;
t153 = -t175 * t259 - t263 * t327;
t154 = -t175 * t263 + t259 * t327;
t298 = pkin(5) * t263 + qJ(6) * t259;
t392 = pkin(5) * t153 - qJ(6) * t154 + t297 * t164 - (qJD(5) * t298 - qJD(6) * t263) * t217 - t384;
t391 = t259 * t99 + t263 * t58;
t64 = t112 * t264 - t108;
t88 = pkin(3) * t190 + t99;
t390 = t259 * t88 + t263 * t64;
t388 = t124 * t30;
t250 = pkin(3) * t260 + pkin(11);
t387 = t250 * t85;
t386 = t263 * t49;
t382 = t300 + t420;
t381 = t420 - t59;
t379 = t122 * t259;
t378 = t124 * t122;
t377 = t124 * t263;
t376 = t289 * t259;
t373 = t189 * t236;
t372 = t190 * t236;
t370 = t217 * t259;
t369 = t217 * t263;
t368 = t236 * t261;
t367 = t236 * t265;
t254 = t257 ^ 2;
t366 = t254 * qJD(1) ^ 2;
t354 = t262 ^ 2 - t266 ^ 2;
t348 = qJD(2) * t265;
t337 = t20 * t342 + t417;
t332 = pkin(3) * t344;
t331 = t262 * t366;
t330 = t259 * t364;
t323 = t250 * t343;
t322 = t254 * t350;
t320 = t254 * t339;
t316 = t120 * t264 - t260 * t130;
t314 = t259 * t164 + t153;
t313 = t164 * t263 + t154;
t305 = t245 + t352;
t304 = 0.2e1 * t320;
t299 = -0.2e1 * pkin(1) * t320;
t72 = pkin(4) * t364 - t316;
t296 = -t289 * t54 - t387;
t294 = t20 * t263 - t21 * t259;
t293 = t20 * t259 + t21 * t263;
t292 = -t259 * t58 + t263 * t99;
t290 = -t259 * t73 + t263 * t94;
t286 = -t120 * t345 - t130 * t344 - t260 * t92 + t264 * t87;
t228 = -pkin(4) - t298;
t285 = -t25 * t423 - t317;
t143 = t157 * t259 + t263 * t364;
t282 = t263 * t17 + t259 * t39 + t94 * t342 - t343 * t73;
t277 = t217 * t342 - t314;
t276 = t217 * t343 + t313;
t275 = t383 * pkin(11);
t274 = t284 * t257;
t18 = -pkin(4) * t325 - t286;
t269 = qJD(5) * t294 + t417;
t251 = -pkin(4) - t397;
t214 = t228 - t397;
t144 = t157 * t263 - t330;
t133 = t217 * t297 - t287;
t111 = -pkin(5) * t216 - t162 * t263 + t177 * t259;
t110 = qJ(6) * t216 + t357;
t77 = pkin(5) * t124 + qJ(6) * t122;
t61 = -qJD(5) * t330 + t157 * t342 + t259 * t97 - t263 * t325;
t60 = qJD(5) * t143 - t259 * t325 - t263 * t97;
t40 = pkin(5) * t143 - qJ(6) * t144 + t72;
t35 = -pkin(5) * t156 - t290;
t34 = qJ(6) * t156 + t291;
t28 = -t122 * t423 - t48;
t27 = -t292 - t415;
t26 = t391 + t414;
t23 = t259 * t64 - t263 * t88 - t415;
t22 = t390 + t414;
t10 = pkin(5) * t61 + qJ(6) * t60 - qJD(6) * t144 + t18;
t8 = -pkin(5) * t98 - t402;
t7 = qJ(6) * t98 + qJD(6) * t156 + t282;
t1 = [0, 0, 0, t262 * t266 * t304, -t354 * t304, t305 * t324, -t305 * t325, 0, -t199 * t258 - t209 * t245 + t262 * t299, -t198 * t258 - t208 * t245 + t266 * t299, t161 * t212 + t167 * t190, -t161 * t211 - t190 * t166 + t167 * t189 - t212 * t329, -t167 * t236 + (-t161 * t266 + (qJD(1) * t212 + t190) * t349) * t257, t166 * t236 + (t329 * t266 + (-qJD(1) * t211 + t189) * t349) * t257 (-t236 * t257 - t322) * t349, -t271 * t236 - t209 * t189 + t200 * t329 + t199 * t211 + t172 * t166 + (-t270 * t266 + (qJD(1) * t311 + t137) * t349) * t257, t279 * t236 + t209 * t190 + t200 * t161 + t199 * t212 + t172 * t167 + (t280 * t266 + (-qJD(1) * t356 - t138) * t349) * t257, t157 * t268 + t288 * t97, -t156 * t268 - t157 * t85 - t288 * t98 + t289 * t97, t97 * t338 + (-t268 * t266 + t288 * t349 + (t157 * t349 - t266 * t97) * qJD(1)) * t257, -t98 * t338 + (t289 * t349 + t85 * t266 + (-t156 * t349 + t266 * t98) * qJD(1)) * t257 (-t274 - t322) * t349, -t286 * t284 - t306 * t364 - t152 * t289 + t160 * t85 + t136 * t156 + t151 * t98 + (qJD(1) * t316 + t58) * t325, t136 * t157 + t151 * t97 + t152 * t288 + t160 * t268 - t239 * t361 + t281 * t284 - t307 * t364 - t325 * t59, -t124 * t60 - t144 * t48, t122 * t60 - t124 * t61 + t143 * t48 - t144 * t49, t124 * t98 + t144 * t85 - t156 * t48 + t423 * t60, -t122 * t98 - t143 * t85 - t156 * t49 + t423 * t61, t156 * t85 - t423 * t98, t18 * t122 + t15 * t143 - t156 * t317 + t24 * t98 + t290 * t85 - t402 * t423 + t72 * t49 + t54 * t61, t18 * t124 + t15 * t144 - t5 * t156 - t25 * t98 + t282 * t423 - t291 * t85 - t72 * t48 - t54 * t60, t10 * t122 + t143 * t9 - t156 * t4 - t20 * t98 + t30 * t61 - t35 * t85 + t40 * t49 + t423 * t8, -t122 * t7 + t124 * t8 - t143 * t2 + t144 * t4 - t20 * t60 - t21 * t61 - t34 * t49 - t35 * t48, -t10 * t124 - t144 * t9 + t156 * t2 + t21 * t98 + t30 * t60 + t34 * t85 + t40 * t48 - t423 * t7, t10 * t30 + t2 * t34 + t20 * t8 + t21 * t7 + t35 * t4 + t40 * t9; 0, 0, 0, -t266 * t331, t354 * t366 (qJD(2) - t245) * t326, t245 * t327 - t239, 0, pkin(1) * t331 + t207 * t245 - t199, pkin(8) * t239 + t204 * t245 + (-t258 * t339 + t366) * t398, t161 * t261 - t190 * t367 (t161 - t373) * t265 + (-t329 + t372) * t261, -t236 * t346 + (t236 * t363 + (qJD(2) * t261 - t190) * t262) * t353, t236 * t347 + (-t266 * t368 + (-t189 + t348) * t262) * t353, t236 * t327, -pkin(2) * t329 - t199 * t265 + t310 * t236 + t207 * t189 + (pkin(9) * t367 + t172 * t261) * qJD(3) + (-t137 * t262 + (-pkin(9) * t349 - t172 * t266) * t261) * t353, -pkin(2) * t161 + t199 * t261 - t355 * t236 - t207 * t190 + (-pkin(9) * t368 + t172 * t265) * qJD(3) + (-t172 * t363 + (-pkin(9) * t348 + t138) * t262) * t353, t217 * t268 + t288 * t359, -t216 * t268 - t217 * t85 - t288 * t358 + t289 * t359 (-t359 * t266 + (qJD(2) * t217 - t288) * t262) * t353 + t359 * t338 (t358 * t266 + (-qJD(2) * t216 - t289) * t262) * t353 - t358 * t338, t274 * t351, t252 * t85 + t136 * t216 + (qJD(2) * t287 - t58) * t327 + t358 * t151 - t412 * t289 + t425 * t284, t136 * t217 + t359 * t151 - t177 * t239 + t252 * t268 + t284 * t438 + t412 * t288 + t327 * t59, -t124 * t276 - t369 * t48, t314 * t124 + t313 * t122 + (t46 - t386 + (-t377 + t379) * qJD(5)) * t217, t124 * t358 - t216 * t48 + t276 * t423 + t369 * t85, -t122 * t358 - t216 * t49 + t277 * t423 - t370 * t85, t216 * t85 - t358 * t423, -t54 * t153 - t287 * t49 - t317 * t216 + t358 * t24 + t384 * t122 + (t54 * qJD(5) * t217 + t162 * t85 - (-qJD(5) * t177 - t433) * t423) * t263 + (t15 * t217 - t54 * t164 - t177 * t85 - (-qJD(5) * t162 + t429) * t423) * t259, t384 * t124 + t15 * t369 - t5 * t216 - t358 * t25 - t276 * t54 + t287 * t48 - t357 * t85 - t411 * t423, -t111 * t85 - t122 * t392 + t133 * t49 - t20 * t358 - t216 * t4 + t277 * t30 + t370 * t9 - t393 * t423, -t110 * t49 - t111 * t48 + t153 * t21 - t154 * t20 - t294 * t164 - t393 * t124 + t394 * t122 + (-qJD(5) * t293 - t2 * t259 + t263 * t4) * t217, t110 * t85 + t124 * t392 + t133 * t48 + t2 * t216 + t21 * t358 + t276 * t30 - t369 * t9 + t394 * t423, t110 * t2 + t111 * t4 + t133 * t9 - t20 * t393 - t21 * t394 - t30 * t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190 * t189, -t189 ^ 2 + t190 ^ 2, t161 + t373, -t329 - t372, t239, -t138 * t236 - t172 * t190 + t270, -t137 * t236 - t172 * t189 - t280, -t427, t422, t419, t403, t239, -t63 * t284 + (t190 * t289 + t239 * t264 + t284 * t345) * pkin(3) + t404, -t64 * t284 + (-t190 * t288 - t239 * t260 + t284 * t344) * pkin(3) + t421, t431, t437, t430, t436, t413, t251 * t49 + t296 * t259 + t300 * t122 - ((-qJD(5) * t250 - t88) * t263 + (t64 - t332) * t259) * t423 + t407, -t251 * t48 + t296 * t263 + t300 * t124 - (-t263 * t332 + t323 + t390) * t423 + t406, t214 * t49 + (-t289 * t30 - t387) * t259 + t382 * t122 - (-t250 * t342 - t259 * t332 + t23) * t423 + t409, t122 * t22 - t124 * t23 + (-t122 * t332 - t289 * t20 + (qJD(5) * t124 - t49) * t250) * t263 + (t124 * t332 + t289 * t21 - t250 * t48 + (t122 * t250 - t21) * qJD(5)) * t259 + t337, t214 * t48 - (-t22 - t323) * t423 - t382 * t124 + (-t332 * t423 + t387 + t435) * t263 + t408, -t20 * t23 - t21 * t22 + t214 * t9 + t250 * t269 + t293 * t332 + t30 * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t427, t422, t419, t403, t239, -t284 * t59 + t404, -t284 * t58 + t421, t431, t437, t430, t436, t413, -pkin(4) * t49 - t59 * t122 + t292 * t423 - t376 * t54 - t275 + t407, pkin(4) * t48 - t391 * t423 - t59 * t124 - t54 * t428 + (-t343 * t423 - t81) * pkin(11) + t406, t122 * t381 + t228 * t49 - t27 * t423 - t30 * t376 - t275 + t409, -t21 * t343 + t122 * t26 - t124 * t27 - t294 * t289 + (-t46 - t386 + (t377 + t379) * qJD(5)) * pkin(11) + t337, t228 * t48 - (-pkin(11) * t343 - t26) * t423 - t381 * t124 + (pkin(11) * t85 + t435) * t263 + t408, pkin(11) * t269 - t20 * t27 - t21 * t26 + t228 * t9 + t30 * t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378, -t122 ^ 2 + t401, t28, -t124 * t423 - t49, t85, -t124 * t54 + t285, t122 * t54 - t24 * t423 - t5, -t122 * t77 + t285 - t388 + 0.2e1 * t399, pkin(5) * t48 - qJ(6) * t49 + (t21 - t25) * t124 + (t20 - t362) * t122, 0.2e1 * t389 - t122 * t30 + t124 * t77 - (0.2e1 * qJD(6) - t24) * t423 + t5, -pkin(5) * t4 + qJ(6) * t2 - t20 * t25 + t21 * t362 - t30 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378 - t85, t28, -t423 ^ 2 - t401, t21 * t423 + t388 + t4;];
tauc_reg  = t1;
