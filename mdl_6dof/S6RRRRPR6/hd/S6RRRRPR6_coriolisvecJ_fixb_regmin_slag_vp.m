% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:12
% EndTime: 2019-03-09 22:23:36
% DurationCPUTime: 9.13s
% Computational Cost: add. (12963->506), mult. (32277->708), div. (0->0), fcn. (23927->10), ass. (0->270)
t281 = sin(qJ(2));
t344 = qJD(1) * qJD(2);
t263 = t281 * t344;
t284 = cos(qJ(3));
t285 = cos(qJ(2));
t333 = t285 * t344;
t280 = sin(qJ(3));
t350 = qJD(3) * t280;
t337 = t281 * t350;
t343 = qJD(2) * qJD(3);
t194 = -qJD(1) * t337 + (t333 + t343) * t284;
t349 = qJD(3) * t284;
t335 = t281 * t349;
t351 = qJD(2) * t285;
t339 = t280 * t351;
t293 = t335 + t339;
t195 = qJD(1) * t293 + t280 * t343;
t355 = qJD(1) * t281;
t340 = t280 * t355;
t345 = t284 * qJD(2);
t228 = t340 - t345;
t353 = qJD(2) * t280;
t230 = t284 * t355 + t353;
t279 = sin(qJ(4));
t283 = cos(qJ(4));
t347 = qJD(4) * t283;
t348 = qJD(4) * t279;
t113 = t283 * t194 - t279 * t195 - t228 * t347 - t230 * t348;
t354 = qJD(1) * t285;
t269 = pkin(7) * t354;
t248 = qJD(2) * pkin(8) + t269;
t241 = -pkin(2) * t285 - pkin(8) * t281 - pkin(1);
t221 = t241 * qJD(1);
t383 = t221 * t280;
t186 = t248 * t284 + t383;
t314 = pkin(2) * t281 - pkin(8) * t285;
t239 = t314 * qJD(2);
t222 = qJD(1) * t239;
t320 = pkin(7) * t263;
t362 = -t284 * t222 - t280 * t320;
t292 = -qJD(3) * t186 - t362;
t107 = pkin(3) * t263 - pkin(9) * t194 + t292;
t299 = t221 * t349 + t280 * t222 - t248 * t350;
t288 = -t284 * t320 + t299;
t119 = -pkin(9) * t195 + t288;
t328 = t283 * t107 - t279 * t119;
t185 = t284 * t221 - t248 * t280;
t158 = -pkin(9) * t230 + t185;
t259 = -qJD(3) + t354;
t149 = -pkin(3) * t259 + t158;
t159 = -pkin(9) * t228 + t186;
t156 = t283 * t159;
t97 = t149 * t279 + t156;
t290 = -qJD(4) * t97 + t328;
t305 = t228 * t279 - t283 * t230;
t16 = pkin(4) * t263 - qJ(5) * t113 + qJD(5) * t305 + t290;
t179 = t283 * t228 + t230 * t279;
t289 = qJD(4) * t305 - t194 * t279 - t283 * t195;
t319 = -t279 * t107 - t283 * t119 - t149 * t347 + t159 * t348;
t18 = qJ(5) * t289 - qJD(5) * t179 - t319;
t276 = sin(pkin(11));
t277 = cos(pkin(11));
t4 = t277 * t16 - t18 * t276;
t57 = t113 * t277 + t276 * t289;
t2 = pkin(5) * t263 - pkin(10) * t57 + t4;
t419 = qJ(5) * t179;
t75 = t97 - t419;
t391 = t277 * t75;
t252 = -qJD(4) + t259;
t418 = qJ(5) * t305;
t154 = t279 * t159;
t96 = t283 * t149 - t154;
t74 = t96 + t418;
t69 = -pkin(4) * t252 + t74;
t30 = t276 * t69 + t391;
t408 = -t179 * t277 + t276 * t305;
t405 = pkin(10) * t408;
t21 = t30 + t405;
t278 = sin(qJ(6));
t346 = qJD(6) * t278;
t20 = t21 * t346;
t282 = cos(qJ(6));
t372 = t282 * t408;
t407 = t179 * t276 + t277 * t305;
t433 = t278 * t407 + t372;
t247 = -qJD(2) * pkin(2) + pkin(7) * t355;
t196 = pkin(3) * t228 + t247;
t145 = pkin(4) * t179 + qJD(5) + t196;
t79 = -pkin(5) * t408 + t145;
t443 = -t278 * t2 - t79 * t433 + t20;
t371 = t284 * t285;
t300 = pkin(3) * t281 - pkin(9) * t371;
t401 = pkin(8) + pkin(9);
t341 = qJD(3) * t401;
t236 = t314 * qJD(1);
t358 = pkin(7) * t340 + t284 * t236;
t442 = qJD(1) * t300 + t284 * t341 + t358;
t215 = t280 * t236;
t373 = t281 * t284;
t374 = t280 * t285;
t439 = t215 + (-pkin(7) * t373 - pkin(9) * t374) * qJD(1) + t280 * t341;
t5 = t276 * t16 + t277 * t18;
t56 = -t113 * t276 + t277 * t289;
t3 = pkin(10) * t56 + t5;
t441 = -t282 * t3 + t443;
t12 = qJD(6) * t372 + t278 * t56 + t282 * t57 + t346 * t407;
t244 = -qJD(6) + t252;
t392 = t244 * t433;
t440 = t12 + t392;
t423 = t278 * t408 - t282 * t407;
t437 = t423 * t433;
t231 = t279 * t280 - t283 * t284;
t296 = t231 * t285;
t402 = qJD(3) + qJD(4);
t435 = -qJD(1) * t296 + t402 * t231;
t232 = t279 * t284 + t280 * t283;
t363 = (-t354 + t402) * t232;
t432 = t423 ^ 2 - t433 ^ 2;
t13 = qJD(6) * t423 + t278 * t57 - t282 * t56;
t393 = t244 * t423;
t429 = -t13 - t393;
t436 = t442 * t283;
t249 = t401 * t280;
t250 = t401 * t284;
t434 = -t249 * t347 - t250 * t348 - t279 * t442 - t439 * t283;
t342 = t282 * t2 - t278 * t3;
t424 = -t79 * t423 + t342;
t430 = pkin(5) * t407;
t420 = pkin(10) * t407;
t359 = -t279 * t249 + t283 * t250;
t427 = -pkin(4) * t355 + qJ(5) * t435 - qJD(4) * t359 - qJD(5) * t232 + t279 * t439 - t436;
t426 = -qJ(5) * t363 - qJD(5) * t231 + t434;
t315 = -t269 + (-t280 * t354 + t350) * pkin(3);
t425 = pkin(4) * t363 + t315;
t422 = pkin(4) * t305;
t417 = t305 * t179;
t416 = t276 * t435 - t277 * t363;
t415 = -t276 * t363 - t277 * t435;
t414 = -t179 ^ 2 + t305 ^ 2;
t413 = -t179 * t252 + t113;
t411 = t179 * t196 + t319;
t410 = t196 * t305 + t290;
t409 = t252 * t305 + t289;
t406 = -0.2e1 * t344;
t396 = -t426 * t276 + t427 * t277;
t395 = t427 * t276 + t426 * t277;
t376 = t277 * t279;
t394 = pkin(3) * qJD(4);
t326 = -t158 * t279 - t156;
t80 = t326 + t419;
t367 = t283 * t158 - t154;
t81 = t367 + t418;
t388 = t276 * t81 - t277 * t80 + (-t276 * t283 - t376) * t394;
t377 = t276 * t279;
t386 = -t276 * t80 - t277 * t81 + (t277 * t283 - t377) * t394;
t207 = t232 * t281;
t227 = t284 * t241;
t399 = pkin(7) * t280;
t184 = -pkin(9) * t373 + t227 + (-pkin(3) - t399) * t285;
t261 = pkin(7) * t371;
t357 = t280 * t241 + t261;
t375 = t280 * t281;
t190 = -pkin(9) * t375 + t357;
t365 = t279 * t184 + t283 * t190;
t403 = t285 * t345 - t337;
t400 = pkin(4) * t276;
t143 = -qJD(2) * t296 - t207 * t402;
t208 = t231 * t281;
t352 = qJD(2) * t281;
t360 = t284 * t239 + t352 * t399;
t138 = t300 * qJD(2) + (-t261 + (pkin(9) * t281 - t241) * t280) * qJD(3) + t360;
t336 = t285 * t350;
t361 = t280 * t239 + t241 * t349;
t142 = -t293 * pkin(9) + (-t281 * t345 - t336) * pkin(7) + t361;
t327 = t283 * t138 - t142 * t279;
t33 = pkin(4) * t352 - qJ(5) * t143 - qJD(4) * t365 + qJD(5) * t208 + t327;
t144 = -t348 * t375 + (t373 * t402 + t339) * t283 + t403 * t279;
t295 = t279 * t138 + t283 * t142 + t184 * t347 - t190 * t348;
t37 = -qJ(5) * t144 - qJD(5) * t207 + t295;
t11 = t276 * t33 + t277 * t37;
t176 = -t231 * t277 - t232 * t276;
t177 = -t231 * t276 + t232 * t277;
t310 = t282 * t176 - t177 * t278;
t398 = qJD(6) * t310 + t278 * t416 + t282 * t415;
t134 = t176 * t278 + t177 * t282;
t397 = qJD(6) * t134 + t278 * t415 - t282 * t416;
t70 = t276 * t75;
t36 = t277 * t74 - t70;
t29 = t277 * t69 - t70;
t19 = -pkin(5) * t252 + t29 + t420;
t390 = t282 * t19;
t389 = t388 + t405;
t387 = t386 - t420;
t385 = -pkin(5) * t416 + t425;
t384 = t194 * t280;
t382 = t228 * t259;
t381 = t230 * t259;
t380 = t247 * t280;
t379 = t247 * t284;
t378 = t259 * t284;
t287 = qJD(1) ^ 2;
t370 = t285 * t287;
t286 = qJD(2) ^ 2;
t369 = t286 * t281;
t368 = t286 * t285;
t325 = t283 * t184 - t190 * t279;
t115 = -pkin(4) * t285 + qJ(5) * t208 + t325;
t121 = -qJ(5) * t207 + t365;
t59 = t276 * t115 + t277 * t121;
t323 = -t283 * t249 - t250 * t279;
t167 = -qJ(5) * t232 + t323;
t168 = -qJ(5) * t231 + t359;
t118 = t276 * t167 + t277 * t168;
t240 = pkin(3) * t375 + t281 * pkin(7);
t274 = t281 ^ 2;
t356 = -t285 ^ 2 + t274;
t197 = pkin(3) * t293 + pkin(7) * t351;
t267 = -pkin(3) * t284 - pkin(2);
t175 = pkin(3) * t195 + pkin(7) * t333;
t331 = qJD(6) * t19 + t3;
t10 = -t276 * t37 + t277 * t33;
t35 = -t276 * t74 - t391;
t329 = pkin(1) * t406;
t58 = t277 * t115 - t121 * t276;
t117 = t277 * t167 - t168 * t276;
t322 = t228 + t345;
t321 = -t230 + t353;
t266 = pkin(3) * t283 + pkin(4);
t210 = -pkin(3) * t377 + t277 * t266;
t318 = pkin(4) * t207 + t240;
t88 = pkin(10) * t176 + t118;
t317 = pkin(5) * t355 + pkin(10) * t415 + qJD(6) * t88 - t396;
t87 = -pkin(10) * t177 + t117;
t316 = pkin(10) * t416 + qJD(6) * t87 + t395;
t313 = pkin(3) * t230 - t422;
t7 = t278 * t19 + t282 * t21;
t161 = -t207 * t276 - t208 * t277;
t44 = -pkin(5) * t285 - pkin(10) * t161 + t58;
t160 = -t207 * t277 + t208 * t276;
t45 = pkin(10) * t160 + t59;
t312 = t278 * t44 + t282 * t45;
t311 = t282 * t160 - t161 * t278;
t109 = t160 * t278 + t161 * t282;
t201 = pkin(5) + t210;
t211 = pkin(3) * t376 + t266 * t276;
t307 = t201 * t282 - t211 * t278;
t306 = t201 * t278 + t211 * t282;
t304 = qJD(1) * t274 - t259 * t285;
t303 = pkin(4) * t231 + t267;
t302 = pkin(4) * t144 + t197;
t85 = -pkin(4) * t289 + t175;
t262 = pkin(4) * t277 + pkin(5);
t298 = t262 * t278 + t282 * t400;
t297 = t262 * t282 - t278 * t400;
t148 = -pkin(5) * t176 + t303;
t126 = -pkin(5) * t160 + t318;
t89 = -t422 - t430;
t86 = t313 - t430;
t84 = t143 * t277 - t144 * t276;
t83 = -t143 * t276 - t144 * t277;
t50 = -pkin(5) * t83 + t302;
t28 = -pkin(5) * t56 + t85;
t25 = qJD(6) * t109 + t278 * t84 - t282 * t83;
t24 = qJD(6) * t311 + t278 * t83 + t282 * t84;
t23 = t36 + t420;
t22 = t35 - t405;
t9 = pkin(10) * t83 + t11;
t8 = pkin(5) * t352 - pkin(10) * t84 + t10;
t6 = -t21 * t278 + t390;
t1 = [0, 0, 0, 0.2e1 * t285 * t263, t356 * t406, t368, -t369, 0, -pkin(7) * t368 + t281 * t329, pkin(7) * t369 + t285 * t329, t194 * t373 + t230 * t403 (-t228 * t284 - t230 * t280) * t351 + (-t384 - t195 * t284 + (t228 * t280 - t230 * t284) * qJD(3)) * t281, t259 * t337 - t194 * t285 + (t230 * t281 + t284 * t304) * qJD(2), t259 * t335 + t195 * t285 + (-t228 * t281 - t280 * t304) * qJD(2) (-t259 - t354) * t352 -(-t241 * t350 + t360) * t259 + (t247 * t349 + pkin(7) * t195 + (qJD(1) * t227 + t185) * qJD(2)) * t281 + ((pkin(7) * t228 + t380) * qJD(2) + (t383 + (pkin(7) * t259 + t248) * t284) * qJD(3) + t362) * t285 (-pkin(7) * t336 + t361) * t259 + t299 * t285 + (pkin(7) * t194 - t247 * t350) * t281 + ((pkin(7) * t230 + t379) * t285 + (-pkin(7) * t378 - qJD(1) * t357 - t186) * t281) * qJD(2), -t113 * t208 - t143 * t305, -t113 * t207 - t143 * t179 + t144 * t305 - t208 * t289, -t113 * t285 - t143 * t252 + (-qJD(1) * t208 - t305) * t352, -t289 * t285 + t144 * t252 + (-qJD(1) * t207 - t179) * t352 (-t252 - t354) * t352, -t327 * t252 - t328 * t285 + t197 * t179 - t240 * t289 + t175 * t207 + t196 * t144 + (t252 * t365 + t285 * t97) * qJD(4) + (qJD(1) * t325 + t96) * t352, t295 * t252 - t319 * t285 - t197 * t305 + t240 * t113 - t175 * t208 + t196 * t143 + (-qJD(1) * t365 - t97) * t352, t10 * t407 + t11 * t408 + t160 * t5 - t161 * t4 - t29 * t84 + t30 * t83 + t56 * t59 - t57 * t58, t29 * t10 + t30 * t11 + t145 * t302 + t318 * t85 + t4 * t58 + t5 * t59, t109 * t12 + t24 * t423, -t109 * t13 + t12 * t311 + t24 * t433 - t25 * t423, -t12 * t285 - t24 * t244 + (qJD(1) * t109 + t423) * t352, t13 * t285 + t244 * t25 + (qJD(1) * t311 + t433) * t352 (-t244 - t354) * t352 -(-t278 * t9 + t282 * t8) * t244 - t342 * t285 - t50 * t433 + t126 * t13 - t28 * t311 + t79 * t25 + (t244 * t312 + t285 * t7) * qJD(6) + ((-t278 * t45 + t282 * t44) * qJD(1) + t6) * t352, t28 * t109 + t126 * t12 - t20 * t285 + t79 * t24 + t50 * t423 + ((-qJD(6) * t45 + t8) * t244 + t2 * t285) * t278 + ((qJD(6) * t44 + t9) * t244 + t331 * t285) * t282 + (-qJD(1) * t312 - t7) * t352; 0, 0, 0, -t281 * t370, t356 * t287, 0, 0, 0, t287 * pkin(1) * t281, pkin(1) * t370, -t230 * t378 + t384 (t194 + t382) * t284 + (-t195 + t381) * t280, -t259 * t349 + (t259 * t371 + t281 * t321) * qJD(1), t259 * t350 + (-t259 * t374 + t281 * t322) * qJD(1), t259 * t355, -pkin(2) * t195 + t358 * t259 + (pkin(8) * t378 + t380) * qJD(3) + ((-pkin(8) * t353 - t185) * t281 + (-pkin(7) * t322 - t380) * t285) * qJD(1), -pkin(2) * t194 - t215 * t259 + (-pkin(8) * t259 * t280 + t379) * qJD(3) + (-t247 * t371 + (-pkin(8) * t345 + t186) * t281 + (t259 * t373 + t285 * t321) * pkin(7)) * qJD(1), t113 * t232 + t305 * t435, -t113 * t231 + t179 * t435 + t232 * t289 + t305 * t363, t435 * t252 + (qJD(2) * t232 + t305) * t355, t363 * t252 + (-qJD(2) * t231 + t179) * t355, t252 * t355, -t267 * t289 + t175 * t231 + (t250 * t347 + (-qJD(4) * t249 - t439) * t279 + t436) * t252 + t363 * t196 + t315 * t179 + (qJD(2) * t323 - t96) * t355, t267 * t113 + t175 * t232 + t434 * t252 - t435 * t196 - t315 * t305 + (-qJD(2) * t359 + t97) * t355, -t117 * t57 + t118 * t56 + t176 * t5 - t177 * t4 - t29 * t415 + t30 * t416 + t395 * t408 + t396 * t407, t4 * t117 + t5 * t118 + t145 * t425 + t396 * t29 + t395 * t30 + t85 * t303, t12 * t134 + t398 * t423, t12 * t310 - t13 * t134 - t397 * t423 + t398 * t433, -t398 * t244 + (qJD(2) * t134 - t423) * t355, t397 * t244 + (qJD(2) * t310 - t433) * t355, t244 * t355, t148 * t13 - t28 * t310 + t397 * t79 - t385 * t433 + (t278 * t316 + t282 * t317) * t244 + ((-t278 * t88 + t282 * t87) * qJD(2) - t6) * t355, t148 * t12 + t28 * t134 + t398 * t79 + t385 * t423 + (-t278 * t317 + t282 * t316) * t244 + (-(t278 * t87 + t282 * t88) * qJD(2) + t7) * t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230 * t228, -t228 ^ 2 + t230 ^ 2, t194 - t382, -t195 - t381, t263, -t186 * t259 - t230 * t247 + t292, -t185 * t259 + t228 * t247 - t288, -t417, t414, t413, t409, t263, t326 * t252 + (-t179 * t230 + t252 * t348 + t263 * t283) * pkin(3) + t410, -t367 * t252 + (t230 * t305 + t252 * t347 - t263 * t279) * pkin(3) + t411, -t210 * t57 + t211 * t56 + (t29 + t386) * t408 + (-t30 + t388) * t407, -t145 * t313 + t4 * t210 + t5 * t211 + t29 * t388 + t30 * t386, -t437, t432, t440, t429, t263, t307 * t263 + t86 * t433 + (t278 * t387 - t282 * t389) * t244 + (t244 * t306 - t7) * qJD(6) + t424, -t306 * t263 - t86 * t423 + (t278 * t389 + t282 * t387) * t244 + (t244 * t307 - t390) * qJD(6) + t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, t414, t413, t409, t263, -t252 * t97 + t410, -t252 * t96 + t411 (t276 * t56 - t277 * t57) * pkin(4) + (-t36 + t29) * t408 + (-t30 - t35) * t407, -t29 * t35 - t30 * t36 + (t145 * t305 + t276 * t5 + t277 * t4) * pkin(4), -t437, t432, t440, t429, t263, t297 * t263 + (t22 * t282 - t23 * t278) * t244 + t89 * t433 + (t244 * t298 - t7) * qJD(6) + t424, -t298 * t263 - (t22 * t278 + t23 * t282) * t244 - t89 * t423 + (t244 * t297 - t390) * qJD(6) + t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t407 ^ 2 - t408 ^ 2, -t29 * t407 - t30 * t408 + t85, 0, 0, 0, 0, 0, t13 - t393, t12 - t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t437, t432, t440, t429, t263 (-qJD(6) - t244) * t7 + t424, -t244 * t6 - t282 * t331 + t443;];
tauc_reg  = t1;