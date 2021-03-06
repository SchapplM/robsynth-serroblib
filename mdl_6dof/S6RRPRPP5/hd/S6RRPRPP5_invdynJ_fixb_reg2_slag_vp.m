% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:24
% EndTime: 2019-03-09 10:06:37
% DurationCPUTime: 7.24s
% Computational Cost: add. (5178->595), mult. (10616->660), div. (0->0), fcn. (6286->6), ass. (0->299)
t230 = sin(qJ(2));
t362 = qJD(1) * t230;
t451 = -t362 - qJD(4);
t232 = cos(qJ(4));
t229 = sin(qJ(4));
t358 = qJD(2) * t229;
t233 = cos(qJ(2));
t361 = qJD(1) * t233;
t139 = t232 * t361 + t358;
t390 = t451 * t139;
t345 = qJD(1) * qJD(2);
t322 = t230 * t345;
t343 = t233 * qJDD(1);
t446 = -t322 + t343;
t63 = qJD(4) * t139 - t232 * qJDD(2) + t229 * t446;
t241 = t63 + t390;
t450 = t63 - t390;
t323 = t229 * t361;
t356 = qJD(2) * t232;
t141 = -t323 + t356;
t284 = t141 * t451;
t64 = -qJD(4) * t323 + t229 * qJDD(2) + (qJD(2) * qJD(4) + t446) * t232;
t449 = t64 + t284;
t448 = t64 - t284;
t329 = t232 * t362;
t321 = t233 * t345;
t344 = t230 * qJDD(1);
t271 = t321 + t344;
t138 = qJDD(4) + t271;
t353 = qJD(4) * t232;
t432 = t138 * t229 - t353 * t451;
t336 = -t329 * t451 + t432;
t270 = t139 * t361 - t336;
t359 = qJD(2) * t141;
t279 = t336 + t359;
t419 = pkin(3) + pkin(7);
t420 = pkin(2) + pkin(8);
t375 = qJ(6) - t420;
t447 = qJ(5) * t138 - qJD(5) * t451;
t399 = t64 * t229;
t445 = t139 * t353 + t399;
t203 = pkin(7) * t362;
t346 = pkin(3) * t362 + qJD(3) + t203;
t357 = qJD(2) * t230;
t402 = t232 * t63;
t444 = t233 * ((t139 * t232 + t141 * t229) * qJD(4) + t399 + t402) - (t139 * t229 - t141 * t232) * t357;
t443 = t229 * t450 - t232 * t448;
t172 = t451 * qJ(5);
t212 = t230 * qJ(3);
t319 = -pkin(1) - t212;
t269 = -t233 * t420 + t319;
t89 = t269 * qJD(1);
t95 = -qJD(2) * t420 + t346;
t44 = t229 * t95 + t232 * t89;
t37 = -t172 + t44;
t405 = t451 * t37;
t354 = qJD(4) * t229;
t184 = pkin(2) * t322;
t397 = qJ(3) * t233;
t294 = pkin(8) * t230 - t397;
t348 = t230 * qJD(3);
t254 = qJD(2) * t294 - t348;
t42 = qJD(1) * t254 + qJDD(1) * t269 + t184;
t182 = pkin(7) * t321;
t200 = pkin(7) * t344;
t320 = qJDD(3) + t182 + t200;
t68 = pkin(3) * t271 - qJDD(2) * t420 + t320;
t318 = t229 * t42 - t232 * t68 + t353 * t89 + t354 * t95;
t125 = t138 * pkin(4);
t430 = t125 - qJDD(5);
t7 = t318 - t430;
t442 = -t7 - t405;
t403 = t451 * t44;
t441 = -t318 - t403;
t440 = qJ(6) * t64 + qJD(6) * t139;
t43 = -t229 * t89 + t232 * t95;
t373 = qJD(5) - t43;
t436 = 0.2e1 * t447;
t113 = t232 * t138;
t435 = -t354 * t451 - t113;
t421 = t141 ^ 2;
t434 = -t451 ^ 2 - t421;
t433 = -qJD(5) * t232 + t346;
t231 = sin(qJ(1));
t234 = cos(qJ(1));
t303 = g(1) * t234 + g(2) * t231;
t422 = t139 ^ 2;
t431 = t421 - t422;
t222 = g(3) * t230;
t377 = t233 * t234;
t382 = t231 * t233;
t429 = -g(1) * t377 - g(2) * t382 - t222;
t428 = t139 * t329 + t445;
t205 = pkin(7) * t361;
t147 = pkin(3) * t361 + t205;
t226 = qJD(2) * qJ(3);
t115 = t226 + t147;
t394 = t138 * t420;
t427 = -t115 * t451 - t394;
t379 = t232 * t234;
t117 = t229 * t231 - t230 * t379;
t383 = t231 * t232;
t119 = t229 * t234 + t230 * t383;
t380 = t232 * t233;
t258 = g(1) * t117 - g(2) * t119 + g(3) * t380 - t318;
t252 = t258 + t430;
t281 = qJ(5) * t141 - t115;
t418 = pkin(4) + pkin(5);
t33 = -t139 * t418 + qJD(6) + t281;
t407 = qJ(6) * t63;
t426 = (qJD(6) + t33) * t141 + t252 - t407;
t381 = t232 * qJ(5);
t425 = t229 * t418 - t381;
t396 = qJ(5) * t229;
t274 = t232 * t418 + t396;
t14 = t230 * (-t358 * t451 - t63) - t233 * (-t359 + t432);
t360 = qJD(2) * t139;
t240 = t230 * (-t356 * t451 - t64) + t233 * (-t360 + t435);
t423 = -0.2e1 * pkin(1);
t417 = pkin(5) * t138;
t237 = qJD(2) ^ 2;
t416 = pkin(7) * t237;
t415 = g(1) * t231;
t223 = g(3) * t233;
t218 = t233 * pkin(2);
t412 = t274 * t451 - t433;
t314 = qJD(4) * t375;
t347 = t232 * qJD(6);
t388 = t229 * t230;
t208 = pkin(2) * t362;
t101 = qJD(1) * t294 + t208;
t60 = -t101 * t229 + t232 * t147;
t411 = t229 * t314 - t347 - (-qJ(6) * t388 - t233 * t418) * qJD(1) + t60;
t349 = t229 * qJD(6);
t61 = t101 * t232 + t147 * t229;
t48 = qJ(5) * t361 + t61;
t410 = qJ(6) * t329 + t232 * t314 + t349 - t48;
t297 = pkin(4) * t232 + t396;
t409 = -t297 * t451 + t433;
t408 = qJ(5) * t64;
t201 = pkin(7) * t343;
t224 = qJDD(2) * qJ(3);
t225 = qJD(2) * qJD(3);
t90 = pkin(7) * t322 - t201 - t224 - t225;
t69 = pkin(3) * t446 - t90;
t11 = pkin(4) * t64 + qJ(5) * t63 - qJD(5) * t141 + t69;
t406 = t11 * t232;
t404 = t451 * t43;
t401 = t232 * t69;
t400 = t420 * t63;
t398 = pkin(7) * qJDD(2);
t395 = qJDD(2) * pkin(2);
t393 = t139 * qJ(5);
t392 = t141 * t139;
t391 = t141 * t420;
t387 = t229 * t233;
t386 = t230 * t231;
t385 = t230 * t232;
t384 = t230 * t234;
t378 = t233 * qJ(6);
t238 = qJD(1) ^ 2;
t376 = t233 * t238;
t30 = t141 * qJ(6) + t43;
t374 = qJD(5) - t30;
t367 = t218 + t212;
t333 = pkin(8) * t233 + t367;
t128 = -pkin(1) - t333;
t166 = t419 * t230;
t71 = t128 * t232 + t166 * t229;
t178 = qJ(3) * t382;
t341 = pkin(4) * t387;
t371 = t231 * t341 + t178;
t180 = qJ(3) * t377;
t370 = t234 * t341 + t180;
t369 = g(1) * t382 - g(2) * t377;
t167 = t419 * t233;
t219 = t234 * pkin(7);
t366 = pkin(3) * t234 + t219;
t365 = pkin(1) * t234 + pkin(7) * t231;
t227 = t230 ^ 2;
t228 = t233 ^ 2;
t364 = t227 - t228;
t363 = t227 + t228;
t355 = qJD(2) * t233;
t352 = qJD(4) * t233;
t351 = qJD(4) * t420;
t350 = qJD(5) * t229;
t339 = t33 * t353;
t338 = t231 * t380;
t337 = t232 * t377;
t66 = qJ(5) * t230 + t71;
t334 = g(1) * t337 + g(2) * t338 + g(3) * t385;
t332 = -g(1) * t384 - g(2) * t386 + t223;
t330 = t229 * t362;
t328 = t229 * t357;
t326 = t232 * t351;
t324 = t451 * t361;
t317 = -qJD(2) * pkin(2) + qJD(3);
t118 = t229 * t384 + t383;
t316 = -pkin(4) * t117 + qJ(5) * t118;
t120 = -t229 * t386 + t379;
t315 = pkin(4) * t119 - qJ(5) * t120;
t70 = -t128 * t229 + t232 * t166;
t312 = pkin(2) * t377 + qJ(3) * t384 + t365;
t311 = -t200 - t332;
t310 = t201 + t429;
t307 = t230 * t321;
t306 = t363 * qJDD(1) * pkin(7);
t305 = -g(1) * t119 - g(2) * t117;
t304 = -g(1) * t120 - g(2) * t118;
t302 = -g(2) * t234 + t415;
t299 = t330 * t451 - t435;
t298 = t420 * t445 - t332;
t296 = -pkin(4) * t229 + t381;
t295 = pkin(5) * t229 - t381;
t34 = pkin(4) * t451 + t373;
t293 = t229 * t34 + t232 * t37;
t292 = -t229 * t43 + t232 * t44;
t154 = t203 + t317;
t160 = -t205 - t226;
t287 = t154 * t233 + t160 * t230;
t286 = g(3) * (pkin(4) * t388 + t333);
t148 = t419 * t355;
t207 = pkin(2) * t357;
t85 = t207 + t254;
t27 = -t128 * t353 + t148 * t232 - t166 * t354 - t229 * t85;
t285 = t451 * t229;
t157 = qJDD(2) * t230 + t233 * t237;
t283 = t319 - t218;
t282 = pkin(4) * t120 + qJ(5) * t119 + t366;
t5 = -pkin(5) * t64 + qJDD(6) - t11;
t278 = t232 * t5 - t33 * t354;
t276 = pkin(3) * t231 + pkin(8) * t377 + t312;
t116 = t283 * qJD(1);
t273 = t116 * t362 + qJDD(3) - t311;
t9 = t229 * t68 + t232 * t42 + t353 * t95 - t354 * t89;
t74 = t138 * t230 - t355 * t451;
t272 = -qJ(3) * t355 - t348;
t26 = -t128 * t354 + t148 * t229 + t166 * t353 + t232 * t85;
t268 = t299 - t360;
t267 = -t229 * t352 - t230 * t356;
t31 = qJ(6) * t139 + t44;
t266 = t274 * t233;
t6 = t9 + t447;
t46 = pkin(4) * t139 - t281;
t265 = -t451 * t46 - t394;
t156 = -pkin(1) - t367;
t264 = t398 + (-qJD(1) * t156 - t116) * qJD(2);
t263 = t303 * t420;
t16 = qJ(5) * t355 + qJD(5) * t230 + t26;
t106 = t207 + t272;
t65 = qJD(1) * t272 + qJDD(1) * t283 + t184;
t262 = qJD(1) * t106 + qJDD(1) * t156 + t416 + t65;
t261 = -t233 * t303 - t222;
t257 = t269 * t415;
t22 = t141 * t285 - t402;
t1 = -qJD(6) * t141 + t407 - t417 + t7;
t21 = t418 * t451 + t374;
t28 = -t172 + t31;
t3 = t6 + t440;
t256 = -t1 * t232 + t3 * t229 + t332 + (t329 + t353) * t28 + (t330 + t354) * t21;
t255 = pkin(4) * t118 + qJ(5) * t117 + t276;
t253 = -t138 + t392;
t251 = g(1) * t118 - g(2) * t120 - g(3) * t387 - t9;
t250 = -t229 * t284 + t402 - t428;
t100 = t320 - t395;
t249 = qJD(2) * t287 + t100 * t230 - t90 * t233;
t248 = -t351 * t451 + t261;
t19 = t139 * t267 + t380 * t64;
t245 = t141 * t46 - t252;
t244 = t251 - t404;
t177 = t230 * t376;
t159 = t364 * t238;
t158 = qJDD(2) * t233 - t230 * t237;
t155 = qJ(3) - t296;
t153 = t375 * t232;
t152 = t375 * t229;
t146 = t419 * t357;
t144 = -qJ(3) * t361 + t208;
t130 = qJDD(1) * t228 - 0.2e1 * t307;
t129 = qJDD(1) * t227 + 0.2e1 * t307;
t124 = -qJ(3) - t425;
t84 = 0.2e1 * t230 * t343 - 0.2e1 * t345 * t364;
t83 = t233 * t297 + t167;
t73 = pkin(4) * t141 + t393;
t72 = -t266 - t167;
t67 = -pkin(4) * t230 - t70;
t50 = -pkin(4) * t361 - t60;
t49 = t232 * t378 + t66;
t47 = -t141 * t418 - t393;
t45 = t229 * t378 - t230 * t418 - t70;
t38 = (qJD(4) * t296 + t350) * t233 + (-t297 - t419) * t357;
t32 = (-t141 * t233 + t388 * t451) * qJD(1) - t435;
t29 = (qJD(4) * t425 - t350) * t233 + (t274 + t419) * t357;
t20 = -pkin(4) * t355 - t27;
t18 = t63 * t387 + (-t232 * t352 + t328) * t141;
t13 = qJ(6) * t267 + t233 * t347 + t16;
t12 = -qJ(6) * t328 + (qJ(6) * t353 - qJD(2) * t418 + t349) * t233 - t27;
t2 = [0, 0, 0, 0, 0, qJDD(1), t302, t303, 0, 0, t129, t84, t157, t130, t158, 0, 0.2e1 * pkin(1) * t446 - pkin(7) * t157 + t369 (t345 * t423 - t398) * t233 + (qJDD(1) * t423 - t302 + t416) * t230, -t303 + 0.2e1 * t306, -g(1) * (-t231 * pkin(1) + t219) - g(2) * t365 + (pkin(7) ^ 2 * t363 + pkin(1) ^ 2) * qJDD(1), 0, -t157, -t158, t129, t84, t130, t306 + t249 - t303, t230 * t264 + t233 * t262 - t369, t264 * t233 + (-t262 + t302) * t230, pkin(7) * t249 - g(1) * t219 - g(2) * t312 + t116 * t106 + t65 * t156 - t283 * t415, t18, t444, t14, t19, t240, t74, t138 * t70 - t139 * t146 + t167 * t64 - t451 * t27 + (-t115 * t356 - t318) * t230 + (qJD(2) * t43 - t115 * t354 + t401) * t233 + t304, -t138 * t71 - t141 * t146 - t167 * t63 + t451 * t26 + (t115 * t358 - t9) * t230 + (-qJD(2) * t44 - t115 * t353 - t229 * t69) * t233 - t305, -t139 * t26 - t141 * t27 + t63 * t70 - t64 * t71 + t292 * t357 + (-t318 * t229 - t232 * t9 + (t229 * t44 + t232 * t43) * qJD(4)) * t233 + t369, -g(1) * t366 - g(2) * t276 - t115 * t146 + t167 * t69 + t26 * t44 + t27 * t43 - t318 * t70 + t71 * t9 - t257, t18, t14, -t444, t74, -t240, t19, -t138 * t67 + t139 * t38 + t451 * t20 + t64 * t83 + (-t356 * t46 - t7) * t230 + (-qJD(2) * t34 - t354 * t46 + t406) * t233 + t304, -t139 * t16 + t141 * t20 - t63 * t67 - t64 * t66 + t293 * t357 + (-t229 * t7 - t232 * t6 + (t229 * t37 - t232 * t34) * qJD(4)) * t233 + t369, t138 * t66 - t141 * t38 - t16 * t451 + t63 * t83 + (-t358 * t46 + t6) * t230 + (qJD(2) * t37 + t11 * t229 + t353 * t46) * t233 + t305, -g(1) * t282 - g(2) * t255 + t11 * t83 + t37 * t16 + t34 * t20 + t46 * t38 + t6 * t66 + t7 * t67 - t257, t18, -t444, -t14, t19, t240, t74, t12 * t451 - t138 * t45 - t139 * t29 - t64 * t72 + (t33 * t356 - t1) * t230 + (-qJD(2) * t21 - t278) * t233 + t304, -t13 * t451 + t138 * t49 + t141 * t29 - t63 * t72 + (t33 * t358 + t3) * t230 + (qJD(2) * t28 - t229 * t5 - t339) * t233 + t305, -t12 * t141 + t13 * t139 + t45 * t63 + t49 * t64 + (-t21 * t229 - t232 * t28) * t357 + (t1 * t229 + t232 * t3 + (t21 * t232 - t229 * t28) * qJD(4)) * t233 - t369, t3 * t49 + t28 * t13 + t1 * t45 + t21 * t12 + t5 * t72 + t33 * t29 - g(1) * (pkin(5) * t120 + t282) - g(2) * (pkin(5) * t118 - qJ(6) * t377 + t255) - (t233 * t375 + t319) * t415; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t159, t344, t177, t343, qJDD(2), pkin(1) * t230 * t238 + t311, pkin(1) * t376 - t310, 0, 0, qJDD(2), -t344, -t343, -t177, t159, t177 (-pkin(2) * t230 + t397) * qJDD(1) + ((-t160 - t226) * t230 + (-t154 + t317) * t233) * qJD(1), -t144 * t361 + t273 - 0.2e1 * t395, 0.2e1 * t224 + 0.2e1 * t225 + (t116 * t233 + t144 * t230) * qJD(1) + t310, -t90 * qJ(3) - t160 * qJD(3) - t100 * pkin(2) - t116 * t144 - g(1) * (-pkin(2) * t384 + t180) - g(2) * (-pkin(2) * t386 + t178) - g(3) * t367 - t287 * qJD(1) * pkin(7), t22, t443, t32, t428, t270, t324, -t43 * t361 + qJ(3) * t64 + t451 * t60 + t346 * t139 + t427 * t232 + (t69 + t248) * t229, t44 * t361 - qJ(3) * t63 + t401 - (t61 + t326) * t451 + t346 * t141 - t427 * t229 - t334, t139 * t61 + t141 * t60 + (-t400 - t441) * t232 + (t43 * t362 - t9 + (t43 - t391) * qJD(4)) * t229 + t298, t69 * qJ(3) - t44 * t61 - t43 * t60 - g(1) * t180 - g(2) * t178 - g(3) * t333 + t263 * t230 + t346 * t115 - (qJD(4) * t292 + t9 * t229 - t232 * t318) * t420, t22, t32, -t443, t324, -t270, t428, t34 * t361 + t155 * t64 - t451 * t50 + t409 * t139 + t265 * t232 + (t11 + t248) * t229, t139 * t48 - t141 * t50 + (-t400 - t442) * t232 + (-t34 * t362 - t6 + (-t34 - t391) * qJD(4)) * t229 + t298, -t37 * t361 - t406 + t155 * t63 - (-t48 - t326) * t451 - t409 * t141 + t265 * t229 + t334, t11 * t155 - t37 * t48 - t34 * t50 - g(1) * (-qJ(5) * t337 + t370) - g(2) * (-qJ(5) * t338 + t371) - t286 + t409 * t46 + (g(3) * t381 + t263) * t230 - (qJD(4) * t293 + t6 * t229 - t7 * t232) * t420, t22, -t443, t141 * t361 - t299, t428, t270, t324, -t339 - t124 * t64 + t138 * t153 + t411 * t451 - t412 * t139 + (t21 * t233 - t33 * t385) * qJD(1) + (-t5 + t261) * t229, -t124 * t63 + t138 * t152 - t410 * t451 + t412 * t141 + (-t233 * t28 - t33 * t388) * qJD(1) + t278 + t334, t139 * t410 - t141 * t411 + t152 * t64 - t153 * t63 + t256, t3 * t152 - t1 * t153 + t5 * t124 - g(1) * t370 - g(2) * t371 - t286 + t412 * t33 + t410 * t28 + (g(3) * qJ(6) - t295 * t303) * t233 + t411 * t21 + (-g(3) * t295 - t303 * t375) * t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t344, qJDD(2) + t177, -t227 * t238 - t237, qJD(2) * t160 + t182 + t273 - t395, 0, 0, 0, 0, 0, 0, -t285 * t451 + t113 - t360, -t279, t250, -qJD(2) * t115 + t441 * t232 + (t9 + t404) * t229 + t332, 0, 0, 0, 0, 0, 0, t268, t250, t279, -qJD(2) * t46 + t442 * t232 + (-t34 * t451 + t6) * t229 + t332, 0, 0, 0, 0, 0, 0, t268, t279, t22 + t428, qJD(2) * t33 + t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, t431, -t241, -t392, -t449, t138, -t115 * t141 + t258 - t403, t115 * t139 + t244, 0, 0, t392, -t241, -t431, t138, t449, -t392, -t139 * t73 + t125 - t245 - t403, pkin(4) * t63 - t408 + (t37 - t44) * t141 + (t34 - t373) * t139, -t139 * t46 + t141 * t73 - t244 + t436, -t7 * pkin(4) - g(1) * t316 - g(2) * t315 + t6 * qJ(5) + t223 * t297 - t34 * t44 + t37 * t373 - t46 * t73, t392, -t431, t241, -t392, -t449, t138 (pkin(5) + t418) * t138 + t139 * t47 - t451 * t31 + t426, t139 * t33 - t141 * t47 + t30 * t451 - t251 + t436 + t440, t408 - t418 * t63 + (-t28 + t31) * t141 + (-t21 + t374) * t139, t3 * qJ(5) - t1 * t418 - t21 * t31 - t33 * t47 - g(1) * (-pkin(5) * t117 + t316) - g(2) * (pkin(5) * t119 + t315) + t374 * t28 + g(3) * t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t241, t434, t245 + t405, 0, 0, 0, 0, 0, 0, t253, t434, t241, t28 * t451 - t417 - t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, -t450, -t421 - t422, -t139 * t28 + t141 * t21 - t429 + t5;];
tau_reg  = t2;
