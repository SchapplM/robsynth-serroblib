% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:12
% EndTime: 2019-03-09 09:05:36
% DurationCPUTime: 16.37s
% Computational Cost: add. (16541->770), mult. (46837->969), div. (0->0), fcn. (37975->12), ass. (0->384)
t262 = sin(qJ(5));
t258 = sin(pkin(11));
t259 = sin(pkin(6));
t263 = sin(qJ(2));
t409 = qJD(1) * qJD(2);
t385 = t263 * t409;
t363 = t259 * t385;
t218 = t258 * t363;
t267 = cos(qJ(2));
t455 = cos(pkin(11));
t216 = -t267 * t258 - t263 * t455;
t307 = t216 * qJDD(1);
t379 = t267 * t455;
t360 = qJD(2) * t379;
t278 = (qJD(1) * t360 - t307) * t259;
t276 = t218 - t278;
t274 = qJDD(5) - t276;
t520 = t262 * t274;
t266 = cos(qJ(5));
t519 = t266 * t274;
t361 = t259 * t379;
t420 = qJD(1) * t263;
t394 = t259 * t420;
t192 = -qJD(1) * t361 + t258 * t394;
t260 = cos(pkin(6));
t421 = qJD(1) * t260;
t239 = qJD(2) + t421;
t442 = t239 * t192;
t518 = t276 - t442;
t517 = t276 + t442;
t435 = t260 * t263;
t245 = pkin(1) * t435;
t437 = t259 * t267;
t470 = pkin(8) + qJ(3);
t173 = (t437 * t470 + t245) * qJD(1);
t164 = t258 * t173;
t246 = t260 * t267 * pkin(1);
t236 = qJD(1) * t246;
t388 = t470 * t263;
t364 = t259 * t388;
t172 = -qJD(1) * t364 + t236;
t109 = t172 * t455 - t164;
t411 = qJD(4) - t109;
t220 = qJDD(1) * t361;
t308 = qJD(2) * t216;
t407 = qJDD(1) * t263;
t383 = t258 * t407;
t130 = -t220 + (-qJD(1) * t308 + t383) * t259;
t380 = qJD(2) * t470;
t317 = qJD(3) * t267 - t263 * t380;
t243 = pkin(8) * t437;
t386 = t260 * t409;
t366 = pkin(1) * t386;
t408 = qJDD(1) * t260;
t400 = pkin(1) * t408;
t395 = qJDD(1) * t243 + t263 * t400 + t267 * t366;
t406 = qJDD(1) * t267;
t102 = (qJ(3) * t406 + qJD(1) * t317) * t259 + t395;
t238 = qJDD(2) + t408;
t233 = t267 * t400;
t328 = -t263 * t366 + t233;
t418 = qJD(3) * t263;
t92 = pkin(2) * t238 + (-qJDD(1) * t388 + (-t267 * t380 - t418) * qJD(1)) * t259 + t328;
t50 = t455 * t102 + t258 * t92;
t500 = t238 * qJ(4) + t239 * qJD(4);
t369 = t50 + t500;
t33 = -pkin(4) * t130 + t369;
t155 = -t266 * t192 + t239 * t262;
t343 = -t262 * t130 - t266 * t238;
t67 = qJD(5) * t155 + t343;
t157 = t192 * t262 + t239 * t266;
t375 = -t266 * t130 + t262 * t238;
t68 = qJD(5) * t157 + t375;
t10 = pkin(5) * t68 + pkin(10) * t67 + t33;
t261 = sin(qJ(6));
t265 = cos(qJ(6));
t200 = t216 * t259;
t296 = -qJD(1) * t200 + qJD(5);
t158 = pkin(2) * t239 + t172;
t95 = t158 * t455 - t164;
t337 = qJD(4) - t95;
t309 = qJD(1) * t216;
t196 = t259 * t309;
t471 = t196 * pkin(4);
t486 = pkin(3) + pkin(9);
t57 = -t239 * t486 + t337 - t471;
t252 = pkin(2) * t267 + pkin(1);
t206 = -qJD(1) * t252 * t259 + qJD(3);
t300 = qJ(4) * t196 + t206;
t75 = t192 * t486 + t300;
t35 = t262 * t57 + t266 * t75;
t24 = pkin(10) * t296 + t35;
t481 = pkin(4) * t192;
t378 = t455 * t173;
t96 = t258 * t158 + t378;
t84 = -t239 * qJ(4) - t96;
t63 = -t84 - t481;
t42 = pkin(5) * t155 - pkin(10) * t157 + t63;
t348 = t24 * t261 - t265 * t42;
t468 = t258 * t102 - t455 * t92;
t396 = qJDD(4) + t468;
t32 = -pkin(4) * t276 - t238 * t486 + t396;
t352 = t252 * qJDD(1);
t405 = pkin(2) * t363 + qJDD(3);
t171 = -t259 * t352 + t405;
t454 = qJ(4) * t276;
t281 = qJD(4) * t196 + t171 + t454;
t39 = t130 * t486 + t281;
t416 = qJD(5) * t266;
t417 = qJD(5) * t262;
t333 = -t262 * t32 - t266 * t39 - t57 * t416 + t417 * t75;
t404 = qJDD(5) - t218;
t498 = t278 + t404;
t5 = pkin(10) * t498 - t333;
t1 = -t348 * qJD(6) + t261 * t10 + t265 * t5;
t154 = qJD(6) + t155;
t516 = t348 * t154 + t1;
t359 = pkin(5) * t266 + pkin(10) * t262;
t515 = qJD(5) * t359 + (-pkin(4) - t359) * t196 + t411;
t514 = qJD(6) * t296 - t67;
t410 = t196 - qJD(5);
t504 = t157 * t410;
t513 = t266 * t504;
t434 = t261 * t262;
t118 = t265 * t192 - t196 * t434;
t512 = t261 * t417 + t118;
t264 = sin(qJ(1));
t268 = cos(qJ(1));
t431 = t263 * t258;
t330 = t379 - t431;
t310 = t260 * t330;
t144 = t216 * t268 - t264 * t310;
t427 = t264 * t267;
t429 = t263 * t268;
t209 = -t260 * t427 - t429;
t511 = pkin(2) * t209 + t144 * pkin(3);
t510 = -t196 * t266 + t416;
t34 = -t262 * t75 + t266 * t57;
t351 = t34 * t410 - t333;
t503 = t260 * t216;
t145 = t264 * t503 + t330 * t268;
t140 = -t264 * t330 + t268 * t503;
t255 = t259 ^ 2;
t509 = 0.2e1 * t255;
t12 = t24 * t265 + t261 * t42;
t2 = -qJD(6) * t12 + t265 * t10 - t261 * t5;
t507 = -t12 * t154 - t2;
t99 = t157 * t261 - t265 * t296;
t506 = t410 * t99;
t371 = t410 * t155;
t505 = -t67 - t371;
t502 = -t471 + t411;
t285 = qJD(5) * t296;
t501 = -t262 * t285 + t519;
t286 = qJD(5) * t410;
t499 = t266 * t286 - t520;
t185 = t196 ^ 2;
t487 = t192 ^ 2;
t496 = -t487 - t185;
t495 = -t487 + t185;
t423 = t243 + t245;
t204 = t423 * qJD(2);
t323 = -g(1) * t145 + g(2) * t140 + g(3) * t200;
t381 = t262 * t39 - t266 * t32;
t8 = -t35 * qJD(5) - t381;
t494 = t35 * t410 - t8;
t349 = -t12 * t265 - t261 * t348;
t6 = -pkin(5) * t498 - t8;
t493 = qJD(5) * t349 + t6;
t101 = t265 * t157 + t261 * t296;
t415 = qJD(6) * t261;
t27 = t157 * t415 - t261 * t274 - t265 * t514;
t413 = qJD(6) * t266;
t306 = -t265 * t413 + t512;
t433 = t261 * t266;
t492 = -t101 * t306 - t27 * t433;
t432 = t262 * t265;
t119 = -t192 * t261 - t196 * t432;
t305 = t261 * t413 + t265 * t417 + t119;
t426 = t265 * t266;
t65 = qJDD(6) + t68;
t491 = -t154 * t305 + t65 * t426;
t141 = t264 * t216 + t268 * t310;
t436 = t259 * t268;
t335 = t141 * t262 + t266 * t436;
t490 = -t140 * t265 + t261 * t335;
t489 = t140 * t261 + t265 * t335;
t170 = pkin(2) * t260 + t246 - t364;
t186 = qJ(3) * t437 + t423;
t345 = -t170 * t455 + t258 * t186;
t76 = -t200 * pkin(4) - t260 * t486 + t345;
t199 = t259 * t431 - t361;
t242 = pkin(2) * t437;
t424 = -t199 * pkin(3) + t242;
t354 = -qJ(4) * t200 + t424;
t478 = pkin(9) * t199;
t321 = t354 - t478;
t485 = pkin(1) * t259;
t93 = -t321 - t485;
t469 = t262 * t76 + t266 * t93;
t194 = t259 * t308;
t419 = qJD(2) * t263;
t393 = t259 * t419;
t195 = -t258 * t393 + t259 * t360;
t235 = pkin(2) * t393;
t331 = -qJ(4) * t195 + qJD(4) * t200 + t235;
t64 = -t194 * t486 + t331;
t237 = qJD(2) * t246;
t159 = t259 * t317 + t237;
t389 = t470 * t259;
t160 = -t259 * t418 + (-t267 * t389 - t245) * qJD(2);
t89 = t159 * t258 - t455 * t160;
t66 = pkin(4) * t195 + t89;
t16 = -qJD(5) * t469 - t262 * t64 + t266 * t66;
t488 = t33 + t323;
t484 = pkin(2) * t258;
t483 = pkin(3) * t130;
t482 = pkin(3) * t238;
t480 = pkin(9) * t141;
t479 = pkin(9) * t144;
t476 = g(1) * t209;
t475 = g(1) * t264;
t472 = g(3) * t267;
t108 = t172 * t258 + t378;
t77 = t108 - t481;
t377 = pkin(2) * t394 + qJ(4) * t192;
t82 = -t196 * t486 + t377;
t44 = t262 * t77 + t266 * t82;
t467 = t101 * t99;
t464 = t154 * t99;
t463 = t261 * t65;
t462 = t261 * t99;
t414 = qJD(6) * t265;
t28 = t157 * t414 + t261 * t514 - t265 * t274;
t461 = t262 * t28;
t460 = t265 * t65;
t459 = t27 * t261;
t458 = t28 * t265;
t339 = pkin(5) * t262 - pkin(10) * t266 + qJ(4);
t213 = t339 + t484;
t251 = -pkin(2) * t455 - pkin(3);
t248 = -pkin(9) + t251;
t161 = t213 * t265 - t248 * t434;
t38 = -pkin(10) * t192 + t44;
t390 = t265 * t416;
t457 = qJD(6) * t161 + t248 * t390 + t261 * t515 - t265 * t38;
t162 = t213 * t261 + t248 * t432;
t391 = t261 * t416;
t456 = -qJD(6) * t162 - t248 * t391 + t261 * t38 + t265 * t515;
t453 = t101 * t154;
t452 = t108 * t239;
t449 = t157 * t155;
t448 = t157 * t192;
t447 = t192 * t155;
t446 = t192 * t196;
t445 = t196 * t239;
t440 = t255 * qJD(1) ^ 2;
t439 = t259 * t263;
t438 = t259 * t264;
t430 = t263 * t264;
t425 = t267 * t268;
t90 = t455 * t159 + t258 * t160;
t113 = t258 * t170 + t455 * t186;
t256 = t263 ^ 2;
t257 = t267 ^ 2;
t422 = t256 - t257;
t412 = qJD(2) - t239;
t403 = t101 * t510 - t27 * t262;
t399 = t267 * t440;
t397 = t260 * t425;
t81 = -t260 * qJD(4) - t90;
t106 = -t260 * qJ(4) - t113;
t387 = pkin(1) * t509;
t384 = t267 * t409;
t382 = g(2) * t436 - g(3) * t260;
t376 = t512 * t154;
t205 = pkin(2) * t435 - t389;
t374 = -t205 * t264 + t268 * t252;
t372 = t154 * t265;
t368 = t239 + t421;
t367 = t238 + t408;
t365 = t263 * t399;
t362 = t263 * t384;
t120 = t144 * t266 + t262 * t438;
t334 = -t141 * t266 + t262 * t436;
t358 = g(1) * t334 + g(2) * t120;
t357 = -g(1) * t141 + g(2) * t144;
t356 = -g(1) * t140 - g(2) * t145;
t355 = g(1) * t268 + g(2) * t264;
t350 = t12 * t261 - t265 * t348;
t41 = -pkin(10) * t200 + t469;
t167 = -t199 * t266 + t260 * t262;
t168 = t199 * t262 + t260 * t266;
t80 = -pkin(4) * t199 - t106;
t52 = pkin(5) * t167 - pkin(10) * t168 + t80;
t18 = t261 * t52 + t265 * t41;
t17 = -t261 * t41 + t265 * t52;
t45 = -t262 * t93 + t266 * t76;
t43 = -t262 * t82 + t266 * t77;
t344 = t130 * t199 - t192 * t194;
t342 = -t195 * t196 + t200 * t276;
t116 = t168 * t265 - t200 * t261;
t115 = t168 * t261 + t200 * t265;
t341 = -t205 * t268 - t252 * t264;
t225 = pkin(2) * t397;
t340 = -pkin(2) * t430 + t141 * pkin(3) + t225;
t58 = pkin(4) * t194 - t81;
t336 = -t154 * t414 - t463;
t15 = t262 * t66 + t266 * t64 + t76 * t416 - t417 * t93;
t23 = -pkin(5) * t296 - t34;
t332 = -pkin(10) * t65 + t154 * t23;
t327 = -g(1) * t438 + t382;
t326 = -t384 - t407;
t325 = g(1) * t120 - g(2) * t334 + g(3) * t167;
t121 = -t144 * t262 + t266 * t438;
t324 = -g(1) * t121 + g(2) * t335 - g(3) * t168;
t322 = g(1) * t144 + g(2) * t141 - g(3) * t199;
t320 = t239 * t89 - t356;
t319 = t145 * pkin(3) - qJ(4) * t144 + t374;
t318 = t325 - t6;
t316 = t130 * t260 - t194 * t239 + t199 * t238;
t315 = t195 * t239 - t200 * t238 - t260 * t276;
t314 = -qJ(4) * t140 + t340;
t313 = t326 * pkin(8);
t304 = -t323 - t50;
t303 = -t322 - t468;
t302 = pkin(3) * t140 + qJ(4) * t141 + t341;
t301 = -t196 * t89 - t259 * t355;
t299 = t336 - t506;
t298 = t130 * t200 - t192 * t195 - t194 * t196 + t199 * t276;
t297 = pkin(4) * t438 + pkin(9) * t145 + t319;
t294 = pkin(10) * qJD(6) * t154 - t318;
t293 = -t28 * t426 + t305 * t99;
t292 = qJ(4) * t145 + t511;
t291 = -qJD(6) * t350 + t1 * t265 - t2 * t261;
t290 = pkin(4) * t436 + pkin(9) * t140 + t302;
t105 = pkin(3) * t192 + t300;
t289 = -t105 * t196 + qJDD(4) - t303;
t288 = t196 * t410;
t287 = t196 * t296;
t280 = (-t352 - t475) * t259 + t382 + t405;
t279 = -t23 * t410 + t291;
t277 = t220 + (qJD(2) * t309 - t383) * t259;
t271 = t262 * t287 + t501;
t270 = -t520 + (-t285 + t287) * t266;
t249 = qJ(4) + t484;
t219 = t238 * t260;
t217 = -t242 - t485;
t211 = -pkin(8) * t439 + t246;
t210 = -t260 * t430 + t425;
t208 = -t260 * t429 - t427;
t207 = -t397 + t430;
t203 = -pkin(8) * t393 + t237;
t202 = t423 * qJD(1);
t201 = -pkin(8) * t394 + t236;
t150 = t259 * t313 + t328;
t149 = -pkin(8) * t363 + t395;
t147 = t196 * t426 + t239 * t261;
t146 = -t196 * t433 + t239 * t265;
t117 = -t354 - t485;
t114 = -pkin(3) * t196 + t377;
t111 = qJD(5) * t168 + t194 * t266;
t110 = t194 * t262 - t199 * t416 + t260 * t417;
t107 = -t260 * pkin(3) + t345;
t91 = pkin(5) * t157 + pkin(10) * t155;
t88 = -pkin(3) * t194 + t331;
t83 = -t239 * pkin(3) + t337;
t70 = t121 * t265 + t145 * t261;
t69 = -t121 * t261 + t145 * t265;
t62 = t266 * t68;
t54 = -qJD(6) * t115 - t110 * t265 + t195 * t261;
t53 = qJD(6) * t116 - t110 * t261 - t195 * t265;
t51 = t281 + t483;
t48 = t396 - t482;
t40 = pkin(5) * t200 - t45;
t37 = pkin(5) * t192 - t43;
t29 = pkin(5) * t111 + pkin(10) * t110 + t58;
t22 = t261 * t91 + t265 * t34;
t21 = -t261 * t34 + t265 * t91;
t14 = -pkin(5) * t195 - t16;
t13 = pkin(10) * t195 + t15;
t4 = -qJD(6) * t18 - t13 * t261 + t265 * t29;
t3 = qJD(6) * t17 + t13 * t265 + t261 * t29;
t7 = [0, 0, 0, 0, 0, qJDD(1), -g(2) * t268 + t475, t355, 0, 0 (qJDD(1) * t256 + 0.2e1 * t362) * t255 (t263 * t406 - t409 * t422) * t509 (qJD(2) * t267 * t368 + t263 * t367) * t259 (qJDD(1) * t257 - 0.2e1 * t362) * t255 (t267 * t367 - t368 * t419) * t259, t219, -g(1) * t208 - g(2) * t210 + t150 * t260 - t204 * t239 + t211 * t238 + (-t385 + t406) * t387, -g(1) * t207 - g(2) * t209 - t149 * t260 - t203 * t239 - t238 * t423 + t326 * t387 ((-qJD(2) * t201 + qJDD(1) * t423 + t149 + (-qJD(2) * t211 + t203) * qJD(1)) * t267 + (-qJD(2) * t202 - qJDD(1) * t211 - t150) * t263 - t355) * t259, t149 * t423 + t202 * t203 + t150 * t211 - t201 * t204 + t255 * qJDD(1) * pkin(1) ^ 2 - g(1) * (-pkin(1) * t264 + pkin(8) * t436) - g(2) * (pkin(1) * t268 + pkin(8) * t438) t342, t298, t315, t344, -t316, t219, t130 * t217 + t171 * t199 + t192 * t235 - t194 * t206 - t238 * t345 - t260 * t468 - t320, -t113 * t238 - t171 * t200 + t195 * t206 - t196 * t235 - t217 * t276 - t239 * t90 - t260 * t50 - t357, -t113 * t130 - t192 * t90 + t194 * t96 - t195 * t95 - t199 * t50 - t200 * t468 - t276 * t345 + t301, -g(1) * t341 - g(2) * t374 + t50 * t113 + t171 * t217 + t206 * t235 + t345 * t468 - t95 * t89 + t96 * t90, t219, -t315, t316, t342, t298, t344, t106 * t130 - t107 * t276 + t192 * t81 - t194 * t84 + t195 * t83 - t199 * t369 - t200 * t48 + t301, t105 * t194 + t107 * t238 - t117 * t130 - t192 * t88 - t199 * t51 + t260 * t48 + t320, -t105 * t195 - t106 * t238 + t117 * t276 + t196 * t88 + t200 * t51 - t239 * t81 + t260 * t369 + t357, -g(1) * t302 - g(2) * t319 + t105 * t88 - t106 * t369 + t48 * t107 + t51 * t117 + t84 * t81 + t83 * t89, -t110 * t157 - t168 * t67, t110 * t155 - t111 * t157 + t167 * t67 - t168 * t68, -t110 * t296 + t157 * t195 + t168 * t274 + t67 * t200, t111 * t155 + t167 * t68, -t111 * t296 - t155 * t195 - t167 * t274 + t68 * t200, -t404 * t200 + qJD(5) * t195 + (t200 * t307 + (-t195 * t216 - t200 * t360) * qJD(1)) * t259, -g(1) * t335 - g(2) * t121 + t63 * t111 + t58 * t155 + t16 * t296 + t33 * t167 + t34 * t195 - t8 * t200 + t274 * t45 + t80 * t68, -t63 * t110 + t15 * t410 + t58 * t157 + t33 * t168 - t35 * t195 - t200 * t333 - t274 * t469 - t80 * t67 + t358, t110 * t34 - t111 * t35 - t15 * t155 - t157 * t16 + t167 * t333 - t168 * t8 + t45 * t67 - t469 * t68 + t356, -g(1) * t290 - g(2) * t297 + t35 * t15 + t34 * t16 + t33 * t80 - t333 * t469 + t8 * t45 + t63 * t58, t101 * t54 - t116 * t27, -t101 * t53 + t115 * t27 - t116 * t28 - t54 * t99, t101 * t111 + t116 * t65 + t154 * t54 - t167 * t27, t115 * t28 + t53 * t99, -t111 * t99 - t115 * t65 - t154 * t53 - t167 * t28, t111 * t154 + t167 * t65, -g(1) * t489 - g(2) * t70 - t111 * t348 + t6 * t115 + t14 * t99 + t4 * t154 + t2 * t167 + t17 * t65 + t23 * t53 + t40 * t28, g(1) * t490 - g(2) * t69 - t1 * t167 + t14 * t101 - t12 * t111 + t6 * t116 - t3 * t154 - t18 * t65 + t23 * t54 - t40 * t27, -t1 * t115 - t101 * t4 - t116 * t2 - t12 * t53 + t17 * t27 - t18 * t28 - t3 * t99 + t348 * t54 - t358, t1 * t18 + t12 * t3 + t2 * t17 - t348 * t4 + t6 * t40 + t23 * t14 - g(1) * (pkin(5) * t335 + pkin(10) * t334 + t290) - g(2) * (pkin(5) * t121 + pkin(10) * t120 + t297); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t365, t422 * t440 (qJD(1) * t267 * t412 + t407) * t259, t365 (-t412 * t420 + t406) * t259, t238, -t476 + g(2) * t207 + t202 * t239 + t233 + (-t386 + t440) * t263 * pkin(1) + (t313 - t472) * t259, pkin(1) * t399 + g(1) * t210 - g(2) * t208 + t201 * t239 + (pkin(8) * t409 + g(3)) * t439 - t395, 0, 0, -t446, t495, -t518, t446, t277 - t445, t238, t452 + t206 * t196 + (-t192 * t394 + t238 * t455) * pkin(2) + t303, t109 * t239 + t192 * t206 + (t196 * t394 - t238 * t258) * pkin(2) + t304 -(-t108 + t96) * t196 + (t109 - t95) * t192 + (-t130 * t258 + t276 * t455) * pkin(2), -g(2) * t225 + t95 * t108 - t96 * t109 + (t50 * t258 - t468 * t455 - t476 + g(2) * t430 + (-t206 * t420 - t472) * t259) * pkin(2), t238, t518, t130 + t445, -t446, t495, t446, -t130 * t249 - t276 * t251 - (-t108 - t84) * t196 + (t83 - t411) * t192, -t452 + t114 * t192 + (-pkin(3) + t251) * t238 + t289, -t105 * t192 - t114 * t196 + t238 * t249 + t239 * t411 - t304 + t500, -g(1) * t292 - g(2) * t314 - g(3) * t354 - t105 * t114 - t83 * t108 + t249 * t369 + t48 * t251 - t411 * t84, t262 * t504 - t266 * t67, -t62 + t513 + (t67 - t371) * t262, t271 + t448, t68 * t262 - t266 * t371, t270 - t447, t296 * t192, t502 * t155 + t34 * t192 + t501 * t248 + t249 * t68 + t488 * t262 - t296 * t43 + t510 * t63, -t35 * t192 - t249 * t67 - t410 * t44 + (t196 * t262 - t417) * t63 + t499 * t248 + t502 * t157 + t488 * t266, t155 * t44 + t157 * t43 + (t196 * t35 + t248 * t67 - t8 + (-t155 * t248 - t35) * qJD(5)) * t266 + (-t196 * t34 - t248 * t68 + t333 + (t157 * t248 + t34) * qJD(5)) * t262 - t322, t33 * t249 - t35 * t44 - t34 * t43 - g(1) * (t292 + t479) - g(2) * (t314 + t480) - g(3) * t321 + t502 * t63 + (-t333 * t262 + t8 * t266 + (-t262 * t34 + t266 * t35) * qJD(5)) * t248, -t101 * t305 - t27 * t426, t293 - t492, t403 + t491, t28 * t433 - t306 * t99, -t461 + (t336 + t506) * t266 + t376, -t154 * t266 * t410 + t65 * t262, -t23 * t118 + t161 * t65 - t37 * t99 + t456 * t154 - t322 * t261 + (t2 + (-t23 * t261 + t248 * t99) * qJD(5) + t323 * t265) * t262 + (t23 * t414 - t248 * t28 + t6 * t261 + t348 * t410) * t266, -t37 * t101 - t23 * t119 - t162 * t65 - t457 * t154 - t322 * t265 + (-t1 + (t101 * t248 - t23 * t265) * qJD(5) - t323 * t261) * t262 + (t12 * t410 - t23 * t415 + t248 * t27 + t6 * t265) * t266, -t348 * t119 + t118 * t12 + t161 * t27 - t162 * t28 - t457 * t99 - t456 * t101 + t350 * t417 + (qJD(6) * t349 - t1 * t261 - t2 * t265 - t323) * t266, t1 * t162 + t2 * t161 - t6 * t266 * t248 - g(1) * (t479 + t511) - g(2) * (t340 + t480) - g(3) * (t424 - t478) + (t248 * t417 - t37) * t23 + t457 * t12 - t456 * t348 + t323 * t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 - t445, -t517, t496, t192 * t96 - t196 * t95 + t280, 0, 0, 0, 0, 0, 0, t496, t277 + t445, t517, t483 + t454 - t192 * t84 - (-qJD(4) - t83) * t196 + t280, 0, 0, 0, 0, 0, 0, t270 + t447, -t519 + t448 + (-t286 + t288) * t262, t262 * t505 - t513 - t62, t192 * t63 + t262 * t494 + t351 * t266 + t327, 0, 0, 0, 0, 0, 0, t266 * t299 + t376 + t461, t403 - t491, t293 + t492, -t118 * t348 - t119 * t12 + t262 * t493 + t279 * t266 + t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518, t238 + t446, -t239 ^ 2 - t185, t239 * t84 + t289 - t482, 0, 0, 0, 0, 0, 0, -t239 * t155 + t271, -t239 * t157 - t266 * t288 + t499, -t505 * t266 + (-t68 - t504) * t262, -t239 * t63 + t351 * t262 - t266 * t494 + t322, 0, 0, 0, 0, 0, 0, -t266 * t28 + (-t146 - t391) * t154 + t299 * t262, t266 * t27 + (t147 - t390) * t154 + (-t101 * t410 + t154 * t415 - t460) * t262, t101 * t146 + t147 * t99 + (t101 * t261 - t265 * t99) * t416 + (-t459 - t458 + (t101 * t265 + t462) * qJD(6)) * t262, -t12 * t147 + t146 * t348 + t279 * t262 - t266 * t493 + t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449, -t155 ^ 2 + t157 ^ 2, t155 * t296 + t192 * t416 - t239 * t417 - t343, -t449, t157 * t296 - t192 * t417 - t239 * t416 - t375, t274, -t63 * t157 + t296 * t35 - t416 * t75 - t417 * t57 + t325 - t381, t63 * t155 - t324 - t351, 0, 0, t101 * t372 - t459 (-t27 - t464) * t265 + (-t28 - t453) * t261, -t101 * t157 + t154 * t372 + t463, t154 * t462 - t458, -t154 ^ 2 * t261 + t157 * t99 + t460, -t154 * t157, -pkin(5) * t28 - t154 * t21 + t157 * t348 + t261 * t332 - t265 * t294 - t35 * t99, pkin(5) * t27 - t101 * t35 + t12 * t157 + t154 * t22 + t261 * t294 + t265 * t332, t101 * t21 + t22 * t99 + ((qJD(6) * t101 - t28) * pkin(10) + t516) * t265 + ((qJD(6) * t99 - t27) * pkin(10) + t507) * t261 + t324, t348 * t21 - t12 * t22 - t23 * t35 + t318 * pkin(5) + (t291 + t324) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t467, t101 ^ 2 - t99 ^ 2, -t27 + t464, -t467, -t28 + t453, t65, -g(1) * t69 - g(2) * t490 + g(3) * t115 - t23 * t101 - t507, g(1) * t70 - g(2) * t489 + g(3) * t116 + t23 * t99 - t516, 0, 0;];
tau_reg  = t7;