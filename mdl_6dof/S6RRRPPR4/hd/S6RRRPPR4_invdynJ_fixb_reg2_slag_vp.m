% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:08
% EndTime: 2019-03-09 15:36:35
% DurationCPUTime: 15.20s
% Computational Cost: add. (14394->771), mult. (32130->958), div. (0->0), fcn. (22842->12), ass. (0->367)
t314 = sin(qJ(6));
t318 = cos(qJ(6));
t315 = sin(qJ(3));
t319 = cos(qJ(3));
t434 = t319 * qJD(2);
t316 = sin(qJ(2));
t439 = qJD(3) * t316;
t554 = qJD(1) * t439 - qJDD(2);
t320 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t410 = t320 * t432;
t431 = t316 * qJDD(1);
t556 = -t410 - t431;
t141 = -qJD(3) * t434 + t315 * t554 + t319 * t556;
t444 = qJD(1) * t320;
t142 = (qJD(2) * (qJD(3) + t444) + t431) * t315 + t554 * t319;
t312 = sin(pkin(10));
t496 = cos(pkin(10));
t79 = -t141 * t312 + t496 * t142;
t80 = -t141 * t496 - t312 * t142;
t445 = qJD(1) * t316;
t420 = t315 * t445;
t234 = t420 - t434;
t418 = t319 * t445;
t443 = qJD(2) * t315;
t236 = t418 + t443;
t150 = t496 * t234 + t236 * t312;
t352 = -t312 * t234 + t236 * t496;
t95 = t150 * t314 + t318 * t352;
t21 = qJD(6) * t95 + t314 * t80 - t318 * t79;
t274 = -qJD(3) + t444;
t262 = qJD(6) + t274;
t506 = t95 * t262;
t570 = t21 - t506;
t553 = -t318 * t150 + t314 * t352;
t520 = t553 ^ 2;
t521 = t95 ^ 2;
t569 = t520 - t521;
t519 = t95 * t553;
t404 = t496 * t315;
t228 = t312 * t319 + t404;
t211 = t228 * qJD(3);
t454 = t228 * t444 - t211;
t403 = t496 * t319;
t440 = qJD(3) * t315;
t212 = qJD(3) * t403 - t312 * t440;
t419 = t315 * t444;
t453 = -t312 * t419 + t403 * t444 - t212;
t308 = g(3) * t320;
t321 = cos(qJ(1));
t317 = sin(qJ(1));
t523 = g(2) * t317;
t388 = g(1) * t321 + t523;
t337 = -t388 * t316 + t308;
t292 = pkin(7) * t431;
t386 = qJDD(2) * pkin(2) - pkin(7) * t410 - t292;
t333 = t386 - t337;
t568 = pkin(8) * qJD(3) * t274 + t333;
t309 = qJ(3) + pkin(10);
t296 = sin(t309);
t457 = t321 * t296;
t297 = cos(t309);
t463 = t317 * t297;
t195 = t320 * t457 - t463;
t458 = t320 * t321;
t196 = t296 * t317 + t297 * t458;
t120 = t195 * t318 - t196 * t314;
t307 = g(3) * t316;
t370 = t296 * t318 - t297 * t314;
t252 = -qJD(2) * pkin(2) + pkin(7) * t445;
t351 = -t234 * pkin(3) - qJD(4) - t252;
t335 = qJ(5) * t352 + t351;
t532 = pkin(4) + pkin(5);
t51 = -t150 * t532 + t335;
t461 = t317 * t320;
t193 = t296 * t461 + t297 * t321;
t194 = t297 * t461 - t457;
t540 = t193 * t318 - t194 * t314;
t567 = -g(1) * t120 - g(2) * t540 - t370 * t307 - t51 * t95;
t517 = qJ(4) + pkin(8);
t249 = t517 * t319;
t471 = t312 * t315;
t168 = t496 * t249 - t517 * t471;
t300 = t320 * qJDD(1);
t542 = -t316 * t432 + t300;
t226 = qJDD(3) - t542;
t566 = t168 * t226 - t296 * t337;
t406 = qJD(3) * t517;
t433 = t319 * qJD(4);
t205 = -t315 * t406 + t433;
t346 = -t315 * qJD(4) - t319 * t406;
t132 = t496 * t205 + t312 * t346;
t391 = pkin(2) * t316 - pkin(8) * t320;
t237 = t391 * qJD(1);
t171 = pkin(7) * t420 + t319 * t237;
t460 = t319 * t320;
t365 = pkin(3) * t316 - qJ(4) * t460;
t137 = qJD(1) * t365 + t171;
t217 = t315 * t237;
t465 = t316 * t319;
t468 = t315 * t320;
t154 = t217 + (-pkin(7) * t465 - qJ(4) * t468) * qJD(1);
t83 = t312 * t137 + t496 * t154;
t73 = qJ(5) * t445 + t83;
t502 = t132 - t73;
t436 = qJD(6) * t318;
t437 = qJD(6) * t314;
t20 = -t150 * t436 - t314 * t79 - t318 * t80 + t352 * t437;
t509 = t262 * t553;
t548 = t20 - t509;
t392 = pkin(2) * t320 + pkin(8) * t316;
t248 = -pkin(1) - t392;
t223 = t248 * qJD(1);
t294 = pkin(7) * t444;
t253 = qJD(2) * pkin(8) + t294;
t160 = t319 * t223 - t253 * t315;
t122 = -qJ(4) * t236 + t160;
t161 = t223 * t315 + t253 * t319;
t123 = -qJ(4) * t234 + t161;
t472 = t312 * t123;
t64 = t496 * t122 - t472;
t456 = qJD(5) - t64;
t503 = (t137 - t346) * t496 + (-t154 + t205) * t312;
t121 = t195 * t314 + t196 * t318;
t108 = -pkin(3) * t274 + t122;
t60 = t108 * t496 - t472;
t359 = qJD(5) - t60;
t551 = pkin(9) * t352;
t37 = t274 * t532 + t359 - t551;
t562 = pkin(9) * t150;
t405 = t496 * t123;
t61 = t312 * t108 + t405;
t57 = -t274 * qJ(5) + t61;
t43 = t57 + t562;
t238 = t391 * qJD(2);
t162 = qJD(1) * t238 + qJDD(1) * t248;
t156 = t319 * t162;
t206 = t542 * pkin(7) + qJDD(2) * pkin(8);
t71 = -qJD(3) * t161 - t315 * t206 + t156;
t42 = t226 * pkin(3) + t141 * qJ(4) - t236 * qJD(4) + t71;
t438 = qJD(3) * t319;
t70 = t315 * t162 + t319 * t206 + t223 * t438 - t253 * t440;
t46 = -qJ(4) * t142 - qJD(4) * t234 + t70;
t514 = t312 * t46 - t496 * t42;
t422 = -qJDD(5) - t514;
t7 = -pkin(9) * t80 - t226 * t532 - t422;
t16 = t312 * t42 + t496 * t46;
t221 = t226 * qJ(5);
t258 = t274 * qJD(5);
t11 = t221 - t258 + t16;
t9 = pkin(9) * t79 + t11;
t362 = -t314 * t7 - t318 * t9 - t37 * t436 + t43 * t437;
t369 = t296 * t314 + t297 * t318;
t372 = t193 * t314 + t194 * t318;
t563 = g(1) * t121 + g(2) * t372 + t369 * t307 + t51 * t553 + t362;
t428 = t316 * t532;
t561 = pkin(9) * t453 + qJD(1) * t428 + t503;
t560 = -pkin(9) * t454 + t502;
t487 = t150 * t274;
t53 = t80 - t487;
t559 = t150 * t352;
t558 = -t551 + t456;
t395 = -t294 + (-t419 + t440) * pkin(3);
t360 = t388 * t320;
t557 = t360 + t307;
t441 = qJD(2) * t320;
t417 = t315 * t441;
t555 = t316 * t438 + t417;
t552 = 0.2e1 * pkin(1);
t491 = t352 ^ 2;
t13 = t314 * t37 + t318 * t43;
t2 = -qJD(6) * t13 - t314 * t9 + t318 * t7;
t550 = t13 * t262 + t2;
t549 = t532 * t79;
t547 = t160 * t274 + t70;
t544 = -t194 * pkin(4) - t193 * qJ(5);
t279 = pkin(7) * t460;
t186 = t315 * t248 + t279;
t543 = qJ(5) * t453 - t228 * qJD(5) + t395;
t541 = -t80 * qJ(5) - t352 * qJD(5);
t227 = -t403 + t471;
t538 = -t150 * t445 + t226 * t227 + t274 * t454;
t66 = pkin(4) * t150 - t335;
t536 = -t66 * t352 - qJDD(5);
t533 = t150 * t453 - t227 * t80 - t228 * t79 + t352 * t454;
t531 = t79 * pkin(4);
t323 = qJD(1) ^ 2;
t530 = pkin(1) * t323;
t529 = pkin(4) * t226;
t528 = pkin(7) * t315;
t525 = g(1) * t317;
t518 = pkin(9) - t517;
t167 = t249 * t312 + t517 * t404;
t134 = -pkin(9) * t228 + t167;
t135 = pkin(9) * t227 + t168;
t68 = t134 * t318 - t135 * t314;
t516 = qJD(6) * t68 + t314 * t561 + t318 * t560;
t69 = t134 * t314 + t135 * t318;
t515 = -qJD(6) * t69 - t314 * t560 + t318 * t561;
t442 = qJD(2) * t316;
t451 = t319 * t238 + t442 * t528;
t84 = -t316 * t433 + t365 * qJD(2) + (-t279 + (qJ(4) * t316 - t248) * t315) * qJD(3) + t451;
t452 = t315 * t238 + t248 * t438;
t99 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t465 + (-qJD(4) * t316 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t320) * t315 + t452;
t41 = t312 * t84 + t496 * t99;
t513 = t454 * t532 - t543;
t12 = -t314 * t43 + t318 * t37;
t512 = t12 * t262;
t63 = t122 * t312 + t405;
t508 = t63 * t352;
t504 = pkin(4) * t445 + t503;
t501 = t132 - t83;
t287 = -pkin(3) * t496 - pkin(4);
t272 = -pkin(5) + t287;
t285 = pkin(3) * t312 + qJ(5);
t187 = t272 * t318 - t285 * t314;
t49 = t63 + t562;
t500 = qJD(6) * t187 - t314 * t49 + t318 * t558;
t188 = t272 * t314 + t285 * t318;
t499 = -qJD(6) * t188 - t314 * t558 - t318 * t49;
t498 = t227 * t436 - t228 * t437 - t314 * t454 - t318 * t453;
t148 = t227 * t314 + t228 * t318;
t497 = qJD(6) * t148 - t314 * t453 + t318 * t454;
t495 = pkin(7) * qJDD(1);
t494 = qJ(5) * t320;
t493 = t141 * t315;
t492 = t142 * t319;
t490 = t352 * t274;
t488 = t150 ^ 2;
t484 = t161 * t274;
t480 = t234 * t274;
t479 = t234 * t315;
t478 = t236 * t234;
t477 = t236 * t274;
t476 = t236 * t319;
t475 = t296 * t316;
t474 = t297 * t316;
t473 = t297 * t320;
t470 = t315 * t316;
t469 = t315 * t317;
t467 = t315 * t321;
t466 = t316 * t317;
t464 = t316 * t321;
t462 = t317 * t319;
t459 = t319 * t321;
t291 = t319 * pkin(3) + pkin(2);
t256 = t320 * t291;
t455 = -pkin(4) * t454 + t543;
t230 = t319 * t248;
t157 = -qJ(4) * t465 + t230 + (-pkin(3) - t528) * t320;
t166 = -qJ(4) * t470 + t186;
t101 = t312 * t157 + t496 * t166;
t276 = pkin(3) * t470;
t450 = qJ(5) * t474 - t276;
t449 = g(1) * t466 - g(2) * t464;
t239 = t316 * pkin(7) + t276;
t448 = t321 * pkin(1) + t317 * pkin(7);
t310 = t316 ^ 2;
t311 = t320 ^ 2;
t447 = t310 - t311;
t446 = t310 + t311;
t427 = t517 * t464;
t426 = t315 * t458;
t425 = t316 * t323 * t320;
t305 = t321 * pkin(7);
t424 = pkin(3) * t467 - t466 * t517 + t305;
t423 = g(1) * t458 + g(2) * t461 + t307;
t175 = pkin(3) * t555 + pkin(7) * t441;
t416 = t320 * t434;
t414 = t274 * t445;
t412 = t518 * t321;
t408 = -pkin(1) - t256;
t402 = -qJ(5) * t296 - t291;
t400 = t262 ^ 2;
t398 = pkin(3) * t469 + t291 * t458 + t448;
t397 = g(3) * (pkin(4) * t473 + t296 * t494 + t256);
t396 = t316 * t410;
t394 = t496 * t441;
t390 = -g(1) * t193 + g(2) * t195;
t389 = g(1) * t194 - g(2) * t196;
t387 = -g(2) * t321 + t525;
t96 = t101 - t494;
t202 = -t312 * t470 + t316 * t403;
t384 = t202 * qJ(5) - t239;
t383 = t312 * t99 - t496 * t84;
t382 = -pkin(3) * t236 - qJ(5) * t150;
t381 = pkin(7) * t234 + t252 * t315;
t380 = pkin(7) * t236 + t252 * t319;
t100 = t157 * t496 - t312 * t166;
t98 = t320 * pkin(4) - t100;
t65 = t320 * pkin(5) - t202 * pkin(9) + t98;
t201 = t228 * t316;
t67 = pkin(9) * t201 + t96;
t30 = -t314 * t67 + t318 * t65;
t31 = t314 * t65 + t318 * t67;
t377 = -pkin(8) * t226 + qJD(3) * t252;
t129 = -t212 * t316 - t312 * t416 - t315 * t394;
t376 = -t129 * t150 + t201 * t79;
t375 = -t488 - t491;
t374 = -t488 + t491;
t373 = -t160 * t319 - t161 * t315;
t128 = t201 * t314 + t202 * t318;
t322 = qJD(2) ^ 2;
t368 = qJDD(2) * t320 - t316 * t322;
t34 = qJ(5) * t442 - qJD(5) * t320 + t41;
t367 = t228 * qJ(5) + t291;
t366 = -g(1) * t195 - g(2) * t193 - g(3) * t475;
t213 = t315 * t461 + t459;
t358 = t79 - t490;
t357 = -t79 - t490;
t355 = t315 * t226 - t274 * t438;
t354 = t319 * t226 + t274 * t440;
t353 = -t150 * t454 + t227 * t79;
t349 = -t132 * t150 + t167 * t80 - t168 * t79 - t423;
t348 = -t366 - t514;
t347 = t196 * pkin(4) + t195 * qJ(5) + t398;
t345 = g(1) * t297 * t464 + g(2) * t316 * t463 - g(3) * t473 - t167 * t226;
t344 = -t142 * pkin(3) - qJDD(4) + t386;
t130 = t211 * t316 + t312 * t417 - t319 * t394;
t342 = -t130 * qJ(5) + t202 * qJD(5) - t175;
t341 = t317 * t408 + t424;
t278 = pkin(3) * t462;
t339 = -pkin(3) * t426 - t195 * pkin(4) + qJ(5) * t196 + t278;
t336 = -t80 - t487;
t334 = -t63 * t274 + t348;
t332 = t129 * t352 + t130 * t150 - t201 * t80 - t202 * t79;
t330 = t129 * t274 + t150 * t442 + t201 * t226 - t320 * t79;
t329 = -pkin(3) * t213 - t193 * pkin(4) + t194 * qJ(5);
t328 = g(1) * t196 + g(2) * t194 + g(3) * t474 - t64 * t274 - t16;
t327 = t344 - t541;
t326 = -t344 + t337;
t325 = t326 + t541;
t222 = -qJDD(6) + t226;
t216 = t319 * t458 + t469;
t215 = -t426 + t462;
t214 = -t317 * t460 + t467;
t185 = -pkin(7) * t468 + t230;
t172 = -pkin(7) * t418 + t217;
t163 = -t226 * t320 - t274 * t442;
t147 = -t318 * t227 + t228 * t314;
t145 = pkin(4) * t227 - t367;
t127 = -t318 * t201 + t202 * t314;
t119 = -t186 * qJD(3) + t451;
t118 = (-t316 * t434 - t320 * t440) * pkin(7) + t452;
t110 = -t227 * t532 + t367;
t109 = pkin(4) * t201 - t384;
t97 = -t201 * t532 + t384;
t72 = pkin(4) * t352 - t382;
t56 = t274 * pkin(4) + t359;
t55 = t228 * t226 + t274 * t453 - t352 * t445;
t54 = -t352 * t532 + t382;
t52 = -pkin(4) * t129 - t342;
t48 = qJD(6) * t128 + t318 * t129 - t314 * t130;
t47 = t314 * t129 + t318 * t130 - t201 * t436 + t202 * t437;
t36 = -pkin(4) * t442 + t383;
t33 = t129 * t532 + t342;
t32 = -t130 * t352 + t202 * t80;
t29 = t228 * t80 - t352 * t453;
t28 = t130 * t274 + t202 * t226 - t320 * t80 + t352 * t442;
t27 = -pkin(9) * t129 + t34;
t26 = t130 * pkin(9) - qJD(2) * t428 + t383;
t19 = -t327 + t531;
t14 = -t422 - t529;
t10 = t327 - t549;
t4 = -qJD(6) * t31 + t318 * t26 - t314 * t27;
t3 = qJD(6) * t30 + t314 * t26 + t318 * t27;
t1 = [0, 0, 0, 0, 0, qJDD(1), t387, t388, 0, 0, qJDD(1) * t310 + 0.2e1 * t396, 0.2e1 * t300 * t316 - 0.2e1 * t432 * t447, qJDD(2) * t316 + t320 * t322, qJDD(1) * t311 - 0.2e1 * t396, t368, 0 (-0.2e1 * pkin(1) * t432 - pkin(7) * qJDD(2)) * t316 + (-pkin(7) * t322 + qJDD(1) * t552 + t387) * t320, -t368 * pkin(7) + t552 * t556 - t449, 0.2e1 * t446 * t495 - t388, -g(1) * (-pkin(1) * t317 + t305) - g(2) * t448 + (pkin(7) ^ 2 * t446 + pkin(1) ^ 2) * qJDD(1), -t141 * t465 + (-t315 * t439 + t416) * t236 (-t234 * t319 - t236 * t315) * t441 + (t493 - t492 + (-t476 + t479) * qJD(3)) * t316 (-t274 * t434 + t141) * t320 + (qJD(2) * t236 + t354) * t316, t142 * t470 + t234 * t555 (t274 * t443 + t142) * t320 + (-qJD(2) * t234 - t355) * t316, t163, -g(1) * t214 - g(2) * t216 - t119 * t274 + t185 * t226 + (qJD(2) * t381 - t71) * t320 + (pkin(7) * t142 + qJD(2) * t160 + t252 * t438 - t315 * t386) * t316, -g(1) * t213 - g(2) * t215 + t118 * t274 - t186 * t226 + (qJD(2) * t380 + t70) * t320 + (-pkin(7) * t141 - qJD(2) * t161 - t252 * t440 - t319 * t386) * t316, -t118 * t234 - t119 * t236 + t185 * t141 - t186 * t142 + t373 * t441 + (-t315 * t70 - t319 * t71 + (t160 * t315 - t161 * t319) * qJD(3)) * t316 + t449, t70 * t186 + t161 * t118 + t71 * t185 + t160 * t119 - g(1) * t305 - g(2) * (t321 * t392 + t448) - t248 * t525 + (t252 * t441 - t316 * t386) * pkin(7), t32, t332, t28, t376, -t330, t163, t100 * t226 + t129 * t351 + t150 * t175 - t201 * t344 + t239 * t79 + t274 * t383 + t320 * t514 + t442 * t60 + t389, -t101 * t226 + t130 * t351 + t16 * t320 + t175 * t352 - t202 * t344 + t239 * t80 + t274 * t41 - t442 * t61 + t390, -t100 * t80 - t101 * t79 + t129 * t61 + t130 * t60 - t150 * t41 - t16 * t201 + t202 * t514 + t352 * t383 + t449, t16 * t101 + t61 * t41 - t514 * t100 - t60 * t383 - t344 * t239 - t351 * t175 - g(1) * t341 - g(2) * (t398 + t427) t32, t28, -t332, t163, t330, t376, t109 * t79 - t129 * t66 + t14 * t320 + t150 * t52 + t19 * t201 - t226 * t98 + t274 * t36 - t442 * t56 + t389, -t11 * t201 + t129 * t57 - t130 * t56 + t14 * t202 - t150 * t34 + t352 * t36 - t79 * t96 + t80 * t98 + t449, -t109 * t80 - t11 * t320 + t130 * t66 - t19 * t202 + t226 * t96 - t274 * t34 - t352 * t52 + t442 * t57 - t390, t11 * t96 + t57 * t34 + t19 * t109 + t66 * t52 + t14 * t98 + t56 * t36 - g(1) * (t341 + t544) - g(2) * (t347 + t427) -t128 * t20 - t47 * t95, t127 * t20 - t128 * t21 + t47 * t553 - t48 * t95, -t128 * t222 - t20 * t320 - t262 * t47 - t442 * t95, t127 * t21 + t48 * t553, t127 * t222 - t21 * t320 - t262 * t48 + t442 * t553, -t222 * t320 - t262 * t442, g(1) * t372 - g(2) * t121 + t10 * t127 - t12 * t442 + t2 * t320 + t97 * t21 - t30 * t222 + t4 * t262 + t33 * t553 + t51 * t48, g(1) * t540 - g(2) * t120 + t10 * t128 + t13 * t442 - t97 * t20 + t31 * t222 - t3 * t262 + t320 * t362 + t33 * t95 - t51 * t47, t12 * t47 + t127 * t362 - t128 * t2 - t13 * t48 + t20 * t30 - t21 * t31 - t3 * t553 - t4 * t95 - t449, -t362 * t31 + t13 * t3 + t2 * t30 + t12 * t4 + t10 * t97 + t51 * t33 - g(1) * (-t194 * pkin(5) + t424 + t544) - g(2) * (t196 * pkin(5) - t316 * t412 + t347) - (pkin(9) * t316 + t408) * t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, t447 * t323, t431, t425, t300, qJDD(2), -t292 - t308 + (t388 + t530) * t316 (-t495 + t530) * t320 + t423, 0, 0, -t274 * t476 - t493 (-t141 + t480) * t319 + (-t142 + t477) * t315 (-t236 * t316 + t274 * t460) * qJD(1) + t355, -t274 * t479 - t492 (t234 * t316 - t274 * t468) * qJD(1) + t354, t414, -pkin(2) * t142 + t171 * t274 + t377 * t315 + (-t160 * t316 - t320 * t381) * qJD(1) + t568 * t319, pkin(2) * t141 - t172 * t274 + t377 * t319 + (t161 * t316 - t320 * t380) * qJD(1) - t568 * t315, t171 * t236 + t172 * t234 + ((qJD(3) * t236 - t142) * pkin(8) + t547) * t319 + (-t71 + t484 + (qJD(3) * t234 - t141) * pkin(8)) * t315 - t423, -t252 * t294 - t160 * t171 - t161 * t172 + t333 * pkin(2) + (qJD(3) * t373 - t71 * t315 + t70 * t319 - t557) * pkin(8), t29, t533, t55, t353, -t538, t414, t150 * t395 - t227 * t344 + t274 * t503 - t291 * t79 + t351 * t454 - t445 * t60 + t345, -t228 * t344 + t274 * t501 - t291 * t80 + t351 * t453 + t352 * t395 + t445 * t61 - t566, t150 * t83 - t16 * t227 + t228 * t514 + t352 * t503 + t453 * t60 + t454 * t61 + t349, t16 * t168 + t514 * t167 + t344 * t291 - g(3) * (t316 * t517 + t256) + t501 * t61 - t503 * t60 - t395 * t351 + t388 * (t291 * t316 - t320 * t517) t29, t55, -t533, t414, t538, t353, t145 * t79 + t150 * t455 + t19 * t227 + t274 * t504 + t445 * t56 - t454 * t66 + t345, -t11 * t227 + t14 * t228 + t150 * t73 + t352 * t504 - t453 * t56 + t454 * t57 + t349, -t145 * t80 - t19 * t228 - t274 * t502 - t352 * t455 - t445 * t57 + t453 * t66 + t566, t11 * t168 + t19 * t145 + t14 * t167 - t397 + t455 * t66 + t502 * t57 + t504 * t56 - t517 * t360 + (-g(3) * t517 + t388 * (pkin(4) * t297 - t402)) * t316, -t148 * t20 + t498 * t95, t147 * t20 - t148 * t21 - t497 * t95 - t498 * t553, -t148 * t222 + t262 * t498 + t445 * t95, t147 * t21 + t497 * t553, t147 * t222 - t262 * t497 - t445 * t553, t262 * t445, t10 * t147 + t110 * t21 + t12 * t445 - t68 * t222 + t515 * t262 - t337 * t369 + t497 * t51 + t513 * t553, t10 * t148 - t110 * t20 - t13 * t445 + t69 * t222 - t516 * t262 - t337 * t370 + t498 * t51 + t513 * t95, -t12 * t498 - t13 * t497 + t147 * t362 - t148 * t2 + t20 * t68 - t21 * t69 - t515 * t95 - t516 * t553 + t423, -t362 * t69 + t2 * t68 + t10 * t110 - t397 + t513 * t51 + t516 * t13 + t515 * t12 + (-g(3) * pkin(5) * t297 + g(1) * t412 + t518 * t523) * t320 + (g(3) * t518 + t388 * (t297 * t532 - t402)) * t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t478, -t234 ^ 2 + t236 ^ 2, -t141 - t480, -t478, -t142 - t477, t226, -t253 * t438 - g(1) * t215 + g(2) * t213 - t484 - t252 * t236 + t156 + (-qJD(3) * t223 - t206 + t307) * t315, g(1) * t216 - g(2) * t214 + g(3) * t465 + t234 * t252 - t547, 0, 0, t559, t374, t53, -t559, t357, t226, t351 * t352 + (-t150 * t236 + t226 * t496) * pkin(3) + t334, -t351 * t150 + (-t226 * t312 - t236 * t352) * pkin(3) + t328, t61 * t352 - t508 + (-t312 * t79 - t496 * t80) * pkin(3) + (t64 - t60) * t150, -g(1) * t278 + t60 * t63 - t61 * t64 + (g(2) * t459 + t16 * t312 + t351 * t236 + t315 * t557 - t514 * t496) * pkin(3), t559, t53, -t374, t226, -t357, -t559, -t72 * t150 + (pkin(4) - t287) * t226 + t334 + t536, -t285 * t79 + t287 * t80 + t352 * t57 - t508 + (-t456 + t56) * t150, -t150 * t66 + t226 * t285 + t352 * t72 + t221 - 0.2e1 * t258 - t328, t11 * t285 + t14 * t287 - t66 * t72 - t56 * t63 - g(1) * t339 - g(2) * t329 - g(3) * (-pkin(4) * t475 + t450) + t456 * t57, -t519, t569, t548, t519, t570, t222, -t187 * t222 + t262 * t499 - t54 * t553 - t2 - t567, t188 * t222 - t262 * t500 - t54 * t95 - t563, t187 * t20 - t188 * t21 + (t12 - t500) * t553 + (-t13 - t499) * t95, -t362 * t188 + t2 * t187 - t51 * t54 - g(1) * (-pkin(5) * t195 + t339) - g(2) * (-t193 * pkin(5) + t329) - g(3) * (-t296 * t428 + t450) + t500 * t13 + t499 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, -t336, t375, t150 * t61 + t352 * t60 + t326, 0, 0, 0, 0, 0, 0, t358, t375, t336, t150 * t57 - t352 * t56 + t325 + t531, 0, 0, 0, 0, 0, 0, -t21 - t506, t20 + t509, t520 + t521, -t12 * t95 - t13 * t553 + t325 + t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226 + t559, t53, -t274 ^ 2 - t491, t274 * t57 - t348 - t529 - t536, 0, 0, 0, 0, 0, 0, -t318 * t222 - t314 * t400 - t352 * t553, t314 * t222 - t318 * t400 - t352 * t95, -t314 * t570 + t548 * t318, -t51 * t352 + t550 * t318 + (-t362 - t512) * t314 + t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t519, -t569, -t548, -t519, -t570, -t222, t550 + t567, t512 + t563, 0, 0;];
tau_reg  = t1;