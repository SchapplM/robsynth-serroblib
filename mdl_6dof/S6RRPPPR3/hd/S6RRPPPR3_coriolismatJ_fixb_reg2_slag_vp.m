% Calculate inertial parameters regressor of coriolis matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPPPR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:46
% EndTime: 2019-03-09 08:15:59
% DurationCPUTime: 9.40s
% Computational Cost: add. (8860->517), mult. (14978->641), div. (0->0), fcn. (15143->6), ass. (0->396)
t342 = sin(pkin(9));
t560 = cos(qJ(6));
t428 = t560 * t342;
t343 = cos(pkin(9));
t345 = sin(qJ(6));
t498 = t345 * t343;
t274 = t428 + t498;
t346 = sin(qJ(2));
t231 = t274 * t346;
t497 = t345 * t346;
t302 = t342 * t497;
t427 = t560 * t343;
t233 = t346 * t427 - t302;
t499 = t345 * t342;
t270 = -t427 + t499;
t43 = t231 * t270 + t233 * t274;
t588 = t43 * qJD(3);
t344 = qJ(3) + pkin(4);
t347 = cos(qJ(2));
t573 = pkin(2) + pkin(3);
t438 = qJ(5) + t573;
t400 = t438 * t347;
t587 = t344 * t346 + t400;
t556 = pkin(8) + t438;
t586 = t438 * t346;
t585 = t573 * t347;
t449 = t274 * qJD(6);
t361 = -t498 / 0.2e1 - t428 / 0.2e1;
t570 = -t274 / 0.2e1;
t351 = (t570 + t361) * t346;
t576 = t351 * qJD(1);
t579 = t576 - t449;
t336 = t342 ^ 2;
t337 = t343 ^ 2;
t306 = t337 + t336;
t450 = t274 * qJD(2);
t303 = t347 * t499;
t234 = t347 * t427 - t303;
t479 = qJD(1) * t234;
t583 = t450 + t479;
t395 = -t427 / 0.2e1;
t569 = t303 / 0.2e1;
t572 = -t270 / 0.2e1;
t168 = t569 + (t395 + t572) * t347;
t582 = -qJD(1) * t168 + t450;
t454 = t270 * qJD(2);
t232 = t274 * t347;
t480 = qJD(1) * t232;
t581 = t454 + t480;
t244 = t232 / 0.2e1;
t359 = t361 * t347;
t166 = t244 + t359;
t580 = -qJD(1) * t166 + t454;
t519 = t233 * t270;
t525 = t231 * t274;
t385 = -t519 + t525;
t578 = qJD(1) * t385;
t577 = qJD(4) * t385;
t470 = qJD(6) * t351;
t564 = t346 / 0.2e1;
t396 = t560 * t564;
t565 = -t346 / 0.2e1;
t390 = t274 * t565 + t343 * t497 / 0.2e1 + t342 * t396;
t471 = qJD(6) * t390;
t490 = -t302 / 0.2e1 + t343 * t396;
t512 = t270 * t346;
t389 = t512 / 0.2e1 + t490;
t157 = t389 * qJD(6);
t575 = t234 ^ 2;
t265 = t270 ^ 2;
t574 = t274 ^ 2;
t243 = -t232 / 0.2e1;
t571 = t270 / 0.2e1;
t332 = t346 * qJ(4);
t568 = t332 / 0.2e1;
t333 = t347 * qJ(4);
t567 = -t333 / 0.2e1;
t566 = -t343 / 0.2e1;
t563 = -t347 / 0.2e1;
t562 = t347 / 0.2e1;
t510 = t274 * t234;
t522 = t232 * t270;
t103 = -t522 / 0.2e1 - t510 / 0.2e1;
t561 = t103 * qJD(5);
t559 = t342 * pkin(5);
t558 = t346 * pkin(7);
t557 = t347 * pkin(7);
t105 = t510 + t522;
t555 = t105 * qJD(5);
t296 = -t332 + t558;
t277 = t343 * t296;
t154 = t277 + t342 * (pkin(1) + t587);
t131 = pkin(8) * t342 * t347 + t154;
t502 = t343 * t347;
t375 = t556 * t502;
t505 = t342 * t296;
t74 = t345 * t131 - t560 * (t343 * pkin(1) - t505 + (t343 * t344 + pkin(5)) * t346 + t375);
t496 = t346 * qJ(3);
t405 = pkin(1) + t496;
t352 = -t505 + t343 * (t346 * pkin(4) + t405);
t75 = t560 * t131 + (t346 * pkin(5) + t352 + t375) * t345;
t22 = t231 * t75 - t233 * t74;
t227 = t347 * t344 - t586;
t298 = -t333 + t557;
t155 = t343 * t227 - t298 * t342;
t503 = t343 * t346;
t130 = pkin(5) * t347 - pkin(8) * t503 + t155;
t430 = t560 * t130;
t156 = t342 * t227 + t298 * t343;
t504 = t342 * t346;
t140 = -pkin(8) * t504 + t156;
t500 = t345 * t140;
t76 = t430 - t500;
t429 = t560 * t140;
t501 = t345 * t130;
t77 = t429 + t501;
t6 = t232 * t77 + t234 * t76 - t22;
t554 = t6 * qJD(1);
t553 = t74 * t274;
t552 = t76 * t274;
t551 = t77 * t270;
t432 = t553 / 0.2e1;
t54 = -t553 / 0.2e1;
t10 = t54 + t432;
t549 = qJD(1) * t10;
t280 = t556 * t343;
t406 = t556 * t342;
t205 = -t345 * t280 - t560 * t406;
t206 = -t560 * t280 + t345 * t406;
t310 = pkin(5) * t343 + t344;
t354 = -t231 * t205 / 0.2e1 - t233 * t206 / 0.2e1 + t310 * t563;
t374 = t77 * t570 + t76 * t571;
t18 = t354 + t374;
t548 = qJD(1) * t18;
t431 = -pkin(7) + t559;
t248 = -t347 * t431 - t333;
t35 = -t248 * t232 - t346 * t74;
t547 = qJD(1) * t35;
t36 = -t234 * t248 - t346 * t75;
t546 = qJD(1) * t36;
t516 = t234 * t270;
t197 = -t516 / 0.2e1;
t409 = t516 / 0.2e1;
t50 = t197 + t409;
t543 = qJD(1) * t50;
t153 = t343 * t400 + t352;
t509 = t298 * t347;
t73 = -t153 * t504 + t154 * t503 + t509;
t542 = qJD(1) * t73;
t387 = t153 * t343 + t154 * t342;
t84 = t387 * t346;
t541 = qJD(1) * t84;
t85 = t387 * t347;
t540 = qJD(1) * t85;
t92 = t231 * t232 + t233 * t234;
t539 = qJD(1) * t92;
t523 = t232 * t233;
t526 = t231 * t234;
t94 = -t523 + t526;
t538 = qJD(1) * t94;
t511 = t274 * t232;
t408 = -t511 / 0.2e1;
t95 = t563 + t408 + t409;
t537 = qJD(1) * t95;
t536 = qJD(3) * t50;
t535 = t10 * qJD(4);
t247 = t346 * t431 + t332;
t11 = t247 * t248 - t74 * t76 + t75 * t77;
t534 = t11 * qJD(1);
t533 = t155 * t342;
t532 = t156 * t343;
t19 = t231 * t248 - t232 * t247 + t346 * t76 - t347 * t74;
t531 = t19 * qJD(1);
t20 = -t248 * t233 + t247 * t234 + t346 * t77 + t347 * t75;
t530 = t20 * qJD(1);
t21 = -t231 * t74 - t233 * t75 - t248 * t347;
t529 = t21 * qJD(1);
t528 = t22 * qJD(1);
t527 = t231 * t206;
t524 = t231 * t346;
t521 = t232 * t347;
t520 = t233 * t205;
t517 = t233 * t346;
t515 = t234 * t347;
t355 = t206 * t565 + t248 * t570 - t310 * t234 / 0.2e1;
t362 = -t500 / 0.2e1 + t430 / 0.2e1;
t24 = -t355 + t362;
t514 = t24 * qJD(1);
t356 = t205 * t564 + t310 * t244 + t248 * t571;
t363 = -t501 / 0.2e1 - t429 / 0.2e1;
t25 = -t356 + t363;
t513 = t25 * qJD(1);
t32 = (-t153 * t346 + t155 * t347) * t343 + (-t154 * t346 + t156 * t347) * t342;
t508 = t32 * qJD(1);
t34 = t153 * t155 + t154 * t156 - t296 * t298;
t507 = t34 * qJD(1);
t341 = t347 ^ 2;
t506 = t341 * t343;
t495 = t347 * qJ(3);
t67 = t153 * t347 + t155 * t346 + (t296 * t347 + t298 * t346) * t342;
t494 = t67 * qJD(1);
t68 = t156 * t346 - t298 * t503 + (t154 - t277) * t347;
t493 = t68 * qJD(1);
t407 = -t337 / 0.2e1 - t336 / 0.2e1;
t357 = t344 * t563 - t407 * t586;
t370 = t155 * t566 - t156 * t342 / 0.2e1;
t70 = t357 + t370;
t492 = t70 * qJD(1);
t93 = t523 + t526;
t491 = t93 * qJD(1);
t340 = t346 ^ 2;
t311 = t340 + t341;
t312 = t341 - t340;
t148 = t521 + t524;
t489 = qJD(1) * t148;
t149 = t521 - t524;
t488 = qJD(1) * t149;
t150 = -t515 + t517;
t487 = qJD(1) * t150;
t151 = -t515 - t517;
t486 = qJD(1) * t151;
t264 = t405 + t585;
t266 = -t573 * t346 + t495;
t180 = t264 * t347 + t266 * t346;
t483 = qJD(1) * t180;
t181 = -t264 * t346 + t266 * t347;
t482 = qJD(1) * t181;
t222 = t296 * t346 + t509;
t481 = qJD(1) * t222;
t478 = qJD(2) * t310;
t477 = qJD(3) * t346;
t273 = t306 * t346;
t476 = qJD(4) * t273;
t475 = qJD(4) * t347;
t474 = qJD(5) * t232;
t473 = qJD(5) * t234;
t472 = qJD(5) * t346;
t469 = qJD(6) * t234;
t468 = qJD(6) * t346;
t100 = t264 * t266;
t467 = t100 * qJD(1);
t466 = t390 * qJD(1);
t167 = -t512 / 0.2e1 + t490;
t464 = t167 * qJD(1);
t463 = t389 * qJD(1);
t171 = t302 / 0.2e1 + (t571 + t395) * t346;
t159 = t171 * qJD(1);
t182 = 0.2e1 * t243;
t461 = t182 * qJD(1);
t391 = -t347 * pkin(2) - t496;
t287 = -pkin(1) + t391;
t297 = t346 * pkin(2) - t495;
t219 = t287 * t347 + t297 * t346;
t460 = t219 * qJD(1);
t220 = -t287 * t346 + t297 * t347;
t459 = t220 * qJD(1);
t437 = -pkin(2) / 0.2e1 - pkin(3) / 0.2e1;
t228 = t495 + (-t573 / 0.2e1 + t437) * t346;
t458 = t228 * qJD(1);
t230 = t347 * t273;
t457 = t230 * qJD(1);
t286 = -0.1e1 / 0.2e1 + t407;
t235 = t286 * t347;
t456 = t235 * qJD(1);
t269 = t306 * t341;
t455 = t269 * qJD(1);
t255 = t270 * qJD(6);
t271 = t311 * t342;
t453 = t271 * qJD(1);
t316 = t340 * t343;
t272 = t316 + t506;
t452 = t272 * qJD(1);
t451 = t273 * qJD(1);
t275 = t312 * t342;
t448 = t275 * qJD(1);
t276 = -t316 + t506;
t447 = t276 * qJD(1);
t446 = t286 * qJD(2);
t445 = t306 * qJD(2);
t444 = t311 * qJD(1);
t443 = t312 * qJD(1);
t299 = t312 * qJD(2);
t328 = t340 * qJD(1);
t327 = t340 * qJD(3);
t442 = t342 * qJD(2);
t441 = t343 * qJD(2);
t330 = t346 * qJD(1);
t329 = t346 * qJD(2);
t440 = t347 * qJD(1);
t331 = t347 * qJD(2);
t439 = t347 * qJD(3);
t436 = pkin(1) * t330;
t435 = pkin(1) * t440;
t434 = pkin(7) * t329;
t433 = pkin(7) * t331;
t426 = t231 * t330;
t425 = t232 * t479;
t424 = t233 * t330;
t423 = t270 * t450;
t422 = t270 * t331;
t421 = t342 * t441;
t420 = t342 * t331;
t419 = t343 * t331;
t418 = t346 * t439;
t417 = t270 * t449;
t416 = t287 * t297 * qJD(1);
t415 = t287 * t330;
t414 = t342 * t328;
t413 = t343 * t328;
t412 = t342 * t440;
t411 = t343 * t440;
t314 = t346 * t331;
t313 = t346 * t440;
t410 = t270 * t440;
t404 = 0.2e1 * t342 * t502;
t226 = t234 * t330;
t402 = qJD(2) * t351 - t226;
t401 = qJD(5) + t478;
t399 = t343 * t313;
t398 = pkin(7) / 0.2e1 - t559 / 0.2e1;
t394 = -t265 / 0.2e1 - t574 / 0.2e1;
t393 = t568 - t558 / 0.2e1;
t392 = t567 + t557 / 0.2e1;
t358 = t205 * t234 / 0.2e1 + t206 * t243 + t75 * t572;
t360 = -t398 * t346 + t568;
t13 = t432 + t358 + t360;
t86 = -t205 * t274 + t206 * t270;
t388 = -qJD(1) * t13 + qJD(2) * t86;
t386 = -t532 + t533;
t229 = t232 ^ 2;
t110 = t229 - t575;
t81 = -t510 + t522;
t384 = qJD(1) * t110 + qJD(2) * t81;
t152 = t265 - t574;
t383 = qJD(1) * t81 + qJD(2) * t152;
t349 = t398 * t347 + t567 + t527 / 0.2e1 - t520 / 0.2e1;
t373 = -t552 / 0.2e1 - t551 / 0.2e1;
t16 = t349 - t373;
t382 = qJD(1) * t16 + t478;
t251 = t306 * t438;
t372 = t153 * t342 / 0.2e1 + t154 * t566;
t79 = -t372 + t393;
t381 = -qJD(1) * t79 + qJD(2) * t251;
t371 = t533 / 0.2e1 - t532 / 0.2e1;
t83 = t371 + t392;
t380 = qJD(1) * t83 + qJD(2) * t344;
t23 = t232 * t75 - t234 * t74;
t379 = -qJD(1) * t23 - qJD(3) * t103;
t109 = -0.1e1 / 0.2e1 + t394;
t378 = qJD(1) * t103 + qJD(2) * t109;
t162 = t229 + t575;
t377 = qJD(1) * t162 + qJD(2) * t105;
t210 = t265 + t574;
t376 = qJD(1) * t105 + qJD(2) * t210;
t199 = t511 / 0.2e1;
t102 = t409 + t199;
t369 = qJD(2) * t102 + t425;
t368 = qJD(1) * t102 + t423;
t104 = t408 + t197;
t367 = -qJD(2) * t104 + t425;
t366 = -qJD(1) * t104 + t423;
t365 = t347 * t472 + t327;
t364 = qJD(2) * t390 - t226 - t469;
t353 = qJD(2) * t391 + t439;
t350 = qJD(2) * t587 - t439 + t472;
t339 = qJ(3) * qJD(3);
t338 = qJD(2) * qJ(3);
t318 = t331 / 0.2e1;
t285 = 0.1e1 / 0.2e1 + t407;
t282 = t298 * qJD(2);
t268 = qJD(6) * t562 + t313;
t246 = t437 * t346 - t565 * t573;
t236 = t306 * t562 + t563;
t183 = t244 + t243;
t169 = t270 * t562 + t347 * t395 + t569;
t165 = t243 + t359;
t158 = t171 * qJD(6);
t114 = t255 + t159;
t108 = 0.1e1 / 0.2e1 + t394;
t106 = -qJD(2) * t171 - t232 * t330;
t99 = t104 * qJD(6);
t97 = t102 * qJD(6);
t96 = t199 + t197 + t563;
t91 = t389 * qJD(2) + (qJD(6) + t330) * t232;
t82 = -t371 + t392;
t80 = t372 + t393;
t78 = t81 * qJD(6);
t69 = t357 - t370;
t42 = t50 * qJD(6);
t27 = t355 + t362;
t26 = t356 + t363;
t17 = t354 - t374;
t15 = t349 + t373;
t14 = t54 - t358 + t360;
t8 = t10 * qJD(6);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t299, 0, -t314, 0, 0, -pkin(1) * t329, -pkin(1) * t331, 0, 0, t314, 0, -t299, 0, 0, -t314, -qJD(2) * t220 + t418, 0, -qJD(2) * t219 + t327 (qJD(2) * t297 - t477) * t287, -t314, t299, 0, t314, 0, 0, qJD(2) * t180 + t327, -qJD(2) * t181 - t418, qJD(4) * t311, qJD(2) * t100 - qJD(4) * t222 + t264 * t477, -t337 * t314, t404 * t329, -t276 * qJD(2), -t336 * t314, t275 * qJD(2), t314, t67 * qJD(2) + t271 * qJD(4) + t343 * t365, -t68 * qJD(2) + t272 * qJD(4) - t342 * t365, qJD(2) * t32 + qJD(3) * t230 + qJD(5) * t269, qJD(2) * t34 + qJD(3) * t84 - qJD(4) * t73 + qJD(5) * t85 (-qJD(2) * t233 - qJD(6) * t232) * t234, qJD(2) * t93 + qJD(6) * t110, qJD(2) * t150 + t232 * t468 (-qJD(2) * t231 + t469) * t232, qJD(2) * t149 + t234 * t468, t314, t19 * qJD(2) + t148 * qJD(4) + t36 * qJD(6) + (qJD(3) * t233 + t473) * t346, -t20 * qJD(2) - t151 * qJD(4) - t35 * qJD(6) + (-qJD(3) * t231 - t474) * t346, qJD(2) * t6 + qJD(3) * t92 + qJD(4) * t94 + qJD(5) * t162, qJD(2) * t11 + qJD(3) * t22 + qJD(4) * t21 + qJD(5) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t443, t331, -t313, -t329, 0, -t433 - t436, t434 - t435, 0, 0, t313, t331, -t443, 0, t329, -t313, -t433 - t459, t353, -t434 - t460, pkin(7) * t353 + t416, -t313, t443, -t329, t313, t331, 0, -qJD(2) * t296 + t483, t282 - t482 (t496 + t585) * qJD(2) - t439, t467 + (-t296 * qJ(3) - t298 * t573) * qJD(2) + t298 * qJD(3) + t246 * qJD(4) (-t337 * t440 - t421) * t346 (qJD(1) * t404 + (t336 - t337) * qJD(2)) * t346, -t420 - t447 (-t336 * t440 + t421) * t346, -t419 + t448, t313, -t296 * t441 + t342 * t350 + t494, t296 * t442 + t343 * t350 - t493, qJD(2) * t386 + t508, t507 + (-t296 * t344 + t386 * t438) * qJD(2) + t82 * qJD(3) + t69 * qJD(4) + t80 * qJD(5), -t583 * t233 + t99, t491 + (t519 + t525) * qJD(2) + t78, -t274 * t331 + t157 + t487, -t581 * t231 + t97, t422 + t488 - t471, t268, t531 + (-t205 * t347 + t231 * t310 - t247 * t270) * qJD(2) + t165 * qJD(3) - t351 * qJD(5) + t27 * qJD(6), -t530 + (-t206 * t347 + t233 * t310 - t247 * t274) * qJD(2) + t169 * qJD(3) + t183 * qJD(4) - t171 * qJD(5) + t26 * qJD(6), t554 + (t520 - t527 + t551 + t552) * qJD(2) + t588 + t555, t534 + (-t205 * t76 + t206 * t77 + t247 * t310) * qJD(2) + t15 * qJD(3) + t17 * qJD(4) + t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t331, t328, -t415 + t433, 0, 0, 0, 0, 0, 0, t328, -t313, -t331, t264 * t330 + t282, 0, 0, 0, 0, 0, 0, t413 - t420, -t414 - t419, t457, qJD(2) * t82 + t541, 0, 0, 0, 0, 0, 0, qJD(2) * t165 + t157 + t424, qJD(2) * t169 - t426 - t471, qJD(2) * t43 + t42 + t539, t15 * qJD(2) + t528 + t561 - t588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t444, qJD(2) * t246 - t481, 0, 0, 0, 0, 0, 0, t453, t452, 0, qJD(2) * t69 + qJD(5) * t236 - t542, 0, 0, 0, 0, 0, 0, t471 + t489, qJD(2) * t183 + t157 - t486, t538, t17 * qJD(2) - qJD(4) * t43 + t96 * qJD(5) + t529 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t411 + t442) * t346 (-t412 + t441) * t346, t455, qJD(2) * t80 + qJD(4) * t236 + t540, 0, 0, 0, 0, 0, 0, -t402, t106, t377, qJD(2) * t14 + qJD(4) * t96 - t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t367, t384, t91, t369, -t364, t318, qJD(2) * t27 + qJD(3) * t389 + qJD(4) * t390 - qJD(6) * t75 + t546, qJD(2) * t26 - qJD(3) * t390 + qJD(4) * t389 + qJD(6) * t74 - t547, t536, t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, -t443, 0, t313, 0, 0, t436, t435, 0, 0, -t313, 0, t443, 0, 0, t313, t459, 0, t460, -t416, t313, -t443, 0, -t313, 0, 0, -t475 - t483, -qJD(4) * t346 + t482, 0, -qJD(4) * t228 - t467, t337 * t313, -0.2e1 * t342 * t399, t447, t336 * t313, -t448, -t313, -t343 * t475 - t494, t342 * t475 + t493, t476 - t508, qJD(3) * t83 + qJD(4) * t70 - qJD(5) * t79 - t507, t233 * t479 + t99, t78 - t491, t158 - t487, t231 * t480 + t97, -t488 - t470, -t268, qJD(3) * t166 - qJD(5) * t390 - qJD(6) * t24 + t270 * t475 - t531, qJD(3) * t168 - qJD(4) * t182 - qJD(5) * t389 - qJD(6) * t25 + t530, -t554 + t555 + t577, qJD(3) * t16 + qJD(4) * t18 - qJD(5) * t13 - t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t339, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, t339, 0, 0, 0, 0, 0, 0, qJD(3) * t343, -qJD(3) * t342, t306 * qJD(5), qJD(3) * t344 + qJD(5) * t251, -t417, t152 * qJD(6), 0, t417, 0, 0, -qJD(3) * t270 - t310 * t449, -qJD(3) * t274 + t255 * t310, qJD(5) * t210, qJD(3) * t310 + qJD(5) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t338, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, t338, 0, 0, 0, 0, 0, 0, t441, -t442, 0, qJD(5) * t285 + t380, 0, 0, 0, 0, 0, 0, -t580, -t582, 0, qJD(5) * t108 + t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t440, -t330, 0, -t458, 0, 0, 0, 0, 0, 0, -t411, t412, t451, t492, 0, 0, 0, 0, 0, 0, t410, -t461, t578, t548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t445, qJD(3) * t285 + t381, 0, 0, 0, 0, 0, 0, -t466, -t463, t376, qJD(3) * t108 + t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t366, t383, t114, t368, -t579, -t440 / 0.2e1, -qJD(6) * t206 - t310 * t450 - t514, qJD(6) * t205 + t310 * t454 - t513, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, 0, -t328, t415, 0, 0, 0, 0, 0, 0, -t328, t313, 0 (-qJD(1) * t264 - qJD(4)) * t346, 0, 0, 0, 0, 0, 0, -t413, t414, -t457, -qJD(2) * t83 - t476 - t541, 0, 0, 0, 0, 0, 0, -qJD(2) * t166 + t158 - t424, -qJD(2) * t168 + t426 - t470, t42 - t539, -qJD(2) * t16 - t528 + t561 - t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t338, 0, 0, 0, 0, 0, 0, -qJD(2), 0, 0, -t338, 0, 0, 0, 0, 0, 0, -t441, t442, 0, qJD(5) * t286 - t380, 0, 0, 0, 0, 0, 0, t580, t582, 0, qJD(5) * t109 - t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t330, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t578; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t446, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t579, t543, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, t329, -t444, qJD(2) * t228 + t477 + t481, 0, 0, 0, 0, 0, 0, t419 - t453, -t420 - t452, -t273 * qJD(2), -qJD(2) * t70 + qJD(3) * t273 - qJD(5) * t235 + t542, 0, 0, 0, 0, 0, 0, -t422 + t470 - t489, qJD(2) * t182 - qJD(6) * t167 + t486, -qJD(2) * t385 - t538, -qJD(2) * t18 + qJD(3) * t385 - qJD(5) * t95 - t529 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, t330, 0, t458, 0, 0, 0, 0, 0, 0, t411, -t412, -t451, -t492, 0, 0, 0, 0, 0, 0, -t410, t461, -t578, -t548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t330, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, t255 - t464, 0, t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t399, t342 * t313, -t455, qJD(2) * t79 + qJD(4) * t235 - t540, 0, 0, 0, 0, 0, 0, t364, t91, -t377, qJD(2) * t13 + qJD(4) * t95 + t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t445, -qJD(3) * t286 - t381, 0, 0, 0, 0, 0, 0, -t449 + t466, t255 + t463, -t376, -qJD(3) * t109 - t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t446, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t456, 0, 0, 0, 0, 0, 0, 0, 0, 0, t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583, t581, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t367, -t384, t106, -t369, t402, t318, qJD(2) * t24 - qJD(3) * t171 - qJD(4) * t351 + t473 - t546, qJD(2) * t25 + qJD(3) * t351 + qJD(4) * t167 - t474 + t547, -t536, -t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t366, -t383, -t159, -t368, t576, t440 / 0.2e1, t274 * t401 + t514, -t270 * t401 + t513, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t576, -t543, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t576, t464, 0, -t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, -t581, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
