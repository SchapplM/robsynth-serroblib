% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% qJDD [7x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% tau_reg [7x45]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [7x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 20:34:24
% EndTime: 2018-11-26 20:36:13
% DurationCPUTime: 49.45s
% Computational Cost: add. (24684->1090), mult. (54706->1597), div. (0->0), fcn. (46069->14), ass. (0->464)
t307 = cos(qJ(4));
t308 = cos(qJ(3));
t309 = cos(qJ(2));
t530 = t308 * t309;
t300 = sin(qJ(4));
t302 = sin(qJ(2));
t555 = t300 * t302;
t251 = t307 * t530 + t555;
t235 = t251 * qJD(1);
t301 = sin(qJ(3));
t515 = qJD(4) * t300;
t471 = t301 * t515;
t498 = t307 * qJD(3);
t473 = t308 * t498;
t361 = -t471 + t473;
t651 = -t235 - t361;
t303 = sin(qJ(1));
t310 = cos(qJ(1));
t547 = t301 * t310;
t246 = t303 * t530 + t547;
t199 = t246 * t307 + t303 * t555;
t529 = t310 * t308;
t543 = t303 * t301;
t245 = t309 * t543 - t529;
t299 = sin(qJ(5));
t306 = cos(qJ(5));
t131 = t199 * t299 + t245 * t306;
t297 = sin(qJ(7));
t304 = cos(qJ(7));
t132 = t199 * t306 - t245 * t299;
t546 = t302 * t307;
t198 = t246 * t300 - t303 * t546;
t298 = sin(qJ(6));
t305 = cos(qJ(6));
t87 = t132 * t305 + t198 * t298;
t650 = -t131 * t297 + t304 * t87;
t520 = qJD(1) * t309;
t479 = t301 * t520;
t514 = qJD(4) * t306;
t537 = t306 * t308;
t524 = (-qJD(5) - t498) * t537 + (t300 * t514 + (qJD(5) * t307 + qJD(3)) * t299) * t301 - t235 * t306 + t299 * t479;
t477 = t308 * t520;
t516 = qJD(3) * t308;
t490 = pkin(2) * t516;
t649 = -pkin(2) * t477 + pkin(3) * t651 - t490;
t481 = -pkin(2) * t307 - pkin(3);
t431 = t301 * t481;
t457 = t301 * t498;
t512 = qJD(4) * t308;
t517 = qJD(3) * t301;
t648 = pkin(3) * t517 + (t300 * t512 + t457) * pkin(2) - t431 * t520;
t521 = qJD(1) * t302;
t234 = t300 * t477 - t307 * t521;
t458 = t300 * t516;
t513 = qJD(4) * t307;
t469 = t301 * t513;
t360 = -t458 - t469;
t616 = t234 - t360;
t292 = t301 * qJD(2);
t377 = -t308 * t521 + t292;
t439 = qJD(3) + t520;
t397 = t300 * t439;
t208 = -t307 * t377 - t397;
t278 = t301 * t521;
t293 = t308 * qJD(2);
t376 = t278 + t293;
t497 = qJD(4) - t376;
t646 = t208 * t497;
t505 = qJD(6) * t300;
t645 = -t301 * t505 + t524;
t538 = t306 * t307;
t380 = t299 * t308 + t301 * t538;
t644 = -qJD(6) * t380 + t616;
t643 = t132 * t298 - t198 * t305;
t492 = t302 * qJDD(1);
t346 = qJD(2) * (-qJD(3) + t520) + t492;
t495 = qJD(1) * qJD(3);
t454 = t302 * t495;
t402 = qJDD(2) + t454;
t161 = (t301 * t346 + t402 * t308) * pkin(2);
t259 = t481 * t308;
t606 = pkin(3) * t307;
t261 = (pkin(2) + t606) * t301;
t188 = t259 * t299 - t261 * t306;
t638 = qJD(5) * t188 + t299 * t649 - t306 * t648;
t206 = -t300 * t377 + t307 * t439;
t540 = t305 * t306;
t109 = -t206 * t540 + t208 * t298;
t508 = qJD(5) * t306;
t462 = t305 * t508;
t642 = -t109 + t462;
t291 = t309 * qJDD(1);
t496 = qJD(1) * qJD(2);
t256 = t302 * t496 - qJDD(3) - t291;
t323 = -t402 * t301 + t346 * t308;
t80 = -qJD(4) * t397 + t323 * t300 + t307 * (-qJD(4) * t377 - t256);
t641 = -qJDD(5) - t80;
t236 = t306 * t497;
t144 = t208 * t299 - t236;
t138 = qJD(6) + t144;
t146 = t306 * t208 + t299 * t497;
t334 = qJD(5) + t206;
t95 = t305 * t146 + t298 * t334;
t42 = t138 * t304 + t297 * t95;
t185 = t305 * t334;
t93 = t146 * t298 - t185;
t528 = qJD(7) - t93;
t640 = t42 * t528;
t44 = -t297 * t138 + t304 * t95;
t639 = t44 * t528;
t581 = -t644 * t298 + t645 * t305;
t617 = qJD(3) * t497;
t635 = t307 * t617;
t506 = qJD(6) * t298;
t461 = t299 * t506;
t634 = -t461 + t642;
t177 = t299 * t377 - t376 * t538;
t509 = qJD(5) * t305;
t465 = t299 * t509;
t502 = qJD(6) * t306;
t539 = t305 * t307;
t563 = t298 * t300;
t526 = -t177 * t305 + t376 * t563 + (qJD(6) - t514) * t539 + (t465 + (-qJD(4) + t502) * t298) * t300;
t533 = t307 * t308;
t550 = t301 * t306;
t523 = t550 * qJD(3) + t380 * qJD(5) + t306 * t479 + (t533 * qJD(3) + t235 - t471) * t299;
t559 = t299 * t307;
t176 = -t306 * t377 - t376 * t559;
t464 = t300 * t508;
t472 = t299 * t513;
t351 = t176 + t464 + t472;
t255 = t376 * pkin(2);
t470 = t255 * t513;
t558 = t300 * t161;
t620 = t470 + t558;
t633 = qJD(2) * qJD(3) - t492;
t244 = -t300 * t309 + t302 * t533;
t551 = t301 * t302;
t487 = t299 * t551;
t197 = t244 * t306 - t487;
t532 = t307 * t309;
t545 = t302 * t308;
t243 = t300 * t545 + t532;
t129 = t197 * t298 - t243 * t305;
t455 = t309 * t496;
t493 = t301 * qJDD(2);
t608 = pkin(2) * t308;
t162 = qJD(2) * t490 + (-t455 - t492) * t608 + (t301 * t454 + t493) * pkin(2);
t445 = t300 * t256 + t307 * t323;
t79 = -t206 * qJD(4) + t445;
t50 = -pkin(3) * t79 + t162;
t354 = -t308 * qJDD(2) + t301 * t633;
t474 = t302 * t516;
t476 = t309 * t292;
t362 = t474 + t476;
t170 = -qJD(1) * t362 + t354;
t169 = -qJDD(4) - t170;
t536 = t307 * t161;
t59 = -pkin(3) * t169 + t255 * t515 - t536;
t404 = t299 * t50 + t306 * t59;
t534 = t307 * t255;
t178 = pkin(3) * t497 - t534;
t254 = t377 * pkin(2);
t338 = pkin(3) * t208 - t254;
t91 = t299 * t178 + t306 * t338;
t20 = -qJD(5) * t91 + t404;
t92 = t306 * t178 - t299 * t338;
t383 = t255 * t563 - t305 * t92;
t554 = t300 * t305;
t8 = qJD(6) * t383 - t161 * t554 - t20 * t298 - t305 * t470;
t253 = t309 * t529 - t543;
t544 = t302 * t310;
t205 = t253 * t307 + t300 * t544;
t252 = -t303 * t308 - t309 * t547;
t137 = t205 * t306 + t252 * t299;
t204 = t253 * t300 - t307 * t544;
t89 = -t137 * t298 + t204 * t305;
t364 = -g(1) * t89 + g(2) * t643 + g(3) * t129 + t8;
t630 = t528 ^ 2;
t631 = pkin(4) * t630 + t364;
t189 = t259 * t306 + t261 * t299;
t579 = qJD(5) * t189 + t299 * t648 + t306 * t649;
t629 = t144 * t334;
t628 = t206 * t439;
t609 = pkin(2) * t300;
t442 = t517 * t609;
t625 = t189 * t506 - t298 * t442 + t305 * t638;
t157 = (-qJD(4) + t293) * t532 + (-t457 + (qJD(2) - t512) * t300) * t302;
t623 = -qJD(5) * t551 + t157;
t330 = qJD(4) * t334;
t327 = t307 * t330;
t622 = t300 * t641 - t327;
t294 = t300 ^ 2;
t621 = t161 * t294 + 0.2e1 * t300 * t470;
t549 = t301 * t307;
t242 = t299 * t549 - t537;
t619 = -t300 * t376 + t515;
t618 = qJD(3) * t439;
t314 = t306 * t641;
t329 = qJD(5) * t334;
t562 = t299 * t300;
t615 = -t300 * t314 + t306 * t327 - t329 * t562;
t435 = t302 * t469;
t614 = -t476 * t609 + (-t302 * t458 - t435) * pkin(2);
t436 = t300 * t479;
t411 = pkin(2) * t436;
t468 = t307 * t512;
t613 = -pkin(2) * t468 + t411 + t442;
t369 = g(1) * t204 + g(2) * t198 + g(3) * t243;
t510 = qJD(5) * t299;
t31 = qJD(5) * t236 - t299 * t169 - t208 * t510 + t306 * t79;
t16 = t95 * qJD(6) + t298 * t31 + t305 * t641;
t612 = t255 * t377;
t428 = g(1) * t310 + g(2) * t303;
t347 = -g(3) * t309 + t428 * t302;
t15 = qJD(6) * t185 - t146 * t506 - t298 * t641 + t305 * t31;
t446 = t306 * t169 + t299 * t79;
t32 = qJD(5) * t146 + t446;
t30 = qJDD(6) + t32;
t501 = qJD(7) * t297;
t4 = -t138 * t501 + t15 * t297 + (qJD(7) * t95 + t30) * t304;
t37 = pkin(4) * t95 + t91;
t38 = -pkin(4) * t138 - t383;
t407 = t297 * t38 + t304 * t37;
t503 = qJD(6) * t305;
t459 = t300 * t503;
t7 = t305 * t20 - t255 * t459 - t298 * t620 - t506 * t92;
t5 = -pkin(4) * t30 + t7;
t451 = t299 * t59 - t306 * t50;
t21 = qJD(5) * t92 + t451;
t6 = t15 * pkin(4) + t21;
t1 = -t407 * qJD(7) - t297 * t6 + t304 * t5;
t611 = t470 + t369;
t610 = g(2) * t87;
t607 = pkin(3) * t306;
t604 = g(2) * t131;
t3 = -t42 * qJD(7) + t304 * t15 - t297 * t30;
t600 = t297 * t3;
t576 = t144 * t305;
t599 = (-pkin(4) * t576 + t92) * t528;
t556 = t300 * t301;
t195 = -t298 * t556 - t305 * t380;
t127 = t195 * t297 - t242 * t304;
t598 = -qJD(7) * t127 + t523 * t297 + t581 * t304;
t597 = t242 * t501 + (qJD(7) * t195 - t523) * t304 + t581 * t297;
t596 = -t189 * t503 + (-t305 * t468 + (t305 * t517 + t308 * t506) * t300) * pkin(2) + t305 * t411 + t638 * t298;
t248 = t298 * t307 - t300 * t540;
t595 = (qJD(7) * t562 + t526) * t304 + (-qJD(7) * t248 + t351) * t297;
t194 = t248 * t304 + t297 * t562;
t594 = qJD(7) * t194 + t526 * t297 - t351 * t304;
t593 = t138 * t93;
t592 = t15 * t298;
t14 = -qJDD(7) + t16;
t591 = t297 * t14;
t590 = t298 * t30;
t589 = t299 * t31;
t588 = t30 * t305;
t587 = t300 * t79;
t586 = t304 * t14;
t585 = t304 * t44;
t584 = t95 * t138;
t504 = qJD(6) * t304;
t541 = t304 * t306;
t565 = t297 * t299;
t583 = -t109 * t304 - t206 * t565 + (qJD(7) + t509) * t541 + (-t298 * t504 + (-qJD(7) * t305 - qJD(5)) * t297) * t299;
t542 = t304 * t305;
t247 = t297 * t306 + t299 * t542;
t561 = t299 * t304;
t582 = t247 * qJD(7) + t334 * t561 + (t540 * qJD(5) - t109 - t461) * t297;
t580 = t645 * t298 + t644 * t305;
t181 = pkin(3) * t377 + t254 * t307;
t186 = t376 * t606 + t255;
t100 = t181 * t306 + t186 * t299;
t578 = t305 * t100 + t254 * t563;
t577 = qJD(5) * t95;
t575 = t146 * t306;
t573 = t206 * t299;
t572 = t206 * t306;
t571 = t208 * t307;
t568 = t254 * t300;
t567 = t255 * t300;
t564 = t298 * t299;
t560 = t299 * t305;
t557 = t300 * t169;
t553 = t300 * t306;
t552 = t301 * t256;
t548 = t301 * t309;
t535 = t307 * t169;
t531 = t308 * t256;
t488 = t255 * t553;
t125 = pkin(3) * t573 + t488;
t527 = t305 * t125 - t298 * t534;
t357 = t298 * t513 + t459;
t466 = t300 * t510;
t525 = t376 * t554 - t305 * t515 + t306 * t357 - t307 * t506 + (t177 - t466) * t298;
t295 = t302 ^ 2;
t522 = -t309 ^ 2 + t295;
t519 = qJD(2) * t302;
t518 = qJD(2) * t309;
t511 = qJD(5) * t298;
t507 = qJD(6) * t297;
t500 = qJD(7) * t304;
t499 = t376 * qJD(2);
t491 = t300 * t608;
t489 = t255 * t562;
t485 = t300 * t551;
t484 = t301 * t554;
t483 = t302 * t550;
t482 = t301 * t544;
t480 = pkin(3) * t305 + pkin(4);
t475 = t302 * t517;
t467 = qJD(5) * t255 * t294;
t453 = -qJD(7) * t38 - t6;
t444 = t138 * t305;
t443 = pkin(2) * t485;
t438 = t298 * t485;
t437 = t302 * t484;
t434 = t299 * t459;
t399 = t305 * t189 - t298 * t491;
t112 = pkin(4) * t242 + t399;
t433 = t581 * pkin(4) + qJD(7) * t112 + t579;
t113 = pkin(4) * t195 + t188;
t432 = pkin(4) * t523 - qJD(7) * t113 + t298 * t411 - t357 * t608 - t625;
t427 = g(1) * t303 - g(2) * t310;
t426 = t138 * t509 + t15;
t425 = t138 * t511 + t16;
t210 = -pkin(2) * t545 - pkin(3) * t244;
t238 = t302 * t431;
t142 = t210 * t299 + t238 * t306;
t363 = t293 * t309 - t475;
t107 = -pkin(2) * t363 - pkin(3) * t157;
t140 = -t362 * pkin(3) + (-t302 * t473 + (t302 * t515 - t307 * t518) * t301) * pkin(2);
t401 = t210 * t306 - t238 * t299;
t35 = qJD(5) * t401 + t107 * t299 + t140 * t306;
t424 = -t142 * t506 + t305 * t35;
t124 = pkin(3) * t572 - t489;
t237 = t480 * t562;
t359 = -t306 * t513 + t466;
t99 = t181 * t299 - t186 * t306;
t422 = pkin(3) * t359 + t526 * pkin(4) + qJD(7) * t237 - t99;
t209 = -pkin(3) * t553 + pkin(4) * t248;
t421 = -qJD(7) * t209 + (t300 * t462 + (-t298 * t505 + t305 * t513) * t299) * pkin(3) - t578 + t351 * pkin(4);
t260 = (pkin(4) * t305 + pkin(3)) * t299;
t420 = -qJD(7) * t260 + (-t298 * t502 - t465) * pkin(3) - t527 + (-t573 - t510) * pkin(4);
t258 = t480 * t306;
t419 = pkin(3) * t508 + t634 * pkin(4) + qJD(7) * t258 + t124;
t418 = t309 * t497;
t108 = -t305 * t208 - t298 * t572;
t415 = -t298 * t508 + t108;
t413 = qJD(4) * t497;
t130 = t197 * t305 + t243 * t298;
t64 = pkin(4) * t130 - t401;
t126 = t305 * t142;
t196 = t244 * t299 + t483;
t74 = -pkin(2) * t438 - pkin(4) * t196 + t126;
t406 = -t297 * t74 - t304 * t64;
t405 = t297 * t64 - t304 * t74;
t403 = qJD(6) * t437;
t84 = t130 * t297 + t196 * t304;
t400 = t254 * t301 + t255 * t308;
t398 = t255 * t439;
t396 = t309 * t439;
t395 = t439 * t208;
t394 = qJD(1) * t439;
t393 = qJD(2) * t439;
t390 = qJD(3) * t376 + t493;
t388 = t497 * t515;
t387 = t307 * t413;
t386 = t300 * t617;
t385 = t500 * t528 - t591;
t384 = t501 * t528 + t586;
t382 = -t138 * t503 - t590;
t381 = -t138 * t506 + t588;
t379 = t302 * t377;
t378 = t377 * t308;
t375 = t309 * t393;
t374 = t309 * t394;
t372 = t308 * t492 - t493;
t136 = t205 * t299 - t252 * t306;
t371 = g(1) * t136 + g(3) * t196 + t604;
t370 = g(1) * t137 + g(2) * t132 + g(3) * t197;
t368 = g(1) * t205 + g(2) * t199 + g(3) * t244;
t367 = pkin(4) * t14;
t366 = t398 - t162;
t365 = t309 * t378;
t355 = (-pkin(2) * t169 - qJD(3) * t254 - t161) * t308;
t353 = -t21 + t371;
t350 = t503 * t299 - t415;
t349 = qJD(5) * t93 + t382;
t348 = -g(1) * t252 + g(2) * t245 + g(3) * t551 + t162;
t345 = t387 - t557;
t344 = (-pkin(4) * t146 - t305 * t91) * t528 + t370;
t343 = -t308 * t618 + t552;
t342 = qJD(5) * t244 + t362;
t340 = -pkin(4) * qJD(6) * t528 + t371;
t337 = -t376 * t497 + t413;
t336 = qJD(1) * t418 + t617;
t36 = qJD(5) * t142 - t107 * t306 + t140 * t299;
t335 = t618 + t374;
t333 = pkin(3) * t334;
t331 = t300 * t334;
t328 = qJD(5) * t333;
t326 = t306 * t329;
t320 = t206 * t334 + t329;
t317 = pkin(3) * t641;
t316 = t299 * t641;
t312 = qJD(1) ^ 2;
t311 = qJD(2) ^ 2;
t250 = t300 * t530 - t546;
t241 = t298 * t553 + t539;
t240 = t297 * t560 - t541;
t224 = t244 * t310;
t223 = t243 * t310;
t222 = t244 * t303;
t221 = t243 * t303;
t220 = t380 * t302;
t219 = t242 * t302;
t203 = t251 * t306 - t299 * t548;
t202 = t251 * t299 + t306 * t548;
t193 = -t298 * t380 + t484;
t192 = t248 * t297 - t300 * t561;
t175 = -t224 * t306 + t299 * t482;
t174 = -t224 * t299 - t306 * t482;
t173 = -t222 * t306 + t303 * t487;
t172 = -t222 * t299 - t303 * t483;
t171 = -t220 * t305 - t438;
t167 = t252 * t538 - t253 * t299;
t166 = t252 * t559 + t253 * t306;
t165 = -t245 * t538 - t246 * t299;
t164 = -t245 * t559 + t246 * t306;
t163 = -t243 * t540 + t244 * t298;
t158 = -t189 * t298 - t305 * t491;
t156 = -t300 * t475 - t309 * t515 - t307 * t519 + (t300 * t518 + t302 * t513) * t308;
t135 = t203 * t305 + t250 * t298;
t128 = t195 * t304 + t242 * t297;
t111 = t175 * t305 - t223 * t298;
t110 = t173 * t305 - t221 * t298;
t105 = -pkin(2) * t437 - t142 * t298;
t104 = t167 * t305 + t252 * t563;
t103 = t165 * t305 - t245 * t563;
t102 = -t204 * t540 + t205 * t298;
t101 = -t198 * t540 + t199 * t298;
t96 = -t125 * t298 - t305 * t534;
t90 = t137 * t305 + t204 * t298;
t85 = t130 * t304 - t196 * t297;
t73 = -t100 * t298 + t254 * t554;
t70 = -t342 * t299 + t623 * t306;
t69 = t623 * t299 + t306 * t342;
t57 = -t144 * t542 - t146 * t297;
t56 = t146 * t304 - t297 * t576;
t53 = -t255 * t554 - t298 * t92;
t41 = -t136 * t297 + t304 * t90;
t40 = -t136 * t304 - t297 * t90;
t33 = t37 * t501;
t27 = -qJD(6) * t129 + t156 * t298 + t305 * t70;
t26 = qJD(6) * t130 - t156 * t305 + t298 * t70;
t22 = -t142 * t503 - t298 * t35 + (-t305 * t435 + (-t305 * t474 + (t302 * t506 - t305 * t518) * t301) * t300) * pkin(2);
t18 = -t297 * t37 + t304 * t38;
t12 = pkin(4) * t27 + t36;
t11 = -pkin(4) * t69 + (-t403 + (-t300 * t362 - t435) * t298) * pkin(2) + t424;
t10 = -qJD(7) * t84 + t27 * t304 - t297 * t69;
t9 = t27 * t297 - t196 * t501 + (qJD(7) * t130 + t69) * t304;
t2 = -t297 * t5 + t304 * t453 + t33;
t13 = [qJDD(1), t427, t428, qJDD(1) * t295 + 0.2e1 * t302 * t455, 0.2e1 * t291 * t302 - 0.2e1 * t496 * t522, qJDD(2) * t302 + t309 * t311, qJDD(2) * t309 - t302 * t311, 0, t427 * t309, -t427 * t302, -qJD(2) * t365 + t323 * t545 + t379 * t517 (-t309 * t499 + (qJD(3) * t377 + t170) * t302) * t308 + (t292 * t518 + ((-qJDD(1) * t308 + t301 * t495) * t302 + (qJD(3) - 0.2e1 * t520) * t293 + t390) * t302) * t301 (-qJD(3) ^ 2 * t301 + qJD(2) * t377 - t531) * t302 + (0.2e1 * qJD(1) * t363 + t372) * t309, -t301 * t375 + t170 * t309 + (t343 + t499) * t302, -t256 * t309 - t302 * t393, -t254 * t519 + g(1) * t246 - g(2) * t253 + t162 * t309 + (-t308 * t375 + (t301 * t618 + t531) * t302) * pkin(2), -t255 * t519 - g(1) * t245 - g(2) * t252 + t161 * t309 + (t618 * t545 + (qJD(2) * t396 - t302 * t256) * t301) * pkin(2), t157 * t208 + t244 * t79, -t156 * t208 - t157 * t206 - t243 * t79 - t244 * t80, t157 * t497 - t244 * t169 - t208 * t362 - t551 * t79, -t156 * t497 + t243 * t169 + t206 * t362 + t551 * t80, -t617 * t545 + (-qJD(2) * t418 + t169 * t302) * t301, g(1) * t199 - g(2) * t205 + t254 * t156 + t162 * t243 + (-t255 * t556 + (-t308 * t206 + t497 * t556) * pkin(2)) * t518 + (-t161 * t556 + t360 * t255 + ((t386 - t80) * t308 + (t206 * qJD(3) + t345) * t301) * pkin(2)) * t302, -g(1) * t198 + g(2) * t204 + t254 * t157 + t162 * t244 + (-t301 * t534 + (-t308 * t208 + t497 * t549) * pkin(2)) * t518 + (-t301 * t536 - t361 * t255 + ((-t79 + t635) * t308 + (t208 * qJD(3) - t388 - t535) * t301) * pkin(2)) * t302, t146 * t70 + t197 * t31, -t144 * t70 - t146 * t69 - t196 * t31 - t197 * t32, t146 * t156 - t197 * t641 + t31 * t243 + t334 * t70, -t144 * t156 + t196 * t641 - t32 * t243 - t334 * t69, t156 * t334 - t243 * t641, g(1) * t132 - g(2) * t137 + t144 * t614 - t91 * t156 - t196 * t620 - t21 * t243 - t32 * t443 - t334 * t36 - t401 * t641 - t567 * t69, -g(1) * t131 + g(2) * t136 + t142 * t641 + t146 * t614 - t92 * t156 - t197 * t620 - t20 * t243 - t31 * t443 - t334 * t35 - t567 * t70, t130 * t15 + t27 * t95, -t129 * t15 - t130 * t16 - t26 * t95 - t27 * t93, t130 * t30 + t138 * t27 + t15 * t196 + t69 * t95, -t129 * t30 - t138 * t26 - t16 * t196 - t69 * t93, t138 * t69 + t196 * t30, g(1) * t87 - g(2) * t90 + t105 * t30 + t129 * t21 + t138 * t22 - t16 * t401 + t196 * t8 + t26 * t91 + t36 * t93 + t53 * t69, -t424 * t138 - t126 * t30 - t7 * t196 + t383 * t69 + t36 * t95 - t401 * t15 + t21 * t130 + t91 * t27 - g(1) * t643 - g(2) * t89 + (t138 * t403 + (t138 * t435 + (t138 * t362 + t30 * t551) * t300) * t298) * pkin(2), t10 * t44 + t3 * t85, -t10 * t42 - t3 * t84 - t4 * t85 - t44 * t9, t10 * t528 - t129 * t3 - t14 * t85 - t26 * t44, t129 * t4 + t14 * t84 + t26 * t42 - t528 * t9, t129 * t14 - t26 * t528 (qJD(7) * t405 - t11 * t297 - t12 * t304) * t528 - t406 * t14 - t2 * t129 + t407 * t26 + t22 * t42 + t105 * t4 + t8 * t84 + t53 * t9 + g(1) * t650 - g(2) * t41 -(qJD(7) * t406 + t11 * t304 - t12 * t297) * t528 - t405 * t14 + t1 * t129 + t18 * t26 + t22 * t44 + t105 * t3 + t8 * t85 + t53 * t10 - g(1) * (t131 * t304 + t297 * t87) - g(2) * t40; 0, 0, 0, -t302 * t312 * t309, t522 * t312, t492, t291, qJDD(2), t347, g(3) * t302 + t309 * t428 -(t308 * t455 + t372) * t301 + qJD(1) * t365 + (t301 * t376 + t378) * qJD(3) (-t292 * t439 - t170) * t301 + (t633 * t308 + (0.2e1 * t475 + (t376 + t278 - t293) * t309) * qJD(1) + t390) * t308, -qJD(1) * t379 - t308 * t374 + t343, t301 * t335 - t376 * t521 + t531, t302 * t394, -pkin(2) * t552 + t254 * t521 + (pkin(2) * t335 + t347) * t308, -pkin(2) * t531 + t255 * t521 + ((-qJD(1) * t396 - t618) * pkin(2) - t347) * t301, t208 * t651 - t79 * t549, t206 * t235 + t208 * t234 + (t206 * t307 + t208 * t300) * t516 + (t587 + t307 * t80 + (-t206 * t300 + t571) * qJD(4)) * t301, -t235 * t497 + (-t79 - t635) * t308 + (t300 * t413 + t395 + t535) * t301, t234 * t497 + (t386 + t80) * t308 + (t345 - t628) * t301, t169 * t308 + t336 * t301, g(1) * t224 + g(2) * t222 - g(3) * t251 - t254 * t234 - t400 * t513 + (t301 * t80 + (t387 + t628) * t308) * pkin(2) + (t355 + (-pkin(2) * t336 + t366) * t301) * t300, -g(1) * t223 - g(2) * t221 + g(3) * t250 - t254 * t235 + t400 * t515 + (t301 * t79 + (-t388 + t395) * t308) * pkin(2) + (t355 + ((-t497 * t520 - t617) * pkin(2) + t366) * t301) * t307, t146 * t524 - t31 * t380, -t144 * t524 + t146 * t523 + t242 * t31 + t32 * t380, -t616 * t146 - t31 * t556 + t524 * t334 + t380 * t641, t616 * t144 - t242 * t641 + t32 * t556 + t523 * t334, -t234 * t334 + t622 * t301 - t331 * t516, -g(1) * t175 - g(2) * t173 - g(3) * t203 + t613 * t144 + t188 * t641 + t21 * t556 + t620 * t242 - t32 * t491 - t579 * t334 + t523 * t567 + t616 * t91, g(1) * t174 + g(2) * t172 + g(3) * t202 + t613 * t146 + t189 * t641 + t20 * t556 - t31 * t491 + t638 * t334 + t620 * t380 - t524 * t567 + t616 * t92, t15 * t195 + t581 * t95, -t15 * t193 - t16 * t195 - t580 * t95 - t581 * t93, t138 * t581 - t15 * t242 + t195 * t30 - t523 * t95, -t138 * t580 + t16 * t242 - t193 * t30 + t523 * t93, -t138 * t523 - t242 * t30, -g(1) * t111 - g(2) * t110 - g(3) * t135 + t138 * t596 + t158 * t30 + t16 * t188 + t193 * t21 - t242 * t8 - t523 * t53 + t579 * t93 + t580 * t91, -t399 * t30 + t7 * t242 + t188 * t15 + t21 * t195 - g(1) * (-t175 * t298 - t223 * t305) - g(2) * (-t173 * t298 - t221 * t305) - g(3) * (-t203 * t298 + t250 * t305) + t579 * t95 + t581 * t91 - t523 * t383 + ((t308 * t459 + (-t436 + t468) * t298) * pkin(2) + t625) * t138, t128 * t3 + t44 * t598, -t127 * t3 - t128 * t4 - t42 * t598 - t44 * t597, -t128 * t14 - t193 * t3 - t44 * t580 + t528 * t598, t127 * t14 + t193 * t4 + t42 * t580 - t528 * t597, t14 * t193 - t528 * t580 -(-t112 * t297 - t113 * t304) * t14 - t2 * t193 + t158 * t4 + t8 * t127 - g(1) * (t111 * t304 - t174 * t297) - g(2) * (t110 * t304 - t172 * t297) - g(3) * (t135 * t304 - t202 * t297) - (t297 * t432 + t304 * t433) * t528 + t597 * t53 + t596 * t42 + t580 * t407 (t112 * t304 - t113 * t297) * t14 + t1 * t193 + t158 * t3 + t8 * t128 - g(1) * (-t111 * t297 - t174 * t304) - g(2) * (-t110 * t297 - t172 * t304) - g(3) * (-t135 * t297 - t202 * t304) - (-t297 * t433 + t304 * t432) * t528 + t598 * t53 + t596 * t44 + t580 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377 * t376, -t376 ^ 2 + t377 ^ 2, t376 * t439 + t323, -t292 * qJD(3) + (-t292 - t377) * t520 + t354, -t256, -t398 + t348, g(1) * t253 + g(2) * t246 + g(3) * t545 + t254 * t439 + t161, -t497 * t571 - t587 (t206 * t497 - t79) * t307 + (t80 + t646) * t300, -t208 * t377 - t307 * t337 + t557, t206 * t377 + t300 * t337 + t535, -t497 * t377, -t255 * t206 - t300 * t612 + t348 * t307, -t255 * t208 - t348 * t300 - t307 * t612, -t31 * t553 + (-t177 + t359) * t146, t144 * t177 + t146 * t176 + (t144 * t306 + t146 * t299) * t513 + (t589 + t306 * t32 + (-t144 * t299 + t575) * qJD(5)) * t300, -t146 * t619 - t177 * t334 + t31 * t307 - t615, t619 * t144 + t176 * t334 - t622 * t299 + t300 * t326 - t32 * t307, -t300 * t330 - t307 * t641 + t331 * t376, t615 * pkin(3) - g(1) * t167 - g(2) * t165 + g(3) * t220 - t144 * t568 + t176 * t567 - t21 * t307 + t621 * t299 + t306 * t467 + t334 * t99 + t619 * t91, g(1) * t166 + g(2) * t164 - g(3) * t219 + t100 * t334 - t146 * t568 + t177 * t567 - t20 * t307 - t299 * t467 + t621 * t306 + t317 * t562 - t328 * t553 - t333 * t472 + t619 * t92, t15 * t248 + t526 * t95, t15 * t241 - t16 * t248 + t525 * t95 - t526 * t93, t138 * t526 - t15 * t562 + t248 * t30 - t351 * t95, t138 * t525 + t16 * t562 + t241 * t30 + t351 * t93, -t138 * t351 - t30 * t562, -t8 * t562 - g(1) * t104 - g(2) * t103 - g(3) * t171 - t138 * t73 - t21 * t241 - t93 * t99 - t525 * t91 - t351 * t53 + ((-t138 * t564 - t306 * t93) * t513 + (t299 * t349 - t306 * t425) * t300) * pkin(3), t7 * t562 + t21 * t248 + t578 * t138 - t99 * t95 - g(1) * (-t167 * t298 + t252 * t554) - g(2) * (-t165 * t298 - t245 * t554) - g(3) * (t220 * t298 - t437) + t526 * t91 - t351 * t383 + ((-t138 * t560 - t306 * t95) * t513 + (-t426 * t306 + (-t381 + t577) * t299) * t300) * pkin(3), t194 * t3 + t44 * t595, -t192 * t3 - t194 * t4 - t42 * t595 - t44 * t594, -t14 * t194 + t241 * t3 + t44 * t525 + t528 * t595, t14 * t192 - t241 * t4 - t42 * t525 - t528 * t594, -t14 * t241 + t525 * t528 -(-t209 * t304 - t237 * t297) * t14 + t2 * t241 + t8 * t192 - t73 * t42 - g(1) * (t104 * t304 - t166 * t297) - g(2) * (t103 * t304 - t164 * t297) - g(3) * (t171 * t304 + t219 * t297) - (t297 * t421 + t304 * t422) * t528 + t594 * t53 - t525 * t407 + (-t42 * t434 + (-t42 * t464 + (-t300 * t4 - t42 * t513) * t299) * t298) * pkin(3) (-t209 * t297 + t237 * t304) * t14 - t1 * t241 + t8 * t194 - t73 * t44 - g(1) * (-t104 * t297 - t166 * t304) - g(2) * (-t103 * t297 - t164 * t304) - g(3) * (-t171 * t297 + t219 * t304) - (-t297 * t422 + t304 * t421) * t528 + t595 * t53 - t525 * t18 + (-t44 * t434 + (-t44 * t464 + (-t3 * t300 - t44 * t513) * t299) * t298) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208 * t206, -t206 ^ 2 + t208 ^ 2, -t206 * t376 + t445, -t80 + t646, -t169, -t254 * t208 + t376 * t534 + t369 + t558, t254 * t206 - t376 * t567 + t368 + t536, t334 * t575 + t589 (t31 - t629) * t306 + (-t146 * t334 - t32) * t299, -t146 * t208 + t306 * t320 - t316, t144 * t208 - t299 * t320 - t314, -t334 * t208, -t124 * t334 + t144 * t534 + t161 * t553 - t206 * t489 + t91 * t208 - t255 * t466 + (t316 - t326) * pkin(3) + t611 * t306, t125 * t334 + t146 * t534 - t206 * t488 + t92 * t208 - t255 * t464 + t306 * t317 + (t328 - t558 - t611) * t299, t15 * t560 + t634 * t95, t108 * t95 + t109 * t93 + (-t298 * t95 - t305 * t93) * t508 + (-t592 - t16 * t305 + (t298 * t93 - t305 * t95) * qJD(6)) * t299, -t15 * t306 + t642 * t138 + (t334 * t95 + t381) * t299, t16 * t306 + t415 * t138 + (-t334 * t93 + t382) * t299, t138 * t299 * t334 - t30 * t306, -g(1) * t102 - g(2) * t101 - g(3) * t163 - t108 * t91 + t124 * t93 - t138 * t96 + (pkin(3) * t349 + t511 * t91 - t8) * t306 + (pkin(3) * t425 + t21 * t298 + t334 * t53 + t503 * t91) * t299, t527 * t138 + t124 * t95 - t91 * t109 - t368 * t305 + (t91 * t509 + t7 + (t577 - t588) * pkin(3) + (pkin(3) * qJD(6) * t138 - t369) * t298) * t306 + (pkin(3) * t426 + t21 * t305 + t334 * t383 - t506 * t91) * t299, t247 * t3 + t44 * t583, -t240 * t3 - t247 * t4 - t42 * t583 - t44 * t582, -t14 * t247 - t3 * t564 - t350 * t44 + t528 * t583, t14 * t240 + t350 * t42 + t4 * t564 - t528 * t582, t14 * t564 - t350 * t528 -(-t258 * t297 - t260 * t304) * t14 + t8 * t240 - t407 * t108 - t96 * t42 - g(1) * (t102 * t304 + t204 * t565) - g(2) * (t101 * t304 + t198 * t565) - g(3) * (t163 * t304 + t243 * t565) - (t297 * t420 + t304 * t419) * t528 + t582 * t53 + (t299 * t407 - t42 * t607) * t503 + (t407 * t508 - t2 * t299 + (-t306 * t4 + t42 * t510) * pkin(3)) * t298 (t258 * t304 - t260 * t297) * t14 + t8 * t247 - t18 * t108 - t96 * t44 - g(1) * (-t102 * t297 + t204 * t561) - g(2) * (-t101 * t297 + t198 * t561) - g(3) * (-t163 * t297 + t243 * t561) - (-t297 * t419 + t304 * t420) * t528 + t583 * t53 + (t18 * t299 - t44 * t607) * t503 + (t18 * t508 + t1 * t299 + (-t3 * t306 + t44 * t510) * pkin(3)) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146 * t144, -t144 ^ 2 + t146 ^ 2, t31 + t629, t146 * t206 - t446, -t641, t146 * t567 + t206 * t92 + t371 - t451, -t144 * t567 - t91 * t206 + t370 - t404, t444 * t95 + t592 (t15 - t593) * t305 + (-t16 - t584) * t298, t138 * t444 - t146 * t95 + t590, -t138 ^ 2 * t298 + t146 * t93 + t588, -t138 * t146, -t146 * t53 + t305 * t353 - t92 * t93, -t146 * t383 - t298 * t353 - t92 * t95, t298 * t3 * t304 + (-t298 * t501 + t304 * t503 - t57) * t44, t42 * t57 + t44 * t56 + (-t297 * t44 - t304 * t42) * t503 + (-t600 - t304 * t4 + (t297 * t42 - t585) * qJD(7)) * t298, -t57 * t528 + (t504 * t528 + t3) * t305 + (-t138 * t44 - t384) * t298, t56 * t528 + (-t507 * t528 - t4) * t305 + (t138 * t42 - t385) * t298, -t138 * t298 * t528 - t14 * t305, t304 * t599 - t53 * t56 + t344 * t297 + (t304 * t340 + t507 * t53 + t2) * t305 + (pkin(4) * t384 + t138 * t407 + t297 * t8 - t91 * t42 + t500 * t53) * t298, -t297 * t599 - t53 * t57 + t344 * t304 + (-t297 * t340 + t504 * t53 - t1) * t305 + (pkin(4) * t385 + t138 * t18 + t8 * t304 - t91 * t44 - t501 * t53) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t93 ^ 2 + t95 ^ 2, t15 + t593, -t16 + t584, t30, -t138 * t383 - t91 * t95 + t364, g(1) * t90 + g(3) * t130 + t138 * t53 + t91 * t93 + t610 - t7, -t528 * t585 - t600 (-t3 + t640) * t304 + (t4 + t639) * t297, -t304 * t630 + t44 * t95 + t591, t297 * t630 - t42 * t95 + t586, t528 * t95, -t367 * t297 + t304 * t631 - t383 * t42 - t407 * t95, -t18 * t95 - t297 * t631 - t367 * t304 - t383 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t42 ^ 2 + t44 ^ 2, t3 + t640, -t4 + t639, -t14, -g(1) * t40 + g(3) * t84 + t18 * t528 - t53 * t44 + t33 + (-t5 + t610) * t297 + (t453 + t604) * t304, g(1) * t41 + g(2) * t650 + g(3) * t85 - t407 * t528 + t53 * t42 - t1;];
tau_reg  = t13;