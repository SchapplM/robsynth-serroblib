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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-03-10 07:43:03
% EndTime: 2019-03-10 07:44:53
% DurationCPUTime: 47.27s
% Computational Cost: add. (24684->1089), mult. (54706->1596), div. (0->0), fcn. (46069->14), ass. (0->462)
t307 = cos(qJ(4));
t308 = cos(qJ(3));
t309 = cos(qJ(2));
t531 = t308 * t309;
t300 = sin(qJ(4));
t302 = sin(qJ(2));
t556 = t300 * t302;
t251 = t307 * t531 + t556;
t235 = t251 * qJD(1);
t301 = sin(qJ(3));
t516 = qJD(4) * t300;
t471 = t301 * t516;
t499 = t307 * qJD(3);
t473 = t308 * t499;
t361 = -t471 + t473;
t652 = -t235 - t361;
t303 = sin(qJ(1));
t310 = cos(qJ(1));
t548 = t301 * t310;
t246 = t303 * t531 + t548;
t199 = t246 * t307 + t303 * t556;
t530 = t310 * t308;
t544 = t303 * t301;
t245 = t309 * t544 - t530;
t299 = sin(qJ(5));
t306 = cos(qJ(5));
t131 = t199 * t299 + t245 * t306;
t297 = sin(qJ(7));
t304 = cos(qJ(7));
t132 = t199 * t306 - t245 * t299;
t547 = t302 * t307;
t198 = t246 * t300 - t303 * t547;
t298 = sin(qJ(6));
t305 = cos(qJ(6));
t87 = t132 * t305 + t198 * t298;
t651 = -t131 * t297 + t304 * t87;
t521 = qJD(1) * t309;
t480 = t301 * t521;
t515 = qJD(4) * t306;
t538 = t306 * t308;
t525 = (-qJD(5) - t499) * t538 + (t300 * t515 + (qJD(5) * t307 + qJD(3)) * t299) * t301 - t235 * t306 + t299 * t480;
t478 = t308 * t521;
t517 = qJD(3) * t308;
t491 = pkin(2) * t517;
t650 = -pkin(2) * t478 + t652 * pkin(3) - t491;
t482 = -pkin(2) * t307 - pkin(3);
t431 = t301 * t482;
t457 = t301 * t499;
t513 = qJD(4) * t308;
t518 = qJD(3) * t301;
t649 = pkin(3) * t518 + (t300 * t513 + t457) * pkin(2) - t431 * t521;
t522 = qJD(1) * t302;
t234 = t300 * t478 - t307 * t522;
t514 = qJD(4) * t307;
t360 = -t300 * t517 - t301 * t514;
t617 = t234 - t360;
t292 = t301 * qJD(2);
t377 = -t308 * t522 + t292;
t439 = qJD(3) + t521;
t397 = t300 * t439;
t208 = -t307 * t377 - t397;
t278 = t301 * t522;
t293 = t308 * qJD(2);
t376 = t278 + t293;
t498 = qJD(4) - t376;
t647 = t208 * t498;
t506 = qJD(6) * t300;
t646 = -t301 * t506 + t525;
t539 = t306 * t307;
t380 = t299 * t308 + t301 * t539;
t645 = -qJD(6) * t380 + t617;
t644 = t132 * t298 - t198 * t305;
t493 = t302 * qJDD(1);
t346 = qJD(2) * (-qJD(3) + t521) + t493;
t496 = qJD(1) * qJD(3);
t454 = t302 * t496;
t402 = qJDD(2) + t454;
t161 = (t301 * t346 + t402 * t308) * pkin(2);
t259 = t482 * t308;
t607 = pkin(3) * t307;
t261 = (pkin(2) + t607) * t301;
t188 = t259 * t299 - t261 * t306;
t639 = qJD(5) * t188 + t650 * t299 - t649 * t306;
t206 = -t300 * t377 + t307 * t439;
t541 = t305 * t306;
t109 = -t206 * t541 + t208 * t298;
t509 = qJD(5) * t306;
t461 = t305 * t509;
t643 = -t109 + t461;
t474 = t302 * t517;
t362 = t309 * t292 + t474;
t291 = t309 * qJDD(1);
t497 = qJD(1) * qJD(2);
t256 = t302 * t497 - qJDD(3) - t291;
t323 = -t402 * t301 + t346 * t308;
t80 = -qJD(4) * t397 + t300 * t323 + t307 * (-qJD(4) * t377 - t256);
t642 = -qJDD(5) - t80;
t236 = t306 * t498;
t144 = t208 * t299 - t236;
t138 = qJD(6) + t144;
t146 = t306 * t208 + t299 * t498;
t334 = qJD(5) + t206;
t95 = t305 * t146 + t298 * t334;
t42 = t138 * t304 + t297 * t95;
t185 = t305 * t334;
t93 = t146 * t298 - t185;
t529 = qJD(7) - t93;
t641 = t42 * t529;
t44 = -t297 * t138 + t304 * t95;
t640 = t44 * t529;
t582 = -t645 * t298 + t646 * t305;
t618 = qJD(3) * t498;
t636 = t307 * t618;
t507 = qJD(6) * t298;
t460 = t299 * t507;
t635 = -t460 + t643;
t177 = t299 * t377 - t376 * t539;
t510 = qJD(5) * t305;
t464 = t299 * t510;
t503 = qJD(6) * t306;
t540 = t305 * t307;
t564 = t298 * t300;
t527 = -t177 * t305 + t376 * t564 + (qJD(6) - t515) * t540 + (t464 + (-qJD(4) + t503) * t298) * t300;
t534 = t307 * t308;
t551 = t301 * t306;
t524 = t551 * qJD(3) + t380 * qJD(5) + t306 * t480 + (t534 * qJD(3) + t235 - t471) * t299;
t560 = t299 * t307;
t176 = -t306 * t377 - t376 * t560;
t463 = t300 * t509;
t472 = t299 * t514;
t351 = t176 + t463 + t472;
t255 = t376 * pkin(2);
t468 = t255 * t514;
t559 = t300 * t161;
t621 = t468 + t559;
t634 = qJD(2) * qJD(3) - t493;
t244 = -t300 * t309 + t302 * t534;
t552 = t301 * t302;
t488 = t299 * t552;
t197 = t244 * t306 - t488;
t533 = t307 * t309;
t546 = t302 * t308;
t243 = t300 * t546 + t533;
t129 = t197 * t298 - t243 * t305;
t455 = t309 * t497;
t494 = t301 * qJDD(2);
t609 = pkin(2) * t308;
t162 = qJD(2) * t491 + (-t455 - t493) * t609 + (t301 * t454 + t494) * pkin(2);
t445 = t300 * t256 + t307 * t323;
t79 = -qJD(4) * t206 + t445;
t50 = -pkin(3) * t79 + t162;
t354 = -t308 * qJDD(2) + t301 * t634;
t170 = -qJD(1) * t362 + t354;
t169 = -qJDD(4) - t170;
t537 = t307 * t161;
t59 = -pkin(3) * t169 + t255 * t516 - t537;
t404 = t299 * t50 + t306 * t59;
t535 = t307 * t255;
t178 = pkin(3) * t498 - t535;
t254 = t377 * pkin(2);
t338 = pkin(3) * t208 - t254;
t91 = t299 * t178 + t306 * t338;
t20 = -qJD(5) * t91 + t404;
t92 = t306 * t178 - t299 * t338;
t383 = t255 * t564 - t305 * t92;
t555 = t300 * t305;
t8 = qJD(6) * t383 - t161 * t555 - t20 * t298 - t305 * t468;
t253 = t309 * t530 - t544;
t545 = t302 * t310;
t205 = t253 * t307 + t300 * t545;
t252 = -t303 * t308 - t309 * t548;
t137 = t205 * t306 + t252 * t299;
t204 = t253 * t300 - t307 * t545;
t89 = -t137 * t298 + t204 * t305;
t364 = -g(1) * t89 + g(2) * t644 + g(3) * t129 + t8;
t631 = t529 ^ 2;
t632 = pkin(4) * t631 + t364;
t189 = t259 * t306 + t261 * t299;
t580 = qJD(5) * t189 + t649 * t299 + t650 * t306;
t630 = t144 * t334;
t629 = t206 * t439;
t610 = pkin(2) * t300;
t441 = t518 * t610;
t626 = t189 * t507 - t298 * t441 + t639 * t305;
t157 = (-qJD(4) + t293) * t533 + (-t457 + (qJD(2) - t513) * t300) * t302;
t624 = -qJD(5) * t552 + t157;
t330 = qJD(4) * t334;
t327 = t307 * t330;
t623 = t300 * t642 - t327;
t294 = t300 ^ 2;
t622 = t161 * t294 + 0.2e1 * t300 * t468;
t550 = t301 * t307;
t242 = t299 * t550 - t538;
t620 = -t300 * t376 + t516;
t619 = qJD(3) * t439;
t314 = t306 * t642;
t329 = qJD(5) * t334;
t563 = t299 * t300;
t616 = -t300 * t314 + t306 * t327 - t329 * t563;
t469 = t302 * t514;
t435 = t301 * t469;
t615 = -pkin(2) * t435 - t362 * t610;
t436 = t300 * t480;
t411 = pkin(2) * t436;
t467 = t307 * t513;
t614 = -pkin(2) * t467 + t411 + t441;
t369 = g(1) * t204 + g(2) * t198 + g(3) * t243;
t511 = qJD(5) * t299;
t31 = qJD(5) * t236 - t299 * t169 - t208 * t511 + t306 * t79;
t16 = qJD(6) * t95 + t298 * t31 + t305 * t642;
t613 = t255 * t377;
t428 = g(1) * t310 + g(2) * t303;
t347 = -g(3) * t309 + t428 * t302;
t15 = qJD(6) * t185 - t146 * t507 - t298 * t642 + t305 * t31;
t446 = t306 * t169 + t299 * t79;
t32 = qJD(5) * t146 + t446;
t30 = qJDD(6) + t32;
t502 = qJD(7) * t297;
t4 = -t138 * t502 + t15 * t297 + t304 * (qJD(7) * t95 + t30);
t37 = pkin(4) * t95 + t91;
t38 = -pkin(4) * t138 - t383;
t407 = t297 * t38 + t304 * t37;
t504 = qJD(6) * t305;
t458 = t300 * t504;
t7 = t305 * t20 - t255 * t458 - t298 * t621 - t507 * t92;
t5 = -pkin(4) * t30 + t7;
t451 = t299 * t59 - t306 * t50;
t21 = qJD(5) * t92 + t451;
t6 = t15 * pkin(4) + t21;
t1 = -t407 * qJD(7) - t297 * t6 + t304 * t5;
t612 = t468 + t369;
t611 = g(2) * t87;
t608 = pkin(3) * t306;
t605 = g(2) * t131;
t3 = -qJD(7) * t42 + t304 * t15 - t297 * t30;
t601 = t297 * t3;
t577 = t144 * t305;
t600 = (-pkin(4) * t577 + t92) * t529;
t557 = t300 * t301;
t195 = -t298 * t557 - t305 * t380;
t127 = t195 * t297 - t242 * t304;
t599 = -qJD(7) * t127 + t524 * t297 + t582 * t304;
t598 = t242 * t502 + (qJD(7) * t195 - t524) * t304 + t582 * t297;
t597 = -t189 * t504 + (-t305 * t467 + (t305 * t518 + t308 * t507) * t300) * pkin(2) + t305 * t411 + t639 * t298;
t248 = t298 * t307 - t300 * t541;
t596 = (qJD(7) * t563 + t527) * t304 + (-qJD(7) * t248 + t351) * t297;
t194 = t248 * t304 + t297 * t563;
t595 = qJD(7) * t194 + t527 * t297 - t351 * t304;
t594 = t138 * t93;
t14 = -qJDD(7) + t16;
t593 = t14 * t297;
t592 = t14 * t304;
t591 = t15 * t298;
t590 = t298 * t30;
t589 = t299 * t31;
t588 = t30 * t305;
t587 = t300 * t79;
t586 = t304 * t44;
t585 = t95 * t138;
t505 = qJD(6) * t304;
t542 = t304 * t306;
t566 = t297 * t299;
t584 = -t109 * t304 - t206 * t566 + (qJD(7) + t510) * t542 + (-t298 * t505 + (-qJD(7) * t305 - qJD(5)) * t297) * t299;
t543 = t304 * t305;
t247 = t297 * t306 + t299 * t543;
t562 = t299 * t304;
t583 = t247 * qJD(7) + t334 * t562 + (t541 * qJD(5) - t109 - t460) * t297;
t581 = t646 * t298 + t645 * t305;
t181 = pkin(3) * t377 + t254 * t307;
t186 = t376 * t607 + t255;
t100 = t181 * t306 + t186 * t299;
t579 = t305 * t100 + t254 * t564;
t578 = qJD(5) * t95;
t576 = t146 * t306;
t574 = t206 * t299;
t573 = t206 * t306;
t572 = t208 * t307;
t569 = t254 * t300;
t568 = t255 * t300;
t565 = t298 * t299;
t561 = t299 * t305;
t558 = t300 * t169;
t554 = t300 * t306;
t553 = t301 * t256;
t549 = t301 * t309;
t536 = t307 * t169;
t532 = t308 * t256;
t489 = t255 * t554;
t125 = pkin(3) * t574 + t489;
t528 = t305 * t125 - t298 * t535;
t357 = t298 * t514 + t458;
t465 = t300 * t511;
t526 = t376 * t555 - t305 * t516 + t306 * t357 - t307 * t507 + (t177 - t465) * t298;
t295 = t302 ^ 2;
t523 = -t309 ^ 2 + t295;
t520 = qJD(2) * t302;
t519 = qJD(2) * t309;
t512 = qJD(5) * t298;
t508 = qJD(6) * t297;
t501 = qJD(7) * t304;
t500 = t376 * qJD(2);
t492 = t300 * t609;
t490 = t255 * t563;
t486 = t300 * t552;
t485 = t301 * t555;
t484 = t302 * t551;
t483 = t301 * t545;
t481 = pkin(3) * t305 + pkin(4);
t475 = t302 * t518;
t466 = qJD(5) * t255 * t294;
t453 = -qJD(7) * t38 - t6;
t444 = t138 * t305;
t443 = pkin(2) * t486;
t438 = t298 * t486;
t437 = t302 * t485;
t434 = t299 * t458;
t399 = t305 * t189 - t298 * t492;
t112 = pkin(4) * t242 + t399;
t433 = t582 * pkin(4) + qJD(7) * t112 + t580;
t113 = pkin(4) * t195 + t188;
t432 = t524 * pkin(4) - qJD(7) * t113 + t298 * t411 - t357 * t609 - t626;
t427 = g(1) * t303 - g(2) * t310;
t426 = t138 * t510 + t15;
t425 = t138 * t512 + t16;
t210 = -pkin(2) * t546 - pkin(3) * t244;
t238 = t302 * t431;
t142 = t210 * t299 + t238 * t306;
t363 = t293 * t309 - t475;
t107 = -pkin(2) * t363 - pkin(3) * t157;
t140 = -t362 * pkin(3) + (-t302 * t473 + (t302 * t516 - t307 * t519) * t301) * pkin(2);
t401 = t210 * t306 - t238 * t299;
t35 = qJD(5) * t401 + t107 * t299 + t140 * t306;
t424 = -t142 * t507 + t305 * t35;
t124 = pkin(3) * t573 - t490;
t237 = t481 * t563;
t359 = -t306 * t514 + t465;
t99 = t181 * t299 - t186 * t306;
t422 = pkin(3) * t359 + t527 * pkin(4) + qJD(7) * t237 - t99;
t209 = -pkin(3) * t554 + pkin(4) * t248;
t421 = -qJD(7) * t209 + (t300 * t461 + (-t298 * t506 + t305 * t514) * t299) * pkin(3) - t579 + t351 * pkin(4);
t260 = (pkin(4) * t305 + pkin(3)) * t299;
t420 = -qJD(7) * t260 + (-t298 * t503 - t464) * pkin(3) - t528 + (-t574 - t511) * pkin(4);
t258 = t481 * t306;
t419 = pkin(3) * t509 + t635 * pkin(4) + qJD(7) * t258 + t124;
t418 = t309 * t498;
t108 = -t305 * t208 - t298 * t573;
t415 = -t298 * t509 + t108;
t413 = qJD(4) * t498;
t130 = t197 * t305 + t243 * t298;
t64 = pkin(4) * t130 - t401;
t126 = t305 * t142;
t196 = t244 * t299 + t484;
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
t390 = qJD(3) * t376 + t494;
t388 = t498 * t516;
t387 = t307 * t413;
t386 = t300 * t618;
t385 = t501 * t529 - t593;
t384 = t502 * t529 + t592;
t382 = -t138 * t504 - t590;
t381 = -t138 * t507 + t588;
t379 = t302 * t377;
t378 = t377 * t308;
t375 = t309 * t393;
t374 = t309 * t394;
t372 = t308 * t493 - t494;
t136 = t205 * t299 - t252 * t306;
t371 = g(1) * t136 + g(3) * t196 + t605;
t370 = g(1) * t137 + g(2) * t132 + g(3) * t197;
t368 = g(1) * t205 + g(2) * t199 + g(3) * t244;
t367 = pkin(4) * t14;
t366 = t398 - t162;
t365 = t309 * t378;
t355 = (-pkin(2) * t169 - qJD(3) * t254 - t161) * t308;
t353 = -t21 + t371;
t350 = t299 * t504 - t415;
t349 = qJD(5) * t93 + t382;
t348 = -g(1) * t252 + g(2) * t245 + g(3) * t552 + t162;
t345 = t387 - t558;
t344 = (-pkin(4) * t146 - t305 * t91) * t529 + t370;
t343 = -t308 * t619 + t553;
t342 = qJD(5) * t244 + t362;
t340 = -pkin(4) * qJD(6) * t529 + t371;
t337 = -t376 * t498 + t413;
t336 = qJD(1) * t418 + t618;
t36 = qJD(5) * t142 - t107 * t306 + t140 * t299;
t335 = t619 + t374;
t333 = pkin(3) * t334;
t331 = t300 * t334;
t328 = qJD(5) * t333;
t326 = t306 * t329;
t320 = t206 * t334 + t329;
t317 = pkin(3) * t642;
t316 = t299 * t642;
t312 = qJD(1) ^ 2;
t311 = qJD(2) ^ 2;
t250 = t300 * t531 - t547;
t241 = t298 * t554 + t540;
t240 = t297 * t561 - t542;
t224 = t244 * t310;
t223 = t243 * t310;
t222 = t244 * t303;
t221 = t243 * t303;
t220 = t380 * t302;
t219 = t242 * t302;
t203 = t251 * t306 - t299 * t549;
t202 = t251 * t299 + t306 * t549;
t193 = -t298 * t380 + t485;
t192 = t248 * t297 - t300 * t562;
t175 = -t224 * t306 + t299 * t483;
t174 = -t224 * t299 - t306 * t483;
t173 = -t222 * t306 + t303 * t488;
t172 = -t222 * t299 - t303 * t484;
t171 = -t220 * t305 - t438;
t167 = t252 * t539 - t253 * t299;
t166 = t252 * t560 + t253 * t306;
t165 = -t245 * t539 - t246 * t299;
t164 = -t245 * t560 + t246 * t306;
t163 = -t243 * t541 + t244 * t298;
t158 = -t189 * t298 - t305 * t492;
t156 = -t300 * t475 - t309 * t516 - t307 * t520 + (t300 * t519 + t469) * t308;
t135 = t203 * t305 + t250 * t298;
t128 = t195 * t304 + t242 * t297;
t111 = t175 * t305 - t223 * t298;
t110 = t173 * t305 - t221 * t298;
t105 = -pkin(2) * t437 - t142 * t298;
t104 = t167 * t305 + t252 * t564;
t103 = t165 * t305 - t245 * t564;
t102 = -t204 * t541 + t205 * t298;
t101 = -t198 * t541 + t199 * t298;
t96 = -t125 * t298 - t305 * t535;
t90 = t137 * t305 + t204 * t298;
t85 = t130 * t304 - t196 * t297;
t73 = -t100 * t298 + t254 * t555;
t70 = -t342 * t299 + t306 * t624;
t69 = t299 * t624 + t306 * t342;
t57 = -t144 * t543 - t146 * t297;
t56 = t146 * t304 - t297 * t577;
t53 = -t255 * t555 - t298 * t92;
t41 = -t136 * t297 + t304 * t90;
t40 = -t136 * t304 - t297 * t90;
t33 = t37 * t502;
t27 = -qJD(6) * t129 + t156 * t298 + t305 * t70;
t26 = qJD(6) * t130 - t156 * t305 + t298 * t70;
t22 = -t142 * t504 - t298 * t35 + (-t305 * t435 + (-t305 * t474 + (t302 * t507 - t305 * t519) * t301) * t300) * pkin(2);
t18 = -t297 * t37 + t304 * t38;
t12 = pkin(4) * t27 + t36;
t11 = -pkin(4) * t69 + (-t403 + (-t300 * t362 - t435) * t298) * pkin(2) + t424;
t10 = -qJD(7) * t84 + t27 * t304 - t297 * t69;
t9 = t27 * t297 - t196 * t502 + (qJD(7) * t130 + t69) * t304;
t2 = -t297 * t5 + t304 * t453 + t33;
t13 = [qJDD(1), t427, t428, qJDD(1) * t295 + 0.2e1 * t302 * t455, 0.2e1 * t291 * t302 - 0.2e1 * t497 * t523, qJDD(2) * t302 + t309 * t311, qJDD(2) * t309 - t302 * t311, 0, t427 * t309, -t427 * t302, -qJD(2) * t365 + t323 * t546 + t379 * t518 (-t309 * t500 + (qJD(3) * t377 + t170) * t302) * t308 + (t292 * t519 + ((-qJDD(1) * t308 + t301 * t496) * t302 + (qJD(3) - 0.2e1 * t521) * t293 + t390) * t302) * t301 (-qJD(3) ^ 2 * t301 + qJD(2) * t377 - t532) * t302 + (0.2e1 * qJD(1) * t363 + t372) * t309, -t301 * t375 + t170 * t309 + (t343 + t500) * t302, -t256 * t309 - t302 * t393, -t254 * t520 + g(1) * t246 - g(2) * t253 + t162 * t309 + (-t308 * t375 + (t301 * t619 + t532) * t302) * pkin(2), -t255 * t520 - g(1) * t245 - g(2) * t252 + t161 * t309 + (t619 * t546 + (qJD(2) * t396 - t302 * t256) * t301) * pkin(2), t157 * t208 + t244 * t79, -t156 * t208 - t157 * t206 - t243 * t79 - t244 * t80, t157 * t498 - t244 * t169 - t208 * t362 - t552 * t79, -t156 * t498 + t243 * t169 + t206 * t362 + t552 * t80, -t618 * t546 + (-qJD(2) * t418 + t169 * t302) * t301, g(1) * t199 - g(2) * t205 + t254 * t156 + t162 * t243 + (-t255 * t557 + (-t308 * t206 + t498 * t557) * pkin(2)) * t519 + (-t161 * t557 + t360 * t255 + ((t386 - t80) * t308 + (t206 * qJD(3) + t345) * t301) * pkin(2)) * t302, -g(1) * t198 + g(2) * t204 + t254 * t157 + t162 * t244 + (-t301 * t535 + (-t308 * t208 + t498 * t550) * pkin(2)) * t519 + (-t301 * t537 - t361 * t255 + ((-t79 + t636) * t308 + (t208 * qJD(3) - t388 - t536) * t301) * pkin(2)) * t302, t146 * t70 + t197 * t31, -t144 * t70 - t146 * t69 - t196 * t31 - t197 * t32, t146 * t156 - t197 * t642 + t31 * t243 + t334 * t70, -t144 * t156 + t196 * t642 - t32 * t243 - t334 * t69, t156 * t334 - t243 * t642, g(1) * t132 - g(2) * t137 + t144 * t615 - t91 * t156 - t196 * t621 - t21 * t243 - t32 * t443 - t334 * t36 - t401 * t642 - t568 * t69, -g(1) * t131 + g(2) * t136 + t142 * t642 + t146 * t615 - t92 * t156 - t197 * t621 - t20 * t243 - t31 * t443 - t334 * t35 - t568 * t70, t130 * t15 + t27 * t95, -t129 * t15 - t130 * t16 - t26 * t95 - t27 * t93, t130 * t30 + t138 * t27 + t15 * t196 + t69 * t95, -t129 * t30 - t138 * t26 - t16 * t196 - t69 * t93, t138 * t69 + t196 * t30, g(1) * t87 - g(2) * t90 + t105 * t30 + t129 * t21 + t138 * t22 - t16 * t401 + t196 * t8 + t26 * t91 + t36 * t93 + t53 * t69, -t424 * t138 - t126 * t30 - t7 * t196 + t383 * t69 + t36 * t95 - t401 * t15 + t21 * t130 + t91 * t27 - g(1) * t644 - g(2) * t89 + (t138 * t403 + (t138 * t435 + (t138 * t362 + t30 * t552) * t300) * t298) * pkin(2), t10 * t44 + t3 * t85, -t10 * t42 - t3 * t84 - t4 * t85 - t44 * t9, t10 * t529 - t129 * t3 - t14 * t85 - t26 * t44, t129 * t4 + t14 * t84 + t26 * t42 - t529 * t9, t129 * t14 - t26 * t529 (qJD(7) * t405 - t11 * t297 - t12 * t304) * t529 - t406 * t14 - t2 * t129 + t407 * t26 + t22 * t42 + t105 * t4 + t8 * t84 + t53 * t9 + g(1) * t651 - g(2) * t41 -(qJD(7) * t406 + t11 * t304 - t12 * t297) * t529 - t405 * t14 + t1 * t129 + t18 * t26 + t22 * t44 + t105 * t3 + t8 * t85 + t53 * t10 - g(1) * (t131 * t304 + t297 * t87) - g(2) * t40; 0, 0, 0, -t302 * t312 * t309, t523 * t312, t493, t291, qJDD(2), t347, g(3) * t302 + t309 * t428 -(t308 * t455 + t372) * t301 + qJD(1) * t365 + (t301 * t376 + t378) * qJD(3) (-t292 * t439 - t170) * t301 + (t634 * t308 + (0.2e1 * t475 + (t376 + t278 - t293) * t309) * qJD(1) + t390) * t308, -qJD(1) * t379 - t308 * t374 + t343, t301 * t335 - t376 * t522 + t532, t302 * t394, -pkin(2) * t553 + t254 * t522 + (pkin(2) * t335 + t347) * t308, -pkin(2) * t532 + t255 * t522 + ((-qJD(1) * t396 - t619) * pkin(2) - t347) * t301, t652 * t208 - t79 * t550, t206 * t235 + t208 * t234 + (t206 * t307 + t208 * t300) * t517 + (t587 + t307 * t80 + (-t206 * t300 + t572) * qJD(4)) * t301, -t235 * t498 + (-t79 - t636) * t308 + (t300 * t413 + t395 + t536) * t301, t234 * t498 + (t386 + t80) * t308 + (t345 - t629) * t301, t169 * t308 + t301 * t336, g(1) * t224 + g(2) * t222 - g(3) * t251 - t254 * t234 - t400 * t514 + (t301 * t80 + (t387 + t629) * t308) * pkin(2) + (t355 + (-pkin(2) * t336 + t366) * t301) * t300, -g(1) * t223 - g(2) * t221 + g(3) * t250 - t254 * t235 + t400 * t516 + (t301 * t79 + (-t388 + t395) * t308) * pkin(2) + (t355 + ((-t498 * t521 - t618) * pkin(2) + t366) * t301) * t307, t146 * t525 - t31 * t380, -t144 * t525 + t146 * t524 + t242 * t31 + t32 * t380, -t146 * t617 - t31 * t557 + t334 * t525 + t380 * t642, t144 * t617 - t242 * t642 + t32 * t557 + t334 * t524, -t234 * t334 + t301 * t623 - t331 * t517, -g(1) * t175 - g(2) * t173 - g(3) * t203 + t144 * t614 + t188 * t642 + t21 * t557 + t242 * t621 - t32 * t492 - t334 * t580 + t524 * t568 + t617 * t91, g(1) * t174 + g(2) * t172 + g(3) * t202 + t614 * t146 + t189 * t642 + t20 * t557 - t31 * t492 + t334 * t639 + t621 * t380 - t525 * t568 + t617 * t92, t15 * t195 + t582 * t95, -t15 * t193 - t16 * t195 - t581 * t95 - t582 * t93, t138 * t582 - t15 * t242 + t195 * t30 - t524 * t95, -t138 * t581 + t16 * t242 - t193 * t30 + t524 * t93, -t138 * t524 - t242 * t30, -g(1) * t111 - g(2) * t110 - g(3) * t135 + t138 * t597 + t158 * t30 + t16 * t188 + t193 * t21 - t242 * t8 - t524 * t53 + t580 * t93 + t581 * t91, -t399 * t30 + t7 * t242 + t188 * t15 + t21 * t195 - g(1) * (-t175 * t298 - t223 * t305) - g(2) * (-t173 * t298 - t221 * t305) - g(3) * (-t203 * t298 + t250 * t305) + t580 * t95 + t582 * t91 - t524 * t383 + ((t308 * t458 + (-t436 + t467) * t298) * pkin(2) + t626) * t138, t128 * t3 + t44 * t599, -t127 * t3 - t128 * t4 - t42 * t599 - t44 * t598, -t128 * t14 - t193 * t3 - t44 * t581 + t529 * t599, t127 * t14 + t193 * t4 + t42 * t581 - t529 * t598, t14 * t193 - t529 * t581 -(-t112 * t297 - t113 * t304) * t14 - t2 * t193 + t158 * t4 + t8 * t127 - g(1) * (t111 * t304 - t174 * t297) - g(2) * (t110 * t304 - t172 * t297) - g(3) * (t135 * t304 - t202 * t297) - (t297 * t432 + t304 * t433) * t529 + t598 * t53 + t597 * t42 + t581 * t407 (t112 * t304 - t113 * t297) * t14 + t1 * t193 + t158 * t3 + t8 * t128 - g(1) * (-t111 * t297 - t174 * t304) - g(2) * (-t110 * t297 - t172 * t304) - g(3) * (-t135 * t297 - t202 * t304) - (-t297 * t433 + t304 * t432) * t529 + t599 * t53 + t597 * t44 + t581 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377 * t376, -t376 ^ 2 + t377 ^ 2, t376 * t439 + t323, -t292 * qJD(3) + (-t292 - t377) * t521 + t354, -t256, -t398 + t348, g(1) * t253 + g(2) * t246 + g(3) * t546 + t254 * t439 + t161, -t498 * t572 - t587 (t206 * t498 - t79) * t307 + (t80 + t647) * t300, -t208 * t377 - t307 * t337 + t558, t206 * t377 + t300 * t337 + t536, -t498 * t377, -t255 * t206 - t300 * t613 + t348 * t307, -t255 * t208 - t348 * t300 - t307 * t613, -t31 * t554 + (-t177 + t359) * t146, t144 * t177 + t146 * t176 + (t144 * t306 + t146 * t299) * t514 + (t589 + t306 * t32 + (-t144 * t299 + t576) * qJD(5)) * t300, -t146 * t620 - t177 * t334 + t31 * t307 - t616, t144 * t620 + t176 * t334 - t299 * t623 + t300 * t326 - t32 * t307, -t300 * t330 - t307 * t642 + t331 * t376, pkin(3) * t616 - g(1) * t167 - g(2) * t165 + g(3) * t220 - t144 * t569 + t176 * t568 - t21 * t307 + t299 * t622 + t306 * t466 + t334 * t99 + t620 * t91, g(1) * t166 + g(2) * t164 - g(3) * t219 + t100 * t334 - t146 * t569 + t177 * t568 - t20 * t307 - t299 * t466 + t306 * t622 + t317 * t563 - t328 * t554 - t333 * t472 + t620 * t92, t15 * t248 + t527 * t95, t15 * t241 - t16 * t248 + t526 * t95 - t527 * t93, t138 * t527 - t15 * t563 + t248 * t30 - t351 * t95, t138 * t526 + t16 * t563 + t241 * t30 + t351 * t93, -t138 * t351 - t30 * t563, -t8 * t563 - g(1) * t104 - g(2) * t103 - g(3) * t171 - t138 * t73 - t21 * t241 - t93 * t99 - t526 * t91 - t351 * t53 + ((-t138 * t565 - t306 * t93) * t514 + (t299 * t349 - t306 * t425) * t300) * pkin(3), t7 * t563 + t21 * t248 + t579 * t138 - t99 * t95 - g(1) * (-t167 * t298 + t252 * t555) - g(2) * (-t165 * t298 - t245 * t555) - g(3) * (t220 * t298 - t437) + t527 * t91 - t351 * t383 + ((-t138 * t561 - t306 * t95) * t514 + (-t426 * t306 + (-t381 + t578) * t299) * t300) * pkin(3), t194 * t3 + t44 * t596, -t192 * t3 - t194 * t4 - t42 * t596 - t44 * t595, -t14 * t194 + t241 * t3 + t44 * t526 + t529 * t596, t14 * t192 - t241 * t4 - t42 * t526 - t529 * t595, -t14 * t241 + t526 * t529 -(-t209 * t304 - t237 * t297) * t14 + t2 * t241 + t8 * t192 - t73 * t42 - g(1) * (t104 * t304 - t166 * t297) - g(2) * (t103 * t304 - t164 * t297) - g(3) * (t171 * t304 + t219 * t297) - (t297 * t421 + t304 * t422) * t529 + t595 * t53 - t526 * t407 + (-t42 * t434 + (-t42 * t463 + (-t300 * t4 - t42 * t514) * t299) * t298) * pkin(3) (-t209 * t297 + t237 * t304) * t14 - t1 * t241 + t8 * t194 - t73 * t44 - g(1) * (-t104 * t297 - t166 * t304) - g(2) * (-t103 * t297 - t164 * t304) - g(3) * (-t171 * t297 + t219 * t304) - (-t297 * t422 + t304 * t421) * t529 + t596 * t53 - t526 * t18 + (-t44 * t434 + (-t44 * t463 + (-t3 * t300 - t44 * t514) * t299) * t298) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208 * t206, -t206 ^ 2 + t208 ^ 2, -t206 * t376 + t445, -t80 + t647, -t169, -t254 * t208 + t376 * t535 + t369 + t559, t254 * t206 - t376 * t568 + t368 + t537, t334 * t576 + t589 (t31 - t630) * t306 + (-t146 * t334 - t32) * t299, -t146 * t208 + t306 * t320 - t316, t144 * t208 - t299 * t320 - t314, -t334 * t208, -t124 * t334 + t144 * t535 + t161 * t554 - t206 * t490 + t91 * t208 - t255 * t465 + (t316 - t326) * pkin(3) + t612 * t306, t125 * t334 + t146 * t535 - t206 * t489 + t92 * t208 - t255 * t463 + t306 * t317 + (t328 - t559 - t612) * t299, t15 * t561 + t635 * t95, t108 * t95 + t109 * t93 + (-t298 * t95 - t305 * t93) * t509 + (-t591 - t16 * t305 + (t298 * t93 - t305 * t95) * qJD(6)) * t299, -t15 * t306 + t643 * t138 + (t334 * t95 + t381) * t299, t16 * t306 + t415 * t138 + (-t334 * t93 + t382) * t299, t138 * t299 * t334 - t30 * t306, -g(1) * t102 - g(2) * t101 - g(3) * t163 - t108 * t91 + t124 * t93 - t138 * t96 + (pkin(3) * t349 + t512 * t91 - t8) * t306 + (pkin(3) * t425 + t21 * t298 + t334 * t53 + t504 * t91) * t299, t528 * t138 + t124 * t95 - t91 * t109 - t368 * t305 + (t91 * t510 + t7 + (t578 - t588) * pkin(3) + (pkin(3) * qJD(6) * t138 - t369) * t298) * t306 + (pkin(3) * t426 + t21 * t305 + t334 * t383 - t507 * t91) * t299, t247 * t3 + t44 * t584, -t240 * t3 - t247 * t4 - t42 * t584 - t44 * t583, -t14 * t247 - t3 * t565 - t350 * t44 + t529 * t584, t14 * t240 + t350 * t42 + t4 * t565 - t529 * t583, t14 * t565 - t350 * t529 -(-t258 * t297 - t260 * t304) * t14 + t8 * t240 - t407 * t108 - t96 * t42 - g(1) * (t102 * t304 + t204 * t566) - g(2) * (t101 * t304 + t198 * t566) - g(3) * (t163 * t304 + t243 * t566) - (t297 * t420 + t304 * t419) * t529 + t583 * t53 + (t299 * t407 - t42 * t608) * t504 + (t407 * t509 - t2 * t299 + (-t306 * t4 + t42 * t511) * pkin(3)) * t298 (t258 * t304 - t260 * t297) * t14 + t8 * t247 - t18 * t108 - t96 * t44 - g(1) * (-t102 * t297 + t204 * t562) - g(2) * (-t101 * t297 + t198 * t562) - g(3) * (-t163 * t297 + t243 * t562) - (-t297 * t419 + t304 * t420) * t529 + t584 * t53 + (t18 * t299 - t44 * t608) * t504 + (t18 * t509 + t1 * t299 + (-t3 * t306 + t44 * t511) * pkin(3)) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146 * t144, -t144 ^ 2 + t146 ^ 2, t31 + t630, t146 * t206 - t446, -t642, t146 * t568 + t206 * t92 + t371 - t451, -t144 * t568 - t91 * t206 + t370 - t404, t444 * t95 + t591 (t15 - t594) * t305 + (-t16 - t585) * t298, t138 * t444 - t146 * t95 + t590, -t138 ^ 2 * t298 + t146 * t93 + t588, -t138 * t146, -t146 * t53 + t305 * t353 - t92 * t93, -t146 * t383 - t298 * t353 - t92 * t95, t298 * t3 * t304 + (-t298 * t502 + t304 * t504 - t57) * t44, t42 * t57 + t44 * t56 + (-t297 * t44 - t304 * t42) * t504 + (-t601 - t304 * t4 + (t297 * t42 - t586) * qJD(7)) * t298, -t57 * t529 + (t505 * t529 + t3) * t305 + (-t138 * t44 - t384) * t298, t56 * t529 + (-t508 * t529 - t4) * t305 + (t138 * t42 - t385) * t298, -t138 * t298 * t529 - t14 * t305, t304 * t600 - t53 * t56 + t344 * t297 + (t304 * t340 + t508 * t53 + t2) * t305 + (pkin(4) * t384 + t138 * t407 + t297 * t8 - t91 * t42 + t501 * t53) * t298, -t297 * t600 - t53 * t57 + t344 * t304 + (-t297 * t340 + t505 * t53 - t1) * t305 + (pkin(4) * t385 + t138 * t18 + t8 * t304 - t91 * t44 - t502 * t53) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t93 ^ 2 + t95 ^ 2, t15 + t594, -t16 + t585, t30, -t138 * t383 - t91 * t95 + t364, g(1) * t90 + g(3) * t130 + t138 * t53 + t91 * t93 + t611 - t7, -t529 * t586 - t601 (-t3 + t641) * t304 + (t4 + t640) * t297, -t304 * t631 + t44 * t95 + t593, t297 * t631 - t42 * t95 + t592, t529 * t95, -t367 * t297 + t304 * t632 - t383 * t42 - t407 * t95, -t18 * t95 - t297 * t632 - t367 * t304 - t383 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t42 ^ 2 + t44 ^ 2, t3 + t641, -t4 + t640, -t14, -g(1) * t40 + g(3) * t84 + t18 * t529 - t53 * t44 + t33 + (-t5 + t611) * t297 + (t453 + t605) * t304, g(1) * t41 + g(2) * t651 + g(3) * t85 - t407 * t529 + t53 * t42 - t1;];
tau_reg  = t13;
