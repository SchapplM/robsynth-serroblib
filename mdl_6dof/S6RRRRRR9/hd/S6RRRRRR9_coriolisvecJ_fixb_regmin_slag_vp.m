% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x38]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:14
% EndTime: 2019-03-10 05:40:01
% DurationCPUTime: 19.80s
% Computational Cost: add. (31731->767), mult. (93132->1090), div. (0->0), fcn. (79137->14), ass. (0->351)
t368 = cos(pkin(7));
t377 = cos(qJ(3));
t378 = cos(qJ(2));
t517 = t377 * t378;
t372 = sin(qJ(3));
t373 = sin(qJ(2));
t522 = t372 * t373;
t405 = -t368 * t522 + t517;
t367 = sin(pkin(6));
t497 = qJD(1) * t367;
t292 = t405 * t497;
t366 = sin(pkin(7));
t494 = qJD(3) * t377;
t470 = t366 * t494;
t587 = -t292 + t470;
t538 = cos(pkin(6));
t476 = pkin(1) * t538;
t441 = t378 * t476;
t350 = qJD(1) * t441;
t431 = t367 * (-pkin(10) * t368 - pkin(9));
t408 = t373 * t431;
t271 = qJD(1) * t408 + t350;
t442 = t373 * t476;
t383 = t378 * t431 - t442;
t272 = t383 * qJD(1);
t552 = pkin(10) * t366;
t395 = (pkin(2) * t373 - t378 * t552) * t367;
t305 = qJD(1) * t395;
t527 = t368 * t372;
t530 = t366 * t372;
t353 = pkin(10) * t530;
t526 = t368 * t377;
t561 = pkin(2) * t526 - t353;
t586 = t561 * qJD(3) - t377 * t271 - t272 * t527 - t305 * t530;
t209 = -t272 * t366 + t368 * t305;
t520 = t373 * t377;
t521 = t372 * t378;
t407 = t368 * t520 + t521;
t291 = t407 * t497;
t585 = pkin(3) * t291 - pkin(11) * t292 + t209 - (pkin(3) * t372 - pkin(11) * t377) * t366 * qJD(3);
t474 = t373 * t497;
t440 = t366 * t474;
t584 = pkin(11) * t440 - t586;
t371 = sin(qJ(4));
t376 = cos(qJ(4));
t317 = -t376 * t368 + t371 * t530;
t506 = -qJD(4) * t317 - t371 * t440 + t376 * t587;
t318 = t368 * t371 + t376 * t530;
t505 = qJD(4) * t318 + t371 * t587 + t376 * t440;
t495 = qJD(3) * t372;
t471 = t366 * t495;
t577 = t291 - t471;
t529 = t366 * t377;
t500 = pkin(2) * t527 + pkin(10) * t529;
t583 = t500 * qJD(3) - t372 * t271;
t473 = t378 * t497;
t437 = t368 * t473;
t329 = t377 * t437;
t455 = t538 * qJD(1);
t412 = t455 + qJD(2);
t396 = t412 * t366;
t439 = t372 * t474;
t235 = -t377 * t396 - t329 + t439;
t233 = qJD(4) + t235;
t503 = t272 * t526 - (-pkin(3) * t474 - t305 * t377) * t366 + t583;
t309 = pkin(11) * t368 + t500;
t310 = (-pkin(3) * t377 - pkin(11) * t372 - pkin(2)) * t366;
t491 = qJD(4) * t376;
t493 = qJD(4) * t371;
t582 = t309 * t493 - t310 * t491 + t371 * t585 + t584 * t376;
t581 = -t309 * t491 - t310 * t493 + t584 * t371 - t376 * t585;
t573 = -t577 * pkin(12) - t582;
t528 = t367 * t378;
t391 = pkin(9) * t528 + t442;
t229 = t391 * qJD(1) + (t396 + t437) * pkin(10);
t386 = pkin(2) * t538 + t408;
t234 = qJD(2) * pkin(2) + qJD(1) * t386 + t350;
t301 = (-pkin(2) * t378 - t373 * t552 - pkin(1)) * t367;
t287 = qJD(1) * t301;
t150 = t377 * t229 + t234 * t527 + t287 * t530;
t580 = t150 - t233 * (pkin(4) * t371 - pkin(12) * t376);
t579 = -t505 * pkin(4) + t506 * pkin(12) - t503;
t406 = t368 * t521 + t520;
t393 = t406 * t367;
t237 = qJD(1) * t393 + t372 * t396;
t370 = sin(qJ(5));
t375 = cos(qJ(5));
t172 = -t235 * t370 * t376 - t375 * t237;
t570 = -t370 * t491 + t172;
t518 = t375 * t376;
t173 = -t235 * t518 + t237 * t370;
t420 = t375 * t491 - t173;
t490 = qJD(5) * t370;
t576 = -t371 * t490 + t420;
t340 = t366 * t473;
t486 = t340 - qJD(3);
t387 = -t368 * t412 + t486;
t195 = t376 * t237 - t371 * t387;
t156 = t195 * t370 - t375 * t233;
t186 = -t234 * t366 + t368 * t287;
t127 = pkin(3) * t235 - pkin(11) * t237 + t186;
t133 = -pkin(11) * t387 + t150;
t61 = t127 * t376 - t371 * t133;
t56 = -pkin(4) * t233 - t61;
t46 = pkin(5) * t156 + t56;
t158 = t195 * t375 + t233 * t370;
t369 = sin(qJ(6));
t374 = cos(qJ(6));
t83 = t374 * t156 + t158 * t369;
t575 = t46 * t83;
t417 = t156 * t369 - t374 * t158;
t574 = t417 * t83;
t489 = qJD(5) * t375;
t571 = t489 * t371 - t570;
t277 = t318 * t370 + t375 * t529;
t513 = -qJD(5) * t277 - t577 * t370 + t506 * t375;
t481 = t370 * t529;
t512 = -qJD(5) * t481 + t318 * t489 + t506 * t370 + t577 * t375;
t286 = t376 * t387;
t193 = t237 * t371 + t286;
t557 = qJD(5) + qJD(6);
t569 = t193 + t557;
t568 = t417 ^ 2 - t83 ^ 2;
t191 = qJD(5) + t193;
t62 = t371 * t127 + t376 * t133;
t57 = pkin(12) * t233 + t62;
t149 = -t372 * t229 + (t234 * t368 + t287 * t366) * t377;
t132 = pkin(3) * t387 - t149;
t68 = t193 * pkin(4) - t195 * pkin(12) + t132;
t31 = -t370 * t57 + t375 * t68;
t25 = -pkin(13) * t158 + t31;
t19 = pkin(5) * t191 + t25;
t32 = t370 * t68 + t375 * t57;
t26 = -pkin(13) * t156 + t32;
t548 = t26 * t374;
t11 = t19 * t369 + t548;
t392 = qJD(3) * t396;
t485 = qJD(1) * qJD(2);
t461 = t367 * t485;
t435 = t378 * t461;
t444 = -t368 * qJD(2) - qJD(3);
t196 = qJD(3) * t329 + t439 * t444 + (t392 + t435) * t377;
t436 = t373 * t461;
t413 = t366 * t436;
t120 = qJD(4) * t195 + t371 * t196 - t376 * t413;
t119 = -qJD(4) * t286 + t376 * t196 - t237 * t493 + t371 * t413;
t555 = qJD(2) * t407 + qJD(3) * t406;
t380 = t555 * t367;
t388 = t372 * t392;
t197 = qJD(1) * t380 + t388;
t54 = t375 * t119 - t195 * t490 + t370 * t197 + t233 * t489;
t274 = t383 * qJD(2);
t252 = qJD(1) * t274;
t306 = qJD(2) * t395;
t296 = qJD(1) * t306;
t198 = -t252 * t366 + t368 * t296;
t105 = pkin(3) * t197 - pkin(11) * t196 + t198;
t422 = qJD(2) * t455;
t410 = pkin(1) * t422;
t344 = t378 * t410;
t399 = qJD(2) * t408;
t251 = qJD(1) * t399 + t344;
t468 = t368 * t494;
t390 = t229 * t495 - t234 * t468 - t377 * t251 - t252 * t527 - t287 * t470 - t296 * t530;
t80 = pkin(11) * t413 - t390;
t22 = t371 * t105 + t127 * t491 - t133 * t493 + t376 * t80;
t17 = pkin(12) * t197 + t22;
t469 = t368 * t495;
t91 = -t229 * t494 - t234 * t469 - t372 * t251 + t252 * t526 - t287 * t471 + t296 * t529;
t81 = -pkin(3) * t413 - t91;
t39 = pkin(4) * t120 - pkin(12) * t119 + t81;
t8 = -qJD(5) * t32 - t17 * t370 + t375 * t39;
t4 = pkin(5) * t120 - pkin(13) * t54 + t8;
t55 = t158 * qJD(5) + t119 * t370 - t375 * t197;
t7 = t375 * t17 + t370 * t39 + t68 * t489 - t490 * t57;
t5 = -pkin(13) * t55 + t7;
t2 = -qJD(6) * t11 - t369 * t5 + t374 * t4;
t567 = t46 * t417 + t2;
t487 = qJD(6) * t374;
t488 = qJD(6) * t369;
t14 = -t156 * t487 - t158 * t488 - t369 * t55 + t374 * t54;
t190 = qJD(6) + t191;
t566 = t190 * t83 + t14;
t384 = qJD(6) * t417 - t369 * t54 - t374 * t55;
t565 = -t190 * t417 + t384;
t545 = t577 * pkin(4) - t581;
t564 = t233 * t371;
t333 = t369 * t375 + t370 * t374;
t314 = t333 * t371;
t483 = pkin(11) * t493;
t563 = t370 * t483 - t580 * t375;
t562 = t579 * t375;
t308 = t353 + (-pkin(2) * t377 - pkin(3)) * t368;
t216 = pkin(4) * t317 - pkin(12) * t318 + t308;
t504 = t376 * t309 + t371 * t310;
t218 = -pkin(12) * t529 + t504;
t509 = t370 * t216 + t375 * t218;
t560 = -t216 * t489 + t218 * t490 + t579 * t370 - t573 * t375;
t559 = t368 * t517 - t522;
t343 = -pkin(4) * t376 - pkin(12) * t371 - pkin(3);
t178 = pkin(3) * t237 + pkin(11) * t235;
t516 = t376 * t149 + t371 * t178;
t79 = pkin(12) * t237 + t516;
t558 = -t343 * t489 + t580 * t370 + t375 * t79;
t457 = t538 * t366;
t254 = (t368 * t528 + t457) * pkin(10) + t391;
t270 = t441 + t386;
t556 = -t372 * t254 + t377 * (t270 * t368 + t301 * t366);
t24 = t26 * t488;
t460 = qJD(6) * t19 + t5;
t1 = t369 * t4 + t460 * t374 - t24;
t379 = qJD(1) ^ 2;
t554 = pkin(12) + pkin(13);
t553 = pkin(5) * t371;
t261 = -t367 * t559 - t377 * t457;
t205 = -t270 * t366 + t368 * t301;
t262 = t372 * t457 + t393;
t145 = pkin(3) * t261 - pkin(11) * t262 + t205;
t456 = t538 * t368;
t316 = t366 * t528 - t456;
t479 = t377 * t254 + t270 * t527 + t301 * t530;
t153 = -pkin(11) * t316 + t479;
t515 = t371 * t145 + t376 * t153;
t71 = pkin(12) * t261 + t515;
t152 = pkin(3) * t316 - t556;
t214 = t262 * t371 + t316 * t376;
t215 = t262 * t376 - t316 * t371;
t92 = pkin(4) * t214 - pkin(12) * t215 + t152;
t551 = t370 * t92 + t375 * t71;
t23 = t376 * t105 - t127 * t493 - t133 * t491 - t371 * t80;
t18 = -pkin(4) * t197 - t23;
t550 = t18 * t370;
t549 = t18 * t375;
t547 = t370 * t54;
t546 = t512 * pkin(5) + t545;
t143 = t371 * t149;
t78 = -pkin(4) * t237 - t178 * t376 + t143;
t544 = t571 * pkin(5) + pkin(11) * t491 - t78;
t137 = pkin(4) * t195 + pkin(12) * t193;
t543 = t370 * t137 + t375 * t61;
t278 = t318 * t375 - t481;
t199 = t374 * t277 + t278 * t369;
t540 = -qJD(6) * t199 - t512 * t369 + t513 * t374;
t200 = -t277 * t369 + t278 * t374;
t539 = qJD(6) * t200 + t513 * t369 + t512 * t374;
t537 = t120 * t376;
t536 = t156 * t191;
t535 = t158 * t191;
t534 = t193 * t233;
t533 = t193 * t370;
t532 = t195 * t233;
t363 = t367 ^ 2;
t531 = t363 * t379;
t525 = t370 * t120;
t524 = t370 * t371;
t523 = t371 * t375;
t519 = t375 * t120;
t332 = t369 * t370 - t374 * t375;
t511 = t172 * t369 - t173 * t374 - t314 * t557 - t332 * t491;
t510 = -t488 * t524 + (t523 * t557 - t570) * t374 + t576 * t369;
t508 = t569 * t332;
t507 = t569 * t333;
t357 = pkin(11) * t518;
t499 = t370 * t343 + t357;
t498 = t373 ^ 2 - t378 ^ 2;
t496 = qJD(2) * t367;
t492 = qJD(4) * t375;
t482 = t378 * t531;
t475 = qJD(5) * t554;
t472 = t373 * t496;
t466 = t191 * t490;
t462 = t363 * t485;
t458 = -t370 * t71 + t375 * t92;
t453 = t145 * t376 - t371 * t153;
t451 = t375 * t216 - t218 * t370;
t210 = -t274 * t366 + t368 * t306;
t450 = -t371 * t309 + t310 * t376;
t447 = t233 * t376;
t446 = t191 * t375;
t445 = 0.2e1 * t462;
t351 = qJD(2) * t441;
t273 = t351 + t399;
t443 = -t254 * t494 - t270 * t469 - t372 * t273 - t301 * t471;
t438 = t366 * t472;
t434 = -t62 + (t490 + t533) * pkin(5);
t117 = pkin(5) * t317 - pkin(13) * t278 + t451;
t433 = t512 * pkin(13) - qJD(6) * t117 + t560;
t134 = -pkin(13) * t277 + t509;
t432 = -t505 * pkin(5) + t513 * pkin(13) + t509 * qJD(5) + qJD(6) * t134 + t573 * t370 + t562;
t430 = qJD(3) * t457;
t429 = t367 * t379 * t538;
t428 = -0.2e1 * pkin(1) * t462;
t281 = -pkin(13) * t524 + t499;
t426 = -pkin(13) * t173 + qJD(6) * t281 - t235 * t553 - t370 * t79 - (-pkin(13) * t518 + t553) * qJD(4) - (-t357 + (pkin(13) * t371 - t343) * t370) * qJD(5) - t563;
t331 = t375 * t343;
t260 = -pkin(13) * t523 + t331 + (-pkin(11) * t370 - pkin(5)) * t376;
t425 = -qJD(6) * t260 - (-t371 * t492 - t376 * t490) * pkin(11) + t558 + t571 * pkin(13);
t348 = t554 * t370;
t424 = pkin(13) * t533 + qJD(6) * t348 + t370 * t475 + t543;
t136 = t375 * t137;
t349 = t554 * t375;
t423 = pkin(5) * t195 + qJD(6) * t349 - t370 * t61 + t136 + (pkin(13) * t193 + t475) * t375;
t217 = pkin(4) * t529 - t450;
t170 = t215 * t375 + t261 * t370;
t30 = pkin(5) * t214 - pkin(13) * t170 + t458;
t169 = t215 * t370 - t261 * t375;
t33 = -pkin(13) * t169 + t551;
t419 = t30 * t374 - t33 * t369;
t418 = t30 * t369 + t33 * t374;
t109 = t374 * t169 + t170 * t369;
t110 = -t169 * t369 + t170 * t374;
t414 = t366 ^ 2 * t436;
t411 = 0.2e1 * t455 + qJD(2);
t207 = t372 * t430 + t380;
t208 = t377 * t430 + (t405 * qJD(2) + qJD(3) * t559) * t367;
t112 = pkin(3) * t207 - pkin(11) * t208 + t210;
t389 = -t254 * t495 + t270 * t468 + t377 * t273 + t274 * t527 + t301 * t470 + t306 * t530;
t98 = pkin(11) * t438 + t389;
t409 = t112 * t376 - t145 * t493 - t153 * t491 - t371 * t98;
t70 = -pkin(4) * t261 - t453;
t400 = t371 * t112 + t145 * t491 - t153 * t493 + t376 * t98;
t28 = pkin(12) * t207 + t400;
t138 = qJD(4) * t215 + t208 * t371 - t376 * t438;
t139 = -qJD(4) * t214 + t208 * t376 + t371 * t438;
t99 = -t274 * t526 + (-pkin(3) * t472 - t306 * t377) * t366 - t443;
t45 = pkin(4) * t138 - pkin(12) * t139 + t99;
t404 = t375 * t28 + t370 * t45 + t92 * t489 - t490 * t71;
t403 = -t191 * t489 - t525;
t402 = -pkin(12) * t120 + t191 * t56;
t401 = -pkin(11) * t197 + t132 * t233;
t29 = -pkin(4) * t207 - t409;
t385 = -qJD(5) * t551 - t28 * t370 + t375 * t45;
t382 = qJD(3) * t387;
t381 = t412 * t391;
t361 = -pkin(5) * t375 - pkin(4);
t338 = (pkin(5) * t370 + pkin(11)) * t371;
t315 = t332 * t371;
t174 = pkin(5) * t277 + t217;
t116 = t120 * t317;
t97 = t120 * t214;
t65 = -qJD(5) * t169 + t139 * t375 + t207 * t370;
t64 = qJD(5) * t170 + t139 * t370 - t207 * t375;
t50 = pkin(5) * t169 + t70;
t21 = qJD(6) * t110 + t369 * t65 + t374 * t64;
t20 = -qJD(6) * t109 - t369 * t64 + t374 * t65;
t13 = pkin(5) * t64 + t29;
t12 = pkin(5) * t55 + t18;
t10 = t19 * t374 - t26 * t369;
t9 = -pkin(13) * t64 + t404;
t6 = pkin(5) * t138 - pkin(13) * t65 + t385;
t3 = [0, 0, 0, t373 * t378 * t445, -t498 * t445, t411 * t378 * t496, -t411 * t472, 0, -qJD(2) * t381 + t373 * t428 - t391 * t422 -(-pkin(9) * t472 + t351) * t412 - (-pkin(9) * t436 + t344) * t538 + t378 * t428, t196 * t262 + t208 * t237, -t196 * t261 - t197 * t262 - t207 * t237 - t208 * t235, -t208 * t387 - t196 * t316 + (qJD(1) * t262 + t237) * t438, t207 * t387 + t197 * t316 + (-qJD(1) * t261 - t235) * t438 (-t340 + (t456 - t316) * qJD(1) - t444) * t438 -((t274 * t368 + t306 * t366) * t377 + t443) * t387 - t91 * t316 + t210 * t235 + t205 * t197 + t198 * t261 + t186 * t207 + (qJD(1) * t556 + t149) * t438, t389 * t387 - t390 * t316 + t210 * t237 + t205 * t196 + t198 * t262 + t186 * t208 + (-qJD(1) * t479 - t150) * t438, t119 * t215 + t139 * t195, -t119 * t214 - t120 * t215 - t138 * t195 - t139 * t193, t119 * t261 + t139 * t233 + t195 * t207 + t197 * t215, -t120 * t261 - t138 * t233 - t193 * t207 - t197 * t214, t197 * t261 + t207 * t233, t152 * t120 + t132 * t138 + t99 * t193 + t197 * t453 + t61 * t207 + t81 * t214 + t23 * t261 + t233 * t409, t152 * t119 + t132 * t139 + t99 * t195 - t197 * t515 - t62 * t207 + t81 * t215 - t22 * t261 - t233 * t400, t158 * t65 + t170 * t54, -t156 * t65 - t158 * t64 - t169 * t54 - t170 * t55, t120 * t170 + t138 * t158 + t191 * t65 + t214 * t54, -t120 * t169 - t138 * t156 - t191 * t64 - t214 * t55, t138 * t191 + t97, t120 * t458 + t31 * t138 + t29 * t156 + t18 * t169 + t191 * t385 + t8 * t214 + t70 * t55 + t56 * t64, -t120 * t551 - t32 * t138 + t29 * t158 + t18 * t170 - t191 * t404 - t7 * t214 + t70 * t54 + t56 * t65, t110 * t14 - t20 * t417, -t109 * t14 + t110 * t384 - t20 * t83 + t21 * t417, t110 * t120 - t138 * t417 + t14 * t214 + t190 * t20, -t109 * t120 - t138 * t83 - t190 * t21 + t214 * t384, t138 * t190 + t97 (-qJD(6) * t418 - t369 * t9 + t374 * t6) * t190 + t419 * t120 + t2 * t214 + t10 * t138 + t13 * t83 - t50 * t384 + t12 * t109 + t46 * t21 -(qJD(6) * t419 + t369 * t6 + t374 * t9) * t190 - t418 * t120 - t1 * t214 - t11 * t138 - t13 * t417 + t50 * t14 + t12 * t110 + t46 * t20; 0, 0, 0, -t373 * t482, t498 * t531, -t378 * t429, t373 * t429, 0, -pkin(9) * t435 + qJD(1) * t381 + (pkin(1) * t531 - t410) * t373, -t344 + (-pkin(9) * t474 + t350) * t455 + pkin(1) * t482 + t350 * qJD(2), t196 * t530 + t237 * t587, t235 * t292 + t237 * t291 + (t196 * t377 - t197 * t372 + (-t235 * t377 - t237 * t372) * qJD(3)) * t366, t372 * t414 + t196 * t368 + t292 * t387 + (-t237 * t474 - t377 * t382) * t366, t377 * t414 - t197 * t368 - t291 * t387 + (t235 * t474 + t372 * t382) * t366 -(t368 * t455 - t486) * t440, t91 * t368 - t366 * pkin(2) * t197 - t198 * t529 - t209 * t235 + (qJD(2) * t561 - t149) * t440 - t577 * t186 + ((t272 * t368 + t305 * t366) * t377 + t583) * t387, t390 * t368 - t209 * t237 - t186 * t292 + (t186 * t494 - pkin(2) * t196 + t198 * t372 + (-qJD(2) * t500 + t150) * t474) * t366 + t586 * t387, t119 * t318 + t195 * t506, -t119 * t317 - t120 * t318 - t193 * t506 - t195 * t505, -t119 * t529 - t195 * t577 + t197 * t318 + t233 * t506, t120 * t529 + t193 * t577 - t197 * t317 - t233 * t505, -t197 * t529 - t233 * t577, t450 * t197 + t308 * t120 + t81 * t317 - t61 * t291 + (-t23 * t377 + t495 * t61) * t366 + t581 * t233 + t503 * t193 + t505 * t132, -t504 * t197 + t308 * t119 + t81 * t318 + t62 * t291 + (t22 * t377 - t495 * t62) * t366 + t582 * t233 + t503 * t195 + t506 * t132, t158 * t513 + t278 * t54, -t156 * t513 - t158 * t512 - t277 * t54 - t278 * t55, t120 * t278 + t158 * t505 + t191 * t513 + t317 * t54, -t120 * t277 - t156 * t505 - t191 * t512 - t317 * t55, t191 * t505 + t116, t451 * t120 + t8 * t317 + t217 * t55 + t18 * t277 + t512 * t56 + t505 * t31 + (-t218 * t489 + (-qJD(5) * t216 - t573) * t370 - t562) * t191 + t545 * t156, -t509 * t120 + t545 * t158 + t18 * t278 + t191 * t560 + t217 * t54 - t7 * t317 - t505 * t32 + t513 * t56, t14 * t200 - t417 * t540, -t14 * t199 + t200 * t384 + t417 * t539 - t540 * t83, t120 * t200 + t14 * t317 + t190 * t540 - t417 * t505, -t120 * t199 - t190 * t539 + t317 * t384 - t505 * t83, t190 * t505 + t116 (t117 * t374 - t134 * t369) * t120 + t2 * t317 - t174 * t384 + t12 * t199 + t546 * t83 + t539 * t46 + (t369 * t433 - t374 * t432) * t190 + t505 * t10 -(t117 * t369 + t134 * t374) * t120 - t1 * t317 + t174 * t14 + t12 * t200 - t546 * t417 + t540 * t46 + (t369 * t432 + t374 * t433) * t190 - t505 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237 * t235, -t235 ^ 2 + t237 ^ 2, -t235 * t387 + t196, -t237 * t387 - t497 * t555 - t388, t413, -t150 * t387 - t186 * t237 + t91, -t149 * t387 + t186 * t235 + t390, t119 * t371 + t195 * t447 (t119 - t534) * t376 + (-t120 - t532) * t371, -t195 * t237 + t371 * t197 + t233 * t447, t193 * t237 + t376 * t197 - t233 * t564, -t233 * t237, -pkin(3) * t120 - t150 * t193 - t61 * t237 - t81 * t376 + (t143 + (-pkin(11) * qJD(4) - t178) * t376) * t233 + t401 * t371, -pkin(3) * t119 - t150 * t195 + t62 * t237 + t81 * t371 + (t483 + t516) * t233 + t401 * t376, t158 * t576 + t54 * t523, t156 * t173 + t158 * t172 + (-t156 * t375 - t158 * t370) * t491 + (-t547 - t375 * t55 + (t156 * t370 - t158 * t375) * qJD(5)) * t371, -t376 * t54 + t420 * t191 + (t158 * t233 - t466 + t519) * t371, t376 * t55 + t570 * t191 + (-t156 * t233 + t403) * t371, t191 * t564 - t537, t331 * t120 - t78 * t156 - t56 * t172 + ((-qJD(5) * t343 + t79) * t370 + t563) * t191 + (t56 * t370 * qJD(4) - t8 + (qJD(4) * t156 + t403) * pkin(11)) * t376 + (pkin(11) * t55 + t233 * t31 + t489 * t56 + t550) * t371, -t499 * t120 - t78 * t158 - t56 * t173 + t558 * t191 + (t56 * t492 + t7 + (qJD(4) * t158 + t466) * pkin(11)) * t376 + (-t56 * t490 + t549 - t233 * t32 + (t191 * t492 + t54) * pkin(11)) * t371, -t14 * t315 - t417 * t511, -t14 * t314 - t315 * t384 + t417 * t510 - t511 * t83, -t120 * t315 - t14 * t376 + t190 * t511 - t417 * t564, -t120 * t314 - t190 * t510 - t376 * t384 - t564 * t83, t190 * t564 - t537 (t260 * t374 - t281 * t369) * t120 - t2 * t376 - t338 * t384 + t12 * t314 + t544 * t83 + t510 * t46 + (t369 * t425 - t374 * t426) * t190 + t10 * t564 -(t260 * t369 + t281 * t374) * t120 + t1 * t376 + t338 * t14 - t12 * t315 - t544 * t417 + t511 * t46 + (t369 * t426 + t374 * t425) * t190 - t11 * t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 * t195, -t193 ^ 2 + t195 ^ 2, t119 + t534, -t120 + t532, t197, -t132 * t195 + t233 * t62 + t23, t132 * t193 + t233 * t61 - t22, t158 * t446 + t547 (t54 - t536) * t375 + (-t55 - t535) * t370, -t158 * t195 + t191 * t446 + t525, -t191 ^ 2 * t370 + t156 * t195 + t519, -t191 * t195, -pkin(4) * t55 - t62 * t156 - t549 - t31 * t195 + (-pkin(12) * t489 - t136) * t191 + (t61 * t191 + t402) * t370, -pkin(4) * t54 - t62 * t158 + t550 + t32 * t195 + (pkin(12) * t490 + t543) * t191 + t402 * t375, t14 * t333 + t417 * t508, -t14 * t332 + t333 * t384 + t417 * t507 + t508 * t83, t120 * t333 - t190 * t508 + t195 * t417, -t120 * t332 - t190 * t507 + t195 * t83, -t190 * t195 (-t348 * t374 - t349 * t369) * t120 - t361 * t384 + t12 * t332 - t10 * t195 + t434 * t83 + t507 * t46 + (t369 * t424 - t374 * t423) * t190 -(-t348 * t369 + t349 * t374) * t120 + t361 * t14 + t12 * t333 + t11 * t195 - t434 * t417 - t508 * t46 + (t369 * t423 + t374 * t424) * t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158 * t156, -t156 ^ 2 + t158 ^ 2, t54 + t536, t535 - t55, t120, -t158 * t56 + t191 * t32 + t8, t156 * t56 + t191 * t31 - t7, -t574, t568, t566, t565, t120 -(-t25 * t369 - t548) * t190 + (t120 * t374 - t158 * t83 - t190 * t488) * pkin(5) + t567, t575 + t24 + (-t190 * t26 - t4) * t369 + (t190 * t25 - t460) * t374 + (-t120 * t369 + t158 * t417 - t190 * t487) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t574, t568, t566, t565, t120, t11 * t190 + t567, t10 * t190 - t1 + t575;];
tauc_reg  = t3;
