% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:28
% EndTime: 2019-03-09 22:10:00
% DurationCPUTime: 17.86s
% Computational Cost: add. (32312->827), mult. (71630->1033), div. (0->0), fcn. (54026->18), ass. (0->401)
t426 = cos(qJ(2));
t611 = cos(qJ(3));
t507 = t611 * t426;
t382 = qJD(1) * t507;
t422 = sin(qJ(3));
t423 = sin(qJ(2));
t526 = qJD(1) * t423;
t506 = t422 * t526;
t326 = -t382 + t506;
t346 = t422 * t426 + t423 * t611;
t328 = t346 * qJD(1);
t243 = pkin(3) * t328 + pkin(9) * t326;
t229 = pkin(2) * t526 + t243;
t428 = -pkin(8) - pkin(7);
t368 = t428 * t426;
t354 = qJD(1) * t368;
t329 = t422 * t354;
t367 = t428 * t423;
t352 = qJD(1) * t367;
t257 = t352 * t611 + t329;
t421 = sin(qJ(4));
t425 = cos(qJ(4));
t175 = t421 * t229 + t425 * t257;
t502 = t611 * qJD(3);
t486 = pkin(2) * t502;
t650 = -t425 * t486 + t175;
t556 = t326 * t421;
t649 = qJ(5) * t556 - t425 * qJD(5);
t174 = t425 * t229 - t257 * t421;
t411 = t425 * qJ(5);
t477 = t328 * pkin(4) + t326 * t411;
t609 = pkin(2) * t422;
t399 = pkin(9) + t609;
t538 = -qJ(5) - t399;
t490 = qJD(4) * t538;
t648 = -t174 - t477 + (-t486 - qJD(5)) * t421 + t425 * t490;
t587 = qJD(2) * pkin(2);
t332 = t352 + t587;
t247 = t332 * t611 + t329;
t176 = t425 * t243 - t247 * t421;
t419 = -qJ(5) - pkin(9);
t496 = qJD(4) * t419;
t647 = -t421 * qJD(5) + t425 * t496 - t176 - t477;
t646 = -t421 * t490 + t649 + t650;
t177 = t421 * t243 + t425 * t247;
t645 = -t421 * t496 + t177 + t649;
t417 = sin(pkin(11));
t418 = cos(pkin(11));
t341 = t417 * t425 + t418 * t421;
t324 = t341 * qJD(4);
t631 = t341 * t326 + t324;
t617 = -t417 * t421 + t418 * t425;
t325 = t617 * qJD(4);
t630 = t617 * t326 + t325;
t517 = qJD(2) + qJD(3);
t487 = t425 * t517;
t272 = t328 * t421 - t487;
t274 = t425 * t328 + t421 * t517;
t195 = t272 * t418 + t417 * t274;
t420 = sin(qJ(6));
t471 = -t272 * t417 + t418 * t274;
t610 = cos(qJ(6));
t110 = t195 * t610 + t420 * t471;
t577 = t110 ^ 2;
t627 = -t420 * t195 + t471 * t610;
t579 = t627 ^ 2;
t644 = -t577 + t579;
t453 = t611 * qJD(2) + t502;
t495 = qJDD(1) * t611;
t521 = qJD(1) * qJD(2);
t501 = t423 * t521;
t529 = qJD(3) * t506 + t422 * t501;
t433 = t423 * t495 + (qJD(1) * t453 + qJDD(1) * t422) * t426 - t529;
t516 = qJDD(2) + qJDD(3);
t524 = qJD(4) * t421;
t161 = -qJD(4) * t487 + t328 * t524 - t421 * t516 - t425 * t433;
t162 = qJD(4) * t274 + t421 * t433 - t425 * t516;
t100 = -t161 * t418 - t162 * t417;
t504 = qJD(6) * t610;
t522 = qJD(6) * t420;
t99 = -t417 * t161 + t162 * t418;
t32 = -t610 * t100 + t195 * t504 + t420 * t99 + t471 * t522;
t318 = qJD(4) + t326;
t300 = qJD(6) + t318;
t576 = t110 * t300;
t643 = -t32 + t576;
t575 = t110 * t627;
t230 = -pkin(3) * t517 - t247;
t185 = t272 * pkin(4) + qJD(5) + t230;
t108 = t195 * pkin(5) + t185;
t413 = qJ(4) + pkin(11);
t407 = qJ(6) + t413;
t391 = sin(t407);
t392 = cos(t407);
t427 = cos(qJ(1));
t416 = qJ(2) + qJ(3);
t410 = cos(t416);
t424 = sin(qJ(1));
t548 = t410 * t424;
t279 = t391 * t427 - t392 * t548;
t547 = t410 * t427;
t281 = t391 * t424 + t392 * t547;
t261 = t517 * t346;
t519 = t423 * qJDD(1);
t475 = t422 * t519 - t426 * t495;
t208 = qJD(1) * t261 + t475;
t205 = qJDD(4) + t208;
t387 = pkin(2) * t501;
t592 = t426 * pkin(2);
t402 = pkin(1) + t592;
t619 = -pkin(9) * t346 - t402;
t120 = t387 - (t382 * t517 - t529) * pkin(9) + t208 * pkin(3) + t619 * qJDD(1);
t115 = t425 * t120;
t500 = t426 * t521;
t267 = qJDD(2) * pkin(2) - t428 * (-t500 - t519);
t518 = t426 * qJDD(1);
t271 = t428 * (-t501 + t518);
t525 = qJD(3) * t422;
t149 = t422 * t267 - t611 * t271 + t332 * t502 + t354 * t525;
t142 = pkin(9) * t516 + t149;
t366 = t402 * qJD(1);
t221 = pkin(3) * t326 - pkin(9) * t328 - t366;
t330 = t611 * t354;
t248 = t422 * t332 - t330;
t231 = pkin(9) * t517 + t248;
t164 = t221 * t421 + t231 * t425;
t53 = -qJD(4) * t164 - t142 * t421 + t115;
t36 = pkin(4) * t205 + qJ(5) * t161 - qJD(5) * t274 + t53;
t523 = qJD(4) * t425;
t52 = t421 * t120 + t425 * t142 + t221 * t523 - t231 * t524;
t40 = -qJ(5) * t162 - qJD(5) * t272 + t52;
t16 = t417 * t36 + t418 * t40;
t10 = -pkin(10) * t99 + t16;
t626 = pkin(10) * t471;
t163 = t425 * t221 - t231 * t421;
t123 = -qJ(5) * t274 + t163;
t107 = pkin(4) * t318 + t123;
t124 = -qJ(5) * t272 + t164;
t116 = t417 * t124;
t72 = t418 * t107 - t116;
t48 = pkin(5) * t318 - t626 + t72;
t637 = pkin(10) * t195;
t545 = t418 * t124;
t73 = t417 * t107 + t545;
t54 = t73 - t637;
t15 = t418 * t36 - t40 * t417;
t9 = pkin(5) * t205 - pkin(10) * t100 + t15;
t3 = t10 * t610 + t420 * t9 + t48 * t504 - t54 * t522;
t409 = sin(t416);
t396 = g(3) * t409;
t642 = g(1) * t281 - g(2) * t279 + t108 * t110 + t392 * t396 - t3;
t586 = t646 * t417 + t418 * t648;
t585 = t417 * t648 - t646 * t418;
t584 = t417 * t645 + t418 * t647;
t583 = t417 * t647 - t418 * t645;
t33 = qJD(6) * t627 + t420 * t100 + t610 * t99;
t578 = t627 * t300;
t641 = -t33 + t578;
t640 = t631 * pkin(10);
t639 = -t328 * pkin(5) - pkin(10) * t630;
t278 = t391 * t548 + t392 * t427;
t280 = -t391 * t547 + t392 * t424;
t23 = t420 * t48 + t54 * t610;
t4 = -qJD(6) * t23 - t420 * t10 + t610 * t9;
t638 = -g(1) * t280 + g(2) * t278 - t108 * t627 + t391 * t396 + t4;
t636 = t639 + t586;
t635 = -t640 + t585;
t634 = t639 + t584;
t633 = -t640 + t583;
t632 = t195 * t471;
t250 = t341 * t610 + t420 * t617;
t537 = -qJD(6) * t250 - t420 * t630 - t610 * t631;
t536 = t341 * t522 + t420 * t631 - t504 * t617 - t610 * t630;
t256 = t352 * t422 - t330;
t482 = pkin(2) * t525 - t256;
t549 = t409 * t427;
t550 = t409 * t424;
t629 = g(1) * t549 + g(2) * t550;
t397 = g(3) * t410;
t479 = g(1) * t427 + g(2) * t424;
t462 = t479 * t409;
t628 = t397 - t462;
t488 = t422 * t517;
t260 = t423 * t488 - t426 * t453;
t505 = t346 * t523;
t460 = -t260 * t421 + t505;
t614 = t397 - t629;
t345 = t422 * t423 - t507;
t246 = pkin(3) * t345 + t619;
t276 = t422 * t367 - t368 * t611;
t268 = t425 * t276;
t189 = t421 * t246 + t268;
t403 = pkin(4) * t524;
t623 = pkin(5) * t631 + t403;
t406 = cos(t413);
t593 = t425 * pkin(4);
t358 = pkin(5) * t406 + t593;
t348 = pkin(3) + t358;
t412 = -pkin(10) + t419;
t622 = t410 * t348 - t409 * t412;
t621 = t611 * t367 + t422 * t368;
t400 = pkin(3) + t593;
t620 = t410 * t400 - t409 * t419;
t480 = t410 * pkin(3) + t409 * pkin(9);
t294 = pkin(4) * t556;
t618 = t294 + t482;
t616 = -t163 * t421 + t164 * t425;
t539 = t425 * t427;
t543 = t421 * t424;
t310 = t410 * t543 + t539;
t540 = t424 * t425;
t542 = t421 * t427;
t312 = -t410 * t542 + t540;
t615 = -g(1) * t312 + g(2) * t310;
t612 = t99 * pkin(5);
t608 = pkin(2) * t423;
t607 = pkin(4) * t417;
t606 = pkin(4) * t421;
t603 = pkin(10) * t341;
t369 = t427 * t402;
t599 = g(2) * t369;
t597 = g(3) * t421;
t596 = g(3) * t426;
t594 = t617 * pkin(5);
t336 = t538 * t421;
t337 = t399 * t425 + t411;
t240 = t418 * t336 - t337 * t417;
t206 = t240 - t603;
t241 = t417 * t336 + t418 * t337;
t335 = t617 * pkin(10);
t207 = t335 + t241;
t136 = t206 * t610 - t420 * t207;
t591 = qJD(6) * t136 + t420 * t636 + t610 * t635;
t137 = t420 * t206 + t207 * t610;
t590 = -qJD(6) * t137 - t420 * t635 + t610 * t636;
t466 = qJ(5) * t260 - qJD(5) * t346;
t514 = t423 * t587;
t184 = pkin(3) * t261 + pkin(9) * t260 + t514;
t508 = qJD(2) * t428;
t353 = t423 * t508;
t355 = t426 * t508;
t200 = qJD(3) * t621 + t611 * t353 + t422 * t355;
t493 = t425 * t184 - t200 * t421;
t59 = pkin(4) * t261 + t466 * t425 + (-t268 + (qJ(5) * t346 - t246) * t421) * qJD(4) + t493;
t510 = t421 * t184 + t425 * t200 + t246 * t523;
t69 = -qJ(5) * t505 + (-qJD(4) * t276 + t466) * t421 + t510;
t27 = t417 * t59 + t418 * t69;
t364 = t419 * t421;
t365 = pkin(9) * t425 + t411;
t269 = t418 * t364 - t365 * t417;
t226 = t269 - t603;
t270 = t417 * t364 + t418 * t365;
t227 = t335 + t270;
t158 = t226 * t610 - t420 * t227;
t589 = qJD(6) * t158 + t420 * t634 + t610 * t633;
t159 = t420 * t226 + t227 * t610;
t588 = -qJD(6) * t159 - t420 * t633 + t610 * t634;
t51 = t52 * t425;
t393 = pkin(4) * t418 + pkin(5);
t308 = t393 * t610 - t420 * t607;
t75 = -t123 * t417 - t545;
t60 = t75 + t637;
t76 = t418 * t123 - t116;
t61 = t76 - t626;
t582 = t308 * qJD(6) - t420 * t60 - t61 * t610;
t309 = t420 * t393 + t607 * t610;
t581 = -t309 * qJD(6) + t420 * t61 - t60 * t610;
t580 = pkin(7) * qJDD(1);
t574 = t161 * t421;
t573 = t162 * t425;
t571 = t163 * t425;
t569 = t471 ^ 2;
t568 = t471 * t318;
t567 = t195 ^ 2;
t566 = t195 * t318;
t565 = t230 * t326;
t563 = t272 * t318;
t562 = t272 * t421;
t561 = t274 * t272;
t560 = t274 * t318;
t559 = t274 * t425;
t558 = t300 * t328;
t557 = t318 * t328;
t555 = t328 * t326;
t554 = t346 * t421;
t553 = t346 * t425;
t188 = t425 * t246 - t276 * t421;
t152 = pkin(4) * t345 - t346 * t411 + t188;
t172 = -qJ(5) * t554 + t189;
t98 = t417 * t152 + t418 * t172;
t533 = t618 + t623;
t204 = -t294 + t248;
t532 = -t204 + t623;
t531 = t403 + t618;
t405 = sin(t413);
t357 = pkin(5) * t405 + t606;
t530 = t357 - t428;
t414 = t423 ^ 2;
t415 = t426 ^ 2;
t528 = t414 - t415;
t527 = t414 + t415;
t515 = t611 * pkin(2);
t513 = pkin(9) * qJD(4) * t318;
t430 = qJD(1) ^ 2;
t511 = t423 * t430 * t426;
t509 = g(1) * t547 + g(2) * t548 + t396;
t499 = -t428 + t606;
t485 = -t611 * t267 - t422 * t271 + t332 * t525 - t354 * t502;
t143 = -pkin(3) * t516 + t485;
t497 = -t143 - t397;
t26 = -t417 * t69 + t418 * t59;
t97 = t418 * t152 - t172 * t417;
t489 = t318 * t425;
t401 = -t515 - pkin(3);
t484 = t423 * t500;
t483 = -g(1) * t550 + g(2) * t549;
t481 = -t204 + t403;
t478 = g(1) * t424 - g(2) * t427;
t223 = pkin(4) * t554 - t621;
t476 = -pkin(9) * t205 + t565;
t473 = t164 * t421 + t571;
t472 = -t205 * t399 + t565;
t469 = t348 * t409 + t410 * t412;
t467 = t400 * t409 + t410 * t419;
t465 = t143 * t421 + t164 * t328 + t230 * t523 + t410 * t597;
t464 = -t163 * t328 + t230 * t524 + t425 * t629;
t463 = -t164 * t556 - t326 * t571 - t509 + t51;
t234 = t617 * t346;
t74 = pkin(5) * t345 - pkin(10) * t234 + t97;
t233 = t341 * t346;
t79 = -pkin(10) * t233 + t98;
t37 = -t420 * t79 + t610 * t74;
t38 = t420 * t74 + t610 * t79;
t461 = -0.2e1 * pkin(1) * t521 - pkin(7) * qJDD(2);
t171 = -t420 * t233 + t234 * t610;
t459 = -t260 * t425 - t346 * t524;
t363 = t401 - t593;
t429 = qJD(2) ^ 2;
t452 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t429 + t478;
t451 = pkin(1) * t430 + t479 - t580;
t450 = -t366 * t326 - t149 + t509;
t449 = t366 * t328 - t485 - t614;
t448 = -qJD(4) * t473 - t53 * t421;
t22 = -t420 * t54 + t48 * t610;
t249 = t341 * t420 - t610 * t617;
t447 = t22 * t536 + t23 * t537 - t3 * t249 - t4 * t250 - t509;
t446 = -t15 * t341 + t16 * t617 - t630 * t72 - t631 * t73 - t509;
t93 = t162 * pkin(4) + qJDD(5) + t143;
t45 = t93 + t612;
t444 = -t108 * t537 - t22 * t328 + t45 * t249 - t392 * t614;
t201 = qJD(3) * t276 + t422 * t353 - t611 * t355;
t443 = t185 * t631 - t328 * t72 - t406 * t614 - t617 * t93;
t442 = t448 + t51;
t135 = pkin(4) * t460 + t201;
t440 = t93 + t614;
t439 = -t108 * t536 + t23 * t328 + t45 * t250 + t391 * t628;
t438 = t185 * t630 + t328 * t73 + t93 * t341 + t405 * t628;
t437 = -g(1) * (-pkin(3) * t549 + pkin(9) * t547) - g(2) * (-pkin(3) * t550 + pkin(9) * t548) - g(3) * t480;
t320 = -qJDD(1) * t402 + t387;
t313 = t410 * t539 + t543;
t311 = -t410 * t540 + t542;
t288 = -t400 - t594;
t286 = t405 * t424 + t406 * t547;
t285 = -t405 * t547 + t406 * t424;
t284 = t405 * t427 - t406 * t548;
t283 = t405 * t548 + t406 * t427;
t277 = t363 - t594;
t210 = -t326 ^ 2 + t328 ^ 2;
t203 = qJDD(6) + t205;
t182 = t326 * t517 + t433;
t173 = pkin(5) * t233 + t223;
t170 = t233 * t610 + t234 * t420;
t160 = pkin(4) * t274 + pkin(5) * t471;
t133 = t205 * t345 + t261 * t318;
t132 = t260 * t617 + t324 * t346;
t131 = t260 * t341 - t325 * t346;
t104 = t205 * t421 - t274 * t328 + t318 * t489;
t103 = -t318 ^ 2 * t421 + t205 * t425 + t272 * t328;
t102 = t318 * t562 - t573;
t101 = t274 * t489 - t574;
t92 = -qJD(4) * t189 + t493;
t91 = -t276 * t524 + t510;
t80 = -t131 * pkin(5) + t135;
t78 = t205 * t341 + t318 * t630 - t328 * t471;
t77 = t195 * t328 + t205 * t617 - t318 * t631;
t56 = qJD(6) * t171 - t131 * t610 - t420 * t132;
t55 = -t420 * t131 + t132 * t610 + t233 * t504 + t234 * t522;
t49 = (-t161 - t563) * t425 + (-t162 - t560) * t421;
t47 = t100 * t341 + t471 * t630;
t46 = t195 * t631 - t617 * t99;
t42 = t110 * t328 - t203 * t249 + t300 * t537;
t41 = t203 * t250 - t300 * t536 - t328 * t627;
t21 = pkin(10) * t131 + t27;
t20 = pkin(5) * t261 + pkin(10) * t132 + t26;
t17 = t100 * t617 - t195 * t630 - t341 * t99 - t471 * t631;
t12 = -t110 * t537 + t249 * t33;
t11 = -t250 * t32 - t536 * t627;
t7 = -qJD(6) * t38 + t20 * t610 - t420 * t21;
t6 = qJD(6) * t37 + t420 * t20 + t21 * t610;
t5 = t110 * t536 + t249 * t32 - t250 * t33 + t537 * t627;
t1 = [0, 0, 0, 0, 0, qJDD(1), t478, t479, 0, 0, qJDD(1) * t414 + 0.2e1 * t484, 0.2e1 * t423 * t518 - 0.2e1 * t521 * t528, qJDD(2) * t423 + t426 * t429, qJDD(1) * t415 - 0.2e1 * t484, qJDD(2) * t426 - t423 * t429, 0, t423 * t461 + t426 * t452, -t423 * t452 + t426 * t461, 0.2e1 * t527 * t580 - t479, -g(1) * (-pkin(1) * t424 + pkin(7) * t427) - g(2) * (pkin(1) * t427 + pkin(7) * t424) + (pkin(7) ^ 2 * t527 + pkin(1) ^ 2) * qJDD(1), -t328 * t260 + t346 * t433, -t346 * t208 + t260 * t326 - t328 * t261 - t345 * t433, -t260 * t517 + t346 * t516, t208 * t345 + t261 * t326, -t261 * t517 - t345 * t516, 0, -t201 * t517 - t402 * t208 - t366 * t261 + t320 * t345 + t326 * t514 + t410 * t478 + t516 * t621, -t200 * t517 + t366 * t260 - t276 * t516 + t320 * t346 + t328 * t514 - t402 * t433 + t483, -t149 * t345 - t200 * t326 + t201 * t328 - t276 * t208 + t247 * t260 - t248 * t261 + t346 * t485 - t433 * t621 - t479, t149 * t276 + t248 * t200 - t485 * t621 - t247 * t201 - t320 * t402 - t366 * t514 - g(1) * (-t402 * t424 - t427 * t428) - g(2) * (-t424 * t428 + t369) -t161 * t553 + t274 * t459 (t272 * t425 + t274 * t421) * t260 + (t574 - t573 + (-t559 + t562) * qJD(4)) * t346, -t161 * t345 + t205 * t553 + t261 * t274 + t318 * t459, t162 * t554 + t272 * t460, -t162 * t345 - t205 * t554 - t261 * t272 - t318 * t460, t133, -g(1) * t311 - g(2) * t313 + t143 * t554 - t162 * t621 + t163 * t261 + t188 * t205 + t201 * t272 + t230 * t460 + t318 * t92 + t345 * t53, -g(1) * t310 - g(2) * t312 + t143 * t553 + t161 * t621 - t164 * t261 - t189 * t205 + t201 * t274 + t230 * t459 - t318 * t91 - t345 * t52, t161 * t188 - t162 * t189 - t272 * t91 - t274 * t92 + t473 * t260 + (-qJD(4) * t616 - t421 * t52 - t425 * t53) * t346 - t483, -t599 - t143 * t621 + t163 * t92 + t164 * t91 + t53 * t188 + t52 * t189 + t230 * t201 + (g(1) * t428 - g(2) * t480) * t427 + (-g(1) * (-t402 - t480) + g(2) * t428) * t424, t100 * t234 - t132 * t471, -t100 * t233 + t131 * t471 + t132 * t195 - t234 * t99, t100 * t345 - t132 * t318 + t205 * t234 + t261 * t471, -t131 * t195 + t233 * t99, t131 * t318 - t195 * t261 - t205 * t233 - t345 * t99, t133, -g(1) * t284 - g(2) * t286 - t131 * t185 + t135 * t195 + t15 * t345 + t205 * t97 + t223 * t99 + t233 * t93 + t26 * t318 + t261 * t72, -g(1) * t283 - g(2) * t285 + t100 * t223 - t132 * t185 + t135 * t471 - t16 * t345 - t205 * t98 + t234 * t93 - t261 * t73 - t27 * t318, -t100 * t97 + t131 * t73 + t132 * t72 - t15 * t234 - t16 * t233 - t195 * t27 - t26 * t471 - t98 * t99 - t483, -t599 + t185 * t135 + t15 * t97 + t16 * t98 + t93 * t223 + t72 * t26 + t73 * t27 + (-g(1) * t499 - g(2) * t620) * t427 + (-g(1) * (-t402 - t620) - g(2) * t499) * t424, -t171 * t32 - t55 * t627, t110 * t55 + t170 * t32 - t171 * t33 - t56 * t627, t171 * t203 + t261 * t627 - t300 * t55 - t32 * t345, t110 * t56 + t170 * t33, -t110 * t261 - t170 * t203 - t300 * t56 - t33 * t345, t203 * t345 + t261 * t300, -g(1) * t279 - g(2) * t281 + t108 * t56 + t110 * t80 + t170 * t45 + t173 * t33 + t203 * t37 + t22 * t261 + t300 * t7 + t345 * t4, -g(1) * t278 - g(2) * t280 - t108 * t55 + t171 * t45 - t173 * t32 - t203 * t38 - t23 * t261 - t3 * t345 - t300 * t6 + t627 * t80, -t110 * t6 - t170 * t3 - t171 * t4 + t22 * t55 - t23 * t56 + t32 * t37 - t33 * t38 - t627 * t7 - t483, -t599 + t108 * t80 + t45 * t173 + t22 * t7 + t23 * t6 + t3 * t38 + t4 * t37 + (-g(1) * t530 - g(2) * t622) * t427 + (-g(1) * (-t402 - t622) - g(2) * t530) * t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t511, t528 * t430, t519, t511, t518, qJDD(2), t423 * t451 - t596, g(3) * t423 + t426 * t451, 0, 0, t555, t210, t182, -t555, -t475, t516, t256 * t517 + (-qJD(3) * t488 - t326 * t526 + t516 * t611) * pkin(2) + t449, t257 * t517 + (-t328 * t526 - t422 * t516 - t502 * t517) * pkin(2) + t450, -t208 * t609 - t433 * t515 + (t248 + t482) * t328 + (-t247 + t257 - t486) * t326, t247 * t256 - t248 * t257 + (-t611 * t485 - t596 + t149 * t422 + (-t247 * t422 + t248 * t611) * qJD(3) + (qJD(1) * t366 + t479) * t423) * pkin(2), t101, t49, t104, t102, t103, -t557, t401 * t162 + t497 * t425 + t472 * t421 + t482 * t272 + (-t399 * t523 - t421 * t486 - t174) * t318 + t464, -t401 * t161 + t472 * t425 - t421 * t462 + t482 * t274 + (t399 * t524 + t650) * t318 + t465, t174 * t274 + t175 * t272 + (-t272 * t486 - t162 * t399 + (t274 * t399 - t163) * qJD(4)) * t425 + (t274 * t486 - t161 * t399 - t53 + (t272 * t399 - t164) * qJD(4)) * t421 + t463, t143 * t401 - t164 * t175 - t163 * t174 - t230 * t256 + (-t596 + t479 * t423 + (t230 * t422 + t611 * t616) * qJD(3)) * pkin(2) + t442 * t399 + t437, t47, t17, t78, t46, t77, -t557, t195 * t531 + t205 * t240 + t318 * t586 + t363 * t99 + t443, t100 * t363 - t205 * t241 - t318 * t585 + t471 * t531 + t438, -t100 * t240 - t195 * t585 - t241 * t99 - t471 * t586 + t446, t16 * t241 + t15 * t240 + t93 * t363 - g(3) * (t620 + t592) + t585 * t73 + t586 * t72 + t531 * t185 + t479 * (t467 + t608) t11, t5, t41, t12, t42, -t558, t110 * t533 + t136 * t203 + t277 * t33 + t300 * t590 + t444, -t137 * t203 - t277 * t32 - t300 * t591 + t533 * t627 + t439, -t110 * t591 + t136 * t32 - t137 * t33 - t590 * t627 + t447, t3 * t137 + t4 * t136 + t45 * t277 - g(3) * (t622 + t592) + t591 * t23 + t590 * t22 + t533 * t108 + t479 * (t469 + t608); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t555, t210, t182, -t555, -t475, t516, t248 * t517 + t449, t247 * t517 + t450, 0, 0, t101, t49, t104, t102, t103, -t557, -pkin(3) * t162 - t176 * t318 - t248 * t272 + t476 * t421 + (t497 - t513) * t425 + t464, pkin(3) * t161 + t177 * t318 - t248 * t274 + t476 * t425 + (-t462 + t513) * t421 + t465, t176 * t274 + t177 * t272 + (-t574 - t573 + (t559 + t562) * qJD(4)) * pkin(9) + t448 + t463, -t143 * pkin(3) + pkin(9) * t442 - t163 * t176 - t164 * t177 - t230 * t248 + t437, t47, t17, t78, t46, t77, -t557, t195 * t481 + t205 * t269 + t318 * t584 - t400 * t99 + t443, -t100 * t400 - t205 * t270 - t318 * t583 + t471 * t481 + t438, -t100 * t269 - t195 * t583 - t270 * t99 - t471 * t584 + t446, -g(3) * t620 + t15 * t269 + t16 * t270 + t185 * t481 - t93 * t400 + t467 * t479 + t583 * t73 + t584 * t72, t11, t5, t41, t12, t42, -t558, t110 * t532 + t158 * t203 + t288 * t33 + t300 * t588 + t444, -t159 * t203 - t288 * t32 - t300 * t589 + t532 * t627 + t439, -t110 * t589 + t158 * t32 - t159 * t33 - t588 * t627 + t447, -g(3) * t622 + t108 * t532 + t4 * t158 + t3 * t159 + t22 * t588 + t23 * t589 + t45 * t288 + t469 * t479; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t561, -t272 ^ 2 + t274 ^ 2, -t161 + t563, -t561, -t162 + t560, t205, -t231 * t523 + t164 * t318 - t230 * t274 + t115 + (-qJD(4) * t221 - t142 + t396) * t421 + t615, g(1) * t313 - g(2) * t311 + t163 * t318 + t230 * t272 + t396 * t425 - t52, 0, 0, t632, -t567 + t569, t100 + t566, -t632, -t99 + t568, t205, t405 * t396 - g(1) * t285 + g(2) * t283 - t185 * t471 - t318 * t75 + (-t195 * t274 + t205 * t418) * pkin(4) + t15, t406 * t396 + g(1) * t286 - g(2) * t284 + t185 * t195 + t318 * t76 + (-t205 * t417 - t274 * t471) * pkin(4) - t16 (-t100 * t418 - t417 * t99) * pkin(4) + (t73 + t75) * t471 + (t76 - t72) * t195, -t72 * t75 - t73 * t76 + (t15 * t418 + t16 * t417 - t185 * t274 + t409 * t597 + t615) * pkin(4), t575, t644, t643, -t575, t641, t203, -t160 * t110 + t308 * t203 + t300 * t581 + t638, -t160 * t627 - t309 * t203 - t300 * t582 + t642, t308 * t32 - t309 * t33 + (t23 - t581) * t627 + (-t582 - t22) * t110, t3 * t309 + t4 * t308 - t108 * t160 - g(1) * (-t357 * t547 + t358 * t424) - g(2) * (-t357 * t548 - t358 * t427) + t357 * t396 + t582 * t23 + t581 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 + t568, t100 - t566, -t567 - t569, t195 * t73 + t471 * t72 + t440, 0, 0, 0, 0, 0, 0, t33 + t578, -t32 - t576, -t577 - t579, t110 * t23 + t22 * t627 + t440 + t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, t644, t643, -t575, t641, t203, t23 * t300 + t638, t22 * t300 + t642, 0, 0;];
tau_reg  = t1;
