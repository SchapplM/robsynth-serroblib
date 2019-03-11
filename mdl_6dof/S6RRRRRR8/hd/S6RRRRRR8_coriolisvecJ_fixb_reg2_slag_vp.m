% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRRR8
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
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:46
% EndTime: 2019-03-10 05:07:58
% DurationCPUTime: 34.13s
% Computational Cost: add. (74608->926), mult. (217814->1283), div. (0->0), fcn. (182502->14), ass. (0->395)
t423 = cos(qJ(2));
t587 = cos(pkin(6));
t534 = pkin(1) * t587;
t404 = t423 * t534;
t394 = qJD(1) * t404;
t421 = sin(qJ(2));
t416 = sin(pkin(6));
t586 = cos(pkin(7));
t464 = t416 * (-pkin(10) * t586 - pkin(9));
t452 = t421 * t464;
t319 = qJD(1) * t452 + t394;
t403 = t421 * t534;
t432 = t423 * t464 - t403;
t320 = t432 * qJD(1);
t415 = sin(pkin(7));
t607 = pkin(10) * t415;
t451 = (pkin(2) * t421 - t423 * t607) * t416;
t349 = qJD(1) * t451;
t610 = cos(qJ(3));
t480 = t586 * t610;
t420 = sin(qJ(3));
t572 = t415 * t420;
t370 = pkin(2) * t480 - pkin(10) * t572;
t511 = t420 * t586;
t622 = t370 * qJD(3) - t610 * t319 - t320 * t511 - t349 * t572;
t443 = -t421 * t511 + t423 * t610;
t551 = qJD(1) * t416;
t338 = t443 * t551;
t518 = qJD(3) * t610;
t489 = t415 * t518;
t621 = t489 - t338;
t251 = -t320 * t415 + t586 * t349;
t441 = t420 * t423 + t421 * t480;
t337 = t441 * t551;
t655 = pkin(3) * t337 - pkin(11) * t338 + t251 - (pkin(3) * t420 - pkin(11) * t610) * t415 * qJD(3);
t522 = t421 * t551;
t494 = t415 * t522;
t654 = pkin(11) * t494 - t622;
t419 = sin(qJ(4));
t422 = cos(qJ(4));
t368 = t419 * t586 + t422 * t572;
t558 = -qJD(4) * t368 - t419 * t621 - t422 * t494;
t367 = t419 * t572 - t422 * t586;
t557 = qJD(4) * t367 + t419 * t494 - t422 * t621;
t549 = qJD(3) * t420;
t519 = t415 * t549;
t623 = t337 - t519;
t529 = t415 * t610;
t372 = pkin(2) * t511 + pkin(10) * t529;
t356 = pkin(11) * t586 + t372;
t357 = (-pkin(3) * t610 - pkin(11) * t420 - pkin(2)) * t415;
t260 = t422 * t356 + t419 * t357;
t567 = qJD(4) * t260 - t654 * t419 + t655 * t422;
t546 = qJD(4) * t422;
t547 = qJD(4) * t419;
t566 = t356 * t547 - t357 * t546 + t655 * t419 + t654 * t422;
t653 = -pkin(4) * t623 + pkin(12) * t557 - t567;
t652 = -pkin(12) * t558 + t566;
t571 = t416 * t423;
t553 = pkin(9) * t571 + t403;
t361 = t553 * qJD(1);
t505 = t587 * qJD(1);
t470 = t505 + qJD(2);
t453 = t415 * t470;
t509 = t423 * t586;
t482 = t416 * t509;
t271 = t361 + (qJD(1) * t482 + t453) * pkin(10);
t435 = pkin(2) * t587 + t452;
t279 = qJD(2) * pkin(2) + qJD(1) * t435 + t394;
t345 = (-pkin(2) * t423 - t421 * t607 - pkin(1)) * t416;
t331 = qJD(1) * t345;
t191 = -t420 * t271 + t279 * t480 + t331 * t529;
t376 = t610 * t453;
t467 = t423 * t480;
t455 = t416 * t467;
t429 = qJD(1) * t455 + t376;
t493 = t420 * t522;
t280 = t493 - t429;
t442 = t420 * t509 + t421 * t610;
t436 = t442 * t416;
t448 = t420 * t453;
t282 = qJD(1) * t436 + t448;
t223 = pkin(3) * t282 + pkin(11) * t280;
t131 = t422 * t191 + t419 * t223;
t612 = -pkin(12) - pkin(11);
t533 = qJD(4) * t612;
t576 = t280 * t419;
t651 = -pkin(12) * t576 + t419 * t533 - t131;
t650 = -t372 * qJD(3) + t420 * t319 - t320 * t480;
t608 = cos(qJ(6));
t516 = qJD(6) * t608;
t536 = t415 * t571;
t340 = qJD(1) * t536 - t470 * t586 - qJD(3);
t235 = -t419 * t282 - t340 * t422;
t236 = t282 * t422 - t340 * t419;
t418 = sin(qJ(5));
t609 = cos(qJ(5));
t169 = -t609 * t235 + t236 * t418;
t638 = t608 * t169;
t649 = t516 + t638;
t527 = t610 * t349;
t563 = (-pkin(3) * t522 - t527) * t415 + t650;
t517 = qJD(5) * t609;
t545 = qJD(5) * t418;
t562 = t367 * t517 + t368 * t545 - t418 * t558 + t557 * t609;
t278 = -t418 * t367 + t368 * t609;
t561 = qJD(5) * t278 - t418 * t557 - t558 * t609;
t525 = t609 * t422;
t570 = t418 * t419;
t379 = -t525 + t570;
t488 = qJD(4) * t525;
t620 = qJD(4) + qJD(5);
t556 = t379 * t280 - t422 * t517 + t570 * t620 - t488;
t569 = t418 * t422;
t380 = t419 * t609 + t569;
t555 = (t280 + t620) * t380;
t417 = sin(qJ(6));
t276 = qJD(4) + t280;
t438 = qJD(5) + t276;
t433 = t608 * t438;
t458 = t418 * t235 + t236 * t609;
t141 = t417 * t458 - t433;
t143 = t417 * t438 + t458 * t608;
t613 = qJD(2) * t441 + qJD(3) * t442;
t428 = t613 * t416;
t440 = qJD(3) * t448;
t237 = qJD(1) * t428 + t440;
t544 = qJD(6) * t417;
t542 = qJD(1) * qJD(2);
t513 = t423 * t542;
t486 = t416 * t513;
t550 = qJD(2) * t421;
t520 = t416 * t550;
t487 = qJD(1) * t520;
t473 = t420 * t487;
t554 = -qJD(3) * t493 - t586 * t473;
t426 = qJD(3) * t429 + t486 * t610 + t554;
t471 = t415 * t487;
t154 = t282 * t547 + t340 * t546 - t419 * t471 - t422 * t426;
t431 = qJD(3) * t376 + t554;
t466 = qJD(3) * t480;
t515 = t610 * qJD(2);
t434 = t423 * (t466 + t515);
t521 = t415 * t550;
t491 = t422 * t521;
t548 = qJD(4) * t236;
t155 = t419 * t431 + (t419 * t434 - t491) * t551 + t548;
t77 = t609 * t154 + t418 * t155 - t235 * t517 + t236 * t545;
t55 = -qJD(6) * t433 - t417 * t237 + t458 * t544 + t608 * t77;
t506 = -t608 * t237 - t417 * t77;
t56 = qJD(6) * t143 + t506;
t636 = qJD(6) + t169;
t646 = t636 * t417;
t648 = -t141 * t649 - t143 * t646 - t417 * t56 - t55 * t608;
t504 = t418 * t154 - t609 * t155;
t78 = qJD(5) * t458 - t504;
t76 = t608 * t78;
t647 = t141 * t458 - t636 * t646 + t76;
t259 = -t419 * t356 + t422 * t357;
t233 = -pkin(4) * t529 - t368 * pkin(12) + t259;
t242 = -pkin(12) * t367 + t260;
t603 = t233 * t517 - t242 * t545 + t418 * t653 - t652 * t609;
t130 = -t191 * t419 + t422 * t223;
t109 = pkin(12) * t280 * t422 + pkin(4) * t282 + t130;
t392 = t612 * t419;
t393 = t612 * t422;
t456 = t392 * t609 + t418 * t393;
t590 = qJD(5) * t456 - t418 * t109 + t533 * t569 + t609 * t651;
t192 = t610 * t271 + (t279 * t586 + t331 * t415) * t420;
t481 = -t192 + (t547 + t576) * pkin(4);
t559 = -pkin(4) * t558 - t563;
t227 = -t279 * t415 + t586 * t331;
t160 = pkin(3) * t280 - pkin(11) * t282 + t227;
t166 = -t340 * pkin(11) + t192;
t104 = t160 * t419 + t166 * t422;
t93 = pkin(12) * t235 + t104;
t539 = t609 * t93;
t103 = t422 * t160 - t166 * t419;
t92 = -pkin(12) * t236 + t103;
t80 = pkin(4) * t276 + t92;
t44 = t418 * t80 + t539;
t42 = pkin(13) * t438 + t44;
t165 = t340 * pkin(3) - t191;
t126 = -t235 * pkin(4) + t165;
t73 = t169 * pkin(5) - pkin(13) * t458 + t126;
t16 = t417 * t73 + t42 * t608;
t595 = t418 * t93;
t43 = t609 * t80 - t595;
t41 = -pkin(5) * t438 - t43;
t390 = qJD(2) * t394;
t446 = qJD(2) * t452;
t299 = qJD(1) * t446 + t390;
t322 = t432 * qJD(2);
t300 = qJD(1) * t322;
t350 = qJD(2) * t451;
t341 = qJD(1) * t350;
t122 = -t271 * t549 + t279 * t466 + t610 * t299 + t300 * t511 + t331 * t489 + t341 * t572;
t120 = pkin(11) * t471 + t122;
t241 = -t415 * t300 + t586 * t341;
t134 = t237 * pkin(3) - pkin(11) * t426 + t241;
t48 = -qJD(4) * t104 - t120 * t419 + t422 * t134;
t28 = pkin(4) * t237 + pkin(12) * t154 + t48;
t47 = t422 * t120 + t419 * t134 + t160 * t546 - t166 * t547;
t36 = -pkin(12) * t155 + t47;
t512 = -t609 * t28 + t418 * t36;
t10 = -qJD(5) * t44 - t512;
t8 = -t237 * pkin(5) - t10;
t645 = t16 * t458 + t41 * t649 + t8 * t417;
t591 = t417 * t78 + t516 * t636;
t643 = -t143 * t458 + t636 * t638 + t591;
t642 = pkin(13) * t623 - t603;
t641 = t169 * t41;
t640 = pkin(5) * t555 + pkin(13) * t556 + t481;
t639 = -pkin(13) * t282 + t590;
t580 = t169 * t458;
t637 = pkin(5) * t561 + pkin(13) * t562 + t559;
t635 = -t169 ^ 2 + t458 ^ 2;
t108 = pkin(5) * t458 + pkin(13) * t169;
t634 = t169 * t438 - t77;
t507 = -t418 * t28 - t609 * t36 - t80 * t517 + t93 * t545;
t633 = t126 * t169 + t507;
t412 = t416 ^ 2;
t631 = -0.2e1 * t412 * t542;
t479 = qJD(3) * t511;
t461 = t271 * t518 + t279 * t479 + t420 * t299 - t300 * t480 + t331 * t519 - t341 * t529;
t121 = -pkin(3) * t471 + t461;
t84 = pkin(4) * t155 + t121;
t18 = pkin(5) * t78 + pkin(13) * t77 + t84;
t7 = pkin(13) * t237 - t507;
t3 = -qJD(6) * t16 + t608 * t18 - t417 * t7;
t629 = -t16 * t636 - t3;
t624 = t418 * t233 + t609 * t242;
t602 = -qJD(5) * t624 + t652 * t418 + t609 * t653;
t318 = t404 + t435;
t245 = -t318 * t415 + t586 * t345;
t508 = t587 * t415;
t468 = t610 * t508;
t568 = t420 * t421;
t309 = t416 * t568 - t455 - t468;
t483 = t420 * t508;
t310 = t483 + t436;
t189 = pkin(3) * t309 - pkin(11) * t310 + t245;
t301 = (t482 + t508) * pkin(10) + t553;
t211 = t610 * t301 + t318 * t511 + t345 * t572;
t366 = -t586 * t587 + t536;
t197 = -pkin(11) * t366 + t211;
t117 = t419 * t189 + t422 * t197;
t254 = t310 * t419 + t366 * t422;
t102 = -pkin(12) * t254 + t117;
t116 = t422 * t189 - t197 * t419;
t255 = t310 * t422 - t366 * t419;
t99 = pkin(4) * t309 - pkin(12) * t255 + t116;
t628 = t609 * t102 + t418 * t99;
t334 = t418 * t392 - t393 * t609;
t588 = qJD(5) * t334 + t109 * t609 + t418 * t651 - t488 * t612;
t627 = -t103 * t276 + t47;
t626 = -t104 * t276 - t48;
t625 = t636 * t458;
t560 = t415 * t527 - t650;
t462 = t417 * t42 - t608 * t73;
t618 = t41 * t544 + t458 * t462 - t608 * t8;
t616 = -t126 * t458 - t512;
t614 = t276 * t458 + t504;
t250 = qJD(3) * t468 + ((t467 - t568) * qJD(3) + t443 * qJD(2)) * t416;
t492 = t415 * t520;
t182 = -qJD(4) * t254 + t250 * t422 + t419 * t492;
t249 = qJD(3) * t483 + t428;
t395 = qJD(2) * t404;
t321 = t395 + t446;
t137 = -t301 * t549 + t318 * t466 + t610 * t321 + t322 * t511 + t345 * t489 + t350 * t572;
t128 = pkin(11) * t492 + t137;
t252 = -t322 * t415 + t586 * t350;
t146 = pkin(3) * t249 - pkin(11) * t250 + t252;
t62 = -qJD(4) * t117 - t128 * t419 + t422 * t146;
t40 = pkin(4) * t249 - pkin(12) * t182 + t62;
t181 = qJD(4) * t255 + t250 * t419 - t416 * t491;
t61 = t422 * t128 + t419 * t146 + t189 * t546 - t197 * t547;
t46 = -pkin(12) * t181 + t61;
t14 = -qJD(5) * t628 + t40 * t609 - t418 * t46;
t2 = -qJD(6) * t462 + t417 * t18 + t608 * t7;
t424 = qJD(1) ^ 2;
t164 = -pkin(13) * t529 + t624;
t277 = t367 * t609 + t368 * t418;
t355 = -pkin(3) * t586 - t370;
t286 = t367 * pkin(4) + t355;
t198 = t277 * pkin(5) - t278 * pkin(13) + t286;
t110 = -t417 * t164 + t198 * t608;
t606 = qJD(6) * t110 + t637 * t417 - t608 * t642;
t111 = t164 * t608 + t417 * t198;
t605 = -qJD(6) * t111 + t417 * t642 + t637 * t608;
t604 = pkin(5) * t623 - t602;
t601 = pkin(4) * qJD(5);
t600 = t462 * t417;
t199 = t254 * t609 + t255 * t418;
t598 = t199 * t78;
t597 = t277 * t78;
t596 = t379 * t78;
t52 = t55 * t417;
t410 = -pkin(4) * t422 - pkin(3);
t297 = pkin(5) * t379 - pkin(13) * t380 + t410;
t239 = t297 * t608 - t417 * t334;
t594 = qJD(6) * t239 + t417 * t640 + t608 * t639;
t240 = t417 * t297 + t334 * t608;
t593 = -qJD(6) * t240 - t417 * t639 + t608 * t640;
t139 = t143 * t516;
t592 = t139 - t52;
t589 = t282 * pkin(5) + t588;
t583 = t141 * t417;
t582 = t143 * t141;
t579 = t235 * t276;
t578 = t236 * t235;
t577 = t236 * t276;
t212 = t237 * t309;
t575 = t282 * t280;
t574 = t380 * t417;
t573 = t412 * t424;
t256 = t417 * t278 + t529 * t608;
t565 = qJD(6) * t256 + t417 * t623 + t562 * t608;
t495 = t417 * t529;
t564 = -qJD(6) * t495 + t278 * t516 - t417 * t562 + t608 * t623;
t552 = t421 ^ 2 - t423 ^ 2;
t1 = t2 * t608;
t54 = t56 * t608;
t538 = t608 * t462;
t537 = t423 * t573;
t531 = t380 * t608;
t408 = pkin(4) * t418 + pkin(13);
t530 = t408 * t608;
t528 = t610 * t237;
t526 = t610 * t350;
t524 = t608 * t143;
t503 = -t282 * t608 + t417 * t556;
t500 = t276 * t419;
t499 = t276 * t422;
t498 = pkin(4) * t517;
t497 = t421 * t537;
t496 = t415 * t528;
t490 = t609 * t608;
t49 = t418 * t92 + t539;
t484 = pkin(4) * t545 - t49;
t478 = t416 * t424 * t587;
t477 = pkin(1) * t631;
t476 = t417 * t282 + t556 * t608;
t472 = t412 * t421 * t513;
t469 = 0.2e1 * t505 + qJD(2);
t58 = pkin(13) * t309 + t628;
t210 = -t420 * t301 + t318 * t480 + t345 * t529;
t196 = t366 * pkin(3) - t210;
t147 = t254 * pkin(4) + t196;
t200 = -t418 * t254 + t255 * t609;
t83 = t199 * pkin(5) - t200 * pkin(13) + t147;
t22 = -t417 * t58 + t608 * t83;
t23 = t417 * t83 + t58 * t608;
t59 = -t418 * t102 + t609 * t99;
t173 = t233 * t609 - t418 * t242;
t158 = t200 * t608 + t309 * t417;
t13 = -t102 * t545 + t418 * t40 + t609 * t46 + t99 * t517;
t454 = -pkin(11) * t237 + t165 * t276;
t351 = -pkin(9) * t487 + t390;
t50 = t609 * t92 - t595;
t94 = pkin(4) * t236 + t108;
t21 = t417 * t94 + t50 * t608;
t449 = -t490 * t601 + t21;
t445 = -t380 * t544 - t476;
t444 = t380 * t516 - t503;
t365 = t553 * qJD(2);
t437 = -t301 * t518 - t318 * t479 - t420 * t321 + t322 * t480 - t345 * t519;
t430 = t1 - t3 * t417 + (-t16 * t417 + t538) * qJD(6);
t129 = (-pkin(3) * t520 - t526) * t415 - t437;
t95 = t181 * pkin(4) + t129;
t425 = t415 * t426;
t411 = t415 ^ 2;
t409 = -pkin(4) * t609 - pkin(5);
t371 = -t421 * t416 * pkin(9) + t404;
t364 = -pkin(9) * t520 + t395;
t359 = -pkin(9) * t522 + t394;
t352 = qJD(1) * t365;
t257 = t278 * t608 - t495;
t163 = pkin(5) * t529 - t173;
t157 = t200 * t417 - t309 * t608;
t138 = t415 * t526 + t437;
t89 = qJD(5) * t200 + t181 * t609 + t418 * t182;
t88 = t418 * t181 - t182 * t609 + t254 * t517 + t255 * t545;
t64 = qJD(6) * t158 - t249 * t608 - t417 * t88;
t63 = t200 * t544 - t249 * t417 - t309 * t516 + t608 * t88;
t57 = -t309 * pkin(5) - t59;
t25 = t417 * t108 + t43 * t608;
t24 = t108 * t608 - t417 * t43;
t20 = -t417 * t50 + t608 * t94;
t19 = t89 * pkin(5) + t88 * pkin(13) + t95;
t12 = -t249 * pkin(5) - t14;
t11 = pkin(13) * t249 + t13;
t5 = -qJD(6) * t23 - t417 * t11 + t19 * t608;
t4 = qJD(6) * t22 + t11 * t608 + t417 * t19;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t472, t552 * t631, t469 * qJD(2) * t571, -0.2e1 * t472, -t469 * t520, 0, -t352 * t587 - t365 * t470 + t421 * t477, -t351 * t587 - t364 * t470 + t423 * t477 (t351 * t423 + t352 * t421 + (-t359 * t423 - t361 * t421) * qJD(2) + (t364 * t423 + t365 * t421 + (-t371 * t423 - t421 * t553) * qJD(2)) * qJD(1)) * t416, t351 * t553 - t352 * t371 - t359 * t365 + t361 * t364, t282 * t250 + t310 * t426, -t310 * t237 - t282 * t249 - t250 * t280 - t309 * t426, -t250 * t340 - t431 * t366 + (t282 * t521 + (t310 * t521 - t366 * t434) * qJD(1)) * t416, t249 * t280 + t212, t237 * t366 + t249 * t340 + (-qJD(1) * t309 - t280) * t492 (-qJD(1) * t366 - t340) * t492, t461 * t366 - t138 * t340 + t227 * t249 + t237 * t245 + t241 * t309 + t252 * t280 + (qJD(1) * t210 + t191) * t492, t122 * t366 + t137 * t340 - t192 * t492 - t211 * t471 + t227 * t250 + t241 * t310 + t245 * t426 + t252 * t282, -t122 * t309 - t137 * t280 - t138 * t282 - t191 * t250 - t192 * t249 - t210 * t426 - t211 * t237 + t310 * t461, t122 * t211 + t137 * t192 + t138 * t191 - t210 * t461 + t227 * t252 + t241 * t245, -t154 * t255 + t182 * t236, t154 * t254 - t155 * t255 - t181 * t236 + t182 * t235, -t154 * t309 + t182 * t276 + t236 * t249 + t237 * t255, t155 * t254 - t181 * t235, -t155 * t309 - t181 * t276 + t235 * t249 - t237 * t254, t249 * t276 + t212, t103 * t249 + t116 * t237 + t121 * t254 - t129 * t235 + t155 * t196 + t165 * t181 + t276 * t62 + t309 * t48, -t104 * t249 - t117 * t237 + t121 * t255 + t129 * t236 - t154 * t196 + t165 * t182 - t276 * t61 - t309 * t47, -t103 * t182 - t104 * t181 + t116 * t154 - t117 * t155 + t235 * t61 - t236 * t62 - t254 * t47 - t255 * t48, t103 * t62 + t104 * t61 + t116 * t48 + t117 * t47 + t121 * t196 + t129 * t165, -t200 * t77 - t458 * t88, t169 * t88 + t199 * t77 - t200 * t78 - t458 * t89, t200 * t237 + t249 * t458 - t77 * t309 - t438 * t88, t169 * t89 + t598, -t169 * t249 - t199 * t237 - t78 * t309 - t438 * t89, t249 * t438 + t212, t10 * t309 + t126 * t89 + t14 * t438 + t147 * t78 + t95 * t169 + t84 * t199 + t59 * t237 + t43 * t249, -t126 * t88 - t13 * t438 - t147 * t77 + t84 * t200 - t237 * t628 - t44 * t249 + t309 * t507 + t458 * t95, -t10 * t200 - t13 * t169 - t14 * t458 + t199 * t507 + t43 * t88 - t44 * t89 + t59 * t77 - t628 * t78, t10 * t59 + t126 * t95 + t13 * t44 + t14 * t43 + t147 * t84 - t507 * t628, -t143 * t63 - t158 * t55, t141 * t63 - t143 * t64 + t157 * t55 - t158 * t56, t143 * t89 + t158 * t78 - t199 * t55 - t63 * t636, t141 * t64 + t157 * t56, -t141 * t89 - t157 * t78 - t199 * t56 - t636 * t64, t636 * t89 + t598, t12 * t141 + t157 * t8 + t199 * t3 + t22 * t78 + t41 * t64 - t462 * t89 + t5 * t636 + t56 * t57, t12 * t143 + t158 * t8 - t16 * t89 - t199 * t2 - t23 * t78 - t4 * t636 - t41 * t63 - t55 * t57, -t141 * t4 - t143 * t5 - t157 * t2 - t158 * t3 - t16 * t64 + t22 * t55 - t23 * t56 - t462 * t63, t12 * t41 + t16 * t4 + t2 * t23 + t22 * t3 - t462 * t5 + t57 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t497, t552 * t573, -t423 * t478, t497, t421 * t478, 0, -pkin(9) * t486 + t361 * t470 + (-qJD(2) * t505 + t573) * t421 * pkin(1), pkin(1) * t537 + t359 * t470 - t351, 0, 0, t282 * t621 + t420 * t425, -t237 * t572 - t280 * t621 + t282 * t623 + t610 * t425, -t282 * t494 - t340 * t621 + t411 * t473 + t426 * t586, -t280 * t623 - t496, -t237 * t586 - t623 * t340 + (t280 * t415 + t411 * t515) * t522 (qJD(2) * t586 + t340) * t494, -t461 * t586 - t227 * t337 - t251 * t280 + t560 * t340 + (t227 * t549 - t610 * t241 - pkin(2) * t237 + (qJD(2) * t370 - t191) * t522) * t415, -pkin(2) * t425 - t122 * t586 + t192 * t494 + t227 * t621 + t241 * t572 - t251 * t282 + t340 * t622 - t372 * t471, t122 * t529 - t191 * t621 + t192 * t623 - t372 * t237 - t280 * t622 + t282 * t560 - t370 * t426 + t461 * t572, -pkin(2) * t241 * t415 + t122 * t372 - t560 * t191 + t192 * t622 - t227 * t251 - t370 * t461, -t154 * t368 - t236 * t557, t154 * t367 - t155 * t368 - t235 * t557 + t236 * t558, t154 * t529 - t236 * t623 + t368 * t237 - t276 * t557, t155 * t367 + t235 * t558, t155 * t529 - t235 * t623 - t367 * t237 + t276 * t558, -t276 * t623 - t496, -t103 * t623 + t121 * t367 + t355 * t155 - t165 * t558 + t235 * t563 + t259 * t237 - t276 * t567 - t48 * t529, t104 * t623 + t121 * t368 - t355 * t154 - t165 * t557 - t236 * t563 - t260 * t237 + t276 * t566 + t47 * t529, t103 * t557 + t104 * t558 + t154 * t259 - t155 * t260 - t235 * t566 + t236 * t567 - t367 * t47 - t368 * t48, -t103 * t567 - t104 * t566 + t121 * t355 - t165 * t563 + t259 * t48 + t260 * t47, -t278 * t77 - t458 * t562, t169 * t562 + t277 * t77 - t278 * t78 - t458 * t561, t278 * t237 - t438 * t562 - t458 * t623 + t529 * t77, t169 * t561 + t597, t169 * t623 - t277 * t237 - t438 * t561 + t529 * t78, -t438 * t337 + (t438 * t549 - t528) * t415, t173 * t237 + t286 * t78 + t84 * t277 - t43 * t337 + (-t10 * t610 + t43 * t549) * t415 + t559 * t169 + t561 * t126 + t602 * t438, -t624 * t237 - t286 * t77 + t84 * t278 + t44 * t337 + (-t44 * t549 - t507 * t610) * t415 + t559 * t458 - t562 * t126 - t603 * t438, -t10 * t278 - t169 * t603 + t173 * t77 + t277 * t507 + t43 * t562 - t44 * t561 - t458 * t602 - t624 * t78, t10 * t173 + t126 * t559 + t286 * t84 + t43 * t602 + t44 * t603 - t507 * t624, -t143 * t565 - t257 * t55, t141 * t565 - t143 * t564 + t256 * t55 - t257 * t56, t143 * t561 + t257 * t78 - t277 * t55 - t565 * t636, t141 * t564 + t256 * t56, -t141 * t561 - t256 * t78 - t277 * t56 - t564 * t636, t561 * t636 + t597, t110 * t78 + t141 * t604 + t163 * t56 + t256 * t8 + t277 * t3 + t41 * t564 - t462 * t561 + t605 * t636, -t111 * t78 + t143 * t604 - t16 * t561 - t163 * t55 - t2 * t277 + t257 * t8 - t41 * t565 - t606 * t636, t110 * t55 - t111 * t56 - t141 * t606 - t143 * t605 - t16 * t564 - t2 * t256 - t257 * t3 - t462 * t565, t110 * t3 + t111 * t2 + t16 * t606 + t163 * t8 + t41 * t604 - t462 * t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, -t280 ^ 2 + t282 ^ 2, -t280 * t340 + t426, -t575, -t282 * t340 - t551 * t613 - t440, t471, -t192 * t340 - t227 * t282 - t461, -t191 * t340 + t227 * t280 - t122, 0, 0, -t154 * t419 + t236 * t499 (-t154 + t579) * t422 + (-t155 - t577) * t419, -t236 * t282 + t237 * t419 + t276 * t499, -t155 * t422 - t235 * t500, -t235 * t282 + t237 * t422 - t276 * t500, -t276 * t282, -pkin(3) * t155 - t103 * t282 - t121 * t422 + t192 * t235 + (-pkin(11) * t546 - t130) * t276 + t454 * t419, pkin(3) * t154 + t104 * t282 + t121 * t419 - t192 * t236 + (pkin(11) * t547 + t131) * t276 + t454 * t422, t130 * t236 - t131 * t235 + ((-t155 + t548) * pkin(11) + t627) * t422 + ((-qJD(4) * t235 - t154) * pkin(11) + t626) * t419, -pkin(3) * t121 - t103 * t130 - t104 * t131 - t165 * t192 + (-t419 * t48 + t422 * t47 + (-t103 * t422 - t104 * t419) * qJD(4)) * pkin(11), -t380 * t77 - t458 * t556, t169 * t556 + t379 * t77 - t380 * t78 - t458 * t555, t380 * t237 - t458 * t282 - t438 * t556, t169 * t555 + t596, t169 * t282 - t379 * t237 - t438 * t555, -t438 * t282, t126 * t555 + t169 * t481 + t237 * t456 - t43 * t282 + t84 * t379 + t410 * t78 - t438 * t588, -t126 * t556 - t334 * t237 + t44 * t282 + t84 * t380 - t410 * t77 - t438 * t590 + t458 * t481, -t10 * t380 - t169 * t590 - t334 * t78 + t379 * t507 + t43 * t556 - t44 * t555 + t456 * t77 + t458 * t588, t10 * t456 + t126 * t481 - t334 * t507 + t410 * t84 - t43 * t588 + t44 * t590, t143 * t445 - t531 * t55, t503 * t143 + t476 * t141 + (-t54 + t52 + (-t524 + t583) * qJD(6)) * t380, t143 * t555 - t55 * t379 + t445 * t636 + t531 * t78, t141 * t444 + t56 * t574, -t141 * t555 - t56 * t379 - t444 * t636 - t574 * t78, t555 * t636 + t596, t141 * t589 + t239 * t78 + t3 * t379 + t41 * t444 - t456 * t56 - t462 * t555 + t574 * t8 + t593 * t636, t143 * t589 - t16 * t555 - t2 * t379 - t240 * t78 + t41 * t445 + t456 * t55 + t531 * t8 - t594 * t636, t239 * t55 - t240 * t56 + t503 * t16 - t476 * t462 - t593 * t143 - t594 * t141 + (-t608 * t3 - t2 * t417 + (-t16 * t608 - t600) * qJD(6)) * t380, t16 * t594 + t2 * t240 + t239 * t3 + t41 * t589 - t456 * t8 - t462 * t593; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t578, -t235 ^ 2 + t236 ^ 2, -t154 - t579, t578, -t282 * t546 + t340 * t547 - t419 * t426 + t422 * t471 + t577, t237, -t165 * t236 - t626, -t165 * t235 - t627, 0, 0, t580, t635, t634, -t580, t614, t237, t49 * t276 + (t49 - t44) * qJD(5) + (-t236 * t169 + t237 * t609 - t438 * t545) * pkin(4) + t616, t50 * t438 + (-t236 * t458 - t418 * t237 - t438 * t517) * pkin(4) + t633, t44 * t458 + t50 * t169 - t43 * t169 - t49 * t458 + (t609 * t77 - t418 * t78 + (-t169 * t609 + t418 * t458) * qJD(5)) * pkin(4), t43 * t49 - t44 * t50 + (t609 * t10 - t126 * t236 - t418 * t507 + (-t418 * t43 + t44 * t609) * qJD(5)) * pkin(4), t169 * t524 + t592, t648, t643, t583 * t636 - t54, t647, -t625, t409 * t56 + (-t408 * t78 + t641) * t417 + t484 * t141 + (-t408 * t516 - t417 * t498 - t20) * t636 + t618, -t78 * t530 - t409 * t55 + t484 * t143 + (t408 * t544 + t449) * t636 + t645, t462 * t638 - t56 * t530 + t20 * t143 + t1 + t449 * t141 + (t408 * t524 + t538) * qJD(6) + (t143 * t498 - t16 * t169 - t408 * t55 - t3 + (t141 * t408 - t16) * qJD(6)) * t417, t462 * t20 - t16 * t21 + t8 * t409 - t41 * t49 + (t16 * t490 + t41 * t418 + t600 * t609) * t601 + t430 * t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t580, t635, t634, -t580, t614, t237, t276 * t44 + t616, t43 * t438 + t633, 0, 0, t143 * t638 + t592, t648, t643, t141 * t646 - t54, t647, -t625, -pkin(5) * t56 - pkin(13) * t591 - t44 * t141 - t24 * t636 + t417 * t641 + t618, pkin(5) * t55 - t44 * t143 + t25 * t636 + (t544 * t636 - t76) * pkin(13) + t645, t25 * t141 + t24 * t143 + t1 + t649 * t462 + (t139 - t54) * pkin(13) + ((qJD(6) * t141 - t55) * pkin(13) + t629) * t417, -t8 * pkin(5) + pkin(13) * t430 - t16 * t25 + t24 * t462 - t41 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, -t141 ^ 2 + t143 ^ 2, t141 * t636 - t55, -t582, -t506 + (-qJD(6) + t636) * t143, t78, -t41 * t143 - t629, t41 * t141 - t462 * t636 - t2, 0, 0;];
tauc_reg  = t6;
