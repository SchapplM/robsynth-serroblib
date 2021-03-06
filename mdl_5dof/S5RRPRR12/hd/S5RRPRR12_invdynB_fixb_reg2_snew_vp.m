% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR12_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:46
% EndTime: 2019-12-31 20:30:59
% DurationCPUTime: 9.86s
% Computational Cost: add. (26545->538), mult. (55549->757), div. (0->0), fcn. (34781->8), ass. (0->380)
t604 = qJD(2) ^ 2;
t598 = sin(qJ(2));
t594 = t598 ^ 2;
t605 = qJD(1) ^ 2;
t679 = t594 * t605;
t571 = t604 + t679;
t602 = cos(qJ(2));
t657 = t602 * t605;
t576 = t598 * t657;
t567 = qJDD(2) - t576;
t658 = t602 * t567;
t520 = -t598 * t571 + t658;
t647 = qJD(1) * qJD(2);
t634 = t602 * t647;
t645 = t598 * qJDD(1);
t556 = 0.2e1 * t634 + t645;
t599 = sin(qJ(1));
t603 = cos(qJ(1));
t478 = t599 * t520 + t603 * t556;
t713 = pkin(5) * t478;
t481 = t603 * t520 - t599 * t556;
t712 = pkin(5) * t481;
t711 = pkin(6) * t520;
t597 = sin(qJ(4));
t601 = cos(qJ(4));
t540 = (t597 * t598 + t601 * t602) * qJD(1);
t557 = t634 + t645;
t636 = t598 * t647;
t643 = t602 * qJDD(1);
t558 = -t636 + t643;
t467 = -t540 * qJD(4) + t601 * t557 - t597 * t558;
t591 = qJD(2) - qJD(4);
t683 = t540 * t591;
t443 = t467 + t683;
t696 = 2 * qJD(3);
t669 = t598 * t567;
t514 = t602 * t571 + t669;
t710 = pkin(1) * t514;
t709 = pkin(6) * t514;
t596 = sin(qJ(5));
t649 = qJD(1) * t602;
t650 = qJD(1) * t598;
t542 = -t597 * t649 + t601 * t650;
t600 = cos(qJ(5));
t504 = t596 * t542 + t600 * t591;
t506 = t600 * t542 - t596 * t591;
t459 = t506 * t504;
t627 = t597 * t557 + t601 * t558;
t466 = -t542 * qJD(4) - t627;
t612 = qJDD(5) - t466;
t699 = -t459 + t612;
t708 = t596 * t699;
t707 = t600 * t699;
t559 = -0.2e1 * t636 + t643;
t659 = t602 * t559;
t672 = t598 * t556;
t496 = -t659 + t672;
t595 = t602 ^ 2;
t565 = (t594 - t595) * t605;
t706 = t599 * t496 + t603 * t565;
t705 = t603 * t496 - t599 * t565;
t678 = t595 * t605;
t573 = -t604 + t678;
t518 = -t602 * t573 + t669;
t642 = t603 * qJDD(1);
t704 = t599 * t518 + t602 * t642;
t703 = t603 * t518 - t599 * t643;
t590 = -qJDD(2) + qJDD(4);
t682 = t542 * t540;
t616 = t590 - t682;
t702 = t597 * t616;
t570 = t603 * g(1) + t599 * g(2);
t544 = -t605 * pkin(1) + qJDD(1) * pkin(6) - t570;
t693 = pkin(2) * t602;
t621 = -qJ(3) * t598 - t693;
t554 = t621 * qJD(1);
t626 = qJD(1) * t554 + t544;
t701 = t598 * t626;
t700 = t601 * t616;
t536 = qJD(5) + t540;
t476 = t536 * t504;
t637 = t504 * qJD(5) - t600 * t467 - t596 * t590;
t406 = t637 + t476;
t631 = t598 * g(3) - t602 * t544;
t613 = qJDD(2) * qJ(3) + qJD(2) * t696 + t554 * t649 - t631;
t698 = t598 * t573 + t658;
t628 = t596 * t467 - t600 * t590;
t401 = (qJD(5) - t536) * t506 + t628;
t501 = t504 ^ 2;
t502 = t506 ^ 2;
t534 = t536 ^ 2;
t538 = t540 ^ 2;
t539 = t542 ^ 2;
t697 = t591 ^ 2;
t695 = pkin(2) + pkin(3);
t566 = qJDD(2) + t576;
t549 = t602 * t566;
t574 = -t604 - t678;
t512 = t598 * t574 + t549;
t694 = pkin(1) * t512;
t670 = t598 * t566;
t517 = t602 * t574 - t670;
t477 = t599 * t517 + t603 * t559;
t692 = pkin(5) * t477;
t651 = t594 + t595;
t561 = t651 * qJDD(1);
t564 = t651 * t605;
t499 = t599 * t561 + t603 * t564;
t691 = pkin(5) * t499;
t690 = pkin(6) * t512;
t689 = t558 * pkin(2);
t688 = t602 * g(3);
t687 = qJ(3) * t602;
t686 = qJDD(2) * pkin(2);
t685 = t536 * t596;
t684 = t536 * t600;
t681 = t591 * t597;
t680 = t591 * t601;
t461 = -t604 * pkin(2) + t613;
t568 = -qJD(2) * pkin(3) - pkin(7) * t650;
t427 = -pkin(3) * t678 - t558 * pkin(7) + qJD(2) * t568 + t461;
t630 = qJDD(3) + t688;
t614 = t604 * qJ(3) - t630;
t606 = -t695 * qJDD(2) + (-t557 + t634) * pkin(7) + (-pkin(3) * t657 + t626) * t598 - t614;
t377 = t597 * t427 - t601 * t606;
t493 = t540 * pkin(4) - t542 * pkin(8);
t362 = -t590 * pkin(4) - t697 * pkin(8) + t542 * t493 + t377;
t677 = t596 * t362;
t413 = t459 + t612;
t676 = t596 * t413;
t569 = t599 * g(1) - t603 * g(2);
t543 = qJDD(1) * pkin(1) + t605 * pkin(6) + t569;
t610 = -pkin(2) * t636 + t543;
t417 = t689 + t557 * qJ(3) + t558 * pkin(3) - pkin(7) * t678 + (qJD(2) * t687 + (t696 + t568) * t598) * qJD(1) + t610;
t675 = t597 * t417;
t490 = -t590 - t682;
t674 = t597 * t490;
t673 = t598 * t543;
t671 = t598 * t559;
t664 = t600 * t362;
t663 = t600 * t413;
t662 = t601 * t417;
t661 = t601 * t490;
t660 = t602 * t543;
t360 = -t443 * pkin(8) + (-t591 * t542 - t466) * pkin(4) + t417;
t378 = t601 * t427 + t597 * t606;
t363 = -t697 * pkin(4) + t590 * pkin(8) - t540 * t493 + t378;
t321 = t596 * t360 + t600 * t363;
t654 = pkin(1) * t559 + pkin(6) * t517;
t653 = pkin(1) * t564 + pkin(6) * t561;
t652 = t564 - t604;
t644 = t599 * qJDD(1);
t641 = t597 * t459;
t640 = t599 * t682;
t639 = t601 * t459;
t638 = t603 * t682;
t632 = pkin(4) * t597 + qJ(3);
t320 = -t600 * t360 + t596 * t363;
t291 = t596 * t320 + t600 * t321;
t523 = t598 * t544 + t688;
t458 = t598 * t523 - t602 * t631;
t509 = -t599 * t569 - t603 * t570;
t625 = pkin(4) * t601 + t695;
t624 = t599 * t576;
t623 = t603 * t576;
t563 = -t599 * t605 + t642;
t622 = -pkin(5) * t563 - t599 * g(3);
t620 = pkin(2) * t598 - t687;
t619 = t557 + t634;
t290 = -t600 * t320 + t596 * t321;
t331 = -t601 * t377 + t597 * t378;
t332 = t597 * t377 + t601 * t378;
t457 = t602 * t523 + t598 * t631;
t618 = t602 * t556 + t671;
t508 = t603 * t569 - t599 * t570;
t611 = (-qJD(4) - t591) * t542 - t627;
t609 = t630 + t701;
t608 = t650 * t696 + t610;
t607 = qJ(3) * t619 + t608;
t572 = t604 - t679;
t562 = t603 * t605 + t644;
t552 = t620 * qJDD(1);
t548 = t651 * t647;
t537 = -pkin(5) * t562 + t603 * g(3);
t530 = -t539 + t697;
t529 = t538 - t697;
t528 = t599 * qJDD(2) + t603 * t548;
t527 = t602 * t557 - t594 * t647;
t526 = -t603 * qJDD(2) + t599 * t548;
t525 = -t598 * t558 - t595 * t647;
t522 = -t539 - t697;
t519 = -t598 * t572 + t549;
t513 = t602 * t572 + t670;
t511 = t619 * t598;
t510 = (t558 - t636) * t602;
t500 = t603 * t561 - t599 * t564;
t497 = pkin(5) * t500;
t494 = t539 - t538;
t489 = t603 * t527 - t624;
t488 = t603 * t525 + t624;
t487 = t599 * t527 + t623;
t486 = t599 * t525 - t623;
t485 = t603 * t519 + t598 * t644;
t484 = t599 * t519 - t598 * t642;
t483 = -t697 - t538;
t480 = t603 * t517 - t599 * t559;
t475 = pkin(5) * t480;
t474 = -t502 + t534;
t473 = t501 - t534;
t472 = (t540 * t601 - t542 * t597) * t591;
t471 = (-t540 * t597 - t542 * t601) * t591;
t470 = -t660 + t709;
t469 = -t673 - t690;
t468 = -t538 - t539;
t465 = t614 + t686 - t701;
t464 = -t631 + t710;
t463 = t523 - t694;
t455 = -t502 + t501;
t454 = t652 * qJ(3) + t609 - t686;
t453 = -t502 - t534;
t452 = t652 * pkin(2) + t613;
t451 = t607 + t689;
t450 = t601 * t529 + t674;
t449 = -t597 * t530 + t700;
t448 = -t597 * t529 + t661;
t447 = -t601 * t530 - t702;
t446 = -t597 * t522 + t661;
t445 = t601 * t522 + t674;
t444 = t467 - t683;
t439 = (qJD(4) - t591) * t542 + t627;
t438 = -t534 - t501;
t437 = t601 * t467 + t542 * t681;
t436 = -t597 * t467 + t542 * t680;
t435 = -t597 * t466 - t540 * t680;
t434 = -t601 * t466 + t540 * t681;
t433 = (t558 + t559) * pkin(2) + t607;
t432 = t689 + (t556 + t619) * qJ(3) + t608;
t431 = t603 * t458 - t599 * t543;
t430 = t599 * t458 + t603 * t543;
t429 = t601 * t483 - t702;
t428 = t597 * t483 + t700;
t426 = t501 + t502;
t421 = -t506 * qJD(5) - t628;
t420 = (-t504 * t600 + t506 * t596) * t536;
t419 = (-t504 * t596 - t506 * t600) * t536;
t418 = -t694 + (-t574 - t604) * qJ(3) + (-qJDD(2) - t566) * pkin(2) + t609;
t416 = -t710 - qJ(3) * t567 + (-t571 + t604) * pkin(2) - t613;
t415 = -t598 * t471 + t602 * t472;
t411 = t602 * t461 - t598 * t465;
t410 = t598 * t461 + t602 * t465;
t409 = -pkin(2) * t672 + t602 * t432 - t709;
t408 = qJ(3) * t659 - t598 * t433 - t690;
t407 = -t598 * t452 + t602 * t454;
t405 = -t476 + t637;
t402 = (-qJD(5) - t536) * t506 - t628;
t400 = -t598 * t448 + t602 * t450;
t399 = -t598 * t447 + t602 * t449;
t398 = -t506 * t685 - t600 * t637;
t397 = t506 * t684 - t596 * t637;
t396 = -t596 * t421 + t504 * t684;
t395 = t600 * t421 + t504 * t685;
t394 = t598 * t445 + t602 * t446;
t393 = -t602 * t445 + t598 * t446;
t392 = t597 * t444 + t601 * t611;
t391 = -t601 * t439 - t597 * t443;
t390 = -t601 * t444 + t597 * t611;
t389 = t597 * t439 - t601 * t443;
t388 = t601 * t420 + t597 * t612;
t387 = -t597 * t420 + t601 * t612;
t386 = t600 * t473 - t676;
t385 = -t596 * t474 + t707;
t384 = t596 * t473 + t663;
t383 = t600 * t474 + t708;
t382 = -t598 * t436 + t602 * t437;
t381 = -t598 * t434 + t602 * t435;
t380 = t598 * t428 + t602 * t429;
t379 = -t602 * t428 + t598 * t429;
t376 = -t596 * t453 - t663;
t375 = t600 * t453 - t676;
t374 = t603 * t411 - t599 * t451;
t373 = t599 * t411 + t603 * t451;
t372 = t600 * t438 - t708;
t371 = t596 * t438 + t707;
t370 = t601 * t398 + t641;
t369 = t601 * t396 - t641;
t368 = -t597 * t398 + t639;
t367 = -t597 * t396 - t639;
t366 = -pkin(1) * t410 - pkin(2) * t465 - qJ(3) * t461;
t365 = t603 * t394 - t443 * t599;
t364 = t599 * t394 + t443 * t603;
t359 = -pkin(7) * t445 + qJ(3) * t443 + t662;
t356 = t603 * t380 - t599 * t439;
t355 = t599 * t380 + t603 * t439;
t354 = -pkin(6) * t410 - t451 * t620;
t353 = -pkin(7) * t428 + qJ(3) * t439 + t675;
t352 = -t401 * t600 - t596 * t405;
t351 = t600 * t402 + t406 * t596;
t350 = -t401 * t596 + t600 * t405;
t349 = t596 * t402 - t406 * t600;
t348 = t598 * t390 + t602 * t392;
t347 = -t598 * t389 + t602 * t391;
t346 = -t602 * t390 + t598 * t392;
t345 = t601 * t386 - t597 * t401;
t344 = t601 * t385 - t597 * t405;
t343 = -t597 * t386 - t601 * t401;
t342 = -t597 * t385 - t601 * t405;
t341 = -pkin(7) * t446 + t443 * t695 - t675;
t340 = -t598 * t387 + t602 * t388;
t339 = t601 * t376 - t597 * t406;
t338 = t597 * t376 + t601 * t406;
t337 = -pkin(7) * t429 + t695 * t439 + t662;
t336 = t601 * t372 - t597 * t402;
t335 = t597 * t372 + t601 * t402;
t334 = t601 * t351 - t597 * t455;
t333 = -t597 * t351 - t601 * t455;
t330 = t603 * t348 - t599 * t468;
t329 = t599 * t348 + t603 * t468;
t328 = t601 * t352 - t597 * t426;
t327 = t597 * t352 + t601 * t426;
t326 = -t598 * t368 + t602 * t370;
t325 = -t598 * t367 + t602 * t369;
t324 = -pkin(8) * t375 + t664;
t323 = -pkin(8) * t371 + t677;
t322 = -pkin(7) * t331 + qJ(3) * t417;
t318 = -pkin(1) * t393 - qJ(3) * t446 + t695 * t445 - t378;
t317 = -pkin(7) * t390 + qJ(3) * t468 - t331;
t316 = -pkin(7) * t332 + t695 * t417;
t315 = -pkin(1) * t379 - qJ(3) * t429 + t695 * t428 - t377;
t314 = -pkin(7) * t392 + t695 * t468 - t332;
t313 = -t598 * t343 + t602 * t345;
t312 = -t598 * t342 + t602 * t344;
t311 = t598 * t338 + t602 * t339;
t310 = -t602 * t338 + t598 * t339;
t309 = -pkin(4) * t375 + t321;
t308 = t598 * t335 + t602 * t336;
t307 = -t602 * t335 + t598 * t336;
t306 = -pkin(4) * t371 + t320;
t305 = -t598 * t333 + t602 * t334;
t304 = t598 * t331 + t602 * t332;
t303 = -t602 * t331 + t598 * t332;
t302 = -pkin(6) * t393 - t598 * t341 + t602 * t359;
t301 = t598 * t327 + t602 * t328;
t300 = -t602 * t327 + t598 * t328;
t299 = -pkin(1) * t346 - qJ(3) * t392 + t695 * t390;
t298 = -pkin(6) * t379 - t598 * t337 + t602 * t353;
t297 = t603 * t304 - t599 * t417;
t296 = t599 * t304 + t603 * t417;
t295 = t603 * t311 - t599 * t375;
t294 = t599 * t311 + t603 * t375;
t293 = t603 * t308 - t599 * t371;
t292 = t599 * t308 + t603 * t371;
t289 = t603 * t301 - t599 * t350;
t288 = t599 * t301 + t603 * t350;
t287 = t601 * t291 + t597 * t362;
t286 = t597 * t291 - t601 * t362;
t285 = -pkin(8) * t350 - t290;
t284 = -pkin(6) * t346 - t598 * t314 + t602 * t317;
t283 = -pkin(7) * t338 + qJ(3) * t375 - t597 * t309 + t601 * t324;
t282 = -pkin(7) * t335 + qJ(3) * t371 - t597 * t306 + t601 * t323;
t281 = -pkin(6) * t303 - t598 * t316 + t602 * t322;
t280 = -pkin(7) * t339 - t601 * t309 - t597 * t324 + t695 * t375;
t279 = -pkin(1) * t303 - qJ(3) * t332 + t695 * t331;
t278 = -pkin(7) * t336 - t601 * t306 - t597 * t323 + t695 * t371;
t277 = -pkin(1) * t310 + pkin(4) * t406 + pkin(8) * t376 - qJ(3) * t339 + t695 * t338 + t677;
t276 = -pkin(7) * t327 + t601 * t285 + t350 * t632;
t275 = -pkin(1) * t307 + pkin(4) * t402 + pkin(8) * t372 - qJ(3) * t336 + t695 * t335 - t664;
t274 = t598 * t286 + t602 * t287;
t273 = -t602 * t286 + t598 * t287;
t272 = -pkin(7) * t328 - t597 * t285 + t350 * t625;
t271 = -pkin(1) * t300 + pkin(4) * t426 + pkin(8) * t352 - qJ(3) * t328 + t695 * t327 + t291;
t270 = t603 * t274 - t599 * t290;
t269 = t599 * t274 + t603 * t290;
t268 = -pkin(6) * t310 - t598 * t280 + t602 * t283;
t267 = -pkin(6) * t307 - t598 * t278 + t602 * t282;
t266 = -pkin(7) * t286 + (-pkin(8) * t601 + t632) * t290;
t265 = -pkin(7) * t287 + (pkin(8) * t597 + t625) * t290;
t264 = -pkin(6) * t300 - t598 * t272 + t602 * t276;
t263 = -pkin(1) * t273 - pkin(4) * t362 + pkin(8) * t291 - qJ(3) * t287 + t695 * t286;
t262 = -pkin(6) * t273 - t598 * t265 + t602 * t266;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t562, -t563, 0, t509, 0, 0, 0, 0, 0, 0, t480, -t481, t500, t431, 0, 0, 0, 0, 0, 0, t480, t500, t481, t374, 0, 0, 0, 0, 0, 0, t356, t365, t330, t297, 0, 0, 0, 0, 0, 0, t293, t295, t289, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t563, -t562, 0, t508, 0, 0, 0, 0, 0, 0, t477, -t478, t499, t430, 0, 0, 0, 0, 0, 0, t477, t499, t478, t373, 0, 0, 0, 0, 0, 0, t355, t364, t329, t296, 0, 0, 0, 0, 0, 0, t292, t294, t288, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t512, -t514, 0, -t457, 0, 0, 0, 0, 0, 0, t512, 0, t514, t410, 0, 0, 0, 0, 0, 0, t379, t393, t346, t303, 0, 0, 0, 0, 0, 0, t307, t310, t300, t273; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t563, 0, -t562, 0, t622, -t537, -t508, -pkin(5) * t508, t489, -t705, t485, t488, -t703, t528, -t599 * t463 + t603 * t469 - t692, -t599 * t464 + t603 * t470 + t713, t603 * t457 - t691, -pkin(5) * t430 - (pkin(1) * t599 - pkin(6) * t603) * t457, t489, t485, t705, t528, t703, t488, t603 * t408 - t599 * t418 - t692, t603 * t407 - t599 * t552 - t691, t603 * t409 - t599 * t416 - t713, -pkin(5) * t373 + t603 * t354 - t599 * t366, t603 * t382 - t640, t603 * t347 - t599 * t494, t603 * t399 - t599 * t444, t603 * t381 + t640, t603 * t400 - t599 * t611, t603 * t415 - t599 * t590, -pkin(5) * t355 + t603 * t298 - t599 * t315, -pkin(5) * t364 + t603 * t302 - t599 * t318, -pkin(5) * t329 + t603 * t284 - t599 * t299, -pkin(5) * t296 - t599 * t279 + t603 * t281, t603 * t326 - t599 * t397, t603 * t305 - t599 * t349, t603 * t312 - t599 * t383, t603 * t325 - t599 * t395, t603 * t313 - t599 * t384, t603 * t340 - t599 * t419, -pkin(5) * t292 + t603 * t267 - t599 * t275, -pkin(5) * t294 + t603 * t268 - t599 * t277, -pkin(5) * t288 + t603 * t264 - t599 * t271, -pkin(5) * t269 + t603 * t262 - t599 * t263; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t562, 0, t563, 0, t537, t622, t509, pkin(5) * t509, t487, -t706, t484, t486, -t704, t526, t603 * t463 + t599 * t469 + t475, t603 * t464 + t599 * t470 - t712, t599 * t457 + t497, pkin(5) * t431 - (-pkin(1) * t603 - pkin(6) * t599) * t457, t487, t484, t706, t526, t704, t486, t599 * t408 + t603 * t418 + t475, t599 * t407 + t603 * t552 + t497, t599 * t409 + t603 * t416 + t712, pkin(5) * t374 + t599 * t354 + t603 * t366, t599 * t382 + t638, t599 * t347 + t603 * t494, t599 * t399 + t603 * t444, t599 * t381 - t638, t599 * t400 + t603 * t611, t599 * t415 + t603 * t590, pkin(5) * t356 + t599 * t298 + t603 * t315, pkin(5) * t365 + t599 * t302 + t603 * t318, pkin(5) * t330 + t599 * t284 + t603 * t299, pkin(5) * t297 + t603 * t279 + t599 * t281, t599 * t326 + t603 * t397, t599 * t305 + t603 * t349, t599 * t312 + t603 * t383, t599 * t325 + t603 * t395, t599 * t313 + t603 * t384, t599 * t340 + t603 * t419, pkin(5) * t293 + t599 * t267 + t603 * t275, pkin(5) * t295 + t599 * t268 + t603 * t277, pkin(5) * t289 + t599 * t264 + t603 * t271, pkin(5) * t270 + t599 * t262 + t603 * t263; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t569, t570, 0, 0, t511, t618, t513, t510, t698, 0, t654 + t660, -pkin(1) * t556 - t673 - t711, t458 + t653, pkin(1) * t543 + pkin(6) * t458, t511, t513, -t618, 0, -t698, t510, qJ(3) * t671 + t602 * t433 + t654, t602 * t452 + t598 * t454 + t653, t711 + t598 * t432 + (pkin(1) + t693) * t556, pkin(6) * t411 + (pkin(1) - t621) * t451, t602 * t436 + t598 * t437, t602 * t389 + t598 * t391, t602 * t447 + t598 * t449, t602 * t434 + t598 * t435, t602 * t448 + t598 * t450, t602 * t471 + t598 * t472, pkin(1) * t439 + pkin(6) * t380 + t602 * t337 + t598 * t353, pkin(1) * t443 + pkin(6) * t394 + t602 * t341 + t598 * t359, pkin(1) * t468 + pkin(6) * t348 + t602 * t314 + t598 * t317, pkin(1) * t417 + pkin(6) * t304 + t602 * t316 + t598 * t322, t602 * t368 + t598 * t370, t602 * t333 + t598 * t334, t602 * t342 + t598 * t344, t602 * t367 + t598 * t369, t602 * t343 + t598 * t345, t602 * t387 + t598 * t388, pkin(1) * t371 + pkin(6) * t308 + t602 * t278 + t598 * t282, pkin(1) * t375 + pkin(6) * t311 + t602 * t280 + t598 * t283, pkin(1) * t350 + pkin(6) * t301 + t602 * t272 + t598 * t276, pkin(1) * t290 + pkin(6) * t274 + t602 * t265 + t598 * t266;];
tauB_reg = t1;
