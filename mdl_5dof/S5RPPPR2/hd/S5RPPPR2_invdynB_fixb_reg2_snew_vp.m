% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:23:13
% EndTime: 2020-01-03 11:23:30
% DurationCPUTime: 17.13s
% Computational Cost: add. (36714->631), mult. (105567->979), div. (0->0), fcn. (74317->10), ass. (0->449)
t651 = sin(pkin(9));
t652 = sin(pkin(8));
t653 = sin(pkin(7));
t710 = qJDD(1) * t653;
t690 = t652 * t710;
t654 = cos(pkin(9));
t655 = cos(pkin(8));
t656 = cos(pkin(7));
t593 = (t651 * t653 * t655 + t654 * t656) * qJD(1);
t713 = t653 * qJD(1);
t595 = -t651 * t656 * qJD(1) + t654 * t655 * t713;
t734 = t595 * t593;
t749 = t690 - t734;
t751 = t651 * t749;
t750 = t654 * t749;
t657 = sin(qJ(5));
t659 = cos(qJ(5));
t693 = t652 * t713;
t554 = t657 * t595 - t659 * t693;
t556 = t659 * t595 + t657 * t693;
t504 = t556 * t554;
t709 = qJDD(1) * t655;
t689 = t653 * t709;
t708 = t656 * qJDD(1);
t714 = t651 * t689 + t654 * t708;
t686 = qJDD(5) + t714;
t745 = -t504 + t686;
t748 = t657 * t745;
t747 = t659 * t745;
t658 = sin(qJ(1));
t660 = cos(qJ(1));
t626 = t658 * g(2) - t660 * g(3);
t661 = qJD(1) ^ 2;
t668 = -t661 * pkin(1) + qJDD(1) * qJ(2) - t626;
t744 = 2 * qJD(2);
t746 = qJD(1) * t744 + t668;
t721 = t656 * t661;
t630 = t651 * t708;
t588 = t654 * t689 - t630;
t499 = -t554 * qJD(5) + t659 * t588 + t657 * t690;
t584 = qJD(5) + t593;
t516 = t584 * t554;
t454 = -t516 + t499;
t680 = t657 * t588 - t659 * t690;
t451 = (qJD(5) - t584) * t556 + t680;
t551 = t554 ^ 2;
t552 = t556 ^ 2;
t582 = t584 ^ 2;
t585 = t593 ^ 2;
t586 = t595 ^ 2;
t743 = 2 * qJD(4);
t742 = pkin(2) * t653;
t741 = pkin(2) * t656;
t740 = pkin(3) * t652;
t739 = pkin(3) * t655;
t738 = pkin(4) * t651;
t737 = t656 * g(1);
t650 = qJDD(1) * pkin(1);
t736 = t584 * t657;
t735 = t584 * t659;
t646 = t653 ^ 2;
t733 = t646 * t661;
t648 = t656 ^ 2;
t639 = t648 * t661;
t560 = -t653 * g(1) + t746 * t656;
t674 = -qJ(3) * t653 - t741;
t526 = t674 * t721 + t560;
t627 = t660 * g(2) + t658 * g(3);
t598 = -t661 * qJ(2) + qJDD(2) + t627 - t650;
t665 = qJDD(1) * t674 + t598;
t681 = t652 * t526 - t655 * t665;
t667 = pkin(3) * t708 - qJ(4) * t639 + qJDD(4) + t681;
t672 = -qJ(4) * t655 + t740;
t597 = t672 * t713;
t679 = ((2 * qJD(3)) + t597) * t655;
t441 = t679 * t713 + t667;
t732 = t651 * t441;
t524 = t734 + t690;
t731 = t651 * t524;
t666 = (qJD(1) * t674 + t744) * qJD(1);
t683 = qJDD(3) + t737;
t522 = t683 + (t666 + t668) * t653;
t730 = t652 * t522;
t701 = t652 * t655 * t661;
t615 = t646 * t701;
t599 = -t615 + t708;
t729 = t652 * t599;
t600 = -t615 - t708;
t728 = t652 * t600;
t727 = t653 * t656;
t726 = t654 * t441;
t725 = t654 * t524;
t724 = t655 * t522;
t723 = t655 * t599;
t722 = t655 * t600;
t537 = t593 * pkin(4) - t595 * pkin(6);
t644 = t652 ^ 2;
t633 = t644 * t733;
t692 = qJD(3) * t713;
t477 = t655 * t526 + (t665 - 0.2e1 * t692) * t652;
t444 = -pkin(3) * t639 - qJ(4) * t708 - t597 * t693 + t477;
t673 = -qJ(4) * t652 - t739;
t664 = ((t656 * t673 - pkin(1)) * t661 + (qJ(2) + t672) * qJDD(1) + t666 - t626) * t653 + t683;
t682 = t651 * t444 - t654 * t664;
t370 = -pkin(4) * t690 - pkin(6) * t633 + (t743 + t537) * t595 + t682;
t720 = t657 * t370;
t484 = t504 + t686;
t719 = t657 * t484;
t718 = t658 * t598;
t717 = t659 * t370;
t716 = t659 * t484;
t715 = t660 * t598;
t396 = t654 * t444 - t593 * t743 + t651 * t664;
t371 = -pkin(4) * t633 + pkin(6) * t690 - t593 * t537 + t396;
t401 = t714 * pkin(4) - t588 * pkin(6) + (t679 + (pkin(4) * t595 + pkin(6) * t593) * t652) * t713 + t667;
t330 = t659 * t371 + t657 * t401;
t711 = qJDD(1) * t652;
t707 = t658 * qJDD(1);
t706 = t660 * qJDD(1);
t647 = t655 ^ 2;
t704 = t647 * t733;
t703 = t651 * t504;
t702 = t652 * t734;
t700 = t652 * t721;
t699 = t653 * t721;
t698 = t654 * t504;
t697 = t655 * t734;
t696 = t655 * t721;
t695 = pkin(4) * t654 + pkin(3);
t694 = qJD(1) * t593 * t652;
t691 = t654 * t709;
t688 = t656 * t707;
t687 = t656 * t706;
t670 = t660 * t661 + t707;
t685 = -pkin(5) * t670 + t660 * g(1);
t684 = -t598 + t650;
t329 = t657 * t371 - t659 * t401;
t559 = t746 * t653 + t737;
t503 = t653 * t559 + t656 * t560;
t570 = -t658 * t626 - t660 * t627;
t645 = t653 * t646;
t678 = t645 * t701;
t677 = t652 * t696;
t676 = t593 * t693;
t675 = t595 * t693;
t671 = t646 * t677;
t300 = -t659 * t329 + t657 * t330;
t301 = t657 * t329 + t659 * t330;
t395 = t595 * t743 + t682;
t336 = -t654 * t395 + t651 * t396;
t337 = t651 * t395 + t654 * t396;
t476 = 0.2e1 * t655 * t692 + t681;
t411 = -t655 * t476 + t652 * t477;
t412 = t652 * t476 + t655 * t477;
t502 = t656 * t559 - t653 * t560;
t571 = t660 * t626 - t658 * t627;
t623 = -t658 * t661 + t706;
t609 = (t646 + t648) * t721;
t567 = -t658 * t609 + t687;
t569 = t660 * t609 + t688;
t530 = t675 - t714;
t638 = t648 * qJDD(1);
t637 = t646 * qJDD(1);
t631 = t644 * t710;
t618 = t639 - t733;
t617 = t639 + t733;
t614 = t652 * t689;
t612 = t638 - t637;
t611 = t638 + t637;
t608 = -t639 - t704;
t607 = t639 - t704;
t606 = (t648 * t653 + t645) * t661;
t605 = -t633 - t639;
t604 = t633 - t639;
t603 = pkin(5) * t623 + t658 * g(1);
t602 = t633 - t704;
t601 = t633 + t704;
t592 = (t700 - t709) * t653;
t591 = (t700 + t709) * t653;
t590 = (t696 - t711) * t653;
t589 = (t696 + t711) * t653;
t581 = t623 * t727;
t580 = t670 * t727;
t577 = (-t644 - t647) * t699;
t576 = (qJDD(1) * t647 + t677) * t653;
t575 = t647 * t699 - t614;
t574 = t644 * t699 + t614;
t573 = -t653 * t677 + t631;
t568 = -t660 * t606 - t653 * t707;
t566 = t658 * t606 - t653 * t706;
t565 = -t586 - t633;
t564 = -t586 + t633;
t563 = t585 - t633;
t562 = -t660 * t611 + t658 * t617;
t561 = t658 * t611 + t660 * t617;
t549 = t656 * t576 + t678;
t548 = t656 * t573 - t678;
t547 = -t652 * t607 + t722;
t546 = -t652 * t608 + t723;
t545 = t655 * t605 - t728;
t544 = t655 * t604 + t729;
t543 = -t655 * t607 - t728;
t542 = t655 * t608 + t729;
t541 = t652 * t605 + t722;
t540 = -t652 * t604 + t723;
t538 = -t586 + t585;
t536 = -t655 * t589 - t652 * t592;
t535 = t655 * t590 - t652 * t591;
t534 = -t652 * t589 + t655 * t592;
t533 = -t652 * t590 - t655 * t591;
t532 = t630 + (-t691 - t694) * t653;
t531 = -t630 + (t691 - t694) * t653;
t529 = -t675 - t714;
t527 = -t633 - t585;
t520 = t654 * t588 - t651 * t675;
t519 = -t651 * t588 - t654 * t675;
t518 = t651 * t714 + t654 * t676;
t517 = -t651 * t676 + t654 * t714;
t515 = -t552 + t582;
t514 = t551 - t582;
t513 = t585 + t586;
t512 = (-t593 * t654 + t595 * t651) * t693;
t511 = (t593 * t651 + t595 * t654) * t693;
t510 = t656 * t547 - t653 * t592;
t509 = t656 * t546 + t653 * t591;
t508 = t656 * t545 - t653 * t590;
t507 = t656 * t544 - t653 * t589;
t506 = t653 * t546 - t656 * t591;
t505 = t653 * t545 + t656 * t590;
t500 = -t552 + t551;
t498 = -t556 * qJD(5) - t680;
t497 = t656 * t536 - t653 * t601;
t496 = t656 * t535 - t653 * t602;
t495 = t653 * t536 + t656 * t601;
t494 = t655 * t512 + t631;
t493 = -t652 * t512 + t614;
t492 = -t552 - t582;
t491 = t654 * t563 - t731;
t490 = -t651 * t565 - t725;
t489 = -t651 * t564 + t750;
t488 = -t651 * t563 - t725;
t487 = t654 * t565 - t731;
t486 = -t654 * t564 - t751;
t482 = -t582 - t551;
t481 = -t660 * t503 - t718;
t480 = t658 * t503 - t715;
t479 = -qJ(3) * t542 + t724;
t478 = -qJ(3) * t541 + t730;
t473 = t655 * t520 + t702;
t472 = t655 * t518 - t702;
t471 = -t652 * t520 + t697;
t470 = -t652 * t518 - t697;
t469 = t654 * t530 - t651 * t532;
t468 = t654 * t529 - t651 * t531;
t467 = t651 * t530 + t654 * t532;
t466 = -t651 * t529 - t654 * t531;
t465 = t551 + t552;
t464 = t654 * t527 - t751;
t463 = t651 * t527 + t750;
t462 = -t660 * t509 - t658 * t542;
t461 = -t660 * t508 - t658 * t541;
t460 = t658 * t509 - t660 * t542;
t459 = t658 * t508 - t660 * t541;
t458 = (-t554 * t659 + t556 * t657) * t584;
t457 = (t554 * t657 + t556 * t659) * t584;
t455 = -t516 - t499;
t452 = (-qJD(5) - t584) * t556 - t680;
t450 = -t660 * t497 - t658 * t534;
t449 = t658 * t497 - t660 * t534;
t448 = t659 * t499 - t556 * t736;
t447 = -t657 * t499 - t556 * t735;
t446 = -t657 * t498 + t554 * t735;
t445 = -t659 * t498 - t554 * t736;
t443 = -pkin(2) * t542 + t477;
t442 = -pkin(2) * t541 + t476;
t439 = t655 * t491 + t530 * t652;
t438 = t655 * t490 + t652 * t531;
t437 = t655 * t489 - t652 * t532;
t436 = -t652 * t491 + t530 * t655;
t435 = t652 * t490 - t655 * t531;
t434 = -t652 * t489 - t655 * t532;
t433 = t656 * t494 - t653 * t511;
t432 = t654 * t458 + t651 * t686;
t431 = -t651 * t458 + t654 * t686;
t430 = t655 * t468 - t652 * t538;
t429 = -t652 * t468 - t655 * t538;
t428 = t659 * t514 - t719;
t427 = -t657 * t515 + t747;
t426 = -t657 * t514 - t716;
t425 = -t659 * t515 - t748;
t424 = t655 * t464 - t652 * t529;
t423 = t652 * t464 + t655 * t529;
t422 = t656 * t473 - t653 * t519;
t421 = t656 * t472 - t653 * t517;
t420 = t655 * t469 - t652 * t513;
t419 = t652 * t469 + t655 * t513;
t418 = -t657 * t492 - t716;
t417 = t659 * t492 - t719;
t416 = -pkin(1) * t505 - pkin(2) * t590 - qJ(3) * t545 + t724;
t415 = -pkin(1) * t506 + pkin(2) * t591 - qJ(3) * t546 - t730;
t414 = t659 * t482 - t748;
t413 = t657 * t482 + t747;
t410 = t654 * t448 + t703;
t409 = t654 * t446 - t703;
t408 = -t651 * t448 + t698;
t407 = -t651 * t446 - t698;
t406 = -qJ(4) * t487 + t726;
t405 = t656 * t439 - t653 * t488;
t404 = t656 * t438 + t653 * t487;
t403 = t656 * t437 - t653 * t486;
t402 = t653 * t438 - t656 * t487;
t398 = -qJ(4) * t463 + t732;
t397 = -qJ(3) * t534 - t411;
t394 = t656 * t412 + t653 * t522;
t393 = t653 * t412 - t656 * t522;
t392 = -t451 * t659 - t657 * t455;
t391 = t659 * t452 - t657 * t454;
t390 = -t451 * t657 + t659 * t455;
t389 = -t657 * t452 - t659 * t454;
t388 = t656 * t430 - t653 * t466;
t387 = t655 * t432 - t652 * t457;
t386 = -t652 * t432 - t655 * t457;
t385 = t656 * t424 + t653 * t463;
t384 = t653 * t424 - t656 * t463;
t383 = t656 * t420 + t653 * t467;
t382 = t653 * t420 - t656 * t467;
t381 = t654 * t428 - t651 * t451;
t380 = t654 * t427 - t651 * t455;
t379 = -t651 * t428 - t654 * t451;
t378 = -t651 * t427 - t654 * t455;
t377 = -qJ(2) * t506 - t653 * t443 + t656 * t479;
t376 = -qJ(2) * t505 - t653 * t442 + t656 * t478;
t375 = t654 * t418 + t454 * t651;
t374 = t651 * t418 - t454 * t654;
t373 = t654 * t414 - t651 * t452;
t372 = t651 * t414 + t654 * t452;
t369 = t654 * t391 - t651 * t500;
t368 = -t651 * t391 - t654 * t500;
t367 = -pkin(1) * t495 - pkin(2) * t601 - qJ(3) * t536 - t412;
t365 = -pkin(3) * t487 + t396;
t364 = t655 * t410 - t652 * t447;
t363 = t655 * t409 - t652 * t445;
t362 = -t652 * t410 - t655 * t447;
t361 = -t652 * t409 - t655 * t445;
t360 = -pkin(3) * t463 + t395;
t359 = t654 * t392 - t651 * t465;
t358 = t651 * t392 + t654 * t465;
t357 = -qJ(2) * t495 + t656 * t397 + t534 * t742;
t356 = -t660 * t404 - t658 * t435;
t355 = t658 * t404 - t660 * t435;
t354 = -pkin(2) * t435 + pkin(3) * t531 - qJ(4) * t490 - t732;
t353 = t656 * t387 - t653 * t431;
t352 = -t660 * t385 - t658 * t423;
t351 = t658 * t385 - t660 * t423;
t350 = t655 * t381 - t652 * t426;
t349 = t655 * t380 - t652 * t425;
t348 = -t652 * t381 - t655 * t426;
t347 = -t652 * t380 - t655 * t425;
t346 = -pkin(2) * t423 - pkin(3) * t529 - qJ(4) * t464 + t726;
t345 = -t660 * t383 - t658 * t419;
t344 = t658 * t383 - t660 * t419;
t343 = -t660 * t394 - t658 * t411;
t342 = t658 * t394 - t660 * t411;
t341 = t655 * t375 + t652 * t417;
t340 = t652 * t375 - t655 * t417;
t339 = -pkin(1) * t393 + pkin(2) * t522 - qJ(3) * t412;
t338 = -pkin(6) * t417 + t717;
t335 = t655 * t373 + t652 * t413;
t334 = t652 * t373 - t655 * t413;
t333 = -pkin(6) * t413 + t720;
t332 = t656 * t364 - t653 * t408;
t331 = t656 * t363 - t653 * t407;
t328 = t655 * t369 - t652 * t389;
t327 = -t652 * t369 - t655 * t389;
t326 = -qJ(4) * t467 - t336;
t325 = t655 * t359 + t652 * t390;
t324 = t652 * t359 - t655 * t390;
t323 = t655 * t337 + t652 * t441;
t322 = t652 * t337 - t655 * t441;
t321 = -qJ(2) * t393 + (-qJ(3) * t656 + t742) * t411;
t320 = -qJ(3) * t435 - t652 * t365 + t655 * t406;
t319 = t656 * t350 - t653 * t379;
t318 = t656 * t349 - t653 * t378;
t317 = -pkin(4) * t417 + t330;
t316 = -pkin(4) * t413 + t329;
t315 = -qJ(3) * t423 - t652 * t360 + t655 * t398;
t314 = t656 * t341 + t653 * t374;
t313 = t653 * t341 - t656 * t374;
t312 = t656 * t335 + t653 * t372;
t311 = t653 * t335 - t656 * t372;
t310 = -pkin(2) * t419 - pkin(3) * t513 - qJ(4) * t469 - t337;
t309 = -qJ(3) * t419 + t655 * t326 + t467 * t740;
t308 = t656 * t328 - t653 * t368;
t307 = -pkin(3) * t374 + pkin(4) * t454 - pkin(6) * t418 - t720;
t306 = -pkin(3) * t372 - pkin(4) * t452 - pkin(6) * t414 + t717;
t305 = t656 * t325 + t653 * t358;
t304 = t653 * t325 - t656 * t358;
t303 = -pkin(1) * t402 + pkin(2) * t487 - qJ(3) * t438 - t655 * t365 - t652 * t406;
t302 = -pkin(1) * t384 + pkin(2) * t463 - qJ(3) * t424 - t655 * t360 - t652 * t398;
t299 = t656 * t323 + t653 * t336;
t298 = t653 * t323 - t656 * t336;
t297 = -pkin(1) * t382 - qJ(3) * t420 - t652 * t326 + (pkin(2) + t739) * t467;
t296 = -t660 * t314 - t658 * t340;
t295 = t658 * t314 - t660 * t340;
t294 = -pkin(2) * t322 + pkin(3) * t441 - qJ(4) * t337;
t293 = -t660 * t312 - t658 * t334;
t292 = t658 * t312 - t660 * t334;
t291 = -qJ(2) * t402 + t656 * t320 - t653 * t354;
t290 = -pkin(6) * t390 - t300;
t289 = t654 * t301 + t651 * t370;
t288 = t651 * t301 - t654 * t370;
t287 = -qJ(2) * t384 + t656 * t315 - t653 * t346;
t286 = -qJ(4) * t374 - t651 * t317 + t654 * t338;
t285 = -qJ(4) * t372 - t651 * t316 + t654 * t333;
t284 = -t660 * t305 - t658 * t324;
t283 = t658 * t305 - t660 * t324;
t282 = -qJ(3) * t322 + t336 * t672;
t281 = -t660 * t299 - t658 * t322;
t280 = t658 * t299 - t660 * t322;
t279 = -pkin(3) * t358 - pkin(4) * t465 - pkin(6) * t392 - t301;
t278 = -qJ(2) * t382 + t656 * t309 - t653 * t310;
t277 = -qJ(4) * t358 + t654 * t290 + t390 * t738;
t276 = -pkin(2) * t340 + pkin(3) * t417 - qJ(4) * t375 - t654 * t317 - t651 * t338;
t275 = -pkin(2) * t334 + pkin(3) * t413 - qJ(4) * t373 - t654 * t316 - t651 * t333;
t274 = t655 * t289 + t652 * t300;
t273 = t652 * t289 - t655 * t300;
t272 = -pkin(2) * t324 - qJ(4) * t359 - t651 * t290 + t390 * t695;
t271 = -pkin(3) * t288 + pkin(4) * t370 - pkin(6) * t301;
t270 = -qJ(3) * t340 + t655 * t286 - t652 * t307;
t269 = -qJ(3) * t334 + t655 * t285 - t652 * t306;
t268 = -pkin(1) * t298 - qJ(3) * t323 + (pkin(2) - t673) * t336;
t267 = -qJ(4) * t288 + (-pkin(6) * t654 + t738) * t300;
t266 = -pkin(1) * t313 + pkin(2) * t374 - qJ(3) * t341 - t652 * t286 - t655 * t307;
t265 = t656 * t274 + t653 * t288;
t264 = t653 * t274 - t656 * t288;
t263 = -qJ(2) * t298 + t656 * t282 - t653 * t294;
t262 = -pkin(1) * t311 + pkin(2) * t372 - qJ(3) * t335 - t652 * t285 - t655 * t306;
t261 = -qJ(3) * t324 + t655 * t277 - t652 * t279;
t260 = -pkin(1) * t304 + pkin(2) * t358 - qJ(3) * t325 - t652 * t277 - t655 * t279;
t259 = -qJ(2) * t313 + t656 * t270 - t653 * t276;
t258 = -qJ(2) * t311 + t656 * t269 - t653 * t275;
t257 = -t660 * t265 - t658 * t273;
t256 = t658 * t265 - t660 * t273;
t255 = -pkin(2) * t273 - qJ(4) * t289 + (pkin(6) * t651 + t695) * t300;
t254 = -qJ(2) * t304 + t656 * t261 - t653 * t272;
t253 = -qJ(3) * t273 + t655 * t267 - t652 * t271;
t252 = -pkin(1) * t264 + pkin(2) * t288 - qJ(3) * t274 - t652 * t267 - t655 * t271;
t251 = -qJ(2) * t264 + t656 * t253 - t653 * t255;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t502, 0, 0, 0, 0, 0, 0, t505, t506, t495, t393, 0, 0, 0, 0, 0, 0, t384, t402, t382, t298, 0, 0, 0, 0, 0, 0, t311, t313, t304, t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t623, -t670, 0, t570, 0, 0, 0, 0, 0, 0, t567, t566, t561, t480, 0, 0, 0, 0, 0, 0, t459, t460, t449, t342, 0, 0, 0, 0, 0, 0, t351, t355, t344, t280, 0, 0, 0, 0, 0, 0, t292, t295, t283, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t670, t623, 0, t571, 0, 0, 0, 0, 0, 0, t569, t568, t562, t481, 0, 0, 0, 0, 0, 0, t461, t462, t450, t343, 0, 0, 0, 0, 0, 0, t352, t356, t345, t281, 0, 0, 0, 0, 0, 0, t293, t296, t284, t257; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t627, t626, 0, 0, t637, 0.2e1 * t653 * t708, 0, t638, 0, 0, -qJ(2) * t609 + t656 * t684, qJ(2) * t606 - t653 * t684, pkin(1) * t617 + qJ(2) * t611 + t503, -pkin(1) * t598 + qJ(2) * t503, t653 * t576 - t671, t653 * t535 + t656 * t602, t653 * t547 + t656 * t592, t653 * t573 + t671, t653 * t544 + t656 * t589, t638, -pkin(1) * t541 + qJ(2) * t508 + t656 * t442 + t653 * t478, -pkin(1) * t542 + qJ(2) * t509 + t656 * t443 + t653 * t479, qJ(2) * t497 + t653 * t397 + (-pkin(1) - t741) * t534, qJ(2) * t394 + (-pkin(1) + t674) * t411, t653 * t473 + t656 * t519, t653 * t430 + t656 * t466, t653 * t437 + t656 * t486, t653 * t472 + t656 * t517, t653 * t439 + t656 * t488, t653 * t494 + t656 * t511, -pkin(1) * t423 + qJ(2) * t385 + t653 * t315 + t656 * t346, -pkin(1) * t435 + qJ(2) * t404 + t653 * t320 + t656 * t354, -pkin(1) * t419 + qJ(2) * t383 + t653 * t309 + t656 * t310, -pkin(1) * t322 + qJ(2) * t299 + t653 * t282 + t656 * t294, t653 * t364 + t656 * t408, t653 * t328 + t656 * t368, t653 * t349 + t656 * t378, t653 * t363 + t656 * t407, t653 * t350 + t656 * t379, t653 * t387 + t656 * t431, -pkin(1) * t334 + qJ(2) * t312 + t653 * t269 + t656 * t275, -pkin(1) * t340 + qJ(2) * t314 + t653 * t270 + t656 * t276, -pkin(1) * t324 + qJ(2) * t305 + t653 * t261 + t656 * t272, -pkin(1) * t273 + qJ(2) * t265 + t653 * t253 + t656 * t255; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t670, 0, t623, 0, t685, -t603, -t571, -pkin(5) * t571, t580, t658 * t612 + t660 * t618, t566, -t580, -t567, 0, -pkin(5) * t569 + t660 * t559 + t653 * t718, -pkin(5) * t568 + t660 * t560 + t656 * t718, -pkin(5) * t562 + t658 * t502, -pkin(5) * t481 - (-pkin(1) * t660 - qJ(2) * t658) * t502, t658 * t549 + t660 * t575, t658 * t496 + t660 * t533, t658 * t510 + t660 * t543, t658 * t548 + t660 * t574, t658 * t507 + t660 * t540, t660 * t577 - t653 * t688, -pkin(5) * t461 + t658 * t376 + t660 * t416, -pkin(5) * t462 + t658 * t377 + t660 * t415, -pkin(5) * t450 + t658 * t357 + t660 * t367, -pkin(5) * t343 + t658 * t321 + t660 * t339, t658 * t422 + t660 * t471, t658 * t388 + t660 * t429, t658 * t403 + t660 * t434, t658 * t421 + t660 * t470, t658 * t405 + t660 * t436, t658 * t433 + t660 * t493, -pkin(5) * t352 + t658 * t287 + t660 * t302, -pkin(5) * t356 + t658 * t291 + t660 * t303, -pkin(5) * t345 + t658 * t278 + t660 * t297, -pkin(5) * t281 + t658 * t263 + t660 * t268, t658 * t332 + t660 * t362, t658 * t308 + t660 * t327, t658 * t318 + t660 * t347, t658 * t331 + t660 * t361, t658 * t319 + t660 * t348, t658 * t353 + t660 * t386, -pkin(5) * t293 + t658 * t258 + t660 * t262, -pkin(5) * t296 + t658 * t259 + t660 * t266, -pkin(5) * t284 + t658 * t254 + t660 * t260, -pkin(5) * t257 + t658 * t251 + t660 * t252; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t623, 0, t670, 0, t603, t685, t570, pkin(5) * t570, -t581, -t660 * t612 + t658 * t618, t568, t581, -t569, 0, pkin(5) * t567 + t658 * t559 - t653 * t715, pkin(5) * t566 + t658 * t560 - t656 * t715, pkin(5) * t561 - t660 * t502, pkin(5) * t480 - (-pkin(1) * t658 + qJ(2) * t660) * t502, -t660 * t549 + t658 * t575, -t660 * t496 + t658 * t533, -t660 * t510 + t658 * t543, -t660 * t548 + t658 * t574, -t660 * t507 + t658 * t540, t658 * t577 + t653 * t687, pkin(5) * t459 - t660 * t376 + t658 * t416, pkin(5) * t460 - t660 * t377 + t658 * t415, pkin(5) * t449 - t660 * t357 + t658 * t367, pkin(5) * t342 - t660 * t321 + t658 * t339, -t660 * t422 + t658 * t471, -t660 * t388 + t658 * t429, -t660 * t403 + t658 * t434, -t660 * t421 + t658 * t470, -t660 * t405 + t658 * t436, -t660 * t433 + t658 * t493, pkin(5) * t351 - t660 * t287 + t658 * t302, pkin(5) * t355 - t660 * t291 + t658 * t303, pkin(5) * t344 - t660 * t278 + t658 * t297, pkin(5) * t280 - t660 * t263 + t658 * t268, -t660 * t332 + t658 * t362, -t660 * t308 + t658 * t327, -t660 * t318 + t658 * t347, -t660 * t331 + t658 * t361, -t660 * t319 + t658 * t348, -t660 * t353 + t658 * t386, pkin(5) * t292 - t660 * t258 + t658 * t262, pkin(5) * t295 - t660 * t259 + t658 * t266, pkin(5) * t283 - t660 * t254 + t658 * t260, pkin(5) * t256 - t660 * t251 + t658 * t252;];
tauB_reg = t1;