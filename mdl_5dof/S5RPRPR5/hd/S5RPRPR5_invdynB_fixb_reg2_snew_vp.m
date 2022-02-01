% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:09
% EndTime: 2022-01-23 09:26:24
% DurationCPUTime: 13.89s
% Computational Cost: add. (66054->629), mult. (172332->978), div. (0->0), fcn. (121642->10), ass. (0->450)
t684 = sin(pkin(8));
t691 = cos(qJ(3));
t688 = sin(qJ(3));
t730 = qJD(1) * t684;
t713 = t688 * t730;
t634 = t691 * t684 * qJDD(1) - qJD(3) * t713;
t683 = sin(pkin(9));
t685 = cos(pkin(9));
t725 = qJDD(1) * t688;
t729 = qJD(1) * t691;
t701 = qJD(3) * t729 + t725;
t696 = t701 * t684;
t580 = t685 * t634 - t683 * t696;
t625 = (-t683 * t691 - t685 * t688) * t730;
t686 = cos(pkin(8));
t728 = t686 * qJD(1);
t670 = -qJD(3) + t728;
t761 = t625 * t670;
t544 = -t580 - t761;
t778 = t580 - t761;
t689 = sin(qJ(1));
t692 = cos(qJ(1));
t661 = t692 * g(1) + t689 * g(2);
t693 = qJD(1) ^ 2;
t642 = -t693 * pkin(1) + qJDD(1) * qJ(2) - t661;
t768 = 2 * qJD(2);
t777 = qJD(1) * t768 + t642;
t724 = t686 * qJDD(1);
t669 = -qJDD(3) + t724;
t712 = t684 * t729;
t627 = -t683 * t713 + t685 * t712;
t760 = t627 * t625;
t700 = -t669 + t760;
t776 = t683 * t700;
t775 = t685 * t700;
t687 = sin(qJ(5));
t659 = -qJDD(5) + t669;
t690 = cos(qJ(5));
t569 = -t690 * t625 + t687 * t627;
t571 = t687 * t625 + t690 * t627;
t762 = t571 * t569;
t697 = -t659 - t762;
t774 = t687 * t697;
t773 = t690 * t697;
t707 = t670 * t713;
t599 = t634 + t707;
t588 = t691 * t599;
t664 = -qJD(5) + t670;
t558 = t569 * t664;
t709 = t683 * t634 + t685 * t696;
t698 = t569 * qJD(5) - t690 * t580 + t687 * t709;
t772 = t558 - t698;
t710 = t687 * t580 + t690 * t709;
t458 = (qJD(5) + t664) * t571 + t710;
t566 = t569 ^ 2;
t567 = t571 ^ 2;
t771 = t625 ^ 2;
t622 = t627 ^ 2;
t658 = t664 ^ 2;
t770 = t670 ^ 2;
t769 = t688 ^ 2;
t767 = 2 * qJD(4);
t766 = pkin(2) * t684;
t765 = pkin(2) * t686;
t764 = t686 * g(3);
t763 = qJDD(1) * pkin(1);
t759 = t664 * t687;
t758 = t664 * t690;
t757 = t670 * t627;
t756 = t670 * t683;
t755 = t670 * t685;
t680 = t684 ^ 2;
t754 = t680 * t693;
t604 = -t684 * g(3) + t686 * t777;
t706 = -pkin(6) * t684 - t765;
t647 = t706 * qJD(1);
t573 = t647 * t728 + t604;
t660 = t689 * g(1) - t692 * g(2);
t699 = -t693 * qJ(2) + qJDD(2) - t660;
t703 = -pkin(1) + t706;
t608 = t703 * qJDD(1) + t699;
t595 = t691 * t608;
t714 = t670 * t730;
t732 = t691 * t693;
t496 = -t669 * pkin(3) - t634 * qJ(4) + t595 + (-pkin(3) * t680 * t732 + qJ(4) * t714 - t573) * t688;
t530 = t691 * t573 + t688 * t608;
t631 = -t670 * pkin(3) - qJ(4) * t712;
t668 = t769 * t754;
t499 = -pkin(3) * t668 - qJ(4) * t696 + t670 * t631 + t530;
t429 = -t685 * t496 + t683 * t499 + t627 * t767;
t395 = t700 * pkin(4) + pkin(7) * t544 - t429;
t430 = t683 * t496 + t685 * t499 + t625 * t767;
t593 = -t670 * pkin(4) - t627 * pkin(7);
t407 = -t771 * pkin(4) - pkin(7) * t709 + t670 * t593 + t430;
t355 = -t690 * t395 + t687 * t407;
t356 = t687 * t395 + t690 * t407;
t324 = -t690 * t355 + t687 * t356;
t753 = t683 * t324;
t727 = t768 + t647;
t513 = t764 + qJDD(4) + pkin(3) * t696 - qJ(4) * t668 + (t642 + (t631 * t691 + t727) * qJD(1)) * t684;
t752 = t683 * t513;
t563 = t669 + t760;
t751 = t683 * t563;
t750 = t684 * t669;
t749 = t684 * t686;
t748 = t685 * t324;
t747 = t685 * t513;
t746 = t685 * t563;
t446 = pkin(4) * t709 - t771 * pkin(7) + t627 * t593 + t513;
t745 = t687 * t446;
t509 = t659 - t762;
t744 = t687 * t509;
t377 = -t685 * t429 + t683 * t430;
t743 = t688 * t377;
t572 = t764 + (t727 * qJD(1) + t642) * t684;
t742 = t688 * t572;
t715 = t688 * t732;
t657 = t680 * t715;
t632 = -t657 + t669;
t741 = t688 * t632;
t633 = -t657 - t669;
t740 = t688 * t633;
t635 = -t699 + t763;
t739 = t689 * t635;
t738 = t690 * t446;
t737 = t690 * t509;
t736 = t691 * t377;
t735 = t691 * t572;
t734 = t691 * t632;
t733 = t691 * t633;
t731 = t692 * t635;
t723 = t689 * qJDD(1);
t722 = t692 * qJDD(1);
t720 = t686 * t762;
t719 = t686 * t760;
t682 = t691 ^ 2;
t718 = t682 * t754;
t717 = t684 * t762;
t716 = t684 * t760;
t711 = t635 + t763;
t325 = t687 * t355 + t690 * t356;
t378 = t683 * t429 + t685 * t430;
t529 = t688 * t573 - t595;
t603 = t684 * t777 + t764;
t552 = t684 * t603 + t686 * t604;
t616 = -t689 * t660 - t692 * t661;
t679 = t684 * t680;
t708 = t679 * t715;
t656 = -t689 * t693 + t722;
t705 = -pkin(5) * t656 - t689 * g(3);
t704 = t686 * t657;
t475 = -t691 * t529 + t688 * t530;
t476 = t688 * t529 + t691 * t530;
t551 = t686 * t603 - t684 * t604;
t615 = t692 * t660 - t689 * t661;
t655 = t692 * t693 + t723;
t540 = t709 + t757;
t681 = t686 ^ 2;
t646 = (t680 + t681) * t686 * t693;
t612 = -t689 * t646 + t686 * t722;
t702 = t692 * t646 + t686 * t723;
t676 = t681 * t693;
t675 = t681 * qJDD(1);
t674 = t680 * qJDD(1);
t653 = t676 - t754;
t652 = t676 + t754;
t651 = t686 * t669;
t650 = t675 - t674;
t649 = t675 + t674;
t645 = (t681 * t684 + t679) * t693;
t644 = t670 * t712;
t641 = t668 - t718;
t640 = t668 + t718;
t639 = t770 - t718;
t638 = -t668 - t770;
t637 = t668 - t770;
t636 = -pkin(5) * t655 + t692 * g(3);
t624 = t656 * t749;
t623 = t655 * t749;
t621 = -t770 - t718;
t613 = t692 * t645 + t684 * t723;
t611 = t689 * t645 - t684 * t722;
t610 = t692 * t649 - t689 * t652;
t609 = t689 * t649 + t692 * t652;
t605 = (-t682 - t769) * t714;
t602 = -t622 + t770;
t601 = -t770 + t771;
t598 = t644 - t696;
t597 = t644 + t696;
t596 = t707 - t634;
t591 = -t688 * t634 + t682 * t714;
t590 = (t769 * t670 * qJD(1) + t691 * t701) * t684;
t589 = -t622 - t770;
t587 = (t725 + (qJD(3) - t670) * t729) * t688 * t684;
t586 = t691 * t638 - t740;
t585 = t691 * t637 + t741;
t584 = -t688 * t639 + t733;
t583 = t688 * t638 + t733;
t582 = -t688 * t637 + t734;
t581 = -t691 * t639 - t740;
t577 = -t622 + t771;
t575 = -t688 * t621 + t734;
t574 = t691 * t621 + t741;
t562 = -t770 - t771;
t561 = t686 * t588 + t708;
t560 = t686 * t587 - t708;
t557 = -t567 + t658;
t556 = t566 - t658;
t555 = (-t625 * t685 - t627 * t683) * t670;
t554 = (-t625 * t683 + t627 * t685) * t670;
t553 = -t622 - t771;
t549 = -t567 - t658;
t548 = t691 * t598 - t688 * t599;
t547 = -t688 * t596 - t691 * t597;
t546 = -t688 * t598 - t588;
t545 = t691 * t596 - t688 * t597;
t539 = t709 - t757;
t538 = t685 * t580 + t627 * t756;
t537 = t683 * t580 - t627 * t755;
t536 = t625 * t755 + t683 * t709;
t535 = t625 * t756 - t685 * t709;
t534 = t686 * t586 - t684 * t598;
t533 = t686 * t585 - t684 * t597;
t532 = t686 * t584 - t684 * t596;
t531 = t684 * t586 + t686 * t598;
t528 = t686 * t575 + t599 * t684;
t527 = t684 * t575 - t599 * t686;
t526 = t685 * t601 + t751;
t525 = -t683 * t602 + t775;
t524 = t683 * t601 - t746;
t523 = t685 * t602 + t776;
t522 = t692 * t552 - t739;
t521 = t689 * t552 + t731;
t520 = -pkin(6) * t583 + t742;
t519 = t686 * t548 - t684 * t641;
t518 = t686 * t547 - t684 * t640;
t517 = t684 * t547 + t686 * t640;
t516 = -t683 * t589 + t746;
t515 = t685 * t589 + t751;
t514 = -pkin(6) * t574 + t735;
t512 = -t567 + t566;
t508 = -t658 - t566;
t507 = t685 * t562 - t776;
t506 = t683 * t562 + t775;
t505 = (t569 * t690 - t571 * t687) * t664;
t504 = (t569 * t687 + t571 * t690) * t664;
t503 = -pkin(2) * t583 + t529;
t502 = t692 * t534 + t689 * t583;
t501 = t689 * t534 - t692 * t583;
t500 = -pkin(2) * t574 + t530;
t498 = -t688 * t554 + t691 * t555;
t497 = -t691 * t554 - t688 * t555;
t492 = t692 * t528 + t689 * t574;
t491 = t689 * t528 - t692 * t574;
t489 = -t571 * qJD(5) - t710;
t488 = -t540 * t685 - t683 * t544;
t487 = -t685 * t539 - t683 * t778;
t486 = -t540 * t683 + t685 * t544;
t485 = -t683 * t539 + t685 * t778;
t484 = t686 * t498 - t750;
t483 = -t566 - t567;
t482 = -t688 * t537 + t691 * t538;
t481 = -t688 * t535 + t691 * t536;
t480 = -t691 * t537 - t688 * t538;
t479 = -t691 * t535 - t688 * t536;
t478 = t692 * t518 + t689 * t545;
t477 = t689 * t518 - t692 * t545;
t474 = t690 * t556 + t744;
t473 = -t687 * t557 + t773;
t472 = t687 * t556 - t737;
t471 = t690 * t557 + t774;
t470 = -t688 * t524 + t691 * t526;
t469 = -t688 * t523 + t691 * t525;
t468 = -t687 * t549 + t737;
t467 = -t691 * t524 - t688 * t526;
t466 = -t691 * t523 - t688 * t525;
t465 = t690 * t549 + t744;
t464 = -t688 * t515 + t691 * t516;
t463 = t691 * t515 + t688 * t516;
t462 = t558 + t698;
t457 = (qJD(5) - t664) * t571 + t710;
t456 = -qJ(4) * t515 + t747;
t455 = -pkin(1) * t531 - pkin(2) * t598 - pkin(6) * t586 + t735;
t454 = t571 * t759 - t690 * t698;
t453 = -t571 * t758 - t687 * t698;
t452 = -t687 * t489 - t569 * t758;
t451 = t690 * t489 - t569 * t759;
t450 = t686 * t482 - t716;
t449 = t686 * t481 + t716;
t448 = -qJ(4) * t506 + t752;
t447 = -pkin(1) * t527 + pkin(2) * t599 - pkin(6) * t575 - t742;
t445 = t690 * t508 - t774;
t444 = t687 * t508 + t773;
t443 = -t688 * t506 + t691 * t507;
t442 = t686 * t476 + t684 * t572;
t441 = t691 * t506 + t688 * t507;
t440 = t684 * t476 - t686 * t572;
t439 = -pkin(6) * t545 - t475;
t438 = -t683 * t504 + t685 * t505;
t437 = t685 * t504 + t683 * t505;
t436 = t686 * t470 - t684 * t540;
t435 = t686 * t469 - t684 * t544;
t434 = -pkin(3) * t778 + qJ(4) * t516 + t752;
t433 = t686 * t464 + t684 * t778;
t432 = t684 * t464 - t686 * t778;
t431 = -pkin(3) * t539 + qJ(4) * t507 - t747;
t428 = t686 * t443 + t684 * t539;
t427 = t684 * t443 - t686 * t539;
t426 = -qJ(2) * t531 - t684 * t503 + t686 * t520;
t425 = -t688 * t486 + t691 * t488;
t424 = -t688 * t485 + t691 * t487;
t423 = t691 * t486 + t688 * t488;
t422 = -t691 * t485 - t688 * t487;
t420 = -qJ(2) * t527 - t684 * t500 + t686 * t514;
t419 = -pkin(1) * t517 - pkin(2) * t640 - pkin(6) * t547 - t476;
t418 = t686 * t424 - t684 * t577;
t417 = -t683 * t472 + t685 * t474;
t416 = -t683 * t471 + t685 * t473;
t415 = t685 * t472 + t683 * t474;
t414 = t685 * t471 + t683 * t473;
t413 = -t683 * t465 + t685 * t468;
t412 = t685 * t465 + t683 * t468;
t411 = t686 * t425 + t684 * t553;
t410 = t684 * t425 - t686 * t553;
t409 = t692 * t442 + t689 * t475;
t408 = t689 * t442 - t692 * t475;
t406 = -t458 * t690 - t687 * t462;
t405 = -t690 * t457 - t687 * t772;
t404 = -t458 * t687 + t690 * t462;
t403 = -t687 * t457 + t690 * t772;
t402 = -pkin(7) * t465 + t738;
t400 = -qJ(2) * t517 + t686 * t439 + t545 * t766;
t399 = -t683 * t453 + t685 * t454;
t398 = -t683 * t451 + t685 * t452;
t397 = t685 * t453 + t683 * t454;
t396 = t685 * t451 + t683 * t452;
t392 = -pkin(1) * t440 + pkin(2) * t572 - pkin(6) * t476;
t391 = -pkin(7) * t444 + t745;
t390 = -pkin(2) * t423 - pkin(3) * t486;
t389 = -t683 * t444 + t685 * t445;
t388 = t685 * t444 + t683 * t445;
t387 = t692 * t433 + t689 * t463;
t386 = t689 * t433 - t692 * t463;
t385 = -t688 * t437 + t691 * t438;
t384 = -t691 * t437 - t688 * t438;
t383 = t686 * t385 - t684 * t659;
t382 = t692 * t428 + t689 * t441;
t381 = t689 * t428 - t692 * t441;
t380 = -qJ(2) * t440 + (-pkin(6) * t686 + t766) * t475;
t379 = -pkin(2) * t463 - pkin(3) * t515 + t430;
t376 = -pkin(4) * t772 + pkin(7) * t468 + t745;
t375 = -pkin(2) * t441 - pkin(3) * t506 + t429;
t374 = -pkin(4) * t457 + pkin(7) * t445 - t738;
t373 = -pkin(6) * t463 - t688 * t434 + t691 * t456;
t372 = t692 * t411 + t689 * t423;
t371 = t689 * t411 - t692 * t423;
t370 = -pkin(3) * t513 + qJ(4) * t378;
t369 = -t688 * t415 + t691 * t417;
t368 = -t688 * t414 + t691 * t416;
t367 = -t691 * t415 - t688 * t417;
t366 = -t691 * t414 - t688 * t416;
t365 = -pkin(6) * t441 - t688 * t431 + t691 * t448;
t364 = -t688 * t412 + t691 * t413;
t363 = t691 * t412 + t688 * t413;
t362 = -qJ(4) * t486 - t377;
t361 = -t683 * t404 + t685 * t406;
t360 = -t683 * t403 + t685 * t405;
t359 = t685 * t404 + t683 * t406;
t358 = t685 * t403 + t683 * t405;
t357 = -pkin(3) * t553 + qJ(4) * t488 + t378;
t353 = -t688 * t397 + t691 * t399;
t352 = -t688 * t396 + t691 * t398;
t351 = -t691 * t397 - t688 * t399;
t350 = -t691 * t396 - t688 * t398;
t349 = -t688 * t388 + t691 * t389;
t348 = t691 * t388 + t688 * t389;
t347 = t686 * t369 - t684 * t458;
t346 = t686 * t368 - t684 * t462;
t345 = t686 * t353 + t717;
t344 = t686 * t352 - t717;
t343 = t686 * t364 + t684 * t772;
t342 = t684 * t364 - t686 * t772;
t341 = -pkin(1) * t432 + pkin(2) * t778 - pkin(6) * t464 - t691 * t434 - t688 * t456;
t340 = -pkin(1) * t427 + pkin(2) * t539 - pkin(6) * t443 - t691 * t431 - t688 * t448;
t339 = t686 * t349 + t684 * t457;
t338 = t684 * t349 - t686 * t457;
t337 = t691 * t378 - t743;
t336 = t688 * t378 + t736;
t335 = t686 * t337 + t684 * t513;
t334 = t684 * t337 - t686 * t513;
t333 = -qJ(4) * t412 - t683 * t376 + t685 * t402;
t332 = -qJ(2) * t432 + t686 * t373 - t684 * t379;
t331 = -qJ(4) * t388 - t683 * t374 + t685 * t391;
t330 = -pkin(3) * t772 + qJ(4) * t413 + t685 * t376 + t683 * t402;
t329 = -t688 * t359 + t691 * t361;
t328 = -t688 * t358 + t691 * t360;
t327 = t691 * t359 + t688 * t361;
t326 = -t691 * t358 - t688 * t360;
t323 = -qJ(2) * t427 + t686 * t365 - t684 * t375;
t322 = t692 * t343 + t689 * t363;
t321 = t689 * t343 - t692 * t363;
t320 = t686 * t328 - t684 * t512;
t319 = -pkin(3) * t457 + qJ(4) * t389 + t685 * t374 + t683 * t391;
t318 = t686 * t329 + t684 * t483;
t317 = t684 * t329 - t686 * t483;
t316 = -pkin(2) * t336 - pkin(3) * t377;
t315 = -pkin(6) * t423 - t688 * t357 + t691 * t362;
t314 = -pkin(4) * t446 + pkin(7) * t325;
t313 = t692 * t339 + t689 * t348;
t312 = t689 * t339 - t692 * t348;
t311 = -pkin(7) * t404 - t324;
t310 = -pkin(2) * t363 - pkin(3) * t412 - pkin(4) * t465 + t356;
t309 = -pkin(4) * t483 + pkin(7) * t406 + t325;
t308 = -pkin(1) * t410 + pkin(2) * t553 - pkin(6) * t425 - t691 * t357 - t688 * t362;
t307 = -pkin(2) * t348 - pkin(3) * t388 - pkin(4) * t444 + t355;
t306 = -pkin(6) * t336 - qJ(4) * t736 - t688 * t370;
t305 = t692 * t335 + t689 * t336;
t304 = t689 * t335 - t692 * t336;
t303 = -qJ(2) * t410 + t686 * t315 - t684 * t390;
t302 = -pkin(2) * t327 - pkin(3) * t359 - pkin(4) * t404;
t301 = t685 * t325 - t753;
t300 = t683 * t325 + t748;
t299 = t692 * t318 + t689 * t327;
t298 = t689 * t318 - t692 * t327;
t297 = -pkin(6) * t363 - t688 * t330 + t691 * t333;
t296 = -pkin(1) * t334 + pkin(2) * t513 - pkin(6) * t337 + qJ(4) * t743 - t691 * t370;
t295 = -pkin(6) * t348 - t688 * t319 + t691 * t331;
t294 = -pkin(1) * t342 + pkin(2) * t772 - pkin(6) * t364 - t691 * t330 - t688 * t333;
t293 = -qJ(4) * t359 - t683 * t309 + t685 * t311;
t292 = -pkin(3) * t483 + qJ(4) * t361 + t685 * t309 + t683 * t311;
t291 = -pkin(1) * t338 + pkin(2) * t457 - pkin(6) * t349 - t691 * t319 - t688 * t331;
t290 = -qJ(2) * t334 + t686 * t306 - t684 * t316;
t289 = -t688 * t300 + t691 * t301;
t288 = t691 * t300 + t688 * t301;
t287 = -qJ(2) * t342 + t686 * t297 - t684 * t310;
t286 = -pkin(7) * t748 - qJ(4) * t300 - t683 * t314;
t285 = t686 * t289 + t684 * t446;
t284 = t684 * t289 - t686 * t446;
t283 = -pkin(3) * t446 - pkin(7) * t753 + qJ(4) * t301 + t685 * t314;
t282 = -qJ(2) * t338 + t686 * t295 - t684 * t307;
t281 = -pkin(2) * t288 - pkin(3) * t300 - pkin(4) * t324;
t280 = -pkin(6) * t327 - t688 * t292 + t691 * t293;
t279 = t692 * t285 + t689 * t288;
t278 = t689 * t285 - t692 * t288;
t277 = -pkin(1) * t317 + pkin(2) * t483 - pkin(6) * t329 - t691 * t292 - t688 * t293;
t276 = -qJ(2) * t317 + t686 * t280 - t684 * t302;
t275 = -pkin(6) * t288 - t688 * t283 + t691 * t286;
t274 = -pkin(1) * t284 + pkin(2) * t446 - pkin(6) * t289 - t691 * t283 - t688 * t286;
t273 = -qJ(2) * t284 + t686 * t275 - t684 * t281;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t655, -t656, 0, t616, 0, 0, 0, 0, 0, 0, -t702, t613, t610, t522, 0, 0, 0, 0, 0, 0, t502, t492, t478, t409, 0, 0, 0, 0, 0, 0, t382, t387, t372, t305, 0, 0, 0, 0, 0, 0, t313, t322, t299, t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t656, -t655, 0, t615, 0, 0, 0, 0, 0, 0, t612, t611, t609, t521, 0, 0, 0, 0, 0, 0, t501, t491, t477, t408, 0, 0, 0, 0, 0, 0, t381, t386, t371, t304, 0, 0, 0, 0, 0, 0, t312, t321, t298, t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t551, 0, 0, 0, 0, 0, 0, t531, t527, t517, t440, 0, 0, 0, 0, 0, 0, t427, t432, t410, t334, 0, 0, 0, 0, 0, 0, t338, t342, t317, t284; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t656, 0, -t655, 0, t705, -t636, -t615, -pkin(5) * t615, t624, t692 * t650 - t689 * t653, t613, -t624, t702, 0, -pkin(5) * t612 - t689 * t603 - t684 * t731, -pkin(5) * t611 - t689 * t604 - t686 * t731, -pkin(5) * t609 + t692 * t551, -pkin(5) * t521 - (pkin(1) * t689 - qJ(2) * t692) * t551, t692 * t561 - t689 * t591, t692 * t519 - t689 * t546, t692 * t532 - t689 * t581, t692 * t560 - t689 * t590, t692 * t533 - t689 * t582, -t689 * t605 - t692 * t750, -pkin(5) * t501 + t692 * t426 - t689 * t455, -pkin(5) * t491 + t692 * t420 - t689 * t447, -pkin(5) * t477 + t692 * t400 - t689 * t419, -pkin(5) * t408 + t692 * t380 - t689 * t392, t692 * t450 - t689 * t480, t692 * t418 - t689 * t422, t692 * t435 - t689 * t466, t692 * t449 - t689 * t479, t692 * t436 - t689 * t467, t692 * t484 - t689 * t497, -pkin(5) * t381 + t692 * t323 - t689 * t340, -pkin(5) * t386 + t692 * t332 - t689 * t341, -pkin(5) * t371 + t692 * t303 - t689 * t308, -pkin(5) * t304 + t692 * t290 - t689 * t296, t692 * t345 - t689 * t351, t692 * t320 - t689 * t326, t692 * t346 - t689 * t366, t692 * t344 - t689 * t350, t692 * t347 - t689 * t367, t692 * t383 - t689 * t384, -pkin(5) * t312 + t692 * t282 - t689 * t291, -pkin(5) * t321 + t692 * t287 - t689 * t294, -pkin(5) * t298 + t692 * t276 - t689 * t277, -pkin(5) * t278 + t692 * t273 - t689 * t274; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t655, 0, t656, 0, t636, t705, t616, pkin(5) * t616, t623, t689 * t650 + t692 * t653, t611, -t623, -t612, 0, -pkin(5) * t702 + t692 * t603 - t684 * t739, pkin(5) * t613 + t692 * t604 - t686 * t739, pkin(5) * t610 + t689 * t551, pkin(5) * t522 - (-pkin(1) * t692 - qJ(2) * t689) * t551, t561 * t689 + t591 * t692, t519 * t689 + t546 * t692, t532 * t689 + t581 * t692, t560 * t689 + t590 * t692, t533 * t689 + t582 * t692, t605 * t692 - t689 * t750, pkin(5) * t502 + t426 * t689 + t455 * t692, pkin(5) * t492 + t420 * t689 + t447 * t692, pkin(5) * t478 + t400 * t689 + t419 * t692, pkin(5) * t409 + t380 * t689 + t392 * t692, t450 * t689 + t480 * t692, t418 * t689 + t422 * t692, t435 * t689 + t466 * t692, t449 * t689 + t479 * t692, t436 * t689 + t467 * t692, t484 * t689 + t497 * t692, pkin(5) * t382 + t323 * t689 + t340 * t692, pkin(5) * t387 + t332 * t689 + t341 * t692, pkin(5) * t372 + t303 * t689 + t308 * t692, pkin(5) * t305 + t290 * t689 + t296 * t692, t345 * t689 + t351 * t692, t320 * t689 + t326 * t692, t346 * t689 + t366 * t692, t344 * t689 + t350 * t692, t347 * t689 + t367 * t692, t383 * t689 + t384 * t692, pkin(5) * t313 + t282 * t689 + t291 * t692, pkin(5) * t322 + t287 * t689 + t294 * t692, pkin(5) * t299 + t276 * t689 + t277 * t692, pkin(5) * t279 + t273 * t689 + t274 * t692; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t660, t661, 0, 0, t674, 0.2e1 * t684 * t724, 0, t675, 0, 0, -qJ(2) * t646 + t711 * t686, qJ(2) * t645 - t711 * t684, pkin(1) * t652 + qJ(2) * t649 + t552, pkin(1) * t635 + qJ(2) * t552, t588 * t684 - t704, t548 * t684 + t641 * t686, t584 * t684 + t596 * t686, t587 * t684 + t704, t585 * t684 + t597 * t686, t651, -pkin(1) * t583 + qJ(2) * t534 + t503 * t686 + t520 * t684, -pkin(1) * t574 + qJ(2) * t528 + t500 * t686 + t514 * t684, qJ(2) * t518 + t684 * t439 + (-pkin(1) - t765) * t545, qJ(2) * t442 + t703 * t475, t482 * t684 + t719, t424 * t684 + t577 * t686, t469 * t684 + t544 * t686, t481 * t684 - t719, t470 * t684 + t540 * t686, t498 * t684 + t651, -pkin(1) * t441 + qJ(2) * t428 + t365 * t684 + t375 * t686, -pkin(1) * t463 + qJ(2) * t433 + t373 * t684 + t379 * t686, -pkin(1) * t423 + qJ(2) * t411 + t315 * t684 + t390 * t686, -pkin(1) * t336 + qJ(2) * t335 + t306 * t684 + t316 * t686, t353 * t684 - t720, t328 * t684 + t512 * t686, t368 * t684 + t462 * t686, t352 * t684 + t720, t369 * t684 + t458 * t686, t385 * t684 + t659 * t686, -pkin(1) * t348 + qJ(2) * t339 + t295 * t684 + t307 * t686, -pkin(1) * t363 + qJ(2) * t343 + t297 * t684 + t310 * t686, -pkin(1) * t327 + qJ(2) * t318 + t280 * t684 + t302 * t686, -pkin(1) * t288 + qJ(2) * t285 + t275 * t684 + t281 * t686;];
tauB_reg = t1;
