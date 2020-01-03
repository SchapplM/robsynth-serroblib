% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRRP8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:20
% EndTime: 2019-12-31 22:02:36
% DurationCPUTime: 10.83s
% Computational Cost: add. (47379->548), mult. (95967->782), div. (0->0), fcn. (67189->8), ass. (0->405)
t676 = sin(qJ(3));
t677 = sin(qJ(2));
t720 = qJD(1) * qJD(2);
t665 = t677 * t720;
t681 = cos(qJ(2));
t717 = t681 * qJDD(1);
t641 = -t665 + t717;
t633 = -qJDD(3) + t641;
t680 = cos(qJ(3));
t725 = qJD(1) * t677;
t636 = t676 * qJD(2) + t680 * t725;
t704 = t680 * qJD(2) - t676 * t725;
t694 = t704 * t636;
t781 = -t633 + t694;
t783 = t676 * t781;
t782 = t680 * t781;
t678 = sin(qJ(1));
t682 = cos(qJ(1));
t651 = t682 * g(1) + t678 * g(2);
t683 = qJD(1) ^ 2;
t627 = -t683 * pkin(1) + qJDD(1) * pkin(6) - t651;
t771 = pkin(2) * t681;
t699 = -pkin(7) * t677 - t771;
t638 = t699 * qJD(1);
t756 = t681 * g(3);
t775 = qJD(2) ^ 2;
t572 = (qJD(1) * t638 + t627) * t677 - qJDD(2) * pkin(2) - t775 * pkin(7) + t756;
t723 = t681 * qJD(1);
t662 = -qJD(3) + t723;
t617 = -t662 * pkin(3) - t636 * pkin(8);
t712 = t681 * t720;
t719 = t677 * qJDD(1);
t640 = t712 + t719;
t706 = -t680 * qJDD(2) + t676 * t640;
t691 = t636 * qJD(3) + t706;
t701 = t704 ^ 2;
t509 = t691 * pkin(3) - t701 * pkin(8) + t636 * t617 + t572;
t629 = -qJDD(4) + t633;
t675 = sin(qJ(4));
t679 = cos(qJ(4));
t598 = t675 * t636 - t679 * t704;
t600 = t679 * t636 + t675 * t704;
t755 = t600 * t598;
t707 = t629 + t755;
t780 = t707 * pkin(4);
t593 = t704 * qJD(3) + t676 * qJDD(2) + t680 * t640;
t708 = t675 * t593 + t679 * t691;
t525 = -t600 * qJD(4) - t708;
t653 = -qJD(4) + t662;
t576 = -t653 * pkin(4) - t600 * qJ(5);
t596 = t598 ^ 2;
t438 = -t525 * pkin(4) - t596 * qJ(5) + t600 * t576 + qJDD(5) + t509;
t747 = t675 * t707;
t736 = t679 * t707;
t526 = -t598 * qJD(4) + t679 * t593 - t675 * t691;
t583 = t598 * t653;
t779 = t583 + t526;
t622 = t704 * t662;
t562 = -t593 - t622;
t597 = t600 ^ 2;
t652 = t653 ^ 2;
t566 = -t597 - t652;
t534 = t629 - t755;
t748 = t675 * t534;
t505 = t679 * t566 + t748;
t737 = t679 * t534;
t506 = -t675 * t566 + t737;
t458 = t680 * t505 + t676 * t506;
t778 = -pkin(2) * t458 - pkin(3) * t505;
t541 = -t652 - t596;
t489 = t675 * t541 - t736;
t490 = t679 * t541 + t747;
t441 = t680 * t489 + t676 * t490;
t777 = -pkin(2) * t441 - pkin(3) * t489;
t496 = (qJD(4) + t653) * t600 + t708;
t557 = (qJD(3) + t662) * t636 + t706;
t632 = t636 ^ 2;
t659 = t662 ^ 2;
t772 = pkin(2) * t677;
t500 = t583 - t526;
t451 = -t496 * t675 + t679 * t500;
t453 = -t496 * t679 - t675 * t500;
t407 = -t676 * t451 + t680 * t453;
t529 = -t596 - t597;
t386 = t681 * t407 + t677 * t529;
t405 = t680 * t451 + t676 * t453;
t355 = t678 * t386 - t682 * t405;
t768 = pkin(5) * t355;
t442 = -t676 * t489 + t680 * t490;
t495 = (qJD(4) - t653) * t600 + t708;
t416 = t681 * t442 + t677 * t495;
t377 = t678 * t416 - t682 * t441;
t767 = pkin(5) * t377;
t459 = -t676 * t505 + t680 * t506;
t422 = t681 * t459 + t677 * t779;
t388 = t678 * t422 - t682 * t458;
t766 = pkin(5) * t388;
t385 = t677 * t407 - t681 * t529;
t765 = pkin(6) * t385;
t415 = t677 * t442 - t681 * t495;
t764 = pkin(6) * t415;
t421 = t677 * t459 - t681 * t779;
t763 = pkin(6) * t421;
t762 = pkin(7) * t405;
t761 = pkin(7) * t441;
t760 = pkin(7) * t458;
t759 = pkin(8) * t451;
t758 = pkin(8) * t489;
t757 = pkin(8) * t505;
t753 = t653 * t675;
t752 = t653 * t679;
t671 = t677 ^ 2;
t751 = t671 * t683;
t650 = t678 * g(1) - t682 * g(2);
t626 = qJDD(1) * pkin(1) + t683 * pkin(6) + t650;
t696 = -t641 + t665;
t697 = t640 + t712;
t556 = pkin(2) * t696 - pkin(7) * t697 - t626;
t616 = -t677 * g(3) + t681 * t627;
t573 = -t775 * pkin(2) + qJDD(2) * pkin(7) + t638 * t723 + t616;
t527 = -t680 * t556 + t676 * t573;
t476 = pkin(3) * t781 + t562 * pkin(8) - t527;
t528 = t676 * t556 + t680 * t573;
t482 = -pkin(3) * t701 - pkin(8) * t691 + t662 * t617 + t528;
t430 = -t679 * t476 + t675 * t482;
t713 = t526 * qJ(5) + t430;
t690 = qJ(5) * t583 - t713;
t724 = qJD(5) * t600;
t401 = t690 - 0.2e1 * t724 - t780;
t750 = t675 * t401;
t749 = t675 * t509;
t431 = t675 * t476 + t679 * t482;
t381 = -t679 * t430 + t675 * t431;
t746 = t676 * t381;
t745 = t676 * t572;
t586 = t633 + t694;
t744 = t676 * t586;
t743 = t676 * t636;
t742 = t677 * t626;
t661 = t681 * t683 * t677;
t648 = -t661 + qJDD(2);
t741 = t677 * t648;
t649 = qJDD(2) + t661;
t740 = t677 * t649;
t739 = t679 * t401;
t738 = t679 * t509;
t735 = t680 * t381;
t734 = t680 * t572;
t733 = t680 * t586;
t732 = t680 * t636;
t731 = t681 * t626;
t730 = t681 * t648;
t729 = -pkin(1) * t405 + pkin(6) * t386;
t728 = -pkin(1) * t441 + pkin(6) * t416;
t727 = -pkin(1) * t458 + pkin(6) * t422;
t672 = t681 ^ 2;
t726 = t671 + t672;
t718 = t678 * qJDD(1);
t716 = t682 * qJDD(1);
t715 = t677 * t755;
t714 = t681 * t755;
t711 = -pkin(3) * t529 + pkin(8) * t453;
t710 = -pkin(3) * t495 + pkin(8) * t490;
t709 = -pkin(3) * t779 + pkin(8) * t506;
t382 = t675 * t430 + t679 * t431;
t615 = t677 * t627 + t756;
t565 = t677 * t615 + t681 * t616;
t608 = -t678 * t650 - t682 * t651;
t703 = t678 * t661;
t702 = t682 * t661;
t560 = -t622 + t593;
t371 = -pkin(2) * t405 - pkin(3) * t451;
t645 = -t678 * t683 + t716;
t698 = -pkin(5) * t645 - t678 * g(3);
t470 = -t680 * t527 + t676 * t528;
t471 = t676 * t527 + t680 * t528;
t564 = t681 * t615 - t677 * t616;
t607 = t682 * t650 - t678 * t651;
t693 = t677 * t694;
t692 = t681 * t694;
t689 = t525 * qJ(5) - 0.2e1 * qJD(5) * t598 + t653 * t576 + t431;
t688 = -pkin(1) * t385 + pkin(2) * t529 - pkin(7) * t407;
t687 = -pkin(1) * t415 + pkin(2) * t495 - pkin(7) * t442;
t686 = -pkin(1) * t421 + pkin(2) * t779 - pkin(7) * t459;
t669 = t672 * t683;
t658 = -t669 - t775;
t657 = t669 - t775;
t656 = -t751 - t775;
t655 = -t751 + t775;
t647 = t669 - t751;
t646 = t669 + t751;
t644 = t682 * t683 + t718;
t643 = t726 * qJDD(1);
t642 = -0.2e1 * t665 + t717;
t639 = 0.2e1 * t712 + t719;
t635 = t681 * t649;
t634 = t726 * t720;
t624 = -pkin(5) * t644 + t682 * g(3);
t621 = -t632 + t659;
t620 = t701 - t659;
t619 = t681 * t640 - t671 * t720;
t618 = -t677 * t641 - t672 * t720;
t614 = -t677 * t656 - t730;
t613 = -t677 * t655 + t635;
t612 = t681 * t658 - t740;
t611 = t681 * t657 - t741;
t610 = t681 * t656 - t741;
t609 = t677 * t658 + t635;
t605 = -t632 + t701;
t604 = t682 * t643 - t678 * t646;
t603 = t678 * t643 + t682 * t646;
t602 = -t632 - t659;
t601 = -t677 * t639 + t681 * t642;
t595 = -t659 - t701;
t592 = 0.2e1 * t724;
t585 = t701 + t632;
t582 = -t597 + t652;
t581 = t596 - t652;
t580 = t682 * t614 + t678 * t639;
t579 = t682 * t612 - t678 * t642;
t578 = t678 * t614 - t682 * t639;
t577 = t678 * t612 + t682 * t642;
t575 = -pkin(6) * t610 - t731;
t574 = -pkin(6) * t609 - t742;
t571 = (-t680 * t704 - t743) * t662;
t570 = (t676 * t704 - t732) * t662;
t568 = -pkin(1) * t610 + t616;
t567 = -pkin(1) * t609 + t615;
t558 = (-qJD(3) + t662) * t636 - t706;
t553 = t680 * t593 + t662 * t743;
t552 = -t676 * t593 + t662 * t732;
t551 = t622 * t680 + t676 * t691;
t550 = -t622 * t676 + t680 * t691;
t549 = t681 * t571 - t677 * t633;
t548 = -t597 + t596;
t547 = t680 * t620 + t744;
t546 = -t676 * t621 + t782;
t545 = -t676 * t620 + t733;
t544 = -t680 * t621 - t783;
t543 = -t676 * t602 + t733;
t542 = t680 * t602 + t744;
t540 = t682 * t565 - t678 * t626;
t539 = t678 * t565 + t682 * t626;
t538 = t680 * t595 - t783;
t537 = t676 * t595 + t782;
t533 = (t598 * t679 - t600 * t675) * t653;
t532 = (t598 * t675 + t600 * t679) * t653;
t531 = t681 * t553 - t693;
t530 = t681 * t551 + t693;
t521 = -t557 * t680 - t676 * t562;
t520 = t680 * t558 - t676 * t560;
t519 = -t557 * t676 + t680 * t562;
t518 = -t676 * t558 - t680 * t560;
t517 = -pkin(7) * t542 + t734;
t516 = t681 * t547 - t677 * t557;
t515 = t681 * t546 - t677 * t562;
t514 = t679 * t581 + t748;
t513 = -t675 * t582 - t736;
t512 = t675 * t581 - t737;
t511 = t679 * t582 - t747;
t510 = -pkin(7) * t537 + t745;
t508 = t681 * t543 + t560 * t677;
t507 = t677 * t543 - t560 * t681;
t503 = t681 * t538 - t677 * t558;
t502 = t677 * t538 + t681 * t558;
t501 = t681 * t520 - t677 * t605;
t494 = t679 * t526 + t600 * t753;
t493 = t675 * t526 - t600 * t752;
t492 = -t675 * t525 - t598 * t752;
t491 = t679 * t525 - t598 * t753;
t487 = t681 * t521 - t677 * t585;
t486 = t677 * t521 + t681 * t585;
t485 = -pkin(2) * t542 + t528;
t484 = -t676 * t532 + t680 * t533;
t483 = -t680 * t532 - t676 * t533;
t481 = -pkin(2) * t537 + t527;
t478 = t681 * t484 - t677 * t629;
t477 = t677 * t484 + t681 * t629;
t473 = t682 * t508 + t678 * t542;
t472 = t678 * t508 - t682 * t542;
t469 = t682 * t503 + t678 * t537;
t468 = t678 * t503 - t682 * t537;
t467 = -pkin(4) * t779 + qJ(5) * t534;
t466 = -t676 * t512 + t680 * t514;
t465 = -t676 * t511 + t680 * t513;
t464 = -t680 * t512 - t676 * t514;
t463 = -t680 * t511 - t676 * t513;
t462 = t681 * t471 + t677 * t572;
t461 = t677 * t471 - t681 * t572;
t460 = t738 - t757;
t456 = t682 * t487 + t678 * t519;
t455 = t678 * t487 - t682 * t519;
t454 = t749 - t758;
t452 = -t679 * t495 - t675 * t779;
t450 = -t675 * t495 + t679 * t779;
t448 = -pkin(1) * t507 + pkin(2) * t560 - pkin(7) * t543 - t745;
t447 = -t676 * t493 + t680 * t494;
t446 = -t676 * t491 + t680 * t492;
t445 = -t680 * t493 - t676 * t494;
t444 = -t680 * t491 - t676 * t492;
t443 = -pkin(1) * t502 - pkin(2) * t558 - pkin(7) * t538 + t734;
t439 = -pkin(7) * t519 - t470;
t437 = t682 * t478 - t678 * t483;
t436 = t678 * t478 + t682 * t483;
t435 = t681 * t447 + t715;
t434 = t681 * t446 - t715;
t433 = t677 * t447 - t714;
t432 = t677 * t446 + t714;
t428 = -qJ(5) * t566 + t438;
t427 = t681 * t466 - t677 * t496;
t426 = t681 * t465 - t677 * t500;
t425 = t677 * t466 + t681 * t496;
t424 = t677 * t465 + t681 * t500;
t423 = t709 + t749;
t419 = -pkin(6) * t507 - t677 * t485 + t681 * t517;
t418 = t710 - t738;
t417 = -pkin(6) * t502 - t677 * t481 + t681 * t510;
t413 = t682 * t462 + t678 * t470;
t412 = t678 * t462 - t682 * t470;
t411 = -pkin(1) * t486 - pkin(2) * t585 - pkin(7) * t521 - t471;
t410 = -pkin(1) * t461 + pkin(2) * t572 - pkin(7) * t471;
t409 = -pkin(4) * t495 + qJ(5) * t541 - t438;
t408 = -pkin(6) * t486 + t681 * t439 + t519 * t772;
t406 = -t676 * t450 + t680 * t452;
t404 = -t680 * t450 - t676 * t452;
t402 = -t596 * pkin(4) + t689;
t400 = t682 * t427 - t678 * t464;
t399 = t682 * t426 - t678 * t463;
t398 = t678 * t427 + t682 * t464;
t397 = t678 * t426 + t682 * t463;
t396 = t681 * t406 - t677 * t548;
t395 = t677 * t406 + t681 * t548;
t394 = -pkin(6) * t461 + (-pkin(7) * t681 + t772) * t470;
t393 = t682 * t435 - t678 * t445;
t392 = t682 * t434 - t678 * t444;
t391 = t678 * t435 + t682 * t445;
t390 = t678 * t434 + t682 * t444;
t389 = t682 * t422 + t678 * t458;
t387 = pkin(5) * t389;
t383 = t679 * t428 - t675 * t467 - t757;
t380 = t431 + t778;
t379 = t592 + (-t500 - t583) * qJ(5) + t780 + t713;
t378 = t682 * t416 + t678 * t441;
t376 = pkin(5) * t378;
t375 = qJ(5) * t736 - t675 * t409 - t758;
t374 = -qJ(5) * t496 + (-t529 - t596) * pkin(4) + t689;
t373 = t430 + t777;
t372 = -pkin(3) * t509 + pkin(8) * t382;
t370 = t675 * t428 + t679 * t467 + t709;
t369 = qJ(5) * t747 + t679 * t409 + t710;
t368 = -t676 * t423 + t680 * t460 - t760;
t367 = -pkin(4) * t438 + qJ(5) * t402;
t366 = -t676 * t418 + t680 * t454 - t761;
t365 = -pkin(4) * t500 + t371;
t364 = (-t566 - t596) * pkin(4) + t689 + t778;
t363 = -t381 - t759;
t362 = t382 + t711;
t361 = t592 - t690 + t777 + 0.2e1 * t780;
t360 = t679 * t402 - t750;
t359 = t675 * t402 + t739;
t358 = t682 * t396 - t678 * t404;
t357 = t678 * t396 + t682 * t404;
t356 = t682 * t386 + t678 * t405;
t354 = pkin(5) * t356;
t353 = -t680 * t423 - t676 * t460 + t686;
t352 = t680 * t382 - t746;
t351 = t676 * t382 + t735;
t350 = -t680 * t418 - t676 * t454 + t687;
t349 = t681 * t352 + t677 * t509;
t348 = t677 * t352 - t681 * t509;
t347 = -t676 * t370 + t680 * t383 - t760;
t346 = -t675 * t374 + t679 * t379 - t759;
t345 = t679 * t374 + t675 * t379 + t711;
t344 = -t676 * t369 + t680 * t375 - t761;
t343 = t681 * t368 - t677 * t380 - t763;
t342 = -pkin(2) * t351 - pkin(3) * t381;
t341 = t681 * t366 - t677 * t373 - t764;
t340 = -t680 * t370 - t676 * t383 + t686;
t339 = -t676 * t359 + t680 * t360;
t338 = t680 * t359 + t676 * t360;
t337 = -t680 * t369 - t676 * t375 + t687;
t336 = -t676 * t362 + t680 * t363 - t762;
t335 = -pkin(8) * t359 - qJ(5) * t739 - t675 * t367;
t334 = t681 * t339 + t677 * t438;
t333 = t677 * t339 - t681 * t438;
t332 = -pkin(3) * t438 + pkin(8) * t360 - qJ(5) * t750 + t679 * t367;
t331 = -pkin(7) * t351 - pkin(8) * t735 - t676 * t372;
t330 = t682 * t349 + t678 * t351;
t329 = t678 * t349 - t682 * t351;
t328 = t681 * t347 - t677 * t364 - t763;
t327 = -t680 * t362 - t676 * t363 + t688;
t326 = t681 * t344 - t677 * t361 - t764;
t325 = -pkin(2) * t338 - pkin(3) * t359 - pkin(4) * t401;
t324 = t681 * t336 - t677 * t371 - t765;
t323 = -t676 * t345 + t680 * t346 - t762;
t322 = -pkin(1) * t348 + pkin(2) * t509 - pkin(7) * t352 + pkin(8) * t746 - t680 * t372;
t321 = t682 * t334 + t678 * t338;
t320 = t678 * t334 - t682 * t338;
t319 = -t680 * t345 - t676 * t346 + t688;
t318 = t681 * t323 - t677 * t365 - t765;
t317 = -pkin(6) * t348 + t681 * t331 - t677 * t342;
t316 = -pkin(7) * t338 - t676 * t332 + t680 * t335;
t315 = -pkin(1) * t333 + pkin(2) * t438 - pkin(7) * t339 - t680 * t332 - t676 * t335;
t314 = -pkin(6) * t333 + t681 * t316 - t677 * t325;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t644, -t645, 0, t608, 0, 0, 0, 0, 0, 0, t579, t580, t604, t540, 0, 0, 0, 0, 0, 0, t469, t473, t456, t413, 0, 0, 0, 0, 0, 0, t378, t389, t356, t330, 0, 0, 0, 0, 0, 0, t378, t389, t356, t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t645, -t644, 0, t607, 0, 0, 0, 0, 0, 0, t577, t578, t603, t539, 0, 0, 0, 0, 0, 0, t468, t472, t455, t412, 0, 0, 0, 0, 0, 0, t377, t388, t355, t329, 0, 0, 0, 0, 0, 0, t377, t388, t355, t320; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t609, t610, 0, -t564, 0, 0, 0, 0, 0, 0, t502, t507, t486, t461, 0, 0, 0, 0, 0, 0, t415, t421, t385, t348, 0, 0, 0, 0, 0, 0, t415, t421, t385, t333; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t645, 0, -t644, 0, t698, -t624, -t607, -pkin(5) * t607, t682 * t619 - t703, t682 * t601 - t678 * t647, t682 * t613 + t677 * t718, t682 * t618 + t703, t682 * t611 + t678 * t717, t678 * qJDD(2) + t682 * t634, -pkin(5) * t577 - t678 * t567 + t682 * t574, -pkin(5) * t578 - t678 * t568 + t682 * t575, -pkin(5) * t603 + t682 * t564, -pkin(5) * t539 - (pkin(1) * t678 - pkin(6) * t682) * t564, t682 * t531 - t678 * t552, t682 * t501 - t678 * t518, t682 * t515 - t678 * t544, t682 * t530 - t678 * t550, t682 * t516 - t678 * t545, t682 * t549 - t678 * t570, -pkin(5) * t468 + t682 * t417 - t678 * t443, -pkin(5) * t472 + t682 * t419 - t678 * t448, -pkin(5) * t455 + t682 * t408 - t678 * t411, -pkin(5) * t412 + t682 * t394 - t678 * t410, t393, t358, t399, t392, t400, t437, t682 * t341 - t678 * t350 - t767, t682 * t343 - t678 * t353 - t766, t682 * t324 - t678 * t327 - t768, -pkin(5) * t329 + t682 * t317 - t678 * t322, t393, t358, t399, t392, t400, t437, t682 * t326 - t678 * t337 - t767, t682 * t328 - t678 * t340 - t766, t682 * t318 - t678 * t319 - t768, -pkin(5) * t320 + t682 * t314 - t678 * t315; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t644, 0, t645, 0, t624, t698, t608, pkin(5) * t608, t678 * t619 + t702, t678 * t601 + t682 * t647, t678 * t613 - t677 * t716, t678 * t618 - t702, t678 * t611 - t681 * t716, -t682 * qJDD(2) + t678 * t634, pkin(5) * t579 + t682 * t567 + t678 * t574, pkin(5) * t580 + t682 * t568 + t678 * t575, pkin(5) * t604 + t678 * t564, pkin(5) * t540 - (-pkin(1) * t682 - pkin(6) * t678) * t564, t678 * t531 + t682 * t552, t678 * t501 + t682 * t518, t678 * t515 + t682 * t544, t678 * t530 + t682 * t550, t678 * t516 + t682 * t545, t678 * t549 + t682 * t570, pkin(5) * t469 + t678 * t417 + t682 * t443, pkin(5) * t473 + t678 * t419 + t682 * t448, pkin(5) * t456 + t678 * t408 + t682 * t411, pkin(5) * t413 + t678 * t394 + t682 * t410, t391, t357, t397, t390, t398, t436, t678 * t341 + t682 * t350 + t376, t678 * t343 + t682 * t353 + t387, t678 * t324 + t682 * t327 + t354, pkin(5) * t330 + t678 * t317 + t682 * t322, t391, t357, t397, t390, t398, t436, t678 * t326 + t682 * t337 + t376, t678 * t328 + t682 * t340 + t387, t678 * t318 + t682 * t319 + t354, pkin(5) * t321 + t678 * t314 + t682 * t315; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t650, t651, 0, 0, t697 * t677, t681 * t639 + t677 * t642, t681 * t655 + t740, -t696 * t681, t677 * t657 + t730, 0, pkin(1) * t642 + pkin(6) * t612 + t731, -pkin(1) * t639 + pkin(6) * t614 - t742, pkin(1) * t646 + pkin(6) * t643 + t565, pkin(1) * t626 + pkin(6) * t565, t677 * t553 + t692, t677 * t520 + t681 * t605, t677 * t546 + t681 * t562, t677 * t551 - t692, t677 * t547 + t681 * t557, t677 * t571 + t681 * t633, -pkin(1) * t537 + pkin(6) * t503 + t681 * t481 + t677 * t510, -pkin(1) * t542 + pkin(6) * t508 + t681 * t485 + t677 * t517, pkin(6) * t487 + t677 * t439 + (-pkin(1) - t771) * t519, pkin(6) * t462 + (-pkin(1) + t699) * t470, t433, t395, t424, t432, t425, t477, t677 * t366 + t681 * t373 + t728, t677 * t368 + t681 * t380 + t727, t677 * t336 + t681 * t371 + t729, -pkin(1) * t351 + pkin(6) * t349 + t677 * t331 + t681 * t342, t433, t395, t424, t432, t425, t477, t677 * t344 + t681 * t361 + t728, t677 * t347 + t681 * t364 + t727, t677 * t323 + t681 * t365 + t729, -pkin(1) * t338 + pkin(6) * t334 + t677 * t316 + t681 * t325;];
tauB_reg = t1;
