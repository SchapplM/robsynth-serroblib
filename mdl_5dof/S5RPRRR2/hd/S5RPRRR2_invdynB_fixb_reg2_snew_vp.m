% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:35
% EndTime: 2019-12-05 18:12:52
% DurationCPUTime: 16.76s
% Computational Cost: add. (113756->633), mult. (284120->979), div. (0->0), fcn. (220478->10), ass. (0->439)
t726 = sin(qJ(1));
t730 = cos(qJ(1));
t703 = t730 * g(1) + t726 * g(2);
t732 = qJD(1) ^ 2;
t801 = -t732 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) - t703;
t723 = sin(qJ(5));
t722 = cos(pkin(9));
t729 = cos(qJ(3));
t721 = sin(pkin(9));
t725 = sin(qJ(3));
t782 = t721 * t725;
t689 = (-t722 * t729 + t782) * qJD(1);
t739 = t721 * t729 + t722 * t725;
t691 = t739 * qJD(1);
t724 = sin(qJ(4));
t728 = cos(qJ(4));
t647 = t728 * t689 + t724 * t691;
t649 = -t724 * t689 + t728 * t691;
t727 = cos(qJ(5));
t602 = t727 * t647 + t723 * t649;
t604 = -t723 * t647 + t727 * t649;
t549 = t604 * t602;
t717 = qJDD(3) + qJDD(4);
t744 = qJDD(5) + t717;
t800 = -t549 + t744;
t807 = t723 * t800;
t610 = t649 * t647;
t798 = -t610 + t717;
t806 = t724 * t798;
t663 = t691 * t689;
t794 = qJDD(3) - t663;
t805 = t725 * t794;
t804 = t727 * t800;
t803 = t728 * t798;
t802 = t729 * t794;
t720 = qJD(3) + qJD(4);
t640 = t720 * t647;
t757 = qJDD(1) * t722;
t687 = qJDD(1) * t782 - t729 * t757;
t761 = t691 * qJD(3);
t658 = -t687 - t761;
t762 = t689 * qJD(3);
t795 = t739 * qJDD(1);
t660 = t795 - t762;
t736 = t647 * qJD(4) - t724 * t658 - t728 * t660;
t561 = -t640 + t736;
t713 = qJD(5) + t720;
t597 = t713 * t602;
t745 = -t728 * t658 + t724 * t660;
t583 = -t649 * qJD(4) - t745;
t735 = t602 * qJD(5) - t723 * t583 + t727 * t736;
t799 = -t597 - t735;
t797 = -t640 - t736;
t718 = t721 ^ 2;
t719 = t722 ^ 2;
t796 = t718 + t719;
t793 = t732 * t796;
t746 = -t727 * t583 - t723 * t736;
t475 = (qJD(5) - t713) * t604 + t746;
t557 = (qJD(4) - t720) * t649 + t745;
t600 = t602 ^ 2;
t601 = t604 ^ 2;
t645 = t647 ^ 2;
t646 = t649 ^ 2;
t685 = t689 ^ 2;
t686 = t691 ^ 2;
t708 = t713 ^ 2;
t716 = t720 ^ 2;
t792 = pkin(2) * t722;
t791 = t722 * g(3);
t790 = qJDD(1) * pkin(1);
t789 = t713 * t723;
t788 = t713 * t727;
t787 = t718 * t732;
t712 = t719 * t732;
t786 = t720 * t724;
t785 = t720 * t728;
t636 = -t791 + (-pkin(6) * qJDD(1) + t732 * t792 - t801) * t721;
t665 = -t721 * g(3) + t722 * t801;
t641 = -pkin(2) * t712 + pkin(6) * t757 + t665;
t588 = -t729 * t636 + t725 * t641;
t589 = t725 * t636 + t729 * t641;
t529 = -t729 * t588 + t725 * t589;
t784 = t721 * t529;
t783 = t721 * t722;
t781 = t722 * t529;
t702 = t726 * g(1) - t730 * g(2);
t741 = -qJDD(2) + t702;
t653 = (pkin(1) + t792) * qJDD(1) + (t796 * pkin(6) + qJ(2)) * t732 + t741;
t740 = qJD(3) * pkin(3) - t691 * pkin(7);
t577 = t658 * pkin(3) + t685 * pkin(7) - t691 * t740 + t653;
t743 = t720 * pkin(4) - t649 * pkin(8);
t486 = t583 * pkin(4) + t645 * pkin(8) - t649 * t743 + t577;
t780 = t723 * t486;
t540 = t549 + t744;
t779 = t723 * t540;
t537 = (-t660 - t762) * pkin(7) + t794 * pkin(3) - t588;
t548 = -t685 * pkin(3) + t658 * pkin(7) - qJD(3) * t740 + t589;
t488 = -t728 * t537 + t724 * t548;
t451 = t798 * pkin(4) + t561 * pkin(8) - t488;
t489 = t724 * t537 + t728 * t548;
t454 = -t645 * pkin(4) + t583 * pkin(8) - t720 * t743 + t489;
t403 = -t727 * t451 + t723 * t454;
t404 = t723 * t451 + t727 * t454;
t373 = -t727 * t403 + t723 * t404;
t778 = t724 * t373;
t777 = t724 * t577;
t607 = t610 + t717;
t776 = t724 * t607;
t436 = -t728 * t488 + t724 * t489;
t775 = t725 * t436;
t774 = t725 * t653;
t655 = qJDD(3) + t663;
t773 = t725 * t655;
t683 = t732 * qJ(2) + t741 + t790;
t772 = t726 * t683;
t771 = t727 * t486;
t770 = t727 * t540;
t769 = t728 * t373;
t768 = t728 * t577;
t767 = t728 * t607;
t766 = t729 * t436;
t765 = t729 * t653;
t764 = t729 * t655;
t763 = t730 * t683;
t756 = t726 * qJDD(1);
t755 = t730 * qJDD(1);
t753 = t726 * t549;
t752 = t726 * t610;
t751 = t726 * t663;
t750 = t730 * t549;
t749 = t730 * t610;
t748 = t730 * t663;
t747 = t683 + t790;
t374 = t723 * t403 + t727 * t404;
t437 = t724 * t488 + t728 * t489;
t530 = t725 * t588 + t729 * t589;
t664 = t721 * t801 + t791;
t619 = t721 * t664 + t722 * t665;
t673 = -t726 * t702 - t730 * t703;
t701 = -t726 * t732 + t755;
t742 = -pkin(5) * t701 - t726 * g(3);
t618 = t722 * t664 - t721 * t665;
t672 = t730 * t702 - t726 * t703;
t700 = t730 * t732 + t756;
t694 = t722 * t793;
t669 = -t726 * t694 + t722 * t755;
t738 = t730 * t694 + t722 * t756;
t731 = qJD(3) ^ 2;
t711 = t719 * qJDD(1);
t710 = t718 * qJDD(1);
t699 = t712 - t787;
t698 = t712 + t787;
t697 = t711 - t710;
t696 = t711 + t710;
t693 = t721 * t793;
t684 = -pkin(5) * t700 + t730 * g(3);
t678 = -t686 - t731;
t677 = -t686 + t731;
t676 = t685 - t731;
t675 = t701 * t783;
t674 = t700 * t783;
t670 = t730 * t693 + t721 * t756;
t668 = t726 * t693 - t721 * t755;
t667 = t730 * t696 - t726 * t698;
t666 = t726 * t696 + t730 * t698;
t662 = -t686 + t685;
t659 = t795 - 0.2e1 * t762;
t657 = t687 + 0.2e1 * t761;
t652 = -t731 - t685;
t644 = (-t689 * t729 + t691 * t725) * qJD(3);
t643 = (-t689 * t725 - t691 * t729) * qJD(3);
t638 = -t646 + t716;
t637 = t645 - t716;
t634 = -t646 - t716;
t631 = -t685 - t686;
t629 = t729 * t660 - t725 * t761;
t628 = t725 * t660 + t729 * t761;
t627 = -t725 * t658 + t729 * t762;
t626 = t729 * t658 + t725 * t762;
t625 = -t725 * t678 - t764;
t624 = -t725 * t677 + t802;
t623 = t729 * t676 - t773;
t622 = t729 * t678 - t773;
t621 = t729 * t677 + t805;
t620 = t725 * t676 + t764;
t616 = -t729 * t657 - t725 * t659;
t615 = -t687 * t729 + t725 * t795;
t614 = -t725 * t657 + t729 * t659;
t613 = -t687 * t725 - t729 * t795;
t612 = t729 * t652 - t805;
t611 = t725 * t652 + t802;
t609 = -t646 + t645;
t605 = -t716 - t645;
t596 = -t721 * t643 + t722 * t644;
t595 = -t601 + t708;
t594 = t600 - t708;
t593 = (-t647 * t728 + t649 * t724) * t720;
t592 = (-t647 * t724 - t649 * t728) * t720;
t591 = t730 * t619 - t772;
t590 = t726 * t619 + t763;
t587 = -t601 - t708;
t585 = -pkin(6) * t622 - t765;
t582 = -pkin(6) * t611 - t774;
t581 = -t645 - t646;
t579 = -t721 * t628 + t722 * t629;
t578 = -t721 * t626 + t722 * t627;
t576 = -t721 * t622 + t722 * t625;
t575 = -t721 * t621 + t722 * t624;
t574 = -t721 * t620 + t722 * t623;
t573 = t722 * t622 + t721 * t625;
t572 = t728 * t637 - t776;
t571 = -t724 * t638 + t803;
t570 = t724 * t637 + t767;
t569 = t728 * t638 + t806;
t568 = -pkin(2) * t659 + pkin(6) * t625 - t774;
t567 = -t724 * t634 - t767;
t566 = t728 * t634 - t776;
t565 = -pkin(2) * t657 + pkin(6) * t612 + t765;
t564 = -t721 * t614 + t722 * t616;
t563 = -t721 * t613 + t722 * t615;
t562 = t722 * t613 + t721 * t615;
t556 = (qJD(4) + t720) * t649 + t745;
t555 = -t721 * t611 + t722 * t612;
t554 = t722 * t611 + t721 * t612;
t553 = -t649 * t786 - t728 * t736;
t552 = t649 * t785 - t724 * t736;
t551 = -t724 * t583 + t647 * t785;
t550 = t728 * t583 + t647 * t786;
t547 = t730 * t576 + t726 * t659;
t546 = t726 * t576 - t730 * t659;
t545 = -t601 + t600;
t544 = t728 * t605 - t806;
t543 = t724 * t605 + t803;
t538 = -t708 - t600;
t536 = (-t602 * t727 + t604 * t723) * t713;
t535 = (-t602 * t723 - t604 * t727) * t713;
t532 = -t725 * t592 + t729 * t593;
t531 = t729 * t592 + t725 * t593;
t528 = t730 * t555 + t726 * t657;
t527 = t726 * t555 - t730 * t657;
t526 = t730 * t563 + t726 * t631;
t525 = t726 * t563 - t730 * t631;
t524 = pkin(2) * t653 + pkin(6) * t530;
t523 = -pkin(1) * t562 - pkin(2) * t613;
t522 = -t600 - t601;
t521 = -pkin(7) * t566 - t768;
t520 = -pkin(6) * t613 - t529;
t519 = -t725 * t570 + t729 * t572;
t518 = -t725 * t569 + t729 * t571;
t517 = t729 * t570 + t725 * t572;
t516 = t729 * t569 + t725 * t571;
t515 = t727 * t594 - t779;
t514 = -t723 * t595 + t804;
t513 = t723 * t594 + t770;
t512 = t727 * t595 + t807;
t511 = -t723 * t587 - t770;
t510 = t727 * t587 - t779;
t509 = -t725 * t566 + t729 * t567;
t508 = t729 * t566 + t725 * t567;
t507 = -pkin(1) * t573 - pkin(2) * t622 + t589;
t506 = -pkin(7) * t543 - t777;
t505 = -pkin(2) * t631 + pkin(6) * t615 + t530;
t504 = -t557 * t728 - t724 * t561;
t503 = -t728 * t556 - t724 * t797;
t502 = -t557 * t724 + t728 * t561;
t501 = -t724 * t556 + t728 * t797;
t499 = -t604 * qJD(5) - t746;
t498 = -pkin(1) * t554 - pkin(2) * t611 + t588;
t497 = -t725 * t552 + t729 * t553;
t496 = -t725 * t550 + t729 * t551;
t495 = t729 * t552 + t725 * t553;
t494 = t729 * t550 + t725 * t551;
t493 = -t725 * t543 + t729 * t544;
t492 = t729 * t543 + t725 * t544;
t491 = t727 * t538 - t807;
t490 = t723 * t538 + t804;
t485 = -qJ(2) * t573 - t721 * t568 + t722 * t585;
t484 = -t724 * t535 + t728 * t536;
t483 = t728 * t535 + t724 * t536;
t482 = -t721 * t531 + t722 * t532;
t481 = t722 * t530 - t784;
t480 = t721 * t530 + t781;
t479 = -t597 + t735;
t474 = (qJD(5) + t713) * t604 + t746;
t473 = -t604 * t789 - t727 * t735;
t472 = t604 * t788 - t723 * t735;
t471 = -t723 * t499 + t602 * t788;
t470 = t727 * t499 + t602 * t789;
t469 = -pkin(3) * t797 + pkin(7) * t567 - t777;
t468 = -qJ(2) * t554 - t721 * t565 + t722 * t582;
t467 = t730 * t481 - t726 * t653;
t466 = t726 * t481 + t730 * t653;
t465 = -pkin(3) * t556 + pkin(7) * t544 + t768;
t464 = -t724 * t513 + t728 * t515;
t463 = -t724 * t512 + t728 * t514;
t462 = t728 * t513 + t724 * t515;
t461 = t728 * t512 + t724 * t514;
t460 = -t721 * t517 + t722 * t519;
t459 = -t721 * t516 + t722 * t518;
t458 = -t724 * t510 + t728 * t511;
t457 = t728 * t510 + t724 * t511;
t456 = -t721 * t508 + t722 * t509;
t455 = t722 * t508 + t721 * t509;
t452 = -pkin(1) * t480 - pkin(2) * t529;
t448 = -pkin(8) * t510 - t771;
t447 = -t725 * t502 + t729 * t504;
t446 = -t725 * t501 + t729 * t503;
t445 = t729 * t502 + t725 * t504;
t444 = t729 * t501 + t725 * t503;
t443 = -t721 * t495 + t722 * t497;
t442 = -t721 * t494 + t722 * t496;
t441 = -t721 * t492 + t722 * t493;
t440 = t722 * t492 + t721 * t493;
t439 = -t724 * t490 + t728 * t491;
t438 = t728 * t490 + t724 * t491;
t435 = -pkin(8) * t490 - t780;
t434 = t730 * t456 + t726 * t797;
t433 = t726 * t456 - t730 * t797;
t432 = -qJ(2) * t562 - t721 * t505 + t722 * t520;
t431 = -t725 * t483 + t729 * t484;
t430 = t729 * t483 + t725 * t484;
t429 = -t475 * t727 - t723 * t479;
t428 = -t727 * t474 - t723 * t799;
t427 = -t475 * t723 + t727 * t479;
t426 = -t723 * t474 + t727 * t799;
t425 = -pkin(6) * t781 - qJ(2) * t480 - t721 * t524;
t424 = -t724 * t472 + t728 * t473;
t423 = -t724 * t470 + t728 * t471;
t422 = t728 * t472 + t724 * t473;
t421 = t728 * t470 + t724 * t471;
t420 = pkin(3) * t577 + pkin(7) * t437;
t419 = t730 * t441 + t726 * t556;
t418 = t726 * t441 - t730 * t556;
t417 = -pkin(6) * t508 - t725 * t469 + t729 * t521;
t416 = -pkin(7) * t502 - t436;
t415 = -pkin(4) * t799 + pkin(8) * t511 - t780;
t414 = -t725 * t462 + t729 * t464;
t413 = -t725 * t461 + t729 * t463;
t412 = t729 * t462 + t725 * t464;
t411 = t729 * t461 + t725 * t463;
t410 = -pkin(6) * t492 - t725 * t465 + t729 * t506;
t409 = -t725 * t457 + t729 * t458;
t408 = t729 * t457 + t725 * t458;
t407 = -pkin(2) * t797 + pkin(6) * t509 + t729 * t469 + t725 * t521;
t406 = -pkin(4) * t474 + pkin(8) * t491 + t771;
t405 = -pkin(3) * t581 + pkin(7) * t504 + t437;
t401 = -pkin(2) * t556 + pkin(6) * t493 + t729 * t465 + t725 * t506;
t400 = -t721 * t445 + t722 * t447;
t399 = -t721 * t444 + t722 * t446;
t398 = t722 * t445 + t721 * t447;
t397 = -pkin(1) * t455 - pkin(2) * t508 - pkin(3) * t566 + t489;
t396 = t730 * t400 + t726 * t581;
t395 = t726 * t400 - t730 * t581;
t394 = -t725 * t438 + t729 * t439;
t393 = t729 * t438 + t725 * t439;
t392 = t729 * t437 - t775;
t391 = t725 * t437 + t766;
t390 = -t721 * t430 + t722 * t431;
t389 = -pkin(1) * t440 - pkin(2) * t492 - pkin(3) * t543 + t488;
t388 = -t724 * t427 + t728 * t429;
t387 = -t724 * t426 + t728 * t428;
t386 = t728 * t427 + t724 * t429;
t385 = t728 * t426 + t724 * t428;
t384 = -t725 * t422 + t729 * t424;
t383 = -t725 * t421 + t729 * t423;
t382 = t729 * t422 + t725 * t424;
t381 = t729 * t421 + t725 * t423;
t380 = -t721 * t412 + t722 * t414;
t379 = -t721 * t411 + t722 * t413;
t378 = -t721 * t408 + t722 * t409;
t377 = t722 * t408 + t721 * t409;
t376 = -pkin(7) * t457 - t724 * t415 + t728 * t448;
t375 = -pkin(1) * t398 - pkin(2) * t445 - pkin(3) * t502;
t372 = -pkin(7) * t438 - t724 * t406 + t728 * t435;
t371 = t730 * t378 + t726 * t799;
t370 = t726 * t378 - t730 * t799;
t369 = -qJ(2) * t455 - t721 * t407 + t722 * t417;
t368 = -pkin(3) * t799 + pkin(7) * t458 + t728 * t415 + t724 * t448;
t367 = pkin(4) * t486 + pkin(8) * t374;
t366 = -pkin(6) * t445 - t725 * t405 + t729 * t416;
t365 = -t721 * t393 + t722 * t394;
t364 = t722 * t393 + t721 * t394;
t363 = -t721 * t391 + t722 * t392;
t362 = t722 * t391 + t721 * t392;
t361 = -pkin(2) * t581 + pkin(6) * t447 + t729 * t405 + t725 * t416;
t360 = -pkin(3) * t474 + pkin(7) * t439 + t728 * t406 + t724 * t435;
t359 = -pkin(6) * t391 - pkin(7) * t766 - t725 * t420;
t358 = t730 * t363 - t726 * t577;
t357 = t726 * t363 + t730 * t577;
t356 = -qJ(2) * t440 - t721 * t401 + t722 * t410;
t355 = pkin(2) * t577 + pkin(6) * t392 - pkin(7) * t775 + t729 * t420;
t354 = -pkin(8) * t427 - t373;
t353 = t730 * t365 + t726 * t474;
t352 = t726 * t365 - t730 * t474;
t351 = -t725 * t386 + t729 * t388;
t350 = -t725 * t385 + t729 * t387;
t349 = t729 * t386 + t725 * t388;
t348 = t729 * t385 + t725 * t387;
t347 = -t721 * t382 + t722 * t384;
t346 = -t721 * t381 + t722 * t383;
t345 = -pkin(4) * t522 + pkin(8) * t429 + t374;
t344 = t728 * t374 - t778;
t343 = t724 * t374 + t769;
t342 = -pkin(1) * t362 - pkin(2) * t391 - pkin(3) * t436;
t341 = -pkin(1) * t377 - pkin(2) * t408 - pkin(3) * t457 - pkin(4) * t510 + t404;
t340 = -pkin(6) * t408 - t725 * t368 + t729 * t376;
t339 = -pkin(2) * t799 + pkin(6) * t409 + t729 * t368 + t725 * t376;
t338 = -pkin(1) * t364 - pkin(2) * t393 - pkin(3) * t438 - pkin(4) * t490 + t403;
t337 = -t721 * t349 + t722 * t351;
t336 = -t721 * t348 + t722 * t350;
t335 = t722 * t349 + t721 * t351;
t334 = t730 * t337 + t726 * t522;
t333 = t726 * t337 - t730 * t522;
t332 = -pkin(6) * t393 - t725 * t360 + t729 * t372;
t331 = -qJ(2) * t398 - t721 * t361 + t722 * t366;
t330 = -pkin(2) * t474 + pkin(6) * t394 + t729 * t360 + t725 * t372;
t329 = -pkin(7) * t386 - t724 * t345 + t728 * t354;
t328 = -pkin(3) * t522 + pkin(7) * t388 + t728 * t345 + t724 * t354;
t327 = -qJ(2) * t362 - t721 * t355 + t722 * t359;
t326 = -t725 * t343 + t729 * t344;
t325 = t729 * t343 + t725 * t344;
t324 = -pkin(7) * t343 - pkin(8) * t769 - t724 * t367;
t323 = pkin(3) * t486 + pkin(7) * t344 - pkin(8) * t778 + t728 * t367;
t322 = -pkin(1) * t335 - pkin(2) * t349 - pkin(3) * t386 - pkin(4) * t427;
t321 = -qJ(2) * t377 - t721 * t339 + t722 * t340;
t320 = -qJ(2) * t364 - t721 * t330 + t722 * t332;
t319 = -t721 * t325 + t722 * t326;
t318 = t722 * t325 + t721 * t326;
t317 = t730 * t319 - t726 * t486;
t316 = t726 * t319 + t730 * t486;
t315 = -pkin(6) * t349 - t725 * t328 + t729 * t329;
t314 = -pkin(2) * t522 + pkin(6) * t351 + t729 * t328 + t725 * t329;
t313 = -pkin(6) * t325 - t725 * t323 + t729 * t324;
t312 = pkin(2) * t486 + pkin(6) * t326 + t729 * t323 + t725 * t324;
t311 = -pkin(1) * t318 - pkin(2) * t325 - pkin(3) * t343 - pkin(4) * t373;
t310 = -qJ(2) * t335 - t721 * t314 + t722 * t315;
t309 = -qJ(2) * t318 - t721 * t312 + t722 * t313;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t700, -t701, 0, t673, 0, 0, 0, 0, 0, 0, -t738, t670, t667, t591, 0, 0, 0, 0, 0, 0, t528, t547, t526, t467, 0, 0, 0, 0, 0, 0, t419, t434, t396, t358, 0, 0, 0, 0, 0, 0, t353, t371, t334, t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t701, -t700, 0, t672, 0, 0, 0, 0, 0, 0, t669, t668, t666, t590, 0, 0, 0, 0, 0, 0, t527, t546, t525, t466, 0, 0, 0, 0, 0, 0, t418, t433, t395, t357, 0, 0, 0, 0, 0, 0, t352, t370, t333, t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t618, 0, 0, 0, 0, 0, 0, t554, t573, t562, t480, 0, 0, 0, 0, 0, 0, t440, t455, t398, t362, 0, 0, 0, 0, 0, 0, t364, t377, t335, t318; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t701, 0, -t700, 0, t742, -t684, -t672, -pkin(5) * t672, t675, t730 * t697 - t726 * t699, t670, -t675, t738, 0, -pkin(5) * t669 - t726 * t664 - t721 * t763, -pkin(5) * t668 - t726 * t665 - t722 * t763, -pkin(5) * t666 + t730 * t618, -pkin(5) * t590 - (pkin(1) * t726 - qJ(2) * t730) * t618, t730 * t579 + t751, t730 * t564 - t726 * t662, t730 * t575 + t726 * t795, t730 * t578 - t751, t730 * t574 - t726 * t687, t726 * qJDD(3) + t730 * t596, -pkin(5) * t527 + t730 * t468 - t726 * t498, -pkin(5) * t546 + t730 * t485 - t726 * t507, -pkin(5) * t525 + t730 * t432 - t726 * t523, -pkin(5) * t466 + t730 * t425 - t726 * t452, t730 * t443 + t752, t730 * t399 - t726 * t609, t730 * t459 - t726 * t561, t730 * t442 - t752, t730 * t460 - t726 * t557, t730 * t482 + t726 * t717, -pkin(5) * t418 + t730 * t356 - t726 * t389, -pkin(5) * t433 + t730 * t369 - t726 * t397, -pkin(5) * t395 + t730 * t331 - t726 * t375, -pkin(5) * t357 + t730 * t327 - t726 * t342, t730 * t347 + t753, t730 * t336 - t726 * t545, t730 * t379 - t726 * t479, t730 * t346 - t753, t730 * t380 - t726 * t475, t730 * t390 + t726 * t744, -pkin(5) * t352 + t730 * t320 - t726 * t338, -pkin(5) * t370 + t730 * t321 - t726 * t341, -pkin(5) * t333 + t730 * t310 - t726 * t322, -pkin(5) * t316 + t730 * t309 - t726 * t311; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t700, 0, t701, 0, t684, t742, t673, pkin(5) * t673, t674, t726 * t697 + t730 * t699, t668, -t674, -t669, 0, -pkin(5) * t738 + t730 * t664 - t721 * t772, pkin(5) * t670 + t730 * t665 - t722 * t772, pkin(5) * t667 + t726 * t618, pkin(5) * t591 - (-pkin(1) * t730 - qJ(2) * t726) * t618, t726 * t579 - t748, t726 * t564 + t730 * t662, t726 * t575 - t730 * t795, t726 * t578 + t748, t726 * t574 + t730 * t687, -t730 * qJDD(3) + t726 * t596, pkin(5) * t528 + t726 * t468 + t730 * t498, pkin(5) * t547 + t726 * t485 + t730 * t507, pkin(5) * t526 + t726 * t432 + t730 * t523, pkin(5) * t467 + t726 * t425 + t730 * t452, t726 * t443 - t749, t726 * t399 + t730 * t609, t726 * t459 + t730 * t561, t726 * t442 + t749, t726 * t460 + t730 * t557, t726 * t482 - t730 * t717, pkin(5) * t419 + t726 * t356 + t730 * t389, pkin(5) * t434 + t726 * t369 + t730 * t397, pkin(5) * t396 + t726 * t331 + t730 * t375, pkin(5) * t358 + t726 * t327 + t730 * t342, t726 * t347 - t750, t726 * t336 + t730 * t545, t726 * t379 + t730 * t479, t726 * t346 + t750, t726 * t380 + t730 * t475, t726 * t390 - t730 * t744, pkin(5) * t353 + t726 * t320 + t730 * t338, pkin(5) * t371 + t726 * t321 + t730 * t341, pkin(5) * t334 + t726 * t310 + t730 * t322, pkin(5) * t317 + t726 * t309 + t730 * t311; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t702, t703, 0, 0, t710, 0.2e1 * t721 * t757, 0, t711, 0, 0, -qJ(2) * t694 + t722 * t747, qJ(2) * t693 - t721 * t747, pkin(1) * t698 + qJ(2) * t696 + t619, pkin(1) * t683 + qJ(2) * t619, t722 * t628 + t721 * t629, t722 * t614 + t721 * t616, t722 * t621 + t721 * t624, t722 * t626 + t721 * t627, t722 * t620 + t721 * t623, t722 * t643 + t721 * t644, -pkin(1) * t657 + qJ(2) * t555 + t722 * t565 + t721 * t582, -pkin(1) * t659 + qJ(2) * t576 + t722 * t568 + t721 * t585, -pkin(1) * t631 + qJ(2) * t563 + t722 * t505 + t721 * t520, pkin(1) * t653 - pkin(6) * t784 + qJ(2) * t481 + t722 * t524, t722 * t495 + t721 * t497, t722 * t444 + t721 * t446, t722 * t516 + t721 * t518, t722 * t494 + t721 * t496, t722 * t517 + t721 * t519, t722 * t531 + t721 * t532, -pkin(1) * t556 + qJ(2) * t441 + t722 * t401 + t721 * t410, -pkin(1) * t797 + qJ(2) * t456 + t722 * t407 + t721 * t417, -pkin(1) * t581 + qJ(2) * t400 + t722 * t361 + t721 * t366, pkin(1) * t577 + qJ(2) * t363 + t722 * t355 + t721 * t359, t722 * t382 + t721 * t384, t722 * t348 + t721 * t350, t722 * t411 + t721 * t413, t722 * t381 + t721 * t383, t722 * t412 + t721 * t414, t722 * t430 + t721 * t431, -pkin(1) * t474 + qJ(2) * t365 + t722 * t330 + t721 * t332, -pkin(1) * t799 + qJ(2) * t378 + t722 * t339 + t721 * t340, -pkin(1) * t522 + qJ(2) * t337 + t722 * t314 + t721 * t315, pkin(1) * t486 + qJ(2) * t319 + t722 * t312 + t721 * t313;];
tauB_reg = t1;
