% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:28
% EndTime: 2019-12-31 20:18:44
% DurationCPUTime: 15.52s
% Computational Cost: add. (90146->679), mult. (211057->1037), div. (0->0), fcn. (154572->10), ass. (0->467)
t759 = cos(qJ(2));
t743 = t759 * qJDD(1);
t755 = sin(qJ(2));
t796 = qJD(1) * qJD(2);
t783 = t755 * t796;
t719 = t743 - t783;
t750 = t759 ^ 2;
t762 = qJD(1) ^ 2;
t803 = qJD(1) * t755;
t771 = qJD(2) * pkin(2) - qJ(3) * t803;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t728 = t756 * g(1) - t760 * g(2);
t772 = qJDD(1) * pkin(1) + t728;
t641 = t719 * pkin(2) - t771 * t803 - qJDD(3) + t772 + (qJ(3) * t750 + pkin(6)) * t762;
t751 = sin(pkin(9));
t752 = cos(pkin(9));
t708 = t752 * t759 * qJD(1) - t751 * t803;
t709 = (t751 * t759 + t752 * t755) * qJD(1);
t680 = t708 * t709;
t836 = qJDD(2) + t680;
t847 = t751 * t836;
t846 = t752 * t836;
t753 = sin(qJ(5));
t754 = sin(qJ(4));
t758 = cos(qJ(4));
t672 = t754 * t708 + t758 * t709;
t748 = qJD(2) + qJD(4);
t757 = cos(qJ(5));
t642 = t753 * t672 - t757 * t748;
t644 = t757 * t672 + t753 * t748;
t605 = t644 * t642;
t782 = t759 * t796;
t795 = t755 * qJDD(1);
t718 = t782 + t795;
t682 = t752 * t718 + t751 * t719;
t778 = t751 * t718 - t752 * t719;
t779 = t754 * t682 + t758 * t778;
t598 = -t672 * qJD(4) - t779;
t770 = qJDD(5) - t598;
t840 = -t605 + t770;
t845 = t753 * t840;
t670 = -t758 * t708 + t754 * t709;
t618 = t672 * t670;
t747 = qJDD(2) + qJDD(4);
t839 = -t618 + t747;
t844 = t754 * t839;
t843 = t757 * t840;
t842 = t758 * t839;
t659 = t748 * t670;
t768 = t670 * qJD(4) - t758 * t682 + t754 * t778;
t841 = t659 + t768;
t706 = t708 ^ 2;
t774 = qJD(2) * pkin(3) - t709 * pkin(7);
t579 = -pkin(3) * t778 + t706 * pkin(7) - t709 * t774 + t641;
t558 = -t642 * qJD(5) + t753 * t747 - t757 * t768;
t667 = qJD(5) + t670;
t611 = t667 * t642;
t535 = -t611 + t558;
t745 = t750 * t762;
t761 = qJD(2) ^ 2;
t734 = -t745 - t761;
t812 = t755 * t762;
t729 = t760 * g(1) + t756 * g(2);
t712 = -t762 * pkin(1) + qJDD(1) * pkin(6) - t729;
t815 = t755 * t712;
t837 = (pkin(2) * t812 + qJ(3) * t796 - g(3)) * t759 + qJDD(2) * pkin(2) - t718 * qJ(3) - t815;
t780 = -t757 * t747 - t753 * t768;
t532 = (qJD(5) - t667) * t644 + t780;
t569 = (qJD(4) - t748) * t672 + t779;
t639 = t642 ^ 2;
t640 = t644 ^ 2;
t666 = t667 ^ 2;
t668 = t670 ^ 2;
t669 = t672 ^ 2;
t707 = t709 ^ 2;
t835 = t748 ^ 2;
t834 = 2 * qJD(3);
t833 = pkin(4) * t754;
t832 = t667 * t753;
t831 = t667 * t757;
t830 = t748 * t754;
t829 = t748 * t758;
t749 = t755 ^ 2;
t828 = t749 * t762;
t696 = -t755 * g(3) + t759 * t712;
t764 = -pkin(2) * t745 + t719 * qJ(3) - qJD(2) * t771 + t696;
t581 = t708 * t834 + t751 * t837 + t752 * t764;
t549 = -t706 * pkin(3) - pkin(7) * t778 - qJD(2) * t774 + t581;
t791 = t709 * t834;
t580 = t751 * t764 - t752 * t837 + t791;
t802 = qJD(2) * t708;
t650 = -t682 + t802;
t763 = pkin(3) * t836 + pkin(7) * t650 - t580;
t481 = t754 * t549 - t758 * t763;
t482 = t758 * t549 + t754 * t763;
t427 = -t758 * t481 + t754 * t482;
t827 = t751 * t427;
t826 = t751 * t641;
t676 = qJDD(2) - t680;
t825 = t751 * t676;
t824 = t752 * t427;
t823 = t752 * t641;
t822 = t752 * t676;
t616 = t670 * pkin(4) - t672 * pkin(8);
t465 = -t747 * pkin(4) - t835 * pkin(8) + t672 * t616 + t481;
t821 = t753 * t465;
t546 = t605 + t770;
t820 = t753 * t546;
t819 = t754 * t579;
t614 = t618 + t747;
t818 = t754 * t614;
t520 = -t752 * t580 + t751 * t581;
t817 = t755 * t520;
t711 = t762 * pkin(6) + t772;
t816 = t755 * t711;
t735 = t759 * t812;
t726 = qJDD(2) + t735;
t814 = t755 * t726;
t727 = qJDD(2) - t735;
t813 = t755 * t727;
t811 = t757 * t465;
t810 = t757 * t546;
t809 = t758 * t579;
t808 = t758 * t614;
t807 = t759 * t520;
t806 = t759 * t711;
t805 = t759 * t727;
t466 = -t835 * pkin(4) + t747 * pkin(8) - t670 * t616 + t482;
t479 = t841 * pkin(8) + (t748 * t672 - t598) * pkin(4) - t579;
t422 = t757 * t466 + t753 * t479;
t804 = t749 + t750;
t801 = qJD(2) * t709;
t800 = qJD(2) * t751;
t799 = qJD(2) * t752;
t794 = t756 * qJDD(1);
t793 = t760 * qJDD(1);
t792 = t760 * qJDD(2);
t790 = t754 * t605;
t789 = t756 * t618;
t788 = t756 * t680;
t787 = t758 * t605;
t786 = t760 * t618;
t785 = t760 * t680;
t784 = -pkin(4) * t758 - pkin(3);
t421 = t753 * t466 - t757 * t479;
t428 = t754 * t481 + t758 * t482;
t521 = t751 * t580 + t752 * t581;
t695 = t759 * g(3) + t815;
t637 = t755 * t695 + t759 * t696;
t688 = -t756 * t728 - t760 * t729;
t777 = t756 * t735;
t776 = t760 * t735;
t723 = -t756 * t762 + t793;
t775 = -pkin(5) * t723 - t756 * g(3);
t386 = -t757 * t421 + t753 * t422;
t387 = t753 * t421 + t757 * t422;
t636 = t759 * t695 - t755 * t696;
t687 = t760 * t728 - t756 * t729;
t648 = -t778 + t801;
t742 = t756 * qJDD(2);
t733 = t745 - t761;
t732 = -t761 - t828;
t731 = t761 - t828;
t725 = t745 - t828;
t724 = t745 + t828;
t722 = t760 * t762 + t794;
t721 = t804 * qJDD(1);
t720 = t743 - 0.2e1 * t783;
t717 = 0.2e1 * t782 + t795;
t715 = t759 * t726;
t714 = t804 * t796;
t705 = -pkin(5) * t722 + t760 * g(3);
t701 = -t707 - t761;
t700 = -t707 + t761;
t699 = t706 - t761;
t698 = t759 * t718 - t749 * t796;
t697 = -t755 * t719 - t750 * t796;
t694 = -t755 * t732 - t805;
t693 = -t755 * t731 + t715;
t692 = t759 * t734 - t814;
t691 = t759 * t733 - t813;
t690 = t759 * t732 - t813;
t689 = t755 * t734 + t715;
t685 = t760 * t721 - t756 * t724;
t684 = t756 * t721 + t760 * t724;
t683 = -t755 * t717 + t759 * t720;
t679 = -t707 + t706;
t674 = -t761 - t706;
t665 = (t708 * t752 + t709 * t751) * qJD(2);
t664 = (t708 * t751 - t709 * t752) * qJD(2);
t663 = t760 * t694 + t756 * t717;
t662 = t760 * t692 - t756 * t720;
t661 = t756 * t694 - t760 * t717;
t660 = t756 * t692 + t760 * t720;
t657 = -t669 + t835;
t656 = t668 - t835;
t655 = -pkin(6) * t690 - t806;
t654 = -pkin(6) * t689 - t816;
t653 = -t669 - t835;
t652 = -pkin(1) * t690 + t696;
t651 = -pkin(1) * t689 + t695;
t649 = t682 + t802;
t646 = t778 + t801;
t645 = -t706 - t707;
t634 = t752 * t682 - t709 * t800;
t633 = t751 * t682 + t709 * t799;
t632 = -t708 * t799 + t751 * t778;
t631 = -t708 * t800 - t752 * t778;
t628 = -t751 * t701 - t822;
t627 = -t751 * t700 + t846;
t626 = t752 * t699 - t825;
t625 = t752 * t701 - t825;
t624 = t752 * t700 + t847;
t623 = t751 * t699 + t822;
t622 = t760 * t637 - t756 * t711;
t621 = t756 * t637 + t760 * t711;
t620 = t752 * t674 - t847;
t619 = t751 * t674 + t846;
t617 = -t669 + t668;
t612 = -t835 - t668;
t610 = -t640 + t666;
t609 = t639 - t666;
t608 = -t755 * t664 + t759 * t665;
t607 = (-t670 * t758 + t672 * t754) * t748;
t606 = (-t670 * t754 - t672 * t758) * t748;
t604 = -t640 + t639;
t603 = t752 * t648 - t751 * t650;
t602 = -t752 * t646 - t751 * t649;
t601 = t751 * t648 + t752 * t650;
t600 = -t751 * t646 + t752 * t649;
t595 = -qJ(3) * t625 - t823;
t594 = -t668 - t669;
t593 = -t755 * t633 + t759 * t634;
t592 = -t755 * t631 + t759 * t632;
t591 = -t755 * t625 + t759 * t628;
t590 = -t755 * t624 + t759 * t627;
t589 = -t755 * t623 + t759 * t626;
t588 = t759 * t625 + t755 * t628;
t587 = -t640 - t666;
t586 = -qJ(3) * t619 - t826;
t585 = t758 * t656 - t818;
t584 = -t754 * t657 + t842;
t583 = t754 * t656 + t808;
t582 = t758 * t657 + t844;
t578 = -t754 * t653 - t808;
t577 = t758 * t653 - t818;
t575 = -t666 - t639;
t574 = t639 + t640;
t573 = -t659 + t768;
t568 = (qJD(4) + t748) * t672 + t779;
t567 = -t672 * t830 - t758 * t768;
t566 = t672 * t829 - t754 * t768;
t565 = -t754 * t598 + t670 * t829;
t564 = t758 * t598 + t670 * t830;
t563 = -t755 * t619 + t759 * t620;
t562 = t759 * t619 + t755 * t620;
t561 = -pkin(2) * t649 + qJ(3) * t628 - t826;
t560 = t758 * t612 - t844;
t559 = t754 * t612 + t842;
t557 = -t644 * qJD(5) - t780;
t556 = -pkin(2) * t646 + qJ(3) * t620 + t823;
t555 = t760 * t591 + t756 * t649;
t554 = t756 * t591 - t760 * t649;
t553 = (-t642 * t757 + t644 * t753) * t667;
t552 = (t642 * t753 + t644 * t757) * t667;
t551 = -t751 * t606 + t752 * t607;
t550 = t752 * t606 + t751 * t607;
t544 = t760 * t563 + t756 * t646;
t543 = -t755 * t601 + t759 * t603;
t542 = -t755 * t600 + t759 * t602;
t541 = t756 * t563 - t760 * t646;
t540 = t759 * t601 + t755 * t603;
t536 = -t611 - t558;
t533 = (-qJD(5) - t667) * t644 - t780;
t531 = t757 * t558 - t644 * t832;
t530 = -t753 * t558 - t644 * t831;
t529 = -t753 * t557 + t642 * t831;
t528 = -t757 * t557 - t642 * t832;
t527 = t760 * t543 + t756 * t645;
t526 = t756 * t543 - t760 * t645;
t525 = -t751 * t583 + t752 * t585;
t524 = -t751 * t582 + t752 * t584;
t523 = t752 * t583 + t751 * t585;
t522 = t752 * t582 + t751 * t584;
t519 = -pkin(7) * t577 - t809;
t518 = -t751 * t577 + t752 * t578;
t517 = t752 * t577 + t751 * t578;
t516 = t758 * t553 + t754 * t770;
t515 = t754 * t553 - t758 * t770;
t514 = t757 * t609 - t820;
t513 = -t753 * t610 + t843;
t512 = -t753 * t609 - t810;
t511 = -t757 * t610 - t845;
t510 = -t569 * t758 - t754 * t573;
t509 = -t758 * t568 + t754 * t841;
t508 = -t569 * t754 + t758 * t573;
t507 = -t754 * t568 - t758 * t841;
t506 = -pkin(7) * t559 - t819;
t505 = -t751 * t566 + t752 * t567;
t504 = -t751 * t564 + t752 * t565;
t503 = t752 * t566 + t751 * t567;
t502 = t752 * t564 + t751 * t565;
t501 = -pkin(1) * t540 - pkin(2) * t601;
t500 = -pkin(1) * t588 - pkin(2) * t625 + t581;
t499 = pkin(2) * t641 + qJ(3) * t521;
t498 = -t753 * t587 - t810;
t497 = t757 * t587 - t820;
t496 = -t751 * t559 + t752 * t560;
t495 = t752 * t559 + t751 * t560;
t494 = t757 * t575 - t845;
t493 = t753 * t575 + t843;
t492 = t758 * t531 + t790;
t491 = t758 * t529 - t790;
t490 = t754 * t531 - t787;
t489 = t754 * t529 + t787;
t488 = -pkin(1) * t562 + t751 * t696 + t752 * t695 + t791 + (t751 * (t719 + t783) - t752 * (-t718 + t782)) * qJ(3) + (-t752 * t726 + t734 * t751 - t619) * pkin(2);
t487 = -qJ(3) * t601 - t520;
t486 = -t755 * t550 + t759 * t551;
t485 = -pkin(6) * t588 - t755 * t561 + t759 * t595;
t484 = pkin(3) * t841 + pkin(7) * t578 - t819;
t483 = -pkin(2) * t645 + qJ(3) * t603 + t521;
t476 = -pkin(3) * t568 + pkin(7) * t560 + t809;
t475 = -pkin(6) * t562 - t755 * t556 + t759 * t586;
t474 = -t532 * t757 - t753 * t536;
t473 = t757 * t533 - t753 * t535;
t472 = -t532 * t753 + t757 * t536;
t471 = -t753 * t533 - t757 * t535;
t470 = -t755 * t523 + t759 * t525;
t469 = -t755 * t522 + t759 * t524;
t468 = t759 * t521 - t817;
t467 = t755 * t521 + t807;
t463 = -t755 * t517 + t759 * t518;
t462 = t759 * t517 + t755 * t518;
t461 = t758 * t514 - t754 * t532;
t460 = t758 * t513 - t754 * t536;
t459 = t754 * t514 + t758 * t532;
t458 = t754 * t513 + t758 * t536;
t457 = t760 * t468 - t756 * t641;
t456 = t756 * t468 + t760 * t641;
t455 = -t751 * t515 + t752 * t516;
t454 = t752 * t515 + t751 * t516;
t453 = t758 * t498 + t535 * t754;
t452 = t754 * t498 - t535 * t758;
t451 = t758 * t494 - t754 * t533;
t450 = t754 * t494 + t758 * t533;
t449 = -t751 * t508 + t752 * t510;
t448 = -t751 * t507 + t752 * t509;
t447 = t752 * t508 + t751 * t510;
t446 = t752 * t507 + t751 * t509;
t445 = t758 * t473 - t754 * t604;
t444 = t754 * t473 + t758 * t604;
t443 = -t755 * t503 + t759 * t505;
t442 = -t755 * t502 + t759 * t504;
t441 = t758 * t474 - t754 * t574;
t440 = t754 * t474 + t758 * t574;
t439 = -t755 * t495 + t759 * t496;
t438 = t759 * t495 + t755 * t496;
t437 = t760 * t463 - t756 * t841;
t436 = t756 * t463 + t760 * t841;
t435 = -t751 * t490 + t752 * t492;
t434 = -t751 * t489 + t752 * t491;
t433 = t752 * t490 + t751 * t492;
t432 = t752 * t489 + t751 * t491;
t431 = -pkin(1) * t467 - pkin(2) * t520;
t430 = t760 * t439 + t756 * t568;
t429 = t756 * t439 - t760 * t568;
t426 = -pkin(8) * t497 + t811;
t425 = -pkin(8) * t493 + t821;
t424 = -qJ(3) * t517 - t751 * t484 + t752 * t519;
t423 = pkin(3) * t579 + pkin(7) * t428;
t420 = -pkin(6) * t540 - t755 * t483 + t759 * t487;
t419 = -qJ(3) * t495 - t751 * t476 + t752 * t506;
t418 = pkin(2) * t841 + qJ(3) * t518 + t752 * t484 + t751 * t519;
t417 = -t751 * t459 + t752 * t461;
t416 = -t751 * t458 + t752 * t460;
t415 = t752 * t459 + t751 * t461;
t414 = t752 * t458 + t751 * t460;
t413 = -pkin(6) * t467 - qJ(3) * t807 - t755 * t499;
t412 = -pkin(7) * t508 - t427;
t411 = -t755 * t454 + t759 * t455;
t410 = -t751 * t452 + t752 * t453;
t409 = t752 * t452 + t751 * t453;
t408 = -t751 * t450 + t752 * t451;
t407 = t752 * t450 + t751 * t451;
t406 = -t755 * t447 + t759 * t449;
t405 = -t755 * t446 + t759 * t448;
t404 = t759 * t447 + t755 * t449;
t403 = -t751 * t444 + t752 * t445;
t402 = t752 * t444 + t751 * t445;
t401 = -pkin(2) * t568 + qJ(3) * t496 + t752 * t476 + t751 * t506;
t400 = -pkin(3) * t594 + pkin(7) * t510 + t428;
t399 = -pkin(4) * t497 + t422;
t398 = -pkin(4) * t493 + t421;
t397 = -t751 * t440 + t752 * t441;
t396 = t752 * t440 + t751 * t441;
t395 = t760 * t406 + t756 * t594;
t394 = t756 * t406 - t760 * t594;
t393 = -pkin(1) * t462 - pkin(2) * t517 - pkin(3) * t577 + t482;
t392 = -t755 * t433 + t759 * t435;
t391 = -t755 * t432 + t759 * t434;
t390 = -pkin(1) * t438 - pkin(2) * t495 - pkin(3) * t559 + t481;
t389 = t752 * t428 - t827;
t388 = t751 * t428 + t824;
t385 = -pkin(1) * t404 - pkin(2) * t447 - pkin(3) * t508;
t384 = -t755 * t415 + t759 * t417;
t383 = -t755 * t414 + t759 * t416;
t382 = -t755 * t409 + t759 * t410;
t381 = t759 * t409 + t755 * t410;
t380 = -t755 * t407 + t759 * t408;
t379 = t759 * t407 + t755 * t408;
t378 = -pkin(8) * t472 - t386;
t377 = -t755 * t402 + t759 * t403;
t376 = t758 * t387 + t754 * t465;
t375 = t754 * t387 - t758 * t465;
t374 = -pkin(6) * t462 - t755 * t418 + t759 * t424;
t373 = -t755 * t396 + t759 * t397;
t372 = t759 * t396 + t755 * t397;
t371 = t760 * t382 + t756 * t497;
t370 = t756 * t382 - t760 * t497;
t369 = t760 * t380 + t756 * t493;
t368 = t756 * t380 - t760 * t493;
t367 = -pkin(7) * t452 - t754 * t399 + t758 * t426;
t366 = -pkin(7) * t450 - t754 * t398 + t758 * t425;
t365 = -qJ(3) * t447 - t751 * t400 + t752 * t412;
t364 = -pkin(6) * t438 - t755 * t401 + t759 * t419;
t363 = -pkin(2) * t594 + qJ(3) * t449 + t752 * t400 + t751 * t412;
t362 = -pkin(3) * t497 + pkin(7) * t453 + t758 * t399 + t754 * t426;
t361 = -pkin(3) * t493 + pkin(7) * t451 + t758 * t398 + t754 * t425;
t360 = t760 * t373 + t756 * t472;
t359 = t756 * t373 - t760 * t472;
t358 = -t755 * t388 + t759 * t389;
t357 = t759 * t388 + t755 * t389;
t356 = -pkin(7) * t824 - qJ(3) * t388 - t751 * t423;
t355 = t760 * t358 - t756 * t579;
t354 = t756 * t358 + t760 * t579;
t353 = pkin(2) * t579 - pkin(7) * t827 + qJ(3) * t389 + t752 * t423;
t352 = -pkin(7) * t440 + t758 * t378 + t472 * t833;
t351 = pkin(7) * t441 + t754 * t378 + t472 * t784;
t350 = -t751 * t375 + t752 * t376;
t349 = t752 * t375 + t751 * t376;
t348 = -pkin(1) * t381 - pkin(2) * t409 - pkin(3) * t452 + pkin(4) * t535 - pkin(8) * t498 - t821;
t347 = -pkin(1) * t379 - pkin(2) * t407 - pkin(3) * t450 - pkin(4) * t533 - pkin(8) * t494 + t811;
t346 = -pkin(1) * t357 - pkin(2) * t388 - pkin(3) * t427;
t345 = -pkin(7) * t375 + (-pkin(8) * t758 + t833) * t386;
t344 = -qJ(3) * t409 - t751 * t362 + t752 * t367;
t343 = -qJ(3) * t407 - t751 * t361 + t752 * t366;
t342 = -pkin(1) * t372 - pkin(2) * t396 - pkin(3) * t440 - pkin(4) * t574 - pkin(8) * t474 - t387;
t341 = -pkin(6) * t404 - t755 * t363 + t759 * t365;
t340 = -pkin(2) * t497 + qJ(3) * t410 + t752 * t362 + t751 * t367;
t339 = -pkin(2) * t493 + qJ(3) * t408 + t752 * t361 + t751 * t366;
t338 = pkin(7) * t376 + (-pkin(8) * t754 + t784) * t386;
t337 = -qJ(3) * t396 - t751 * t351 + t752 * t352;
t336 = -pkin(2) * t472 + qJ(3) * t397 + t752 * t351 + t751 * t352;
t335 = -t755 * t349 + t759 * t350;
t334 = t759 * t349 + t755 * t350;
t333 = -pkin(6) * t357 - t755 * t353 + t759 * t356;
t332 = t760 * t335 + t756 * t386;
t331 = t756 * t335 - t760 * t386;
t330 = -pkin(6) * t381 - t755 * t340 + t759 * t344;
t329 = -pkin(6) * t379 - t755 * t339 + t759 * t343;
t328 = -qJ(3) * t349 - t751 * t338 + t752 * t345;
t327 = -pkin(1) * t334 - pkin(2) * t349 - pkin(3) * t375 + pkin(4) * t465 - pkin(8) * t387;
t326 = -pkin(6) * t372 - t755 * t336 + t759 * t337;
t325 = -pkin(2) * t386 + qJ(3) * t350 + t752 * t338 + t751 * t345;
t324 = -pkin(6) * t334 - t755 * t325 + t759 * t328;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t722, -t723, 0, t688, 0, 0, 0, 0, 0, 0, t662, t663, t685, t622, 0, 0, 0, 0, 0, 0, t544, t555, t527, t457, 0, 0, 0, 0, 0, 0, t430, t437, t395, t355, 0, 0, 0, 0, 0, 0, t369, t371, t360, t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t723, -t722, 0, t687, 0, 0, 0, 0, 0, 0, t660, t661, t684, t621, 0, 0, 0, 0, 0, 0, t541, t554, t526, t456, 0, 0, 0, 0, 0, 0, t429, t436, t394, t354, 0, 0, 0, 0, 0, 0, t368, t370, t359, t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t689, t690, 0, -t636, 0, 0, 0, 0, 0, 0, t562, t588, t540, t467, 0, 0, 0, 0, 0, 0, t438, t462, t404, t357, 0, 0, 0, 0, 0, 0, t379, t381, t372, t334; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t723, 0, -t722, 0, t775, -t705, -t687, -pkin(5) * t687, t760 * t698 - t777, t760 * t683 - t756 * t725, t760 * t693 + t755 * t794, t760 * t697 + t777, t760 * t691 + t743 * t756, t760 * t714 + t742, -pkin(5) * t660 - t756 * t651 + t760 * t654, -pkin(5) * t661 - t756 * t652 + t760 * t655, -pkin(5) * t684 + t760 * t636, -pkin(5) * t621 - (pkin(1) * t756 - pkin(6) * t760) * t636, t760 * t593 - t788, t760 * t542 - t756 * t679, t760 * t590 - t756 * t650, t760 * t592 + t788, t760 * t589 + t648 * t756, t760 * t608 + t742, -pkin(5) * t541 + t760 * t475 - t756 * t488, -pkin(5) * t554 + t760 * t485 - t756 * t500, -pkin(5) * t526 + t760 * t420 - t756 * t501, -pkin(5) * t456 + t760 * t413 - t756 * t431, t760 * t443 + t789, t760 * t405 - t756 * t617, t760 * t469 - t756 * t573, t760 * t442 - t789, t760 * t470 - t756 * t569, t760 * t486 + t756 * t747, -pkin(5) * t429 + t760 * t364 - t756 * t390, -pkin(5) * t436 + t760 * t374 - t756 * t393, -pkin(5) * t394 + t760 * t341 - t756 * t385, -pkin(5) * t354 + t760 * t333 - t756 * t346, t760 * t392 - t756 * t530, t760 * t377 - t756 * t471, t760 * t383 - t756 * t511, t760 * t391 - t756 * t528, t760 * t384 - t756 * t512, t760 * t411 - t756 * t552, -pkin(5) * t368 + t760 * t329 - t756 * t347, -pkin(5) * t370 + t760 * t330 - t756 * t348, -pkin(5) * t359 + t760 * t326 - t756 * t342, -pkin(5) * t331 + t760 * t324 - t756 * t327; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t722, 0, t723, 0, t705, t775, t688, pkin(5) * t688, t756 * t698 + t776, t756 * t683 + t760 * t725, t756 * t693 - t755 * t793, t756 * t697 - t776, t756 * t691 - t759 * t793, t756 * t714 - t792, pkin(5) * t662 + t760 * t651 + t756 * t654, pkin(5) * t663 + t760 * t652 + t756 * t655, pkin(5) * t685 + t756 * t636, pkin(5) * t622 - (-pkin(1) * t760 - pkin(6) * t756) * t636, t756 * t593 + t785, t756 * t542 + t760 * t679, t756 * t590 + t760 * t650, t756 * t592 - t785, t756 * t589 - t648 * t760, t756 * t608 - t792, pkin(5) * t544 + t756 * t475 + t760 * t488, pkin(5) * t555 + t756 * t485 + t760 * t500, pkin(5) * t527 + t756 * t420 + t760 * t501, pkin(5) * t457 + t756 * t413 + t760 * t431, t756 * t443 - t786, t756 * t405 + t760 * t617, t756 * t469 + t760 * t573, t756 * t442 + t786, t756 * t470 + t760 * t569, t756 * t486 - t760 * t747, pkin(5) * t430 + t756 * t364 + t760 * t390, pkin(5) * t437 + t756 * t374 + t760 * t393, pkin(5) * t395 + t756 * t341 + t760 * t385, pkin(5) * t355 + t756 * t333 + t760 * t346, t756 * t392 + t760 * t530, t756 * t377 + t760 * t471, t756 * t383 + t760 * t511, t756 * t391 + t760 * t528, t756 * t384 + t760 * t512, t756 * t411 + t760 * t552, pkin(5) * t369 + t756 * t329 + t760 * t347, pkin(5) * t371 + t756 * t330 + t760 * t348, pkin(5) * t360 + t756 * t326 + t760 * t342, pkin(5) * t332 + t756 * t324 + t760 * t327; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t728, t729, 0, 0, (t718 + t782) * t755, t759 * t717 + t755 * t720, t759 * t731 + t814, (t719 - t783) * t759, t755 * t733 + t805, 0, pkin(1) * t720 + pkin(6) * t692 + t806, -pkin(1) * t717 + pkin(6) * t694 - t816, pkin(1) * t724 + pkin(6) * t721 + t637, pkin(1) * t711 + pkin(6) * t637, t759 * t633 + t755 * t634, t759 * t600 + t755 * t602, t759 * t624 + t755 * t627, t759 * t631 + t755 * t632, t759 * t623 + t755 * t626, t759 * t664 + t755 * t665, -pkin(1) * t646 + pkin(6) * t563 + t759 * t556 + t755 * t586, -pkin(1) * t649 + pkin(6) * t591 + t759 * t561 + t755 * t595, -pkin(1) * t645 + pkin(6) * t543 + t759 * t483 + t755 * t487, pkin(1) * t641 + pkin(6) * t468 - qJ(3) * t817 + t759 * t499, t759 * t503 + t755 * t505, t759 * t446 + t755 * t448, t759 * t522 + t755 * t524, t759 * t502 + t755 * t504, t759 * t523 + t755 * t525, t759 * t550 + t755 * t551, -pkin(1) * t568 + pkin(6) * t439 + t759 * t401 + t755 * t419, pkin(1) * t841 + pkin(6) * t463 + t759 * t418 + t755 * t424, -pkin(1) * t594 + pkin(6) * t406 + t759 * t363 + t755 * t365, pkin(1) * t579 + pkin(6) * t358 + t759 * t353 + t755 * t356, t759 * t433 + t755 * t435, t759 * t402 + t755 * t403, t759 * t414 + t755 * t416, t759 * t432 + t755 * t434, t759 * t415 + t755 * t417, t759 * t454 + t755 * t455, -pkin(1) * t493 + pkin(6) * t380 + t759 * t339 + t755 * t343, -pkin(1) * t497 + pkin(6) * t382 + t759 * t340 + t755 * t344, -pkin(1) * t472 + pkin(6) * t373 + t759 * t336 + t755 * t337, -pkin(1) * t386 + pkin(6) * t335 + t759 * t325 + t755 * t328;];
tauB_reg = t1;
