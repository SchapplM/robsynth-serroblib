% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:33
% EndTime: 2019-12-31 21:14:48
% DurationCPUTime: 15.63s
% Computational Cost: add. (95710->678), mult. (216607->1037), div. (0->0), fcn. (158388->10), ass. (0->467)
t790 = cos(qJ(2));
t773 = t790 * qJDD(1);
t786 = sin(qJ(2));
t826 = qJD(1) * qJD(2);
t814 = t786 * t826;
t748 = t773 - t814;
t781 = t790 ^ 2;
t793 = qJD(1) ^ 2;
t787 = sin(qJ(1));
t791 = cos(qJ(1));
t757 = t787 * g(1) - t791 * g(2);
t800 = qJDD(1) * pkin(1) + t757;
t829 = qJD(1) * t786;
t801 = qJD(2) * pkin(2) - pkin(7) * t829;
t679 = t748 * pkin(2) - t801 * t829 + t800 + (pkin(7) * t781 + pkin(6)) * t793;
t782 = sin(pkin(9));
t785 = sin(qJ(3));
t789 = cos(qJ(3));
t737 = t789 * t790 * qJD(1) - t785 * t829;
t738 = (t785 * t790 + t786 * t789) * qJD(1);
t783 = cos(pkin(9));
t700 = -t783 * t737 + t782 * t738;
t702 = t782 * t737 + t783 * t738;
t643 = t702 * t700;
t778 = qJDD(2) + qJDD(3);
t869 = -t643 + t778;
t876 = t782 * t869;
t875 = t783 * t869;
t784 = sin(qJ(5));
t779 = qJD(2) + qJD(3);
t788 = cos(qJ(5));
t671 = t784 * t702 - t788 * t779;
t673 = t788 * t702 + t784 * t779;
t628 = t673 * t671;
t813 = t790 * t826;
t825 = t786 * qJDD(1);
t747 = t813 + t825;
t808 = t785 * t747 - t789 * t748;
t676 = -t738 * qJD(3) - t808;
t677 = t737 * qJD(3) + t789 * t747 + t785 * t748;
t809 = -t783 * t676 + t782 * t677;
t802 = qJDD(5) + t809;
t870 = -t628 + t802;
t874 = t784 * t870;
t708 = t737 * t738;
t867 = t708 + t778;
t873 = t785 * t867;
t872 = t788 * t870;
t871 = t789 * t867;
t735 = t737 ^ 2;
t803 = t779 * pkin(3) - t738 * qJ(4);
t589 = t676 * pkin(3) + t735 * qJ(4) - t738 * t803 - qJDD(4) + t679;
t860 = t779 * t702;
t595 = t809 + t860;
t626 = t782 * t676 + t783 * t677;
t692 = t779 * t700;
t599 = t626 - t692;
t584 = -t671 * qJD(5) + t788 * t626 + t784 * t778;
t696 = qJD(5) + t700;
t634 = t696 * t671;
t559 = -t634 + t584;
t729 = t779 * t737;
t655 = -t677 + t729;
t868 = t677 + t729;
t775 = t781 * t793;
t792 = qJD(2) ^ 2;
t763 = -t775 - t792;
t840 = t786 * t793;
t758 = t791 * g(1) + t787 * g(2);
t740 = -t793 * pkin(1) + qJDD(1) * pkin(6) - t758;
t843 = t786 * t740;
t866 = (pkin(2) * t840 + pkin(7) * t826 - g(3)) * t790 + qJDD(2) * pkin(2) - t747 * pkin(7) - t843;
t810 = t784 * t626 - t788 * t778;
t556 = (qJD(5) - t696) * t673 + t810;
t651 = (qJD(3) - t779) * t738 + t808;
t669 = t671 ^ 2;
t670 = t673 ^ 2;
t695 = t696 ^ 2;
t697 = t700 ^ 2;
t698 = t702 ^ 2;
t736 = t738 ^ 2;
t865 = t779 ^ 2;
t864 = 2 * qJD(4);
t863 = pkin(4) * t782;
t862 = t696 * t784;
t861 = t696 * t788;
t859 = t779 * t782;
t858 = t779 * t783;
t857 = t779 * t785;
t856 = t779 * t789;
t780 = t786 ^ 2;
t855 = t780 * t793;
t854 = t782 * t589;
t637 = t643 + t778;
t853 = t782 * t637;
t852 = t783 * t589;
t851 = t783 * t637;
t639 = t700 * pkin(4) - t702 * pkin(8);
t723 = -t786 * g(3) + t790 * t740;
t795 = -pkin(2) * t775 + t748 * pkin(7) - qJD(2) * t801 + t723;
t622 = t785 * t866 + t789 * t795;
t571 = -t735 * pkin(3) + t676 * qJ(4) - t779 * t803 + t622;
t621 = t785 * t795 - t789 * t866;
t794 = pkin(3) * t867 + t655 * qJ(4) - t621;
t812 = t782 * t571 - t783 * t794;
t472 = -t778 * pkin(4) - t865 * pkin(8) + (t864 + t639) * t702 + t812;
t850 = t784 * t472;
t573 = t628 + t802;
t849 = t784 * t573;
t501 = t702 * t864 + t812;
t502 = -0.2e1 * qJD(4) * t700 + t783 * t571 + t782 * t794;
t448 = -t783 * t501 + t782 * t502;
t848 = t785 * t448;
t847 = t785 * t679;
t705 = -t708 + t778;
t846 = t785 * t705;
t564 = -t789 * t621 + t785 * t622;
t845 = t786 * t564;
t739 = t793 * pkin(6) + t800;
t844 = t786 * t739;
t765 = t790 * t840;
t755 = qJDD(2) + t765;
t842 = t786 * t755;
t756 = qJDD(2) - t765;
t841 = t786 * t756;
t839 = t787 * t778;
t838 = t788 * t472;
t837 = t788 * t573;
t836 = t789 * t448;
t835 = t789 * t679;
t834 = t789 * t705;
t833 = t790 * t564;
t832 = t790 * t739;
t831 = t790 * t756;
t473 = -pkin(4) * t865 + t778 * pkin(8) - t700 * t639 + t502;
t506 = pkin(4) * t595 - pkin(8) * t599 - t589;
t445 = t788 * t473 + t784 * t506;
t830 = t780 + t781;
t824 = t787 * qJDD(1);
t823 = t791 * qJDD(1);
t822 = t782 * t628;
t821 = t783 * t628;
t820 = t787 * t643;
t819 = t787 * t708;
t818 = t791 * t643;
t817 = t791 * t708;
t816 = -pkin(4) * t783 - pkin(3);
t444 = t784 * t473 - t788 * t506;
t449 = t782 * t501 + t783 * t502;
t565 = t785 * t621 + t789 * t622;
t722 = t790 * g(3) + t843;
t666 = t786 * t722 + t790 * t723;
t714 = -t787 * t757 - t791 * t758;
t807 = t787 * t765;
t806 = t791 * t765;
t752 = -t787 * t793 + t823;
t804 = -pkin(5) * t752 - t787 * g(3);
t408 = -t788 * t444 + t784 * t445;
t409 = t784 * t444 + t788 * t445;
t665 = t790 * t722 - t786 * t723;
t713 = t791 * t757 - t787 * t758;
t596 = t809 - t860;
t768 = t791 * t778;
t762 = t775 - t792;
t761 = -t792 - t855;
t760 = t792 - t855;
t754 = t775 - t855;
t753 = t775 + t855;
t751 = t791 * t793 + t824;
t750 = t830 * qJDD(1);
t749 = t773 - 0.2e1 * t814;
t746 = 0.2e1 * t813 + t825;
t744 = t790 * t755;
t743 = t830 * t826;
t734 = -pkin(5) * t751 + t791 * g(3);
t727 = -t736 + t865;
t726 = t735 - t865;
t725 = t790 * t747 - t780 * t826;
t724 = -t786 * t748 - t781 * t826;
t721 = -t736 - t865;
t720 = -t786 * t761 - t831;
t719 = -t786 * t760 + t744;
t718 = t790 * t763 - t842;
t717 = t790 * t762 - t841;
t716 = t790 * t761 - t841;
t715 = t786 * t763 + t744;
t711 = t791 * t750 - t787 * t753;
t710 = t787 * t750 + t791 * t753;
t709 = -t786 * t746 + t790 * t749;
t707 = -t736 + t735;
t703 = -t865 - t735;
t690 = t791 * t720 + t787 * t746;
t689 = t791 * t718 - t787 * t749;
t688 = t787 * t720 - t791 * t746;
t687 = t787 * t718 + t791 * t749;
t686 = -t698 + t865;
t685 = t697 - t865;
t684 = (t737 * t789 + t738 * t785) * t779;
t683 = (t737 * t785 - t738 * t789) * t779;
t682 = -pkin(6) * t716 - t832;
t681 = -pkin(6) * t715 - t844;
t680 = -t698 - t865;
t678 = -t735 - t736;
t675 = -pkin(1) * t716 + t723;
t674 = -pkin(1) * t715 + t722;
t661 = t789 * t726 - t846;
t660 = -t785 * t727 + t871;
t659 = t785 * t726 + t834;
t658 = t789 * t727 + t873;
t657 = -t785 * t721 - t834;
t656 = t789 * t721 - t846;
t650 = (qJD(3) + t779) * t738 + t808;
t649 = t789 * t677 - t738 * t857;
t648 = t785 * t677 + t738 * t856;
t647 = -t785 * t676 - t737 * t856;
t646 = t789 * t676 - t737 * t857;
t645 = t791 * t666 - t787 * t739;
t644 = t787 * t666 + t791 * t739;
t642 = t789 * t703 - t873;
t641 = t785 * t703 + t871;
t640 = -t698 + t697;
t635 = -t865 - t697;
t633 = -t670 + t695;
t632 = t669 - t695;
t631 = (-t700 * t783 + t702 * t782) * t779;
t630 = (-t700 * t782 - t702 * t783) * t779;
t629 = -t786 * t683 + t790 * t684;
t627 = -t670 + t669;
t620 = -t697 - t698;
t618 = -pkin(7) * t656 - t835;
t617 = -t670 - t695;
t616 = -pkin(7) * t641 - t847;
t615 = -t786 * t659 + t790 * t661;
t614 = -t786 * t658 + t790 * t660;
t613 = t783 * t685 - t853;
t612 = -t782 * t686 + t875;
t611 = t782 * t685 + t851;
t610 = t783 * t686 + t876;
t609 = -t782 * t680 - t851;
t608 = t783 * t680 - t853;
t607 = -t695 - t669;
t606 = -t786 * t656 + t790 * t657;
t605 = t790 * t656 + t786 * t657;
t604 = -t651 * t789 - t785 * t655;
t603 = -t789 * t650 - t785 * t868;
t602 = -t651 * t785 + t789 * t655;
t601 = -t785 * t650 + t789 * t868;
t600 = -t626 - t692;
t594 = t669 + t670;
t593 = t783 * t626 - t702 * t859;
t592 = t782 * t626 + t702 * t858;
t591 = t700 * t858 + t782 * t809;
t590 = t700 * t859 - t783 * t809;
t588 = -t786 * t648 + t790 * t649;
t587 = -t786 * t646 + t790 * t647;
t586 = -t786 * t641 + t790 * t642;
t585 = t790 * t641 + t786 * t642;
t583 = -t673 * qJD(5) - t810;
t582 = t783 * t635 - t876;
t581 = t782 * t635 + t875;
t580 = (-t671 * t788 + t673 * t784) * t696;
t579 = (t671 * t784 + t673 * t788) * t696;
t578 = -pkin(2) * t868 + pkin(7) * t657 - t847;
t577 = -t785 * t630 + t789 * t631;
t576 = t789 * t630 + t785 * t631;
t575 = -pkin(2) * t650 + pkin(7) * t642 + t835;
t569 = t791 * t606 + t787 * t868;
t568 = t787 * t606 - t791 * t868;
t563 = t791 * t586 + t787 * t650;
t562 = t787 * t586 - t791 * t650;
t560 = -t634 - t584;
t557 = (-qJD(5) - t696) * t673 - t810;
t555 = t788 * t584 - t673 * t862;
t554 = -t784 * t584 - t673 * t861;
t553 = -t784 * t583 + t671 * t861;
t552 = -t788 * t583 - t671 * t862;
t551 = -t785 * t611 + t789 * t613;
t550 = -t785 * t610 + t789 * t612;
t549 = t789 * t611 + t785 * t613;
t548 = t789 * t610 + t785 * t612;
t547 = -t785 * t608 + t789 * t609;
t546 = t789 * t608 + t785 * t609;
t545 = pkin(2) * t679 + pkin(7) * t565;
t544 = -t786 * t602 + t790 * t604;
t543 = -t786 * t601 + t790 * t603;
t542 = t790 * t602 + t786 * t604;
t541 = -qJ(4) * t608 - t852;
t540 = t783 * t580 + t782 * t802;
t539 = t782 * t580 - t783 * t802;
t538 = -t596 * t783 - t782 * t600;
t537 = -t783 * t595 - t782 * t599;
t536 = -t596 * t782 + t783 * t600;
t535 = -t782 * t595 + t783 * t599;
t534 = t788 * t632 - t849;
t533 = -t784 * t633 + t872;
t532 = -t784 * t632 - t837;
t531 = -t788 * t633 - t874;
t530 = -t785 * t592 + t789 * t593;
t529 = -t785 * t590 + t789 * t591;
t528 = t789 * t592 + t785 * t593;
t527 = t789 * t590 + t785 * t591;
t526 = -pkin(1) * t605 - pkin(2) * t656 + t622;
t525 = -qJ(4) * t581 - t854;
t524 = -t784 * t617 - t837;
t523 = t788 * t617 - t849;
t522 = t791 * t544 + t787 * t678;
t521 = t787 * t544 - t791 * t678;
t520 = -t785 * t581 + t789 * t582;
t519 = t789 * t581 + t785 * t582;
t518 = t788 * t607 - t874;
t517 = t784 * t607 + t872;
t516 = -pkin(1) * t585 + t785 * t723 + t789 * t722 + (t785 * (t748 + t814) - t789 * (-t747 + t813)) * pkin(7) + (-t789 * t755 + t763 * t785 - t641) * pkin(2);
t515 = t783 * t555 + t822;
t514 = t783 * t553 - t822;
t513 = t782 * t555 - t821;
t512 = t782 * t553 + t821;
t511 = -pkin(7) * t602 - t564;
t510 = -t786 * t576 + t790 * t577;
t509 = -pkin(2) * t678 + pkin(7) * t604 + t565;
t508 = -pkin(1) * t542 - pkin(2) * t602;
t507 = -pkin(3) * t599 + qJ(4) * t609 - t854;
t503 = -pkin(6) * t605 - t786 * t578 + t790 * t618;
t500 = -pkin(3) * t595 + qJ(4) * t582 + t852;
t499 = t790 * t565 - t845;
t498 = t786 * t565 + t833;
t496 = -pkin(6) * t585 - t786 * t575 + t790 * t616;
t495 = -t556 * t788 - t784 * t560;
t494 = t788 * t557 - t784 * t559;
t493 = -t556 * t784 + t788 * t560;
t492 = -t784 * t557 - t788 * t559;
t491 = t791 * t499 - t787 * t679;
t490 = t787 * t499 + t791 * t679;
t489 = -t786 * t549 + t790 * t551;
t488 = -t786 * t548 + t790 * t550;
t487 = t783 * t534 - t782 * t556;
t486 = t783 * t533 - t782 * t560;
t485 = t782 * t534 + t783 * t556;
t484 = t782 * t533 + t783 * t560;
t483 = -t786 * t546 + t790 * t547;
t482 = t790 * t546 + t786 * t547;
t481 = -t785 * t539 + t789 * t540;
t480 = t789 * t539 + t785 * t540;
t479 = t783 * t524 + t559 * t782;
t478 = t782 * t524 - t559 * t783;
t477 = -t785 * t536 + t789 * t538;
t476 = -t785 * t535 + t789 * t537;
t475 = t789 * t536 + t785 * t538;
t474 = t789 * t535 + t785 * t537;
t470 = t783 * t518 - t782 * t557;
t469 = t782 * t518 + t783 * t557;
t468 = t783 * t494 - t782 * t627;
t467 = t782 * t494 + t783 * t627;
t466 = -t786 * t528 + t790 * t530;
t465 = -t786 * t527 + t790 * t529;
t464 = t783 * t495 - t782 * t594;
t463 = t782 * t495 + t783 * t594;
t462 = -t786 * t519 + t790 * t520;
t461 = t790 * t519 + t786 * t520;
t460 = t791 * t483 + t599 * t787;
t459 = t787 * t483 - t599 * t791;
t458 = -pkin(1) * t498 - pkin(2) * t564;
t457 = -t785 * t513 + t789 * t515;
t456 = -t785 * t512 + t789 * t514;
t455 = t789 * t513 + t785 * t515;
t454 = t789 * t512 + t785 * t514;
t453 = t791 * t462 + t787 * t595;
t452 = t787 * t462 - t791 * t595;
t451 = -pkin(8) * t523 + t838;
t450 = -pkin(8) * t517 + t850;
t447 = -pkin(6) * t498 - pkin(7) * t833 - t786 * t545;
t446 = -pkin(7) * t546 - t785 * t507 + t789 * t541;
t443 = -pkin(6) * t542 - t786 * t509 + t790 * t511;
t442 = pkin(3) * t589 + qJ(4) * t449;
t441 = -t785 * t485 + t789 * t487;
t440 = -t785 * t484 + t789 * t486;
t439 = t789 * t485 + t785 * t487;
t438 = t789 * t484 + t785 * t486;
t437 = -pkin(7) * t519 - t785 * t500 + t789 * t525;
t436 = -pkin(2) * t599 + pkin(7) * t547 + t789 * t507 + t785 * t541;
t435 = -t786 * t480 + t790 * t481;
t434 = -t785 * t478 + t789 * t479;
t433 = t789 * t478 + t785 * t479;
t432 = -t786 * t475 + t790 * t477;
t431 = -t786 * t474 + t790 * t476;
t430 = t790 * t475 + t786 * t477;
t429 = -t785 * t469 + t789 * t470;
t428 = t789 * t469 + t785 * t470;
t427 = -t785 * t467 + t789 * t468;
t426 = t789 * t467 + t785 * t468;
t425 = -qJ(4) * t536 - t448;
t424 = -pkin(2) * t595 + pkin(7) * t520 + t789 * t500 + t785 * t525;
t423 = t791 * t432 + t787 * t620;
t422 = t787 * t432 - t791 * t620;
t421 = -t785 * t463 + t789 * t464;
t420 = t789 * t463 + t785 * t464;
t419 = -pkin(3) * t620 + qJ(4) * t538 + t449;
t418 = -pkin(4) * t523 + t445;
t417 = -pkin(4) * t517 + t444;
t416 = -pkin(1) * t482 - pkin(2) * t546 - pkin(3) * t608 + t502;
t415 = -t786 * t455 + t790 * t457;
t414 = -t786 * t454 + t790 * t456;
t413 = -pkin(1) * t461 - pkin(2) * t519 - pkin(3) * t581 + t501;
t412 = t789 * t449 - t848;
t411 = t785 * t449 + t836;
t410 = -pkin(1) * t430 - pkin(2) * t475 - pkin(3) * t536;
t407 = -t786 * t439 + t790 * t441;
t406 = -t786 * t438 + t790 * t440;
t405 = -t786 * t433 + t790 * t434;
t404 = t790 * t433 + t786 * t434;
t403 = -t786 * t428 + t790 * t429;
t402 = t790 * t428 + t786 * t429;
t401 = -t786 * t426 + t790 * t427;
t400 = -pkin(8) * t493 - t408;
t399 = -t786 * t420 + t790 * t421;
t398 = t790 * t420 + t786 * t421;
t397 = t783 * t409 + t782 * t472;
t396 = t782 * t409 - t783 * t472;
t395 = -pkin(6) * t482 - t786 * t436 + t790 * t446;
t394 = t791 * t405 + t787 * t523;
t393 = t787 * t405 - t791 * t523;
t392 = t791 * t403 + t787 * t517;
t391 = t787 * t403 - t791 * t517;
t390 = -qJ(4) * t478 - t782 * t418 + t783 * t451;
t389 = -qJ(4) * t469 - t782 * t417 + t783 * t450;
t388 = -pkin(6) * t461 - t786 * t424 + t790 * t437;
t387 = t791 * t399 + t787 * t493;
t386 = t787 * t399 - t791 * t493;
t385 = -pkin(7) * t475 - t785 * t419 + t789 * t425;
t384 = -pkin(3) * t523 + qJ(4) * t479 + t783 * t418 + t782 * t451;
t383 = -pkin(3) * t517 + qJ(4) * t470 + t783 * t417 + t782 * t450;
t382 = -pkin(2) * t620 + pkin(7) * t477 + t789 * t419 + t785 * t425;
t381 = -qJ(4) * t463 + t783 * t400 + t493 * t863;
t380 = -t786 * t411 + t790 * t412;
t379 = t790 * t411 + t786 * t412;
t378 = -pkin(7) * t411 - qJ(4) * t836 - t785 * t442;
t377 = t791 * t380 - t787 * t589;
t376 = t787 * t380 + t791 * t589;
t375 = pkin(2) * t589 + pkin(7) * t412 - qJ(4) * t848 + t789 * t442;
t374 = qJ(4) * t464 + t782 * t400 + t816 * t493;
t373 = -pkin(1) * t404 - pkin(2) * t433 - pkin(3) * t478 + pkin(4) * t559 - pkin(8) * t524 - t850;
t372 = -t785 * t396 + t789 * t397;
t371 = t789 * t396 + t785 * t397;
t370 = -pkin(1) * t402 - pkin(2) * t428 - pkin(3) * t469 - pkin(4) * t557 - pkin(8) * t518 + t838;
t369 = -pkin(1) * t379 - pkin(2) * t411 - pkin(3) * t448;
t368 = -qJ(4) * t396 + (-pkin(8) * t783 + t863) * t408;
t367 = -pkin(7) * t433 - t785 * t384 + t789 * t390;
t366 = -pkin(1) * t398 - pkin(2) * t420 - pkin(3) * t463 - pkin(4) * t594 - pkin(8) * t495 - t409;
t365 = -pkin(7) * t428 - t785 * t383 + t789 * t389;
t364 = -pkin(2) * t523 + pkin(7) * t434 + t789 * t384 + t785 * t390;
t363 = -pkin(6) * t430 - t786 * t382 + t790 * t385;
t362 = -pkin(2) * t517 + pkin(7) * t429 + t789 * t383 + t785 * t389;
t361 = qJ(4) * t397 + (-pkin(8) * t782 + t816) * t408;
t360 = -pkin(7) * t420 - t785 * t374 + t789 * t381;
t359 = -pkin(2) * t493 + pkin(7) * t421 + t789 * t374 + t785 * t381;
t358 = -t786 * t371 + t790 * t372;
t357 = t790 * t371 + t786 * t372;
t356 = -pkin(6) * t379 - t786 * t375 + t790 * t378;
t355 = t791 * t358 + t787 * t408;
t354 = t787 * t358 - t791 * t408;
t353 = -pkin(6) * t404 - t786 * t364 + t790 * t367;
t352 = -pkin(6) * t402 - t786 * t362 + t790 * t365;
t351 = -pkin(7) * t371 - t785 * t361 + t789 * t368;
t350 = -pkin(1) * t357 - pkin(2) * t371 - pkin(3) * t396 + pkin(4) * t472 - pkin(8) * t409;
t349 = -pkin(6) * t398 - t786 * t359 + t790 * t360;
t348 = -pkin(2) * t408 + pkin(7) * t372 + t789 * t361 + t785 * t368;
t347 = -pkin(6) * t357 - t786 * t348 + t790 * t351;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t751, -t752, 0, t714, 0, 0, 0, 0, 0, 0, t689, t690, t711, t645, 0, 0, 0, 0, 0, 0, t563, t569, t522, t491, 0, 0, 0, 0, 0, 0, t453, t460, t423, t377, 0, 0, 0, 0, 0, 0, t392, t394, t387, t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t752, -t751, 0, t713, 0, 0, 0, 0, 0, 0, t687, t688, t710, t644, 0, 0, 0, 0, 0, 0, t562, t568, t521, t490, 0, 0, 0, 0, 0, 0, t452, t459, t422, t376, 0, 0, 0, 0, 0, 0, t391, t393, t386, t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t715, t716, 0, -t665, 0, 0, 0, 0, 0, 0, t585, t605, t542, t498, 0, 0, 0, 0, 0, 0, t461, t482, t430, t379, 0, 0, 0, 0, 0, 0, t402, t404, t398, t357; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t752, 0, -t751, 0, t804, -t734, -t713, -pkin(5) * t713, t791 * t725 - t807, t791 * t709 - t787 * t754, t791 * t719 + t786 * t824, t791 * t724 + t807, t791 * t717 + t773 * t787, t787 * qJDD(2) + t791 * t743, -pkin(5) * t687 - t787 * t674 + t791 * t681, -pkin(5) * t688 - t787 * t675 + t791 * t682, -pkin(5) * t710 + t791 * t665, -pkin(5) * t644 - (pkin(1) * t787 - pkin(6) * t791) * t665, t791 * t588 - t819, t791 * t543 - t787 * t707, t791 * t614 - t787 * t655, t791 * t587 + t819, t791 * t615 - t787 * t651, t791 * t629 + t839, -pkin(5) * t562 + t791 * t496 - t787 * t516, -pkin(5) * t568 + t791 * t503 - t787 * t526, -pkin(5) * t521 + t791 * t443 - t787 * t508, -pkin(5) * t490 + t791 * t447 - t787 * t458, t791 * t466 + t820, t791 * t431 - t787 * t640, t791 * t488 - t787 * t600, t791 * t465 - t820, t791 * t489 - t787 * t596, t791 * t510 + t839, -pkin(5) * t452 + t791 * t388 - t787 * t413, -pkin(5) * t459 + t791 * t395 - t787 * t416, -pkin(5) * t422 + t791 * t363 - t787 * t410, -pkin(5) * t376 + t791 * t356 - t787 * t369, t791 * t415 - t787 * t554, t791 * t401 - t787 * t492, t791 * t406 - t787 * t531, t791 * t414 - t787 * t552, t791 * t407 - t787 * t532, t791 * t435 - t787 * t579, -pkin(5) * t391 + t791 * t352 - t787 * t370, -pkin(5) * t393 + t791 * t353 - t787 * t373, -pkin(5) * t386 + t791 * t349 - t787 * t366, -pkin(5) * t354 + t791 * t347 - t787 * t350; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t751, 0, t752, 0, t734, t804, t714, pkin(5) * t714, t787 * t725 + t806, t787 * t709 + t791 * t754, t787 * t719 - t786 * t823, t787 * t724 - t806, t787 * t717 - t790 * t823, -t791 * qJDD(2) + t787 * t743, pkin(5) * t689 + t791 * t674 + t787 * t681, pkin(5) * t690 + t791 * t675 + t787 * t682, pkin(5) * t711 + t787 * t665, pkin(5) * t645 - (-pkin(1) * t791 - pkin(6) * t787) * t665, t787 * t588 + t817, t787 * t543 + t791 * t707, t787 * t614 + t791 * t655, t787 * t587 - t817, t787 * t615 + t791 * t651, t787 * t629 - t768, pkin(5) * t563 + t787 * t496 + t791 * t516, pkin(5) * t569 + t787 * t503 + t791 * t526, pkin(5) * t522 + t787 * t443 + t791 * t508, pkin(5) * t491 + t787 * t447 + t791 * t458, t787 * t466 - t818, t787 * t431 + t791 * t640, t787 * t488 + t791 * t600, t787 * t465 + t818, t787 * t489 + t791 * t596, t787 * t510 - t768, pkin(5) * t453 + t787 * t388 + t791 * t413, pkin(5) * t460 + t787 * t395 + t791 * t416, pkin(5) * t423 + t787 * t363 + t791 * t410, pkin(5) * t377 + t787 * t356 + t791 * t369, t787 * t415 + t791 * t554, t787 * t401 + t791 * t492, t787 * t406 + t791 * t531, t787 * t414 + t791 * t552, t787 * t407 + t791 * t532, t787 * t435 + t791 * t579, pkin(5) * t392 + t787 * t352 + t791 * t370, pkin(5) * t394 + t787 * t353 + t791 * t373, pkin(5) * t387 + t787 * t349 + t791 * t366, pkin(5) * t355 + t787 * t347 + t791 * t350; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t757, t758, 0, 0, (t747 + t813) * t786, t790 * t746 + t786 * t749, t790 * t760 + t842, (t748 - t814) * t790, t786 * t762 + t831, 0, pkin(1) * t749 + pkin(6) * t718 + t832, -pkin(1) * t746 + pkin(6) * t720 - t844, pkin(1) * t753 + pkin(6) * t750 + t666, pkin(1) * t739 + pkin(6) * t666, t790 * t648 + t786 * t649, t790 * t601 + t786 * t603, t790 * t658 + t786 * t660, t790 * t646 + t786 * t647, t790 * t659 + t786 * t661, t790 * t683 + t786 * t684, -pkin(1) * t650 + pkin(6) * t586 + t790 * t575 + t786 * t616, -pkin(1) * t868 + pkin(6) * t606 + t790 * t578 + t786 * t618, -pkin(1) * t678 + pkin(6) * t544 + t790 * t509 + t786 * t511, pkin(1) * t679 + pkin(6) * t499 - pkin(7) * t845 + t790 * t545, t790 * t528 + t786 * t530, t790 * t474 + t786 * t476, t790 * t548 + t786 * t550, t790 * t527 + t786 * t529, t790 * t549 + t786 * t551, t790 * t576 + t786 * t577, -pkin(1) * t595 + pkin(6) * t462 + t790 * t424 + t786 * t437, -pkin(1) * t599 + pkin(6) * t483 + t790 * t436 + t786 * t446, -pkin(1) * t620 + pkin(6) * t432 + t790 * t382 + t786 * t385, pkin(1) * t589 + pkin(6) * t380 + t790 * t375 + t786 * t378, t790 * t455 + t786 * t457, t790 * t426 + t786 * t427, t790 * t438 + t786 * t440, t790 * t454 + t786 * t456, t790 * t439 + t786 * t441, t790 * t480 + t786 * t481, -pkin(1) * t517 + pkin(6) * t403 + t790 * t362 + t786 * t365, -pkin(1) * t523 + pkin(6) * t405 + t790 * t364 + t786 * t367, -pkin(1) * t493 + pkin(6) * t399 + t790 * t359 + t786 * t360, -pkin(1) * t408 + pkin(6) * t358 + t790 * t348 + t786 * t351;];
tauB_reg = t1;
