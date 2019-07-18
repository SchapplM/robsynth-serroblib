% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRRPR7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:02:07
% EndTime: 2019-05-05 23:02:26
% DurationCPUTime: 17.63s
% Computational Cost: add. (120088->730), mult. (261495->1085), div. (0->0), fcn. (187932->10), ass. (0->487)
t830 = sin(pkin(10));
t833 = sin(qJ(4));
t837 = cos(qJ(4));
t838 = cos(qJ(3));
t883 = qJD(1) * t838;
t834 = sin(qJ(3));
t884 = qJD(1) * t834;
t784 = -t833 * t883 - t837 * t884;
t785 = -t833 * t884 + t837 * t883;
t831 = cos(pkin(10));
t741 = -t831 * t784 + t785 * t830;
t743 = t830 * t784 + t831 * t785;
t686 = t743 * t741;
t825 = qJDD(3) + qJDD(4);
t934 = -t686 + t825;
t941 = t830 * t934;
t940 = t831 * t934;
t832 = sin(qJ(6));
t826 = qJD(3) + qJD(4);
t836 = cos(qJ(6));
t717 = t743 * t832 - t836 * t826;
t719 = t743 * t836 + t826 * t832;
t665 = t719 * t717;
t880 = qJD(1) * qJD(3);
t864 = t838 * t880;
t878 = qJDD(1) * t834;
t793 = -t864 - t878;
t865 = t834 * t880;
t876 = qJDD(1) * t838;
t794 = -t865 + t876;
t860 = -t837 * t793 + t833 * t794;
t720 = -qJD(4) * t785 - t860;
t721 = t784 * qJD(4) + t833 * t793 + t837 * t794;
t861 = -t831 * t720 + t721 * t830;
t851 = qJDD(6) + t861;
t935 = -t665 + t851;
t939 = t832 * t935;
t749 = t784 * t785;
t933 = t749 + t825;
t938 = t833 * t933;
t937 = t836 * t935;
t936 = t837 * t933;
t907 = t743 * t826;
t631 = t861 + t907;
t663 = t830 * t720 + t831 * t721;
t908 = t741 * t826;
t853 = t663 - t908;
t619 = -t717 * qJD(6) + t836 * t663 + t832 * t825;
t737 = qJD(6) + t741;
t675 = t737 * t717;
t596 = -t675 + t619;
t903 = t784 * t826;
t696 = -t721 + t903;
t841 = qJD(1) ^ 2;
t887 = t838 * t841;
t835 = sin(qJ(1));
t839 = cos(qJ(1));
t803 = t835 * g(1) - t839 * g(2);
t854 = qJDD(2) - t803;
t845 = -t841 * qJ(2) + t854;
t928 = pkin(7) + pkin(1);
t771 = -qJDD(1) * t928 + t845;
t888 = t838 * t771;
t708 = qJDD(3) * pkin(3) - t794 * pkin(8) + t888 + (-pkin(3) * t887 - pkin(8) * t880 + g(3)) * t834;
t752 = -t838 * g(3) + t834 * t771;
t847 = qJD(3) * pkin(3) - pkin(8) * t883;
t828 = t834 ^ 2;
t892 = t828 * t841;
t709 = -pkin(3) * t892 + t793 * pkin(8) - qJD(3) * t847 + t752;
t658 = t833 * t708 + t837 * t709;
t782 = t784 ^ 2;
t852 = pkin(4) * t826 - qJ(5) * t785;
t604 = -t782 * pkin(4) + t720 * qJ(5) - t826 * t852 + t658;
t850 = t837 * t708 - t833 * t709;
t842 = pkin(4) * t933 + qJ(5) * t696 + t850;
t539 = -0.2e1 * qJD(5) * t741 + t831 * t604 + t830 * t842;
t862 = t832 * t663 - t836 * t825;
t593 = (qJD(6) - t737) * t719 + t862;
t804 = t839 * g(1) + t835 * g(2);
t827 = qJDD(1) * qJ(2);
t849 = t804 - t827;
t932 = -t793 * pkin(3) - (pkin(8) * t828 + t928) * t841 + t847 * t883 - t849;
t931 = t720 * pkin(4) + t782 * qJ(5) - t785 * t852 - qJDD(5);
t715 = t717 ^ 2;
t716 = t719 ^ 2;
t736 = t737 ^ 2;
t738 = t741 ^ 2;
t739 = t743 ^ 2;
t783 = t785 ^ 2;
t930 = t826 ^ 2;
t929 = 0.2e1 * qJD(5);
t927 = pkin(5) * t830;
t925 = qJDD(1) * pkin(1);
t863 = t830 * t604 - t831 * t842;
t538 = t743 * t929 + t863;
t484 = -t538 * t831 + t539 * t830;
t924 = t484 * t833;
t923 = t484 * t837;
t680 = pkin(5) * t741 - pkin(9) * t743;
t506 = -t825 * pkin(5) - t930 * pkin(9) + (t929 + t680) * t743 + t863;
t922 = t506 * t832;
t921 = t506 * t836;
t591 = t658 * t833 + t837 * t850;
t920 = t591 * t834;
t919 = t591 * t838;
t609 = t665 + t851;
t918 = t609 * t832;
t917 = t609 * t836;
t879 = qJD(2) * qJD(1);
t874 = -0.2e1 * t879;
t710 = t874 - t932;
t622 = t710 + t931;
t916 = t622 * t830;
t915 = t622 * t831;
t678 = t686 + t825;
t914 = t678 * t830;
t913 = t678 * t831;
t912 = t710 * t833;
t911 = t710 * t837;
t910 = t737 * t832;
t909 = t737 * t836;
t746 = -t749 + t825;
t906 = t746 * t833;
t905 = t746 * t837;
t829 = t838 ^ 2;
t885 = t828 + t829;
t796 = t885 * qJDD(1);
t902 = t796 * t835;
t901 = t796 * t839;
t867 = t834 * t887;
t801 = qJDD(3) + t867;
t900 = t801 * t834;
t899 = t801 * t838;
t802 = qJDD(3) - t867;
t898 = t802 * t834;
t897 = t802 * t838;
t896 = t826 * t830;
t895 = t826 * t831;
t894 = t826 * t833;
t893 = t826 * t837;
t891 = t829 * t841;
t766 = t841 * t928 + t849 + t874;
t890 = t834 * t766;
t889 = t838 * t766;
t507 = -pkin(5) * t930 + pkin(9) * t825 - t680 * t741 + t539;
t822 = 0.2e1 * t879;
t544 = pkin(5) * t631 - t853 * pkin(9) + t822 - t931 + t932;
t482 = t836 * t507 + t832 * t544;
t877 = qJDD(1) * t835;
t875 = qJDD(1) * t839;
t873 = t830 * t665;
t872 = t831 * t665;
t871 = t835 * t686;
t870 = t839 * t686;
t869 = t835 * t749;
t868 = t839 * t749;
t866 = -pkin(5) * t831 - pkin(4);
t481 = t507 * t832 - t836 * t544;
t445 = t481 * t832 + t836 * t482;
t485 = t538 * t830 + t831 * t539;
t592 = t837 * t658 - t833 * t850;
t772 = -t841 * pkin(1) + t822 - t849;
t777 = -t845 + t925;
t731 = t839 * t772 - t777 * t835;
t756 = -t803 * t835 - t839 * t804;
t858 = t835 * t867;
t857 = t839 * t867;
t797 = -t835 * t841 + t875;
t856 = pkin(6) * t797 + g(3) * t835;
t798 = t839 * t841 + t877;
t855 = -pkin(6) * t798 + g(3) * t839;
t751 = t834 * g(3) + t888;
t444 = -t481 * t836 + t482 * t832;
t700 = t838 * t751 + t834 * t752;
t701 = -t751 * t834 + t752 * t838;
t728 = t772 * t835 + t777 * t839;
t755 = t803 * t839 - t804 * t835;
t848 = t721 + t903;
t846 = -t861 + t907;
t844 = (-qJD(4) + t826) * t785 - t860;
t840 = qJD(3) ^ 2;
t816 = t839 * t825;
t814 = t835 * t825;
t809 = -t840 - t891;
t808 = t840 - t891;
t807 = -t840 - t892;
t806 = -t840 + t892;
t800 = (-t828 + t829) * t841;
t799 = t885 * t841;
t795 = -0.2e1 * t865 + t876;
t792 = 0.2e1 * t864 + t878;
t790 = t885 * t880;
t770 = -t783 + t930;
t769 = t782 - t930;
t768 = -t794 * t834 - t829 * t880;
t767 = -t793 * t838 - t828 * t880;
t763 = -t783 - t930;
t762 = -t809 * t834 - t899;
t761 = t807 * t838 - t898;
t760 = t809 * t838 - t900;
t759 = -t808 * t838 - t898;
t758 = t807 * t834 + t897;
t757 = -t806 * t834 - t899;
t754 = -t799 * t839 - t902;
t753 = -t799 * t835 + t901;
t750 = t792 * t834 - t795 * t838;
t748 = t783 - t782;
t744 = -t930 - t782;
t733 = t760 * t835 + t795 * t839;
t732 = t758 * t835 + t792 * t839;
t730 = -t760 * t839 + t795 * t835;
t729 = -t758 * t839 + t792 * t835;
t727 = -t739 + t930;
t726 = t738 - t930;
t725 = (t784 * t837 + t785 * t833) * t826;
t724 = (t784 * t833 - t785 * t837) * t826;
t723 = -t739 - t930;
t722 = -t782 - t783;
t705 = t769 * t837 - t906;
t704 = -t770 * t833 + t936;
t703 = t769 * t833 + t905;
t702 = t770 * t837 + t938;
t699 = -t763 * t833 - t905;
t698 = t763 * t837 - t906;
t692 = (qJD(4) + t826) * t785 + t860;
t691 = t721 * t837 - t785 * t894;
t690 = t721 * t833 + t785 * t893;
t689 = -t720 * t833 - t784 * t893;
t688 = t720 * t837 - t784 * t894;
t687 = -pkin(2) * t799 - t701;
t685 = t744 * t837 - t938;
t684 = t744 * t833 + t936;
t683 = t739 - t738;
t682 = pkin(2) * t760 - qJ(2) * t762 - t752;
t681 = pkin(2) * t758 - qJ(2) * t761 + t751;
t676 = -t930 - t738;
t674 = pkin(2) * t792 - t761 * t928 - t889;
t673 = pkin(2) * t795 - t762 * t928 + t890;
t672 = -t716 + t736;
t671 = t715 - t736;
t670 = (-t741 * t831 + t743 * t830) * t826;
t669 = (-t741 * t830 - t743 * t831) * t826;
t668 = t700 * t835 - t766 * t839;
t667 = -t700 * t839 - t766 * t835;
t666 = -t724 * t838 - t725 * t834;
t664 = -t716 + t715;
t659 = -t738 - t739;
t656 = -pkin(8) * t698 - t911;
t654 = -t716 - t736;
t653 = pkin(2) * t700 - qJ(2) * t701;
t652 = t726 * t831 - t914;
t651 = -t727 * t830 + t940;
t650 = -t703 * t838 - t705 * t834;
t649 = -t702 * t838 - t704 * t834;
t648 = t726 * t830 + t913;
t647 = t727 * t831 + t941;
t646 = -t723 * t830 - t913;
t645 = t723 * t831 - t914;
t644 = -pkin(8) * t684 - t912;
t643 = -t736 - t715;
t642 = -t698 * t834 + t699 * t838;
t641 = t698 * t838 + t699 * t834;
t640 = -t696 * t833 + t837 * t844;
t639 = -t692 * t837 - t833 * t848;
t638 = t696 * t837 + t833 * t844;
t637 = -t692 * t833 + t837 * t848;
t635 = -t663 - t908;
t630 = t715 + t716;
t629 = t663 * t831 - t743 * t896;
t628 = t663 * t830 + t743 * t895;
t627 = t741 * t895 + t830 * t861;
t626 = t741 * t896 - t831 * t861;
t625 = -pkin(2) * t766 - t701 * t928;
t624 = -t690 * t838 - t691 * t834;
t623 = -t688 * t838 - t689 * t834;
t621 = -t684 * t834 + t685 * t838;
t620 = t684 * t838 + t685 * t834;
t618 = -qJD(6) * t719 - t862;
t617 = t676 * t831 - t941;
t616 = t676 * t830 + t940;
t615 = (-t717 * t836 + t719 * t832) * t737;
t614 = (-t717 * t832 - t719 * t836) * t737;
t613 = -t669 * t833 + t670 * t837;
t612 = t669 * t837 + t670 * t833;
t611 = -pkin(3) * t848 + pkin(8) * t699 - t912;
t607 = -pkin(3) * t692 + pkin(8) * t685 + t911;
t606 = t641 * t835 + t839 * t848;
t605 = -t641 * t839 + t835 * t848;
t600 = t620 * t835 + t692 * t839;
t599 = -t620 * t839 + t692 * t835;
t597 = -t675 - t619;
t594 = (-qJD(6) - t737) * t719 - t862;
t590 = t619 * t836 - t719 * t910;
t589 = t619 * t832 + t719 * t909;
t588 = -t618 * t832 + t717 * t909;
t587 = t618 * t836 + t717 * t910;
t586 = -t648 * t833 + t652 * t837;
t585 = -t647 * t833 + t651 * t837;
t584 = t648 * t837 + t652 * t833;
t583 = t647 * t837 + t651 * t833;
t582 = -t645 * t833 + t646 * t837;
t581 = t645 * t837 + t646 * t833;
t580 = -t638 * t834 + t640 * t838;
t579 = t638 * t838 + t640 * t834;
t578 = -t637 * t838 - t639 * t834;
t577 = t615 * t831 + t830 * t851;
t576 = t615 * t830 - t831 * t851;
t575 = -t635 * t830 + t831 * t846;
t574 = -t631 * t831 - t830 * t853;
t573 = t635 * t831 + t830 * t846;
t572 = -t631 * t830 + t831 * t853;
t571 = t671 * t836 - t918;
t570 = -t672 * t832 + t937;
t569 = t671 * t832 + t917;
t568 = t672 * t836 + t939;
t567 = -qJ(5) * t645 - t915;
t566 = -t628 * t833 + t629 * t837;
t565 = -t626 * t833 + t627 * t837;
t564 = t628 * t837 + t629 * t833;
t563 = t626 * t837 + t627 * t833;
t562 = pkin(3) * t710 + pkin(8) * t592;
t561 = -qJ(5) * t616 - t916;
t560 = -t654 * t832 - t917;
t559 = t654 * t836 - t918;
t558 = t579 * t835 + t722 * t839;
t557 = -t579 * t839 + t722 * t835;
t556 = -t616 * t833 + t617 * t837;
t555 = t616 * t837 + t617 * t833;
t554 = t643 * t836 - t939;
t553 = t643 * t832 + t937;
t552 = t590 * t831 + t873;
t551 = t588 * t831 - t873;
t550 = t590 * t830 - t872;
t549 = t588 * t830 + t872;
t548 = -t612 * t838 - t613 * t834;
t547 = -pkin(8) * t638 - t591;
t546 = -pkin(3) * t722 + pkin(8) * t640 + t592;
t545 = -pkin(4) * t853 + qJ(5) * t646 - t916;
t541 = pkin(2) * t641 + pkin(3) * t698 - qJ(2) * t642 - t658;
t540 = -pkin(4) * t631 + qJ(5) * t617 + t915;
t536 = -t593 * t836 - t597 * t832;
t535 = t594 * t836 - t596 * t832;
t534 = -t593 * t832 + t597 * t836;
t533 = t594 * t832 + t596 * t836;
t532 = t592 * t838 - t920;
t531 = t592 * t834 + t919;
t530 = pkin(2) * t620 + pkin(3) * t684 - qJ(2) * t621 + t850;
t529 = -t584 * t838 - t586 * t834;
t528 = -t583 * t838 - t585 * t834;
t527 = t571 * t831 - t593 * t830;
t526 = t570 * t831 - t597 * t830;
t525 = t571 * t830 + t593 * t831;
t524 = t570 * t830 + t597 * t831;
t523 = -t581 * t834 + t582 * t838;
t522 = t581 * t838 + t582 * t834;
t521 = t531 * t835 - t710 * t839;
t520 = -t531 * t839 - t710 * t835;
t519 = -t576 * t833 + t577 * t837;
t518 = t576 * t837 + t577 * t833;
t517 = t560 * t831 + t596 * t830;
t516 = t560 * t830 - t596 * t831;
t515 = -t573 * t833 + t575 * t837;
t514 = -t572 * t833 + t574 * t837;
t513 = t573 * t837 + t575 * t833;
t512 = t572 * t837 + t574 * t833;
t511 = t554 * t831 - t594 * t830;
t510 = t554 * t830 + t594 * t831;
t509 = t535 * t831 - t664 * t830;
t508 = t535 * t830 + t664 * t831;
t504 = -t564 * t838 - t566 * t834;
t503 = -t563 * t838 - t565 * t834;
t502 = t536 * t831 - t630 * t830;
t501 = t536 * t830 + t630 * t831;
t500 = pkin(2) * t848 - t838 * t611 - t642 * t928 - t834 * t656;
t499 = -t555 * t834 + t556 * t838;
t498 = t555 * t838 + t556 * t834;
t497 = t522 * t835 + t839 * t853;
t496 = -t522 * t839 + t835 * t853;
t495 = pkin(2) * t579 + pkin(3) * t638 - qJ(2) * t580;
t494 = -t550 * t833 + t552 * t837;
t493 = -t549 * t833 + t551 * t837;
t492 = t550 * t837 + t552 * t833;
t491 = t549 * t837 + t551 * t833;
t490 = pkin(2) * t692 - t838 * t607 - t621 * t928 - t834 * t644;
t489 = t498 * t835 + t631 * t839;
t488 = -t498 * t839 + t631 * t835;
t487 = -pkin(9) * t559 + t921;
t486 = -pkin(9) * t553 + t922;
t483 = -pkin(8) * t581 - t545 * t833 + t567 * t837;
t479 = -t525 * t833 + t527 * t837;
t478 = -t524 * t833 + t526 * t837;
t477 = t525 * t837 + t527 * t833;
t476 = t524 * t837 + t526 * t833;
t475 = pkin(4) * t622 + qJ(5) * t485;
t474 = -pkin(8) * t555 - t540 * t833 + t561 * t837;
t473 = -pkin(3) * t853 + pkin(8) * t582 + t545 * t837 + t567 * t833;
t472 = -t518 * t838 - t519 * t834;
t471 = -t516 * t833 + t517 * t837;
t470 = t516 * t837 + t517 * t833;
t469 = -t513 * t834 + t515 * t838;
t468 = t513 * t838 + t515 * t834;
t467 = -t512 * t838 - t514 * t834;
t466 = -t510 * t833 + t511 * t837;
t465 = t510 * t837 + t511 * t833;
t464 = -t508 * t833 + t509 * t837;
t463 = t508 * t837 + t509 * t833;
t462 = pkin(2) * t531 + pkin(3) * t591 - qJ(2) * t532;
t461 = -qJ(5) * t573 - t484;
t460 = -pkin(3) * t631 + pkin(8) * t556 + t540 * t837 + t561 * t833;
t459 = t468 * t835 + t659 * t839;
t458 = -t468 * t839 + t659 * t835;
t457 = -t501 * t833 + t502 * t837;
t456 = t501 * t837 + t502 * t833;
t455 = -pkin(4) * t659 + qJ(5) * t575 + t485;
t454 = -pkin(5) * t559 + t482;
t453 = -pkin(5) * t553 + t481;
t452 = pkin(2) * t722 - t838 * t546 - t834 * t547 - t580 * t928;
t451 = -t492 * t838 - t494 * t834;
t450 = -t491 * t838 - t493 * t834;
t449 = -pkin(2) * t710 + pkin(8) * t920 - t532 * t928 - t838 * t562;
t448 = t485 * t837 - t924;
t447 = t485 * t833 + t923;
t446 = pkin(2) * t522 + pkin(3) * t581 + pkin(4) * t645 - qJ(2) * t523 - t539;
t443 = -t477 * t838 - t479 * t834;
t442 = -t476 * t838 - t478 * t834;
t441 = -t470 * t834 + t471 * t838;
t440 = t470 * t838 + t471 * t834;
t439 = pkin(2) * t498 + pkin(3) * t555 + pkin(4) * t616 - qJ(2) * t499 - t538;
t438 = -t465 * t834 + t466 * t838;
t437 = t465 * t838 + t466 * t834;
t436 = -t463 * t838 - t464 * t834;
t435 = -pkin(9) * t534 - t444;
t434 = -t456 * t834 + t457 * t838;
t433 = t456 * t838 + t457 * t834;
t432 = t440 * t835 + t559 * t839;
t431 = -t440 * t839 + t559 * t835;
t430 = t445 * t831 + t506 * t830;
t429 = t445 * t830 - t506 * t831;
t428 = t437 * t835 + t553 * t839;
t427 = -t437 * t839 + t553 * t835;
t426 = -qJ(5) * t516 - t454 * t830 + t487 * t831;
t425 = -qJ(5) * t510 - t453 * t830 + t486 * t831;
t424 = t433 * t835 + t534 * t839;
t423 = -t433 * t839 + t534 * t835;
t422 = -pkin(8) * t513 - t455 * t833 + t461 * t837;
t421 = pkin(2) * t468 + pkin(3) * t513 + pkin(4) * t573 - qJ(2) * t469;
t420 = -pkin(4) * t559 + qJ(5) * t517 + t454 * t831 + t487 * t830;
t419 = -pkin(4) * t553 + qJ(5) * t511 + t453 * t831 + t486 * t830;
t418 = -pkin(3) * t659 + pkin(8) * t515 + t455 * t837 + t461 * t833;
t417 = pkin(2) * t853 - t838 * t473 - t834 * t483 - t523 * t928;
t416 = -qJ(5) * t501 + t435 * t831 + t534 * t927;
t415 = -t447 * t834 + t448 * t838;
t414 = t447 * t838 + t448 * t834;
t413 = pkin(2) * t631 - t838 * t460 - t834 * t474 - t499 * t928;
t412 = -pkin(8) * t447 - qJ(5) * t923 - t475 * t833;
t411 = t414 * t835 - t622 * t839;
t410 = -t414 * t839 - t622 * t835;
t409 = pkin(3) * t622 + pkin(8) * t448 - qJ(5) * t924 + t475 * t837;
t408 = qJ(5) * t502 + t830 * t435 + t534 * t866;
t407 = -t429 * t833 + t430 * t837;
t406 = t429 * t837 + t430 * t833;
t405 = -qJ(5) * t429 + (-pkin(9) * t831 + t927) * t444;
t404 = -pkin(8) * t470 - t420 * t833 + t426 * t837;
t403 = -pkin(8) * t465 - t419 * t833 + t425 * t837;
t402 = pkin(2) * t440 + pkin(3) * t470 + pkin(4) * t516 - pkin(5) * t596 + pkin(9) * t560 - qJ(2) * t441 + t922;
t401 = -pkin(3) * t559 + pkin(8) * t471 + t420 * t837 + t426 * t833;
t400 = pkin(2) * t437 + pkin(3) * t465 + pkin(4) * t510 + pkin(5) * t594 + pkin(9) * t554 - qJ(2) * t438 - t921;
t399 = -pkin(3) * t553 + pkin(8) * t466 + t419 * t837 + t425 * t833;
t398 = pkin(2) * t659 - t838 * t418 - t834 * t422 - t469 * t928;
t397 = qJ(5) * t430 + (-pkin(9) * t830 + t866) * t444;
t396 = -pkin(8) * t456 - t408 * t833 + t416 * t837;
t395 = -pkin(3) * t534 + pkin(8) * t457 + t408 * t837 + t416 * t833;
t394 = pkin(2) * t433 + pkin(3) * t456 + pkin(4) * t501 + pkin(5) * t630 + pkin(9) * t536 - qJ(2) * t434 + t445;
t393 = pkin(2) * t414 + pkin(3) * t447 + pkin(4) * t484 - qJ(2) * t415;
t392 = -t406 * t834 + t407 * t838;
t391 = t406 * t838 + t407 * t834;
t390 = t391 * t835 + t444 * t839;
t389 = -t391 * t839 + t444 * t835;
t388 = -pkin(2) * t622 - t838 * t409 - t834 * t412 - t415 * t928;
t387 = pkin(2) * t559 - t838 * t401 - t834 * t404 - t441 * t928;
t386 = pkin(2) * t553 - t838 * t399 - t834 * t403 - t438 * t928;
t385 = -pkin(8) * t406 - t397 * t833 + t405 * t837;
t384 = -pkin(3) * t444 + pkin(8) * t407 + t397 * t837 + t405 * t833;
t383 = pkin(2) * t534 - t838 * t395 - t834 * t396 - t434 * t928;
t382 = pkin(2) * t391 + pkin(3) * t406 + pkin(4) * t429 - pkin(5) * t506 + pkin(9) * t445 - qJ(2) * t392;
t381 = pkin(2) * t444 - t838 * t384 - t834 * t385 - t392 * t928;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t798, -t797, 0, t756, 0, 0, 0, 0, 0, 0, 0, t798, t797, t731, 0, 0, 0, 0, 0, 0, t732, t733, t754, t668, 0, 0, 0, 0, 0, 0, t600, t606, t558, t521, 0, 0, 0, 0, 0, 0, t489, t497, t459, t411, 0, 0, 0, 0, 0, 0, t428, t432, t424, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t797, -t798, 0, t755, 0, 0, 0, 0, 0, 0, 0, -t797, t798, t728, 0, 0, 0, 0, 0, 0, t729, t730, t753, t667, 0, 0, 0, 0, 0, 0, t599, t605, t557, t520, 0, 0, 0, 0, 0, 0, t488, t496, t458, t410, 0, 0, 0, 0, 0, 0, t427, t431, t423, t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t761, t762, 0, t701, 0, 0, 0, 0, 0, 0, t621, t642, t580, t532, 0, 0, 0, 0, 0, 0, t499, t523, t469, t415, 0, 0, 0, 0, 0, 0, t438, t441, t434, t392; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t797, 0, -t798, 0, -t856, -t855, -t755, -pkin(6) * t755, 0, -t797, t798, 0, 0, 0, -t728, t856, t855, -pkin(6) * t728 + (-pkin(1) * t835 + qJ(2) * t839) * g(3), -t768 * t835 + t857, -t750 * t835 + t800 * t839, -t759 * t835 + t838 * t875, -t767 * t835 - t857, -t757 * t835 - t834 * t875, qJDD(3) * t839 - t790 * t835, -pkin(6) * t729 - t674 * t835 + t681 * t839, -pkin(6) * t730 - t673 * t835 + t682 * t839, -pkin(2) * t901 - pkin(6) * t753 - t687 * t835, -pkin(6) * t667 - t625 * t835 + t653 * t839, -t624 * t835 - t868, -t578 * t835 + t748 * t839, -t649 * t835 - t696 * t839, -t623 * t835 + t868, -t650 * t835 + t839 * t844, -t666 * t835 + t816, -pkin(6) * t599 - t490 * t835 + t530 * t839, -pkin(6) * t605 - t500 * t835 + t541 * t839, -pkin(6) * t557 - t452 * t835 + t495 * t839, -pkin(6) * t520 - t449 * t835 + t462 * t839, -t504 * t835 + t870, -t467 * t835 + t683 * t839, -t528 * t835 - t635 * t839, -t503 * t835 - t870, -t529 * t835 + t839 * t846, -t548 * t835 + t816, -pkin(6) * t488 - t413 * t835 + t439 * t839, -pkin(6) * t496 - t417 * t835 + t446 * t839, -pkin(6) * t458 - t398 * t835 + t421 * t839, -pkin(6) * t410 - t388 * t835 + t393 * t839, -t451 * t835 + t589 * t839, -t436 * t835 + t533 * t839, -t442 * t835 + t568 * t839, -t450 * t835 + t587 * t839, -t443 * t835 + t569 * t839, -t472 * t835 + t614 * t839, -pkin(6) * t427 - t386 * t835 + t400 * t839, -pkin(6) * t431 - t387 * t835 + t402 * t839, -pkin(6) * t423 - t383 * t835 + t394 * t839, -pkin(6) * t389 - t381 * t835 + t382 * t839; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t798, 0, t797, 0, t855, -t856, t756, pkin(6) * t756, 0, -t798, -t797, 0, 0, 0, t731, -t855, t856, pkin(6) * t731 + (pkin(1) * t839 + qJ(2) * t835) * g(3), t768 * t839 + t858, t750 * t839 + t800 * t835, t759 * t839 + t835 * t876, t767 * t839 - t858, t757 * t839 - t834 * t877, qJDD(3) * t835 + t790 * t839, pkin(6) * t732 + t674 * t839 + t681 * t835, pkin(6) * t733 + t673 * t839 + t682 * t835, -pkin(2) * t902 + pkin(6) * t754 + t687 * t839, pkin(6) * t668 + t625 * t839 + t653 * t835, t624 * t839 - t869, t578 * t839 + t748 * t835, t649 * t839 - t696 * t835, t623 * t839 + t869, t650 * t839 + t835 * t844, t666 * t839 + t814, pkin(6) * t600 + t490 * t839 + t530 * t835, pkin(6) * t606 + t500 * t839 + t541 * t835, pkin(6) * t558 + t452 * t839 + t495 * t835, pkin(6) * t521 + t449 * t839 + t462 * t835, t504 * t839 + t871, t467 * t839 + t683 * t835, t528 * t839 - t635 * t835, t503 * t839 - t871, t529 * t839 + t835 * t846, t548 * t839 + t814, pkin(6) * t489 + t413 * t839 + t439 * t835, pkin(6) * t497 + t417 * t839 + t446 * t835, pkin(6) * t459 + t398 * t839 + t421 * t835, pkin(6) * t411 + t388 * t839 + t393 * t835, t451 * t839 + t589 * t835, t436 * t839 + t533 * t835, t442 * t839 + t568 * t835, t450 * t839 + t587 * t835, t443 * t839 + t569 * t835, t472 * t839 + t614 * t835, pkin(6) * t428 + t386 * t839 + t400 * t835, pkin(6) * t432 + t387 * t839 + t402 * t835, pkin(6) * t424 + t383 * t839 + t394 * t835, pkin(6) * t390 + t381 * t839 + t382 * t835; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t803, t804, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t854 - 0.2e1 * t925, -t804 + t822 + 0.2e1 * t827, pkin(1) * t777 + qJ(2) * t772, (t794 - t865) * t838, -t792 * t838 - t795 * t834, -t808 * t834 + t897, (-t793 + t864) * t834, t806 * t838 - t900, 0, qJ(2) * t792 - t758 * t928 - t890, qJ(2) * t795 - t760 * t928 - t889, -qJ(2) * t799 + t796 * t928 - t700, -qJ(2) * t766 - t700 * t928, -t690 * t834 + t691 * t838, -t637 * t834 + t639 * t838, -t702 * t834 + t704 * t838, -t688 * t834 + t689 * t838, -t703 * t834 + t705 * t838, -t724 * t834 + t725 * t838, qJ(2) * t692 - t834 * t607 - t620 * t928 + t838 * t644, qJ(2) * t848 - t834 * t611 - t641 * t928 + t838 * t656, qJ(2) * t722 - t834 * t546 + t838 * t547 - t579 * t928, -pkin(8) * t919 - qJ(2) * t710 - t531 * t928 - t834 * t562, -t564 * t834 + t566 * t838, -t512 * t834 + t514 * t838, -t583 * t834 + t585 * t838, -t563 * t834 + t565 * t838, -t584 * t834 + t586 * t838, -t612 * t834 + t613 * t838, qJ(2) * t631 - t834 * t460 + t838 * t474 - t498 * t928, qJ(2) * t853 - t834 * t473 + t838 * t483 - t522 * t928, qJ(2) * t659 - t834 * t418 + t838 * t422 - t468 * t928, -qJ(2) * t622 - t834 * t409 + t838 * t412 - t414 * t928, -t492 * t834 + t494 * t838, -t463 * t834 + t464 * t838, -t476 * t834 + t478 * t838, -t491 * t834 + t493 * t838, -t477 * t834 + t479 * t838, -t518 * t834 + t519 * t838, qJ(2) * t553 - t834 * t399 + t838 * t403 - t437 * t928, qJ(2) * t559 - t834 * t401 + t838 * t404 - t440 * t928, qJ(2) * t534 - t834 * t395 + t838 * t396 - t433 * t928, qJ(2) * t444 - t834 * t384 + t838 * t385 - t391 * t928;];
tauB_reg  = t1;