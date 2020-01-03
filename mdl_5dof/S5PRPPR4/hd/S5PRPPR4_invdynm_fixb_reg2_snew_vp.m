% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRPPR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:04
% EndTime: 2019-12-31 17:37:11
% DurationCPUTime: 7.60s
% Computational Cost: add. (17404->491), mult. (39688->604), div. (0->0), fcn. (26045->8), ass. (0->321)
t787 = sin(pkin(8));
t782 = t787 ^ 2;
t789 = cos(pkin(8));
t783 = t789 ^ 2;
t796 = qJD(2) ^ 2;
t882 = (t782 + t783) * t796;
t739 = t787 * t882;
t792 = sin(qJ(2));
t794 = cos(qJ(2));
t836 = t794 * qJDD(2);
t707 = -t792 * t739 + t787 * t836;
t837 = t792 * qJDD(2);
t710 = t794 * t739 + t787 * t837;
t788 = sin(pkin(7));
t790 = cos(pkin(7));
t891 = t790 * t707 - t788 * t710;
t905 = qJ(1) * t891;
t791 = sin(qJ(5));
t793 = cos(qJ(5));
t806 = t787 * t791 + t789 * t793;
t731 = t806 * qJD(2);
t841 = qJD(2) * t789;
t842 = qJD(2) * t787;
t733 = -t791 * t841 + t793 * t842;
t863 = t733 * t731;
t900 = qJDD(5) - t863;
t904 = t791 * t900;
t903 = t793 * t900;
t754 = t790 * g(1) + t788 * g(2);
t824 = t788 * g(1) - t790 * g(2);
t800 = t794 * t754 - t792 * t824;
t822 = t792 * t754 + t794 * t824;
t823 = -t792 * t822 - t794 * t800;
t642 = t792 * t800 - t794 * t822;
t854 = t790 * t642;
t902 = -t788 * t823 + t854;
t860 = t788 * t642;
t901 = t790 * t823 + t860;
t899 = pkin(1) * t707;
t898 = pkin(5) * t707;
t897 = pkin(5) * t710;
t649 = t788 * t707 + t790 * t710;
t896 = qJ(1) * t649;
t752 = t794 * t796 + t837;
t785 = g(3) - qJDD(1);
t722 = pkin(5) * t752 - t794 * t785;
t753 = -t792 * t796 + t836;
t810 = -pkin(5) * t753 - t792 * t785;
t885 = t790 * t752 + t788 * t753;
t895 = qJ(1) * t885 + t790 * t722 - t788 * t810;
t700 = -t788 * t752 + t790 * t753;
t894 = -qJ(1) * t700 + t788 * t722 + t790 * t810;
t775 = t782 * qJDD(2);
t777 = t783 * qJDD(2);
t748 = t777 - t775;
t779 = t782 * t796;
t862 = t783 * t796;
t750 = -t779 + t862;
t695 = t792 * t748 + t794 * t750;
t698 = t794 * t748 - t792 * t750;
t893 = t790 * t695 + t788 * t698;
t892 = t788 * t695 - t790 * t698;
t740 = t789 * t882;
t826 = t789 * t836;
t709 = -t792 * t740 + t826;
t712 = t794 * t740 + t789 * t837;
t647 = t790 * t709 - t788 * t712;
t890 = t788 * t709 + t790 * t712;
t689 = -t796 * pkin(2) + qJDD(2) * qJ(3) - t800;
t811 = -pkin(3) * t789 - qJ(4) * t787;
t745 = t811 * qJD(2);
t805 = t689 + ((2 * qJD(3)) + t745) * qJD(2);
t768 = t789 * t785;
t835 = qJDD(4) + t768;
t855 = t789 * t796;
t609 = (-pkin(4) * t855 - pkin(6) * qJDD(2) + t805) * t787 + t835;
t838 = qJD(2) * qJD(3);
t828 = t789 * t838;
t765 = 0.2e1 * t828;
t846 = t789 * t689 - t787 * t785;
t829 = -t745 * t841 - t846;
t626 = t765 - t829;
t778 = t789 * qJDD(2);
t615 = -pkin(4) * t862 - pkin(6) * t778 + t626;
t549 = -t793 * t609 + t791 * t615;
t551 = t791 * t609 + t793 * t615;
t519 = -t793 * t549 + t791 * t551;
t666 = t768 + (t689 + 0.2e1 * t838) * t787;
t667 = t765 + t846;
t598 = t787 * t666 + t789 * t667;
t776 = t787 * qJDD(2);
t884 = -pkin(2) * t776 + qJ(3) * t739;
t883 = t806 * qJDD(2);
t520 = t791 * t549 + t793 * t551;
t874 = pkin(3) + pkin(4);
t881 = qJ(4) * t520 - t519 * t874;
t727 = t731 ^ 2;
t795 = qJD(5) ^ 2;
t679 = -t795 - t727;
t616 = t791 * t679 + t903;
t617 = t793 * t679 - t904;
t880 = qJ(4) * t617 - t616 * t874 + t549;
t730 = t776 * t793 - t791 * t778;
t618 = -t793 * t730 - t791 * t883;
t620 = t791 * t730 - t793 * t883;
t879 = qJ(4) * t620 - t618 * t874;
t728 = t733 ^ 2;
t718 = -t728 - t795;
t681 = qJDD(5) + t863;
t850 = t791 * t681;
t635 = t793 * t718 - t850;
t848 = t793 * t681;
t638 = -t791 * t718 - t848;
t878 = qJ(4) * t638 - t635 * t874 + t551;
t877 = -t790 * t754 - t788 * t824;
t875 = -t788 * t754 + t790 * t824;
t747 = t777 + t775;
t749 = t779 + t862;
t694 = t792 * t747 + t794 * t749;
t873 = pkin(5) * t694;
t872 = pkin(5) * t709;
t871 = pkin(6) * t519;
t870 = pkin(6) * t520;
t697 = t794 * t747 - t792 * t749;
t869 = qJ(1) * (t790 * t694 + t788 * t697);
t868 = qJ(1) * t647;
t683 = -qJDD(2) * pkin(2) - t796 * qJ(3) + qJDD(3) - t822;
t673 = t787 * t683;
t861 = t787 * t789;
t856 = t788 * t785;
t674 = t789 * t683;
t853 = t790 * t785;
t769 = qJ(4) * t776;
t772 = pkin(3) * t778;
t801 = -t683 + t772;
t827 = qJD(4) * t842;
t657 = -t769 - t801 - 0.2e1 * t827;
t624 = -pkin(4) * t778 + pkin(6) * t882 + t657;
t851 = t791 * t624;
t849 = t792 * t683;
t622 = t793 * t624;
t847 = t794 * t683;
t845 = -pkin(2) * t683 + qJ(3) * t598;
t844 = pkin(2) * t749 + qJ(3) * t747;
t843 = pkin(2) * t778 - qJ(3) * t740;
t840 = t731 * qJD(5);
t839 = t733 * qJD(5);
t834 = pkin(3) * t776;
t760 = t787 * t855;
t833 = t792 * t863;
t832 = t794 * t863;
t831 = t673 + t884;
t830 = -t674 + t843;
t759 = t787 * t778;
t825 = -pkin(6) * t635 - t622;
t668 = -t727 - t728;
t803 = -pkin(6) * t620 - t520;
t506 = t668 * t874 + t803;
t804 = -pkin(6) * t618 - t519;
t514 = qJ(4) * t668 + t804;
t565 = t787 * t618 + t789 * t620;
t821 = pkin(2) * t668 + qJ(3) * t565 + t789 * t506 + t787 * t514;
t502 = t787 * t519 + t789 * t520;
t510 = -t624 * t874 - t870;
t516 = -qJ(4) * t624 - t871;
t820 = -pkin(2) * t624 + qJ(3) * t502 + t789 * t510 + t787 * t516;
t684 = t883 + 0.2e1 * t839;
t808 = -pkin(6) * t617 - t622;
t526 = t684 * t874 + t808;
t809 = -pkin(6) * t616 - t851;
t540 = qJ(4) * t684 + t809;
t561 = t787 * t616 + t789 * t617;
t819 = pkin(2) * t684 + qJ(3) * t561 + t789 * t526 + t787 * t540;
t686 = -0.2e1 * t840 + t730;
t807 = -pkin(6) * t638 + t851;
t532 = t686 * t874 + t807;
t548 = qJ(4) * t686 + t825;
t584 = t787 * t635 + t789 * t638;
t818 = pkin(2) * t686 + qJ(3) * t584 + t789 * t532 + t787 * t548;
t625 = t787 * t805 + t835;
t613 = qJ(4) * t749 + t625;
t614 = pkin(3) * t749 + t626;
t817 = t787 * t613 + t789 * t614 + t844;
t763 = 0.2e1 * t827;
t658 = t763 + 0.2e1 * t769 + t801;
t816 = pkin(3) * t759 + t787 * t658 - t884;
t815 = t844 + t598;
t659 = -t683 + t763 + t769 + 0.2e1 * t772;
t814 = qJ(4) * t759 + t789 * t659 + t843;
t812 = -pkin(3) * t625 + qJ(4) * t626;
t597 = t789 * t666 - t787 * t667;
t714 = t752 * t861;
t715 = -t760 * t792 + t787 * t826;
t661 = t790 * t714 + t788 * t715;
t664 = t788 * t714 - t790 * t715;
t575 = t787 * t625 + t789 * t626;
t799 = qJ(3) * t575 + (-pkin(2) + t811) * t657;
t770 = qJ(4) * t778;
t758 = -0.2e1 * t759;
t757 = 0.2e1 * t759;
t742 = -t770 + t834;
t717 = -t728 + t795;
t716 = t727 - t795;
t703 = pkin(1) * t709;
t702 = pkin(5) * t712;
t692 = pkin(1) * t694;
t691 = pkin(5) * t697;
t688 = t728 - t727;
t687 = -t840 + t730;
t685 = -t883 - t839;
t672 = -pkin(1) * t752 + t800;
t671 = pkin(1) * t753 + t822;
t670 = (-t731 * t793 + t733 * t791) * qJD(5);
t669 = (-t731 * t791 - t733 * t793) * qJD(5);
t651 = t793 * t687 - t791 * t839;
t648 = t791 * t687 + t793 * t839;
t645 = -t791 * t685 + t793 * t840;
t644 = -t793 * t685 - t791 * t840;
t639 = qJ(1) * t890;
t637 = -t791 * t717 + t903;
t636 = t793 * t716 - t850;
t634 = t793 * t717 + t904;
t633 = t791 * t716 + t848;
t632 = pkin(1) * t642;
t631 = -pkin(3) * t775 + t789 * t658;
t630 = qJ(4) * t777 - t787 * t659;
t628 = pkin(1) * t785 + pkin(5) * t823;
t627 = qJ(1) * (-t788 * t694 + t790 * t697);
t621 = -t793 * t684 - t791 * t686;
t619 = -t791 * t684 + t793 * t686;
t612 = -0.2e1 * t828 + (-pkin(3) * t782 + qJ(4) * t861) * t796 + t829;
t608 = qJ(4) * t862 + (-pkin(3) * t855 + t805) * t787 + t835;
t606 = t703 + t830;
t605 = t831 - t899;
t602 = t787 * t669 + t789 * t670;
t601 = -t789 * t669 + t787 * t670;
t600 = -t792 * qJDD(5) + t794 * t602;
t599 = t794 * qJDD(5) + t792 * t602;
t594 = t816 + t899;
t593 = t703 + t814;
t592 = -t792 * t667 + t789 * t847 + t898;
t591 = -t792 * t666 + t787 * t847 - t872;
t590 = t794 * t667 + t789 * t849 + t897;
t589 = t794 * t666 + t787 * t849 - t702;
t588 = t787 * t648 + t789 * t651;
t587 = -t787 * t644 + t789 * t645;
t586 = -t789 * t648 + t787 * t651;
t585 = t789 * t644 + t787 * t645;
t583 = t787 * t634 + t789 * t637;
t582 = t787 * t633 + t789 * t636;
t581 = -t789 * t635 + t787 * t638;
t580 = -t789 * t634 + t787 * t637;
t579 = -t789 * t633 + t787 * t636;
t577 = t794 * t597 - t873;
t576 = t792 * t597 + t691;
t574 = -t789 * t625 + t787 * t626;
t572 = t794 * t598 + t849;
t571 = t792 * t598 - t847;
t570 = t794 * t583 - t792 * t730;
t569 = t794 * t582 + t792 * t883;
t568 = t792 * t583 + t794 * t730;
t567 = t792 * t582 - t794 * t883;
t566 = t787 * t619 + t789 * t621;
t564 = -t789 * t619 + t787 * t621;
t563 = -t789 * t618 + t787 * t620;
t560 = -t789 * t616 + t787 * t617;
t558 = t794 * t588 - t833;
t557 = t794 * t587 + t833;
t556 = t792 * t588 + t832;
t555 = t792 * t587 - t832;
t554 = t692 + t815;
t553 = t789 * t613 - t787 * t614;
t552 = t794 * t584 - t792 * t686;
t550 = t792 * t584 + t794 * t686;
t544 = -t792 * t612 + t794 * t631 - t898;
t543 = t794 * t612 + t792 * t631 - t897;
t542 = -t792 * t608 + t794 * t630 - t872;
t541 = t794 * t608 + t792 * t630 - t702;
t538 = t794 * t566 - t792 * t688;
t537 = t792 * t566 + t794 * t688;
t536 = t794 * t561 - t792 * t684;
t535 = t792 * t561 + t794 * t684;
t534 = t794 * t575 + t792 * t657;
t533 = t792 * t575 - t794 * t657;
t530 = t794 * t565 - t792 * t668;
t529 = t792 * t565 + t794 * t668;
t528 = t794 * t553 - t792 * t742 - t873;
t527 = t792 * t553 + t794 * t742 + t691;
t525 = t692 + t817;
t523 = -qJ(3) * t574 + (pkin(3) * t787 - qJ(4) * t789) * t657;
t522 = -pkin(2) * t574 - t812;
t521 = pkin(1) * t571 + t845;
t518 = -pkin(5) * t571 - (pkin(2) * t792 - qJ(3) * t794) * t597;
t517 = -pkin(2) * t563 - t879;
t512 = pkin(5) * t572 - (-pkin(2) * t794 - qJ(3) * t792 - pkin(1)) * t597;
t511 = -pkin(2) * t581 - t878;
t508 = -qJ(3) * t581 - t787 * t532 + t789 * t548;
t507 = pkin(1) * t533 + t799;
t504 = -pkin(2) * t560 - t880;
t503 = -qJ(3) * t560 - t787 * t526 + t789 * t540;
t501 = -t789 * t519 + t787 * t520;
t499 = pkin(1) * t550 + t818;
t498 = -pkin(5) * t533 - t792 * t522 + t794 * t523;
t497 = t794 * t502 + t792 * t624;
t496 = t792 * t502 - t794 * t624;
t495 = pkin(1) * t535 + t819;
t494 = -pkin(1) * t574 + pkin(5) * t534 + t794 * t522 + t792 * t523;
t493 = -qJ(3) * t563 - t787 * t506 + t789 * t514;
t492 = -pkin(5) * t550 + t794 * t508 - t792 * t511;
t491 = -pkin(1) * t581 + pkin(5) * t552 + t792 * t508 + t794 * t511;
t490 = -pkin(5) * t535 + t794 * t503 - t792 * t504;
t489 = pkin(1) * t529 + t821;
t488 = -pkin(1) * t560 + pkin(5) * t536 + t792 * t503 + t794 * t504;
t487 = -qJ(3) * t501 - t787 * t510 + t789 * t516;
t486 = -pkin(2) * t501 - t881;
t485 = -pkin(5) * t529 + t794 * t493 - t792 * t517;
t484 = -pkin(1) * t563 + pkin(5) * t530 + t792 * t493 + t794 * t517;
t483 = pkin(1) * t496 + t820;
t482 = -pkin(5) * t496 - t792 * t486 + t794 * t487;
t481 = -pkin(1) * t501 + pkin(5) * t497 + t794 * t486 + t792 * t487;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t856, -t853, -t875, -qJ(1) * t875, 0, 0, t700, 0, -t885, 0, t894, t895, t902, pkin(5) * t854 + qJ(1) * t902 - t788 * t628, -t664, -t892, t649, t664, t890, 0, -t788 * t589 + t790 * t591 - t868, -t788 * t590 + t790 * t592 + t905, -t788 * t576 + t790 * t577 - t869, t790 * t518 - t788 * t512 - qJ(1) * (t790 * t571 + t788 * t572), -t664, t649, t892, 0, -t890, t664, -t788 * t541 + t790 * t542 - t868, -t788 * t527 + t790 * t528 - t869, -t788 * t543 + t790 * t544 - t905, t790 * t498 - t788 * t494 - qJ(1) * (t790 * t533 + t788 * t534), -t788 * t556 + t790 * t558, -t788 * t537 + t790 * t538, -t788 * t568 + t790 * t570, -t788 * t555 + t790 * t557, -t788 * t567 + t790 * t569, -t788 * t599 + t790 * t600, t790 * t490 - t788 * t488 - qJ(1) * (t790 * t535 + t788 * t536), t790 * t492 - t788 * t491 - qJ(1) * (t790 * t550 + t788 * t552), t790 * t485 - t788 * t484 - qJ(1) * (t790 * t529 + t788 * t530), t790 * t482 - t788 * t481 - qJ(1) * (t790 * t496 + t788 * t497); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t853, -t856, t877, qJ(1) * t877, 0, 0, t885, 0, t700, 0, -t895, t894, t901, pkin(5) * t860 + qJ(1) * t901 + t790 * t628, t661, t893, -t891, -t661, -t647, 0, t790 * t589 + t788 * t591 - t639, t790 * t590 + t788 * t592 + t896, t790 * t576 + t788 * t577 + t627, t788 * t518 + t790 * t512 + qJ(1) * (-t788 * t571 + t790 * t572), t661, -t891, -t893, 0, t647, -t661, t790 * t541 + t788 * t542 - t639, t790 * t527 + t788 * t528 + t627, t790 * t543 + t788 * t544 - t896, t788 * t498 + t790 * t494 + qJ(1) * (-t788 * t533 + t790 * t534), t790 * t556 + t788 * t558, t790 * t537 + t788 * t538, t790 * t568 + t788 * t570, t790 * t555 + t788 * t557, t790 * t567 + t788 * t569, t790 * t599 + t788 * t600, t788 * t490 + t790 * t488 + qJ(1) * (-t788 * t535 + t790 * t536), t788 * t492 + t790 * t491 + qJ(1) * (-t788 * t550 + t790 * t552), t788 * t485 + t790 * t484 + qJ(1) * (-t788 * t529 + t790 * t530), t788 * t482 + t790 * t481 + qJ(1) * (-t788 * t496 + t790 * t497); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t824, t754, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t671, t672, 0, -t632, t775, t757, 0, t777, 0, 0, t606, t605, t554, t521, t775, 0, t758, 0, 0, t777, t593, t525, t594, t507, t586, t564, t580, t585, t579, t601, t495, t499, t489, t483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t785, -t824, 0, 0, 0, t753, 0, -t752, 0, t810, t722, t642, pkin(5) * t642, t715, t698, t710, -t715, t712, 0, t591, t592, t577, t518, t715, t710, -t698, 0, -t712, -t715, t542, t528, t544, t498, t558, t538, t570, t557, t569, t600, t490, t492, t485, t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t785, 0, -t754, 0, 0, 0, t752, 0, t753, 0, -t722, t810, t823, t628, t714, t695, -t707, -t714, -t709, 0, t589, t590, t576, t512, t714, -t707, -t695, 0, t709, -t714, t541, t527, t543, t494, t556, t537, t568, t555, t567, t599, t488, t491, t484, t481; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t824, t754, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t671, t672, 0, -t632, t775, t757, 0, t777, 0, 0, t606, t605, t554, t521, t775, 0, t758, 0, 0, t777, t593, t525, t594, t507, t586, t564, t580, t585, t579, t601, t495, t499, t489, t483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t796, 0, 0, -t785, -t822, 0, t759, t748, t739, -t759, t740, 0, t673, t674, t597, qJ(3) * t597, t759, t739, -t748, 0, -t740, -t759, t630, t553, t631, t523, t588, t566, t583, t587, t582, t602, t503, t508, t493, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t796, 0, qJDD(2), 0, t785, 0, -t800, 0, t760, t750, -t776, -t760, -t778, 0, t666, t667, 0, pkin(2) * t597, t760, -t776, -t750, 0, t778, -t760, t608, t742, t612, t522, t863, t688, t730, -t863, -t883, qJDD(5), t504, t511, t517, t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t822, t800, 0, 0, t775, t757, 0, t777, 0, 0, t830, t831, t815, t845, t775, 0, t758, 0, 0, t777, t814, t817, t816, t799, t586, t564, t580, t585, t579, t601, t819, t818, t821, t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t776, t778, t760, 0, t862, 0, 0, t683, t666, 0, t776, t760, -t778, 0, -t862, 0, t770, t613, t658, -qJ(4) * t657, t651, t621, t637, t645, t636, t670, t540, t548, t514, t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t776, -t779, t778, -t760, 0, -t683, 0, t667, 0, 0, -t779, -t776, 0, t760, t778, t659, t614, t834, -pkin(3) * t657, -t648, -t619, -t634, t644, -t633, -t669, t526, t532, t506, t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t760, -t750, t776, t760, t778, 0, -t666, -t667, 0, 0, -t760, t776, t750, 0, -t778, t760, -t608, -t742, -t612, t812, -t863, -t688, -t730, t863, t883, -qJDD(5), t880, t878, t879, t881; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t776, t760, -t778, 0, -t862, 0, 0, t625, -t657, 0, t651, t621, t637, t645, t636, t670, t809, t825, t804, -t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t760, t776, t750, 0, -t778, t760, -t625, 0, t626, 0, -t863, -t688, -t730, t863, t883, -qJDD(5), -pkin(4) * t616 + t549, -pkin(4) * t635 + t551, -pkin(4) * t618, -pkin(4) * t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t779, t776, 0, -t760, -t778, t657, -t626, 0, 0, t648, t619, t634, -t644, t633, t669, -pkin(4) * t684 - t808, -pkin(4) * t686 - t807, -pkin(4) * t668 - t803, pkin(4) * t624 + t870; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, -t684, t900, t840, t716, -t840, 0, -t624, t549, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t839, t686, t717, t685, t681, -t839, t624, 0, t551, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t863, t688, t730, -t863, -t883, qJDD(5), -t549, -t551, 0, 0;];
m_new_reg = t1;