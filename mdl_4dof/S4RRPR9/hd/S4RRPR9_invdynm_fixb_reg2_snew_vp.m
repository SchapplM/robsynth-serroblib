% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPR9_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:15
% EndTime: 2019-12-31 17:10:21
% DurationCPUTime: 6.58s
% Computational Cost: add. (29663->486), mult. (65221->661), div. (0->0), fcn. (43424->8), ass. (0->327)
t823 = sin(qJ(4));
t821 = sin(pkin(7));
t822 = cos(pkin(7));
t824 = sin(qJ(2));
t864 = qJD(1) * t824;
t783 = -t822 * qJD(2) + t821 * t864;
t784 = t821 * qJD(2) + t822 * t864;
t826 = cos(qJ(4));
t742 = t826 * t783 + t823 * t784;
t744 = -t823 * t783 + t826 * t784;
t696 = t744 * t742;
t861 = qJD(1) * qJD(2);
t812 = t824 * t861;
t827 = cos(qJ(2));
t858 = t827 * qJDD(1);
t792 = -t812 + t858;
t785 = -qJDD(4) + t792;
t895 = -t696 - t785;
t899 = t823 * t895;
t898 = t826 * t895;
t886 = t784 * t783;
t838 = -t792 - t886;
t897 = t821 * t838;
t896 = t822 * t838;
t813 = t827 * t861;
t860 = t824 * qJDD(1);
t791 = t813 + t860;
t766 = t821 * qJDD(2) + t822 * t791;
t840 = t822 * qJDD(2) - t821 * t791;
t682 = -t742 * qJD(4) + t826 * t766 + t823 * t840;
t863 = t827 * qJD(1);
t809 = -qJD(4) + t863;
t731 = t742 * t809;
t894 = t731 + t682;
t773 = t784 * t863;
t734 = t840 - t773;
t851 = t823 * t766 - t826 * t840;
t655 = (qJD(4) + t809) * t744 + t851;
t740 = t742 ^ 2;
t741 = t744 ^ 2;
t893 = t783 ^ 2;
t782 = t784 ^ 2;
t807 = t809 ^ 2;
t892 = qJD(2) ^ 2;
t891 = pkin(2) * t824;
t890 = pkin(2) * t827;
t825 = sin(qJ(1));
t828 = cos(qJ(1));
t801 = t825 * g(1) - t828 * g(2);
t829 = qJD(1) ^ 2;
t778 = qJDD(1) * pkin(1) + t829 * pkin(5) + t801;
t842 = t791 + t813;
t709 = -t842 * qJ(3) + (-t792 + t812) * pkin(2) - t778;
t802 = t828 * g(1) + t825 * g(2);
t779 = -t829 * pkin(1) + qJDD(1) * pkin(5) - t802;
t762 = -t824 * g(3) + t827 * t779;
t843 = -qJ(3) * t824 - t890;
t789 = t843 * qJD(1);
t718 = -t892 * pkin(2) + qJDD(2) * qJ(3) + t789 * t863 + t762;
t662 = 0.2e1 * qJD(3) * t784 - t822 * t709 + t821 * t718;
t852 = t783 * t863;
t845 = -t766 + t852;
t630 = pkin(3) * t838 + t845 * pkin(6) - t662;
t663 = -0.2e1 * qJD(3) * t783 + t821 * t709 + t822 * t718;
t767 = -pkin(3) * t863 - t784 * pkin(6);
t635 = -t893 * pkin(3) + pkin(6) * t840 + t767 * t863 + t663;
t597 = -t826 * t630 + t823 * t635;
t598 = t823 * t630 + t826 * t635;
t568 = -t826 * t597 + t823 * t598;
t889 = pkin(3) * t568;
t659 = -t731 + t682;
t624 = -t655 * t823 - t826 * t659;
t888 = pkin(3) * t624;
t887 = t827 * g(3);
t885 = t809 * t744;
t884 = t809 * t823;
t883 = t809 * t826;
t818 = t824 ^ 2;
t882 = t818 * t829;
t881 = t821 * t568;
t717 = t887 + qJDD(3) - t892 * qJ(3) - qJDD(2) * pkin(2) + (qJD(1) * t789 + t779) * t824;
t880 = t821 * t717;
t737 = t792 - t886;
t879 = t821 * t737;
t878 = t822 * t568;
t877 = t822 * t717;
t876 = t822 * t737;
t666 = -pkin(3) * t840 - t893 * pkin(6) + t784 * t767 + t717;
t875 = t823 * t666;
t683 = -t696 + t785;
t874 = t823 * t683;
t873 = t824 * t778;
t872 = t824 * t792;
t808 = t827 * t829 * t824;
t799 = qJDD(2) + t808;
t871 = t824 * t799;
t800 = qJDD(2) - t808;
t870 = t824 * t800;
t869 = t826 * t666;
t868 = t826 * t683;
t867 = t827 * t778;
t866 = t827 * t800;
t819 = t827 ^ 2;
t865 = t818 + t819;
t859 = t825 * qJDD(1);
t857 = t828 * qJDD(1);
t856 = t824 * t696;
t855 = t824 * t886;
t854 = t827 * t696;
t853 = t827 * t886;
t569 = t823 * t597 + t826 * t598;
t628 = t821 * t662 + t822 * t663;
t761 = t824 * t779 + t887;
t712 = t824 * t761 + t827 * t762;
t850 = -t825 * t801 - t828 * t802;
t849 = t825 * t808;
t848 = t828 * t808;
t847 = t821 * t852;
t796 = -t825 * t829 + t857;
t846 = -pkin(4) * t796 - t825 * g(3);
t844 = -pkin(2) * t717 + qJ(3) * t628;
t627 = -t822 * t662 + t821 * t663;
t711 = t827 * t761 - t824 * t762;
t841 = t828 * t801 - t825 * t802;
t693 = -t807 - t740;
t638 = t823 * t693 + t898;
t839 = pkin(3) * t638 - t597;
t815 = t819 * t829;
t745 = -t815 - t893;
t691 = t822 * t745 - t897;
t733 = t773 + t840;
t837 = pkin(2) * t733 + qJ(3) * t691 - t877;
t771 = -t782 - t815;
t700 = -t821 * t771 + t876;
t736 = t766 + t852;
t836 = -pkin(2) * t736 + qJ(3) * t700 + t880;
t714 = -t741 - t807;
t640 = t826 * t714 + t874;
t835 = pkin(3) * t640 - t598;
t680 = t822 * t734 - t821 * t845;
t729 = t782 + t893;
t834 = pkin(2) * t729 + qJ(3) * t680 + t628;
t626 = -t655 * t826 + t823 * t659;
t674 = -t740 - t741;
t559 = -pkin(3) * t674 + pkin(6) * t626 + t569;
t561 = -pkin(6) * t624 - t568;
t582 = -t821 * t624 + t822 * t626;
t833 = -pkin(2) * t674 + qJ(3) * t582 + t822 * t559 + t821 * t561;
t639 = t826 * t693 - t899;
t654 = (qJD(4) - t809) * t744 + t851;
t589 = -pkin(3) * t654 + pkin(6) * t639 - t869;
t606 = -t821 * t638 + t822 * t639;
t610 = -pkin(6) * t638 + t875;
t832 = -pkin(2) * t654 + qJ(3) * t606 + t822 * t589 + t821 * t610;
t641 = -t823 * t714 + t868;
t599 = -pkin(3) * t894 + pkin(6) * t641 + t875;
t612 = -t821 * t640 + t822 * t641;
t622 = -pkin(6) * t640 + t869;
t831 = -pkin(2) * t894 + qJ(3) * t612 + t822 * t599 + t821 * t622;
t557 = t822 * t569 - t881;
t564 = -pkin(3) * t666 + pkin(6) * t569;
t830 = -pkin(2) * t666 - pkin(6) * t881 + qJ(3) * t557 + t822 * t564;
t806 = -t815 - t892;
t805 = t815 - t892;
t804 = -t882 - t892;
t803 = -t882 + t892;
t798 = -t815 + t882;
t797 = t815 + t882;
t795 = t828 * t829 + t859;
t794 = t865 * qJDD(1);
t793 = -0.2e1 * t812 + t858;
t790 = 0.2e1 * t813 + t860;
t787 = t827 * t799;
t786 = t865 * t861;
t780 = t827 * t792;
t774 = -pkin(4) * t795 + t828 * g(3);
t770 = -t782 + t815;
t769 = -t815 + t893;
t768 = t822 * t773;
t765 = t827 * t791 - t818 * t861;
t764 = -t819 * t861 - t872;
t760 = -t824 * t804 - t866;
t759 = -t824 * t803 + t787;
t758 = t827 * t806 - t871;
t757 = t827 * t805 - t870;
t756 = t827 * t804 - t870;
t755 = t827 * t803 + t871;
t754 = t824 * t806 + t787;
t753 = t824 * t805 + t866;
t752 = t842 * t824;
t751 = -t824 * t813 + t780;
t748 = t782 - t893;
t747 = -t824 * t790 + t827 * t793;
t746 = t827 * t790 + t824 * t793;
t728 = (t783 * t822 - t784 * t821) * t863;
t727 = -t768 - t847;
t726 = -t741 + t807;
t725 = t740 - t807;
t724 = -pkin(5) * t756 - t867;
t723 = -pkin(5) * t754 - t873;
t722 = t822 * t766 + t821 * t773;
t721 = t821 * t766 - t768;
t720 = -t821 * t840 - t822 * t852;
t719 = t822 * t840 - t847;
t716 = -pkin(1) * t756 + t762;
t715 = -pkin(1) * t754 + t761;
t708 = pkin(1) * t793 + pkin(5) * t758 + t867;
t707 = -pkin(1) * t790 + pkin(5) * t760 - t873;
t704 = t827 * t728 - t872;
t703 = t824 * t728 + t780;
t702 = t822 * t769 + t879;
t701 = -t821 * t770 + t896;
t699 = t821 * t769 - t876;
t698 = t822 * t770 + t897;
t697 = t822 * t771 + t879;
t695 = t741 - t740;
t694 = pkin(1) * t778 + pkin(5) * t712;
t692 = pkin(1) * t797 + pkin(5) * t794 + t712;
t690 = t821 * t745 + t896;
t689 = t827 * t722 + t855;
t688 = t827 * t720 - t855;
t687 = t824 * t722 - t853;
t686 = t824 * t720 + t853;
t681 = -t744 * qJD(4) - t851;
t679 = t822 * t733 - t821 * t736;
t678 = t821 * t734 + t822 * t845;
t677 = t821 * t733 + t822 * t736;
t676 = (t742 * t826 - t744 * t823) * t809;
t675 = (t742 * t823 + t744 * t826) * t809;
t673 = t827 * t702 + t734 * t824;
t672 = t827 * t701 - t824 * t845;
t671 = t827 * t700 + t824 * t736;
t670 = t824 * t702 - t734 * t827;
t669 = t824 * t701 + t827 * t845;
t668 = t824 * t700 - t827 * t736;
t667 = -qJ(3) * t697 + t877;
t665 = t827 * t679 + t824 * t748;
t664 = t824 * t679 - t827 * t748;
t661 = t827 * t691 - t824 * t733;
t660 = t824 * t691 + t827 * t733;
t652 = t826 * t725 + t874;
t651 = -t823 * t726 + t898;
t650 = t823 * t725 - t868;
t649 = t826 * t726 + t899;
t648 = t826 * t682 + t744 * t884;
t647 = t823 * t682 - t744 * t883;
t646 = -t823 * t681 - t742 * t883;
t645 = t826 * t681 - t742 * t884;
t644 = -qJ(3) * t690 + t880;
t643 = t827 * t680 - t824 * t729;
t642 = t824 * t680 + t827 * t729;
t637 = -t821 * t675 + t822 * t676;
t636 = t822 * t675 + t821 * t676;
t634 = -pkin(2) * t697 + t663;
t633 = t827 * t637 - t824 * t785;
t632 = t824 * t637 + t827 * t785;
t631 = -pkin(2) * t690 + t662;
t625 = -t826 * t654 - t823 * t894;
t623 = -t823 * t654 + t826 * t894;
t621 = -t821 * t650 + t822 * t652;
t620 = -t821 * t649 + t822 * t651;
t619 = t822 * t650 + t821 * t652;
t618 = t822 * t649 + t821 * t651;
t617 = -t821 * t647 + t822 * t648;
t616 = -t821 * t645 + t822 * t646;
t615 = t822 * t647 + t821 * t648;
t614 = t822 * t645 + t821 * t646;
t613 = -pkin(1) * t668 - t836;
t611 = t822 * t640 + t821 * t641;
t609 = -pkin(1) * t660 - t837;
t608 = t827 * t628 + t824 * t717;
t607 = t824 * t628 - t827 * t717;
t605 = t822 * t638 + t821 * t639;
t604 = -qJ(3) * t678 - t627;
t603 = t827 * t617 + t856;
t602 = t827 * t616 - t856;
t601 = t824 * t617 - t854;
t600 = t824 * t616 + t854;
t596 = t827 * t621 - t824 * t655;
t595 = t827 * t620 + t824 * t659;
t594 = t824 * t621 + t827 * t655;
t593 = t824 * t620 - t827 * t659;
t591 = t827 * t612 + t824 * t894;
t590 = t824 * t612 - t827 * t894;
t588 = -pkin(5) * t668 - t824 * t634 + t827 * t667;
t587 = t827 * t606 + t824 * t654;
t586 = t824 * t606 - t827 * t654;
t585 = -pkin(5) * t660 - t824 * t631 + t827 * t644;
t584 = -pkin(1) * t697 + pkin(5) * t671 + t827 * t634 + t824 * t667;
t583 = -pkin(1) * t642 - t834;
t581 = -t821 * t623 + t822 * t625;
t580 = t822 * t624 + t821 * t626;
t579 = t822 * t623 + t821 * t625;
t578 = -pkin(1) * t690 + pkin(5) * t661 + t827 * t631 + t824 * t644;
t577 = -pkin(5) * t642 + t827 * t604 + t678 * t891;
t576 = t827 * t581 + t824 * t695;
t575 = t824 * t581 - t827 * t695;
t574 = -pkin(1) * t607 - t844;
t573 = t827 * t582 + t824 * t674;
t572 = t824 * t582 - t827 * t674;
t571 = pkin(5) * t643 + t824 * t604 + (-pkin(1) - t890) * t678;
t570 = -pkin(2) * t580 - t888;
t567 = -pkin(2) * t611 - t835;
t566 = -pkin(5) * t607 + (-qJ(3) * t827 + t891) * t627;
t565 = -pkin(2) * t605 - t839;
t563 = -qJ(3) * t611 - t821 * t599 + t822 * t622;
t562 = -qJ(3) * t605 - t821 * t589 + t822 * t610;
t560 = pkin(5) * t608 + (-pkin(1) + t843) * t627;
t558 = -pkin(1) * t590 - t831;
t556 = t821 * t569 + t878;
t555 = -pkin(1) * t586 - t832;
t554 = t827 * t557 + t824 * t666;
t553 = t824 * t557 - t827 * t666;
t552 = -pkin(5) * t590 + t827 * t563 - t824 * t567;
t551 = -pkin(5) * t586 + t827 * t562 - t824 * t565;
t550 = -pkin(2) * t556 - t889;
t549 = -pkin(1) * t611 + pkin(5) * t591 + t824 * t563 + t827 * t567;
t548 = -pkin(1) * t605 + pkin(5) * t587 + t824 * t562 + t827 * t565;
t547 = -qJ(3) * t580 - t821 * t559 + t822 * t561;
t546 = -pkin(6) * t878 - qJ(3) * t556 - t821 * t564;
t545 = -pkin(1) * t572 - t833;
t544 = -pkin(5) * t572 + t827 * t547 - t824 * t570;
t543 = -pkin(1) * t580 + pkin(5) * t573 + t824 * t547 + t827 * t570;
t542 = -pkin(1) * t553 - t830;
t541 = -pkin(5) * t553 + t827 * t546 - t824 * t550;
t540 = -pkin(1) * t556 + pkin(5) * t554 + t824 * t546 + t827 * t550;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t796, 0, -t795, 0, t846, -t774, -t841, -pkin(4) * t841, t828 * t765 - t849, t828 * t747 + t825 * t798, t828 * t759 + t824 * t859, t828 * t764 + t849, t828 * t757 + t825 * t858, t825 * qJDD(2) + t828 * t786, t828 * t723 - t825 * t715 - pkin(4) * (t825 * t758 + t828 * t793), t828 * t724 - t825 * t716 - pkin(4) * (t825 * t760 - t828 * t790), t828 * t711 - pkin(4) * (t825 * t794 + t828 * t797), -pkin(4) * (t825 * t712 + t828 * t778) - (t825 * pkin(1) - t828 * pkin(5)) * t711, t828 * t689 + t825 * t721, t828 * t665 + t825 * t677, t828 * t672 + t825 * t698, t828 * t688 + t825 * t719, t828 * t673 + t825 * t699, t828 * t704 - t825 * t727, t828 * t585 - t825 * t609 - pkin(4) * (t825 * t661 - t828 * t690), t828 * t588 - t825 * t613 - pkin(4) * (t825 * t671 - t828 * t697), t828 * t577 - t825 * t583 - pkin(4) * (t825 * t643 - t828 * t678), t828 * t566 - t825 * t574 - pkin(4) * (t825 * t608 - t828 * t627), t828 * t603 + t825 * t615, t828 * t576 + t825 * t579, t828 * t595 + t825 * t618, t828 * t602 + t825 * t614, t828 * t596 + t825 * t619, t828 * t633 + t825 * t636, t828 * t551 - t825 * t555 - pkin(4) * (t825 * t587 - t828 * t605), t828 * t552 - t825 * t558 - pkin(4) * (t825 * t591 - t828 * t611), t828 * t544 - t825 * t545 - pkin(4) * (t825 * t573 - t828 * t580), t828 * t541 - t825 * t542 - pkin(4) * (t825 * t554 - t828 * t556); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t795, 0, t796, 0, t774, t846, t850, pkin(4) * t850, t825 * t765 + t848, t825 * t747 - t828 * t798, t825 * t759 - t824 * t857, t825 * t764 - t848, t825 * t757 - t827 * t857, -t828 * qJDD(2) + t825 * t786, t825 * t723 + t828 * t715 + pkin(4) * (t828 * t758 - t825 * t793), t825 * t724 + t828 * t716 + pkin(4) * (t828 * t760 + t825 * t790), t825 * t711 + pkin(4) * (t828 * t794 - t825 * t797), pkin(4) * (t828 * t712 - t825 * t778) - (-t828 * pkin(1) - t825 * pkin(5)) * t711, t825 * t689 - t828 * t721, t825 * t665 - t828 * t677, t825 * t672 - t828 * t698, t825 * t688 - t828 * t719, t825 * t673 - t828 * t699, t825 * t704 + t828 * t727, t825 * t585 + t828 * t609 + pkin(4) * (t828 * t661 + t825 * t690), t825 * t588 + t828 * t613 + pkin(4) * (t828 * t671 + t825 * t697), t825 * t577 + t828 * t583 + pkin(4) * (t828 * t643 + t825 * t678), t825 * t566 + t828 * t574 + pkin(4) * (t828 * t608 + t825 * t627), t825 * t603 - t828 * t615, t825 * t576 - t828 * t579, t825 * t595 - t828 * t618, t825 * t602 - t828 * t614, t825 * t596 - t828 * t619, t825 * t633 - t828 * t636, t825 * t551 + t828 * t555 + pkin(4) * (t828 * t587 + t825 * t605), t825 * t552 + t828 * t558 + pkin(4) * (t828 * t591 + t825 * t611), t825 * t544 + t828 * t545 + pkin(4) * (t828 * t573 + t825 * t580), t825 * t541 + t828 * t542 + pkin(4) * (t828 * t554 + t825 * t556); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t801, t802, 0, 0, t752, t746, t755, t751, t753, 0, t708, t707, t692, t694, t687, t664, t669, t686, t670, t703, t578, t584, t571, t560, t601, t575, t593, t600, t594, t632, t548, t549, t543, t540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t829, 0, 0, -g(3), -t801, 0, t765, t747, t759, t764, t757, t786, t723, t724, t711, pkin(5) * t711, t689, t665, t672, t688, t673, t704, t585, t588, t577, t566, t603, t576, t595, t602, t596, t633, t551, t552, t544, t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t829, 0, qJDD(1), 0, g(3), 0, -t802, 0, t808, -t798, -t860, -t808, -t858, -qJDD(2), t715, t716, 0, pkin(1) * t711, -t721, -t677, -t698, -t719, -t699, t727, t609, t613, t583, t574, -t615, -t579, -t618, -t614, -t619, -t636, t555, t558, t545, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t801, t802, 0, 0, t752, t746, t755, t751, t753, 0, t708, t707, t692, t694, t687, t664, t669, t686, t670, t703, t578, t584, t571, t560, t601, t575, t593, t600, t594, t632, t548, t549, t543, t540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t791, t793, t799, -t813, t805, t813, 0, -t778, t761, 0, t722, t679, t701, t720, t702, t728, t644, t667, t604, -qJ(3) * t627, t617, t581, t620, t616, t621, t637, t562, t563, t547, t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, t790, t803, t792, t800, -t812, t778, 0, t762, 0, -t886, -t748, t845, t886, -t734, t792, t631, t634, -pkin(2) * t678, -pkin(2) * t627, -t696, -t695, -t659, t696, t655, t785, t565, t567, t570, t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t808, t798, t860, t808, t858, qJDD(2), -t761, -t762, 0, 0, t721, t677, t698, t719, t699, -t727, t837, t836, t834, t844, t615, t579, t618, t614, t619, t636, t832, t831, t833, t830; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t766, t733, t838, -t852, t769, t852, 0, t717, t662, 0, t648, t625, t651, t646, t652, t676, t610, t622, t561, -pkin(6) * t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, t736, t770, t840, -t737, t773, -t717, 0, t663, 0, t647, t623, t649, t645, t650, t675, t589, t599, t559, t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t886, t748, -t845, -t886, t734, -t792, -t662, -t663, 0, 0, t696, t695, t659, -t696, -t655, -t785, t839, t835, t888, t889; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t682, -t654, t895, -t731, t725, t731, 0, t666, t597, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t885, t894, t726, t681, -t683, t885, -t666, 0, t598, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t696, t695, t659, -t696, -t655, -t785, -t597, -t598, 0, 0;];
m_new_reg = t1;
