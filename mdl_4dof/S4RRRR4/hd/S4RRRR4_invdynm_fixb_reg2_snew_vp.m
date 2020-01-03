% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRRR4
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRRR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:27
% EndTime: 2019-12-31 17:26:34
% DurationCPUTime: 6.75s
% Computational Cost: add. (31484->484), mult. (65512->667), div. (0->0), fcn. (44767->8), ass. (0->331)
t857 = sin(qJ(4));
t858 = sin(qJ(3));
t859 = sin(qJ(2));
t862 = cos(qJ(3));
t863 = cos(qJ(2));
t819 = (t858 * t863 + t859 * t862) * qJD(1);
t861 = cos(qJ(4));
t903 = qJD(2) + qJD(3);
t789 = t819 * t857 - t861 * t903;
t791 = t861 * t819 + t857 * t903;
t761 = t791 * t789;
t905 = qJD(1) * qJD(2);
t848 = t863 * t905;
t904 = t859 * qJDD(1);
t827 = t848 + t904;
t851 = t863 * qJDD(1);
t893 = t859 * t905;
t828 = t851 - t893;
t891 = t858 * t827 - t862 * t828;
t765 = -t819 * qJD(3) - t891;
t762 = qJDD(4) - t765;
t930 = -t761 + t762;
t935 = t857 * t930;
t906 = qJD(1) * t859;
t817 = -t862 * t863 * qJD(1) + t858 * t906;
t783 = t819 * t817;
t854 = qJDD(2) + qJDD(3);
t929 = -t783 + t854;
t934 = t858 * t929;
t933 = t861 * t930;
t932 = t862 * t929;
t766 = -qJD(3) * t817 + t827 * t862 + t828 * t858;
t810 = t903 * t817;
t931 = t766 - t810;
t902 = t903 ^ 2;
t856 = t863 ^ 2;
t866 = qJD(1) ^ 2;
t860 = sin(qJ(1));
t864 = cos(qJ(1));
t837 = t860 * g(1) - t864 * g(2);
t873 = qJDD(1) * pkin(1) + t837;
t874 = qJD(2) * pkin(2) - pkin(6) * t906;
t768 = t828 * pkin(2) + (pkin(6) * t856 + pkin(5)) * t866 - t874 * t906 + t873;
t889 = t903 * t819;
t666 = -t931 * pkin(7) + (-t765 + t889) * pkin(3) - t768;
t838 = g(1) * t864 + g(2) * t860;
t869 = -pkin(1) * t866 + qJDD(1) * pkin(5) - t838;
t804 = -t859 * g(3) + t863 * t869;
t852 = t856 * t866;
t760 = -pkin(2) * t852 + t828 * pkin(6) - qJD(2) * t874 + t804;
t868 = t859 * t869;
t909 = t859 * t866;
t867 = -t868 - t827 * pkin(6) + qJDD(2) * pkin(2) + (pkin(2) * t909 + pkin(6) * t905 - g(3)) * t863;
t716 = t862 * t760 + t858 * t867;
t781 = pkin(3) * t817 - pkin(7) * t819;
t674 = -t902 * pkin(3) + t854 * pkin(7) - t817 * t781 + t716;
t627 = -t861 * t666 + t674 * t857;
t628 = t666 * t857 + t674 * t861;
t596 = t857 * t627 + t861 * t628;
t786 = t789 ^ 2;
t787 = t791 ^ 2;
t813 = qJD(4) + t817;
t812 = t813 ^ 2;
t815 = t817 ^ 2;
t816 = t819 ^ 2;
t715 = t760 * t858 - t862 * t867;
t657 = -t715 * t862 + t716 * t858;
t928 = pkin(2) * t657;
t740 = qJD(2) * t819 - t891;
t743 = t766 + t810;
t691 = t740 * t858 - t743 * t862;
t927 = pkin(2) * t691;
t926 = pkin(3) * t858;
t925 = t657 * t859;
t924 = t657 * t863;
t718 = t761 + t762;
t923 = t718 * t857;
t922 = t718 * t861;
t921 = t768 * t858;
t920 = t768 * t862;
t779 = t783 + t854;
t919 = t779 * t858;
t918 = t779 * t862;
t917 = t813 * t857;
t916 = t813 * t861;
t820 = t866 * pkin(5) + t873;
t915 = t820 * t859;
t914 = t820 * t863;
t844 = t863 * t909;
t835 = qJDD(2) + t844;
t913 = t835 * t859;
t836 = qJDD(2) - t844;
t912 = t836 * t859;
t911 = t836 * t863;
t855 = t859 ^ 2;
t910 = t855 * t866;
t673 = -t854 * pkin(3) - t902 * pkin(7) + t781 * t819 + t715;
t670 = t857 * t673;
t671 = t861 * t673;
t908 = -pkin(3) * t673 + pkin(7) * t596;
t907 = t855 + t856;
t901 = t858 * t761;
t900 = t862 * t761;
t899 = t860 * t783;
t898 = t864 * t783;
t750 = -t787 - t812;
t677 = -t750 * t857 - t922;
t878 = -t766 * t861 - t854 * t857;
t712 = (qJD(4) + t813) * t789 + t878;
t897 = pkin(3) * t712 + pkin(7) * t677 + t670;
t737 = -t812 - t786;
t669 = t737 * t861 - t935;
t892 = -t766 * t857 + t861 * t854;
t725 = -qJD(4) * t791 + t892;
t776 = t813 * t791;
t708 = t725 - t776;
t896 = pkin(3) * t708 + pkin(7) * t669 - t671;
t895 = -pkin(3) * t862 - pkin(2);
t658 = t715 * t858 + t862 * t716;
t803 = t863 * g(3) + t868;
t759 = t803 * t859 + t863 * t804;
t890 = -t837 * t860 - t864 * t838;
t888 = t860 * t844;
t887 = t864 * t844;
t709 = (-qJD(4) + t813) * t791 + t892;
t726 = -qJD(4) * t789 - t878;
t775 = t813 * t789;
t711 = t726 + t775;
t656 = t709 * t861 + t711 * t857;
t728 = t786 + t787;
t886 = pkin(3) * t728 + pkin(7) * t656 + t596;
t592 = t596 * t858 - t673 * t862;
t885 = pkin(2) * t592 + t908;
t802 = -t816 - t902;
t744 = t802 * t862 - t919;
t884 = pkin(2) * t744 - t716;
t832 = qJDD(1) * t864 - t860 * t866;
t883 = -pkin(4) * t832 - g(3) * t860;
t882 = t858 * t810;
t881 = t858 * t889;
t880 = t862 * t810;
t879 = t862 * t889;
t595 = -t627 * t861 + t628 * t857;
t758 = t803 * t863 - t804 * t859;
t877 = t837 * t864 - t838 * t860;
t637 = t669 * t858 + t708 * t862;
t876 = pkin(2) * t637 + t896;
t639 = t677 * t858 + t712 * t862;
t875 = pkin(2) * t639 + t897;
t777 = -t902 - t815;
t729 = t777 * t858 + t932;
t872 = pkin(2) * t729 - t715;
t631 = t656 * t858 + t728 * t862;
t871 = pkin(2) * t631 + t886;
t865 = qJD(2) ^ 2;
t842 = -t852 - t865;
t841 = t852 - t865;
t840 = -t865 - t910;
t839 = t865 - t910;
t834 = -t852 + t910;
t833 = t852 + t910;
t831 = qJDD(1) * t860 + t864 * t866;
t830 = t907 * qJDD(1);
t829 = t851 - 0.2e1 * t893;
t826 = 0.2e1 * t848 + t904;
t824 = t863 * t835;
t823 = t907 * t905;
t814 = -pkin(4) * t831 + g(3) * t864;
t808 = -t816 + t902;
t807 = t815 - t902;
t806 = t827 * t863 - t855 * t905;
t805 = -t828 * t859 - t856 * t905;
t801 = -t840 * t859 - t911;
t800 = -t839 * t859 + t824;
t799 = t842 * t863 - t913;
t798 = t841 * t863 - t912;
t797 = t840 * t863 - t912;
t796 = t839 * t863 + t913;
t795 = t842 * t859 + t824;
t794 = t841 * t859 + t911;
t793 = (t827 + t848) * t859;
t792 = (t828 - t893) * t863;
t785 = -t826 * t859 + t829 * t863;
t784 = t826 * t863 + t829 * t859;
t782 = t816 - t815;
t774 = -t787 + t812;
t773 = t786 - t812;
t772 = -t880 + t881;
t771 = -t882 - t879;
t770 = -pkin(5) * t797 - t914;
t769 = -pkin(5) * t795 - t915;
t767 = -t815 - t816;
t764 = -pkin(1) * t797 + t804;
t763 = -pkin(1) * t795 + t803;
t756 = t787 - t786;
t752 = pkin(1) * t829 + pkin(5) * t799 + t914;
t751 = -pkin(1) * t826 + pkin(5) * t801 - t915;
t749 = t807 * t862 - t919;
t748 = -t808 * t858 + t932;
t747 = t807 * t858 + t918;
t746 = t808 * t862 + t934;
t745 = -t802 * t858 - t918;
t738 = (0.2e1 * qJD(3) + qJD(2)) * t819 + t891;
t736 = pkin(1) * t820 + pkin(5) * t759;
t735 = t862 * t766 - t881;
t734 = t858 * t766 + t879;
t733 = -t858 * t765 + t880;
t732 = t862 * t765 + t882;
t731 = pkin(1) * t833 + pkin(5) * t830 + t759;
t730 = t777 * t862 - t934;
t723 = (-t789 * t861 + t791 * t857) * t813;
t722 = (-t789 * t857 - t791 * t861) * t813;
t721 = -t771 * t859 + t772 * t863;
t720 = t771 * t863 + t772 * t859;
t713 = -pkin(6) * t744 - t920;
t710 = t726 - t775;
t707 = -t725 - t776;
t706 = -pkin(6) * t729 - t921;
t703 = -t747 * t859 + t749 * t863;
t702 = -t746 * t859 + t748 * t863;
t701 = t747 * t863 + t749 * t859;
t700 = t746 * t863 + t748 * t859;
t699 = t726 * t861 - t791 * t917;
t698 = t726 * t857 + t791 * t916;
t697 = -t725 * t857 + t789 * t916;
t696 = -t861 * t725 - t789 * t917;
t695 = -t744 * t859 + t745 * t863;
t694 = t744 * t863 + t745 * t859;
t693 = t740 * t862 + t743 * t858;
t692 = -t738 * t862 - t858 * t931;
t690 = -t738 * t858 + t862 * t931;
t689 = t723 * t862 + t762 * t858;
t688 = t723 * t858 - t762 * t862;
t687 = t773 * t861 - t923;
t686 = -t774 * t857 + t933;
t685 = t773 * t857 + t922;
t684 = t774 * t861 + t935;
t683 = -t734 * t859 + t735 * t863;
t682 = -t732 * t859 + t733 * t863;
t681 = t734 * t863 + t735 * t859;
t680 = t732 * t863 + t733 * t859;
t679 = -t729 * t859 + t730 * t863;
t678 = t729 * t863 + t730 * t859;
t676 = t750 * t861 - t923;
t668 = t737 * t857 + t933;
t665 = -pkin(2) * t931 + pkin(6) * t745 - t921;
t663 = t699 * t862 + t901;
t662 = t697 * t862 - t901;
t661 = t699 * t858 - t900;
t660 = t697 * t858 + t900;
t659 = -pkin(2) * t738 + pkin(6) * t730 + t920;
t655 = t708 * t861 - t710 * t857;
t654 = t709 * t857 - t711 * t861;
t653 = t708 * t857 + t710 * t861;
t651 = pkin(2) * t768 + pkin(6) * t658;
t650 = -t691 * t859 + t693 * t863;
t649 = -t690 * t859 + t692 * t863;
t648 = t691 * t863 + t693 * t859;
t647 = t690 * t863 + t692 * t859;
t646 = t687 * t862 - t707 * t858;
t645 = t686 * t862 + t711 * t858;
t644 = t687 * t858 + t707 * t862;
t643 = t686 * t858 - t711 * t862;
t642 = -t859 * t688 + t863 * t689;
t641 = t688 * t863 + t689 * t859;
t640 = t677 * t862 - t712 * t858;
t638 = t669 * t862 - t708 * t858;
t636 = -pkin(1) * t694 - t884;
t635 = t655 * t862 + t756 * t858;
t634 = t655 * t858 - t756 * t862;
t633 = -pkin(7) * t676 + t671;
t632 = t656 * t862 - t728 * t858;
t630 = -pkin(1) * t678 - t872;
t629 = -pkin(7) * t668 + t670;
t624 = -pkin(6) * t691 - t657;
t623 = -t661 * t859 + t663 * t863;
t622 = -t660 * t859 + t662 * t863;
t621 = t661 * t863 + t663 * t859;
t620 = t660 * t863 + t662 * t859;
t619 = -pkin(2) * t767 + pkin(6) * t693 + t658;
t618 = -pkin(1) * t648 - t927;
t617 = -pkin(5) * t694 - t665 * t859 + t713 * t863;
t616 = t658 * t863 - t925;
t615 = t658 * t859 + t924;
t614 = -pkin(5) * t678 - t659 * t859 + t706 * t863;
t613 = -pkin(1) * t931 + pkin(5) * t695 + t665 * t863 + t713 * t859;
t612 = -pkin(3) * t676 + t628;
t611 = -pkin(3) * t668 + t627;
t610 = -pkin(1) * t738 + pkin(5) * t679 + t659 * t863 + t706 * t859;
t609 = -t644 * t859 + t646 * t863;
t608 = -t643 * t859 + t645 * t863;
t607 = t644 * t863 + t646 * t859;
t606 = t643 * t863 + t645 * t859;
t605 = -t639 * t859 + t640 * t863;
t604 = t639 * t863 + t640 * t859;
t603 = -t637 * t859 + t638 * t863;
t602 = t637 * t863 + t638 * t859;
t601 = -t634 * t859 + t635 * t863;
t600 = t634 * t863 + t635 * t859;
t599 = -t631 * t859 + t632 * t863;
t598 = t631 * t863 + t632 * t859;
t597 = -pkin(1) * t615 - t928;
t593 = t596 * t862 + t673 * t858;
t591 = -pkin(7) * t654 - t595;
t590 = -pkin(5) * t615 - pkin(6) * t924 - t651 * t859;
t589 = pkin(1) * t768 + pkin(5) * t616 - pkin(6) * t925 + t651 * t863;
t588 = -pkin(5) * t648 - t619 * t859 + t624 * t863;
t587 = -pkin(1) * t767 + pkin(5) * t650 + t619 * t863 + t624 * t859;
t586 = -pkin(6) * t639 - t612 * t858 + t633 * t862;
t585 = -pkin(6) * t637 - t611 * t858 + t629 * t862;
t584 = -pkin(2) * t676 + pkin(6) * t640 + t612 * t862 + t633 * t858;
t583 = -pkin(1) * t604 - t875;
t582 = -pkin(2) * t668 + pkin(6) * t638 + t611 * t862 + t629 * t858;
t581 = -pkin(1) * t602 - t876;
t580 = -pkin(6) * t631 + t591 * t862 + t654 * t926;
t579 = -t592 * t859 + t593 * t863;
t578 = t592 * t863 + t593 * t859;
t577 = pkin(6) * t632 + t591 * t858 + t654 * t895;
t576 = -pkin(1) * t598 - t871;
t575 = -pkin(6) * t592 + (-pkin(7) * t862 + t926) * t595;
t574 = pkin(6) * t593 + (-pkin(7) * t858 + t895) * t595;
t573 = -pkin(5) * t604 - t584 * t859 + t586 * t863;
t572 = -pkin(5) * t602 - t582 * t859 + t585 * t863;
t571 = -pkin(1) * t676 + pkin(5) * t605 + t584 * t863 + t586 * t859;
t570 = -pkin(1) * t668 + pkin(5) * t603 + t582 * t863 + t585 * t859;
t569 = -pkin(1) * t578 - t885;
t568 = -pkin(5) * t598 - t577 * t859 + t580 * t863;
t567 = -pkin(1) * t654 + pkin(5) * t599 + t577 * t863 + t580 * t859;
t566 = -pkin(5) * t578 - t574 * t859 + t575 * t863;
t565 = -pkin(1) * t595 + pkin(5) * t579 + t574 * t863 + t575 * t859;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t832, 0, -t831, 0, t883, -t814, -t877, -pkin(4) * t877, t806 * t864 - t888, t785 * t864 + t834 * t860, t800 * t864 + t860 * t904, t805 * t864 + t888, t798 * t864 + t851 * t860, qJDD(2) * t860 + t823 * t864, t864 * t769 - t860 * t763 - pkin(4) * (t799 * t860 + t829 * t864), t864 * t770 - t860 * t764 - pkin(4) * (t801 * t860 - t826 * t864), t864 * t758 - pkin(4) * (t830 * t860 + t833 * t864), -pkin(4) * (t759 * t860 + t820 * t864) - (pkin(1) * t860 - pkin(5) * t864) * t758, t683 * t864 + t899, t649 * t864 + t782 * t860, t702 * t864 + t743 * t860, t682 * t864 - t899, t703 * t864 + t740 * t860, t721 * t864 + t854 * t860, t864 * t614 - t860 * t630 - pkin(4) * (t679 * t860 - t738 * t864), t864 * t617 - t860 * t636 - pkin(4) * (t695 * t860 - t864 * t931), t864 * t588 - t860 * t618 - pkin(4) * (t650 * t860 - t767 * t864), t864 * t590 - t860 * t597 - pkin(4) * (t616 * t860 + t768 * t864), t623 * t864 + t698 * t860, t601 * t864 + t653 * t860, t608 * t864 + t684 * t860, t622 * t864 - t696 * t860, t609 * t864 + t685 * t860, t642 * t864 + t722 * t860, t864 * t572 - t860 * t581 - pkin(4) * (t603 * t860 - t668 * t864), t864 * t573 - t860 * t583 - pkin(4) * (t605 * t860 - t676 * t864), t864 * t568 - t860 * t576 - pkin(4) * (t599 * t860 - t654 * t864), t864 * t566 - t860 * t569 - pkin(4) * (t579 * t860 - t595 * t864); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t831, 0, t832, 0, t814, t883, t890, pkin(4) * t890, t806 * t860 + t887, t785 * t860 - t834 * t864, t800 * t860 - t864 * t904, t805 * t860 - t887, t798 * t860 - t851 * t864, -qJDD(2) * t864 + t823 * t860, t860 * t769 + t864 * t763 + pkin(4) * (t799 * t864 - t829 * t860), t860 * t770 + t864 * t764 + pkin(4) * (t801 * t864 + t826 * t860), t860 * t758 + pkin(4) * (t830 * t864 - t833 * t860), pkin(4) * (t759 * t864 - t820 * t860) - (-pkin(1) * t864 - pkin(5) * t860) * t758, t683 * t860 - t898, t649 * t860 - t782 * t864, t702 * t860 - t743 * t864, t682 * t860 + t898, t703 * t860 - t740 * t864, t721 * t860 - t854 * t864, t860 * t614 + t864 * t630 + pkin(4) * (t679 * t864 + t738 * t860), t860 * t617 + t864 * t636 + pkin(4) * (t695 * t864 + t860 * t931), t860 * t588 + t864 * t618 + pkin(4) * (t650 * t864 + t767 * t860), t860 * t590 + t864 * t597 + pkin(4) * (t616 * t864 - t768 * t860), t623 * t860 - t698 * t864, t601 * t860 - t653 * t864, t608 * t860 - t684 * t864, t622 * t860 + t696 * t864, t609 * t860 - t685 * t864, t642 * t860 - t722 * t864, t860 * t572 + t864 * t581 + pkin(4) * (t603 * t864 + t668 * t860), t860 * t573 + t864 * t583 + pkin(4) * (t605 * t864 + t676 * t860), t860 * t568 + t864 * t576 + pkin(4) * (t599 * t864 + t654 * t860), t860 * t566 + t864 * t569 + pkin(4) * (t579 * t864 + t595 * t860); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t837, t838, 0, 0, t793, t784, t796, t792, t794, 0, t752, t751, t731, t736, t681, t647, t700, t680, t701, t720, t610, t613, t587, t589, t621, t600, t606, t620, t607, t641, t570, t571, t567, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t866, 0, 0, -g(3), -t837, 0, t806, t785, t800, t805, t798, t823, t769, t770, t758, pkin(5) * t758, t683, t649, t702, t682, t703, t721, t614, t617, t588, t590, t623, t601, t608, t622, t609, t642, t572, t573, t568, t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t866, 0, qJDD(1), 0, g(3), 0, -t838, 0, t844, -t834, -t904, -t844, -t851, -qJDD(2), t763, t764, 0, pkin(1) * t758, -t783, -t782, -t743, t783, -t740, -t854, t630, t636, t618, t597, -t698, -t653, -t684, t696, -t685, -t722, t581, t583, t576, t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t837, t838, 0, 0, t793, t784, t796, t792, t794, 0, t752, t751, t731, t736, t681, t647, t700, t680, t701, t720, t610, t613, t587, t589, t621, t600, t606, t620, t607, t641, t570, t571, t567, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t827, t829, t835, -t848, t841, t848, 0, -t820, t803, 0, t735, t692, t748, t733, t749, t772, t706, t713, t624, -pkin(6) * t657, t663, t635, t645, t662, t646, t689, t585, t586, t580, t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t893, t826, t839, t828, t836, -t893, t820, 0, t804, 0, t734, t690, t746, t732, t747, t771, t659, t665, t619, t651, t661, t634, t643, t660, t644, t688, t582, t584, t577, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t844, t834, t904, t844, t851, qJDD(2), -t803, -t804, 0, 0, t783, t782, t743, -t783, t740, t854, t872, t884, t927, t928, t698, t653, t684, -t696, t685, t722, t876, t875, t871, t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t766, -t738, t929, t810, t807, -t810, 0, -t768, t715, 0, t699, t655, t686, t697, t687, t723, t629, t633, t591, -pkin(7) * t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t889, t931, t808, t765, t779, -t889, t768, 0, t716, 0, -t761, -t756, -t711, t761, t707, -t762, t611, t612, -pkin(3) * t654, -pkin(3) * t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, t782, t743, -t783, t740, t854, -t715, -t716, 0, 0, t698, t653, t684, -t696, t685, t722, t896, t897, t886, t908; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t726, t708, t930, t775, t773, -t775, 0, t673, t627, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t776, t710, t774, t725, t718, -t776, -t673, 0, t628, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t761, t756, t711, -t761, -t707, t762, -t627, -t628, 0, 0;];
m_new_reg = t1;
