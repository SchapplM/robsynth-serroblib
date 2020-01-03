% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPPR9_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:00
% EndTime: 2019-12-31 19:42:08
% DurationCPUTime: 8.53s
% Computational Cost: add. (17472->574), mult. (37632->589), div. (0->0), fcn. (19042->6), ass. (0->335)
t1002 = qJD(2) ^ 2;
t862 = sin(qJ(2));
t855 = t862 ^ 2;
t867 = qJD(1) ^ 2;
t951 = t855 * t867;
t829 = t951 + t1002;
t865 = cos(qJ(2));
t1020 = t865 * t867;
t834 = t862 * t1020;
t824 = -t834 + qJDD(2);
t946 = t865 * t824;
t774 = -t829 * t862 + t946;
t929 = qJD(1) * qJD(2);
t843 = t865 * t929;
t926 = t862 * qJDD(1);
t811 = 0.2e1 * t843 + t926;
t863 = sin(qJ(1));
t866 = cos(qJ(1));
t728 = pkin(5) * (t774 * t866 - t811 * t863);
t993 = pkin(5) * (t774 * t863 + t811 * t866);
t823 = t834 + qJDD(2);
t801 = t865 * t823;
t856 = t865 ^ 2;
t950 = t856 * t867;
t831 = t950 + t1002;
t771 = -t831 * t862 + t801;
t998 = pkin(1) * t771;
t989 = pkin(6) * t771;
t959 = t823 * t862;
t775 = t831 * t865 + t959;
t842 = t862 * t929;
t925 = t865 * qJDD(1);
t814 = -0.2e1 * t842 + t925;
t729 = pkin(5) * (t775 * t866 + t814 * t863);
t992 = pkin(5) * (t775 * t863 - t814 * t866);
t944 = pkin(1) * t811 + pkin(6) * t774;
t958 = t824 * t862;
t766 = t829 * t865 + t958;
t1032 = pkin(1) * t766;
t1031 = pkin(3) * t823;
t990 = pkin(6) * t766;
t861 = sin(qJ(5));
t864 = cos(qJ(5));
t938 = qJD(1) * t865;
t804 = -qJD(2) * t864 + t861 * t938;
t805 = qJD(2) * t861 + t864 * t938;
t759 = t804 * t805;
t812 = t843 + t926;
t796 = qJDD(5) + t812;
t1016 = -t759 + t796;
t1030 = t1016 * t861;
t1029 = t1016 * t864;
t1015 = pkin(2) * t829 + qJ(3) * t824;
t927 = qJDD(2) * qJ(3);
t852 = t862 * g(3);
t940 = t1002 * pkin(2) + t852;
t886 = -0.2e1 * qJD(3) * qJD(2) - t927 + t940;
t826 = g(1) * t866 + g(2) * t863;
t793 = -pkin(1) * t867 + qJDD(1) * pkin(6) - t826;
t899 = -pkin(2) * t865 - qJ(3) * t862;
t809 = t899 * qJD(1);
t907 = qJD(1) * t809 + t793;
t718 = t865 * t907 - t886;
t1028 = t718 + t1015;
t988 = t865 * g(3);
t891 = -qJDD(2) * pkin(2) - t1002 * qJ(3) + qJDD(3) + t988;
t721 = t862 * t907 + t891;
t983 = qJ(3) * t831;
t1027 = pkin(2) * t823 - t721 - t983;
t943 = pkin(1) * t814 - pkin(6) * t775;
t963 = t814 * t865;
t969 = t811 * t862;
t754 = -t963 + t969;
t821 = (t855 - t856) * t867;
t1026 = t754 * t863 + t821 * t866;
t1025 = t754 * t866 - t821 * t863;
t832 = t950 - t1002;
t777 = -t832 * t865 + t958;
t1022 = t777 * t863 + t866 * t925;
t1021 = t777 * t866 - t863 * t925;
t939 = t855 + t856;
t817 = t939 * qJDD(1);
t820 = t939 * t867;
t755 = pkin(5) * (t817 * t866 - t820 * t863);
t813 = -t842 + t925;
t1011 = t813 * pkin(3) - qJ(4) * t950 + qJDD(4);
t825 = g(1) * t863 - t866 * g(2);
t792 = qJDD(1) * pkin(1) + pkin(6) * t867 + t825;
t878 = -pkin(2) * t842 + t792;
t875 = t878 + t1011;
t1001 = 0.2e1 * qJD(3);
t935 = t862 * qJD(1);
t822 = -qJD(2) * pkin(3) - qJ(4) * t935;
t934 = t1001 + t822;
t908 = t934 * t862;
t987 = pkin(4) + qJ(3);
t641 = (pkin(2) + pkin(7)) * t813 + t987 * t812 + (t908 + (-pkin(7) * t862 + t865 * t987) * qJD(2)) * qJD(1) + t875;
t948 = t862 * t793;
t876 = -t812 * qJ(4) - t1031 + t891 + t948;
t903 = pkin(4) * t862 + pkin(7) * t865;
t936 = qJD(2) * t865;
t917 = qJ(4) * t936;
t933 = -(2 * qJD(4)) + t809;
t666 = -t1002 * pkin(4) - qJDD(2) * pkin(7) + (t917 + (-qJD(1) * t903 + t933) * t862) * qJD(1) + t876;
t616 = -t864 * t641 + t861 * t666;
t617 = t641 * t861 + t666 * t864;
t604 = t861 * t616 + t617 * t864;
t1018 = -t864 * t616 + t617 * t861;
t748 = qJD(5) * t804 - t861 * qJDD(2) - t813 * t864;
t836 = qJD(5) + t935;
t788 = t836 * t804;
t1017 = t748 + t788;
t942 = pkin(1) * t820 + pkin(6) * t817;
t1014 = pkin(3) * t950 + qJ(4) * t813;
t769 = t832 * t862 + t946;
t1013 = pkin(3) * t814 + qJ(4) * t831;
t898 = t812 + t843;
t906 = t935 * t1001;
t871 = qJ(3) * t898 + t878 + t906;
t687 = (t813 + t814) * pkin(2) + t871;
t892 = qJDD(2) * t864 - t813 * t861;
t712 = (qJD(5) - t836) * t805 - t892;
t1000 = pkin(2) + pkin(3);
t877 = qJD(2) * t934 - t1014 - t940;
t893 = qJD(1) * t933 + t793;
t873 = t865 * t893 + t877;
t665 = -pkin(7) * t1002 + qJDD(2) * t987 - t903 * t1020 + t873;
t945 = pkin(4) * t665 - pkin(7) * t604;
t1010 = qJ(3) * t665 - t1000 * t604 + t945;
t714 = t748 - t788;
t657 = t712 * t864 + t861 * t714;
t1003 = t805 ^ 2;
t795 = t804 ^ 2;
t735 = -t795 - t1003;
t887 = -pkin(4) * t735 + pkin(7) * t657 + t604;
t1009 = qJ(3) * t735 - t1000 * t657 - t887;
t677 = t873 + t927;
t678 = (t862 * t933 + t917) * qJD(1) + t876;
t1008 = qJ(3) * t677 - t1000 * t678;
t833 = t836 ^ 2;
t749 = -t833 - t795;
t682 = t749 * t864 - t1030;
t710 = (-qJD(5) - t836) * t805 + t892;
t947 = t864 * t665;
t919 = pkin(4) * t710 - pkin(7) * t682 + t947;
t1007 = qJ(3) * t710 - t1000 * t682 + t919;
t916 = -t833 - t1003;
t737 = t759 + t796;
t972 = t737 * t864;
t689 = -t861 * t916 - t972;
t660 = t861 * t665;
t920 = -pkin(4) * t1017 + pkin(7) * t689 + t660;
t1006 = qJ(3) * t1017 - t1000 * t689 - t920;
t1005 = -t1000 * t823 + t983;
t923 = 0.2e1 * qJD(1) * qJD(4);
t1004 = -qJD(2) * t822 + t865 * t923 + t1014;
t999 = -pkin(3) - pkin(7);
t997 = pkin(2) * t813;
t994 = pkin(3) * t829;
t991 = pkin(5) * (t817 * t863 + t820 * t866);
t985 = qJ(3) * t814;
t984 = qJ(3) * t820;
t982 = qJ(4) * t665;
t981 = qJ(4) * t677;
t980 = qJ(4) * t678;
t978 = qJ(4) * t823;
t977 = qJ(4) * t824;
t976 = qJ(4) * t829;
t971 = t792 * t862;
t970 = t792 * t865;
t967 = t811 * t865;
t965 = t814 * t862;
t954 = t836 * t805;
t953 = t836 * t861;
t952 = t836 * t864;
t949 = t861 * t737;
t931 = pkin(2) - t999;
t930 = qJ(4) * qJDD(1);
t924 = t833 - t1003;
t922 = t862 * t759;
t921 = t865 * t759;
t681 = t861 * t749 + t1029;
t918 = -pkin(4) * t681 + t616;
t602 = pkin(4) * t1018;
t913 = -qJ(4) * t604 + t602;
t655 = t712 * t861 - t714 * t864;
t654 = pkin(4) * t655;
t912 = -qJ(4) * t657 + t654;
t911 = -qJ(4) * t1017 - t947;
t780 = t948 + t988;
t781 = t793 * t865 - t852;
t717 = t780 * t862 + t865 * t781;
t909 = -t825 * t863 - t866 * t826;
t905 = t863 * t834;
t904 = t866 * t834;
t819 = qJDD(1) * t866 - t863 * t867;
t902 = -pkin(5) * t819 - g(3) * t863;
t900 = -pkin(2) * t721 + qJ(3) * t718;
t807 = -pkin(2) * t926 + qJ(3) * t925;
t896 = -qJ(4) * t710 - t660;
t716 = t780 * t865 - t781 * t862;
t751 = t965 + t967;
t894 = t825 * t866 - t826 * t863;
t888 = -qJ(4) * t682 - t918;
t885 = -qJ(4) * t735 + t1018;
t688 = t864 * t916 - t949;
t884 = -pkin(4) * t688 + t617;
t881 = qJD(5) * t805 - t892;
t879 = -qJ(4) * t689 - t884;
t671 = -t1028 - t1032;
t872 = t997 + (t811 + t898) * qJ(3) + t878;
t870 = t862 * t923 + (-t862 * t809 - t917) * qJD(1) - t876;
t869 = qJ(4) * t926 + t870;
t868 = t677 + t994;
t670 = t997 + qJ(3) * t812 + (qJ(3) * t936 + t908) * qJD(1) + t875;
t846 = pkin(3) * t926;
t830 = t951 - t1002;
t818 = qJDD(1) * t863 + t866 * t867;
t800 = t939 * t929;
t790 = -pkin(5) * t818 + g(3) * t866;
t789 = t807 - t846;
t786 = t795 - t833;
t785 = qJDD(2) * t863 + t800 * t866;
t784 = t812 * t865 - t855 * t929;
t783 = -qJDD(2) * t866 + t800 * t863;
t782 = -t813 * t862 - t856 * t929;
t773 = t830 * t862 + t801;
t767 = t898 * t862;
t765 = -t830 * t865 + t959;
t764 = (t813 - t842) * t865;
t762 = -t978 - t985;
t757 = -t795 + t1003;
t747 = t784 * t866 - t905;
t746 = t782 * t866 + t905;
t745 = t784 * t863 + t904;
t744 = t782 * t863 - t904;
t742 = t773 * t866 + t863 * t926;
t740 = t773 * t863 - t866 * t926;
t738 = t1000 * t811 - t977;
t727 = -t970 + t990;
t726 = -t971 - t989;
t725 = (t804 * t864 - t805 * t861) * t836;
t724 = (t804 * t861 + t805 * t864) * t836;
t720 = t781 + t1032;
t719 = t780 - t998;
t707 = t943 + t970;
t706 = -t944 - t971;
t703 = t721 + t984;
t702 = pkin(2) * t820 + t718;
t701 = t748 * t864 + t805 * t953;
t700 = t748 * t861 - t805 * t952;
t699 = t804 * t952 + t861 * t881;
t698 = t804 * t953 - t864 * t881;
t697 = t871 + t997;
t696 = t725 * t862 + t796 * t865;
t695 = -t725 * t865 + t796 * t862;
t694 = t786 * t864 - t949;
t693 = -t861 * t924 + t1029;
t692 = t786 * t861 + t972;
t691 = t864 * t924 + t1030;
t690 = pkin(1) * t792 + pkin(6) * t717;
t686 = t872 + t906;
t683 = t717 + t942;
t676 = -t998 - t1027;
t675 = t701 * t862 + t921;
t674 = -t699 * t862 - t921;
t673 = -t701 * t865 + t922;
t672 = t699 * t865 - t922;
t669 = t869 - t984;
t668 = t718 * t865 + t721 * t862;
t667 = t718 * t862 - t721 * t865;
t659 = -t1000 * t820 + (-t907 + t930) * t865 + t886 + t1004;
t658 = -t1017 * t861 - t710 * t864;
t656 = t1017 * t864 - t861 * t710;
t652 = t908 * qJD(1) + t1011 + t872 + t976;
t651 = -pkin(2) * t969 + t686 * t865 - t990;
t650 = qJ(3) * t963 - t687 * t862 - t989;
t649 = -t702 * t862 + t703 * t865;
t648 = t693 * t862 + t714 * t865;
t647 = t694 * t862 + t712 * t865;
t646 = -t693 * t865 + t714 * t862;
t645 = -t694 * t865 + t712 * t862;
t644 = -t822 * t935 - t1011 - t1013 - t687;
t643 = pkin(2) * t967 + t686 * t862 + t944;
t642 = qJ(3) * t965 + t687 * t865 + t943;
t640 = t1017 * t865 + t689 * t862;
t639 = t1017 * t862 - t689 * t865;
t637 = t671 - t994 + t1004;
t636 = t682 * t862 + t710 * t865;
t635 = -t682 * t865 + t710 * t862;
t634 = -t1005 + t870 + t998;
t633 = t702 * t865 + t703 * t862 + t942;
t632 = t658 * t862 + t757 * t865;
t631 = -t658 * t865 + t757 * t862;
t630 = t657 * t862 + t735 * t865;
t629 = -t657 * t865 + t735 * t862;
t628 = t677 * t865 + t678 * t862;
t627 = t677 * t862 - t678 * t865;
t626 = qJ(3) * t670 - t980;
t625 = -t644 * t862 + t762 * t865 + t989;
t624 = t652 * t865 - t738 * t862 - t990;
t623 = t644 * t865 + t762 * t862 - t943;
t622 = t652 * t862 + t738 * t865 + t944;
t621 = -pkin(1) * t667 - t900;
t620 = -t659 * t862 + t669 * t865;
t619 = -pkin(6) * t667 + (-pkin(2) * t862 + qJ(3) * t865) * t697;
t618 = t659 * t865 + t669 * t862 - t942;
t612 = t1000 * t670 - t981;
t611 = pkin(6) * t668 + (pkin(1) - t899) * t697;
t610 = qJ(3) * t655 + t912;
t609 = t688 * t931 + t911;
t608 = t681 * t931 + t896;
t607 = -pkin(1) * t627 - t1008;
t606 = qJ(3) * t688 + t879;
t605 = qJ(3) * t681 + t888;
t600 = -pkin(1) * t639 - t1006;
t599 = -pkin(1) * t635 - t1007;
t598 = t604 * t862 + t665 * t865;
t597 = -t604 * t865 + t665 * t862;
t596 = -pkin(6) * t627 - t612 * t862 + t626 * t865;
t595 = pkin(1) * t670 + pkin(6) * t628 + t612 * t865 + t626 * t862;
t594 = t655 * t931 + t885;
t593 = -pkin(6) * t639 + t606 * t865 - t609 * t862;
t592 = -pkin(6) * t635 + t605 * t865 - t608 * t862;
t591 = pkin(1) * t688 + pkin(6) * t640 + t606 * t862 + t609 * t865;
t590 = -pkin(1) * t629 - t1009;
t589 = pkin(1) * t681 + pkin(6) * t636 + t605 * t862 + t608 * t865;
t588 = -pkin(6) * t629 - t594 * t862 + t610 * t865;
t587 = qJ(3) * t1018 + t913;
t586 = pkin(1) * t655 + pkin(6) * t630 + t594 * t865 + t610 * t862;
t585 = t1018 * t931 - t982;
t584 = -pkin(1) * t597 - t1010;
t583 = -pkin(6) * t597 - t585 * t862 + t587 * t865;
t582 = pkin(1) * t1018 + pkin(6) * t598 + t585 * t865 + t587 * t862;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t819, 0, -t818, 0, t902, -t790, -t894, -pkin(5) * t894, t747, -t1025, t742, t746, -t1021, t785, -t719 * t863 + t726 * t866 + t992, -t863 * t720 + t866 * t727 + t993, t716 * t866 - t991, -pkin(5) * (t717 * t863 + t792 * t866) - (pkin(1) * t863 - pkin(6) * t866) * t716, t747, t742, t1025, t785, t1021, t746, t650 * t866 - t676 * t863 + t992, t649 * t866 + t807 * t863 - t991, t651 * t866 - t671 * t863 - t993, t866 * t619 - t863 * t621 - pkin(5) * (t668 * t863 + t697 * t866), t746, -t1025, -t1021, t747, t742, t785, t624 * t866 - t637 * t863 - t993, t866 * t625 - t863 * t634 - t992, t866 * t620 - t863 * t789 + t991, t866 * t596 - t863 * t607 - pkin(5) * (t628 * t863 + t670 * t866), t675 * t866 - t700 * t863, t632 * t866 - t656 * t863, t648 * t866 - t691 * t863, t674 * t866 + t698 * t863, t647 * t866 - t692 * t863, t696 * t866 - t724 * t863, t866 * t592 - t863 * t599 - pkin(5) * (t636 * t863 + t681 * t866), t866 * t593 - t863 * t600 - pkin(5) * (t640 * t863 + t688 * t866), t866 * t588 - t863 * t590 - pkin(5) * (t630 * t863 + t655 * t866), t866 * t583 - t863 * t584 - pkin(5) * (t1018 * t866 + t598 * t863); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t818, 0, t819, 0, t790, t902, t909, pkin(5) * t909, t745, -t1026, t740, t744, -t1022, t783, t719 * t866 + t726 * t863 - t729, t866 * t720 + t863 * t727 - t728, t716 * t863 + t755, pkin(5) * (t717 * t866 - t792 * t863) - (-pkin(1) * t866 - pkin(6) * t863) * t716, t745, t740, t1026, t783, t1022, t744, t650 * t863 + t676 * t866 - t729, t649 * t863 - t807 * t866 + t755, t651 * t863 + t671 * t866 + t728, t863 * t619 + t866 * t621 + pkin(5) * (t668 * t866 - t697 * t863), t744, -t1026, -t1022, t745, t740, t783, t624 * t863 + t637 * t866 + t728, t863 * t625 + t866 * t634 + t729, t863 * t620 + t866 * t789 - t755, t863 * t596 + t866 * t607 + pkin(5) * (t628 * t866 - t670 * t863), t675 * t863 + t700 * t866, t632 * t863 + t656 * t866, t648 * t863 + t691 * t866, t674 * t863 - t698 * t866, t647 * t863 + t692 * t866, t696 * t863 + t724 * t866, t863 * t592 + t866 * t599 + pkin(5) * (t636 * t866 - t681 * t863), t863 * t593 + t866 * t600 + pkin(5) * (t640 * t866 - t688 * t863), t863 * t588 + t866 * t590 + pkin(5) * (t630 * t866 - t655 * t863), t863 * t583 + t866 * t584 + pkin(5) * (-t1018 * t863 + t598 * t866); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t825, t826, 0, 0, t767, t751, t765, t764, t769, 0, t707, t706, t683, t690, t767, t765, -t751, 0, -t769, t764, t642, t633, t643, t611, t764, t751, t769, t767, t765, 0, t622, t623, t618, t595, t673, t631, t646, t672, t645, t695, t589, t591, t586, t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t867, 0, 0, -g(3), -t825, 0, t784, -t754, t773, t782, -t777, t800, t726, t727, t716, pkin(6) * t716, t784, t773, t754, t800, t777, t782, t650, t649, t651, t619, t782, -t754, -t777, t784, t773, t800, t624, t625, t620, t596, t675, t632, t648, t674, t647, t696, t592, t593, t588, t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t867, 0, qJDD(1), 0, g(3), 0, -t826, 0, t834, -t821, -t926, -t834, -t925, -qJDD(2), t719, t720, 0, pkin(1) * t716, t834, -t926, t821, -qJDD(2), t925, -t834, t676, -t807, t671, t621, -t834, -t821, -t925, t834, -t926, -qJDD(2), t637, t634, t789, t607, t700, t656, t691, -t698, t692, t724, t599, t600, t590, t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t825, t826, 0, 0, t767, t751, t765, t764, t769, 0, t707, t706, t683, t690, t767, t765, -t751, 0, -t769, t764, t642, t633, t643, t611, t764, t751, t769, t767, t765, 0, t622, t623, t618, t595, t673, t631, t646, t672, t645, t695, t589, t591, t586, t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, t814, t823, -t843, t832, t843, 0, -t792, t780, 0, t812, t823, -t814, t843, -t832, -t843, t985, t703, t686, qJ(3) * t697, -t843, t814, t832, t812, t823, t843, t652, t762, t669, t626, t759, t757, t714, -t759, t712, t796, t605, t606, t610, t587; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t842, t811, -t830, t813, t824, -t842, t792, 0, t781, 0, t842, -t830, -t811, -t842, -t824, t813, t687, t702, pkin(2) * t811, pkin(2) * t697, t813, t811, t824, t842, -t830, -t842, t738, t644, t659, t612, -t701, -t658, -t693, t699, -t694, -t725, t608, t609, t594, t585; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t834, t821, t926, t834, t925, qJDD(2), -t780, -t781, 0, 0, -t834, t926, -t821, qJDD(2), -t925, t834, t1027, t807, t1028, t900, t834, t821, t925, -t834, t926, qJDD(2), t868 + t1015, t1005 + t678, -t789, t1008, -t700, -t656, -t691, t698, -t692, -t724, t1007, t1006, t1009, t1010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, t823, -t814, t843, -t832, -t843, 0, t721, t697, 0, -t843, t814, t832, t812, t823, t843, t670 + t976, -t978, t869, -t980, t759, t757, t714, -t759, t712, t796, t888, t879, t912, t913; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t834, t926, -t821, qJDD(2), -t925, t834, -t721, 0, t718, 0, t834, t821, t925, -t834, t926, qJDD(2), t868, t678 - t1031, t846, -pkin(3) * t678, -t700, -t656, -t691, t698, -t692, -t724, -pkin(3) * t682 + t919, -pkin(3) * t689 - t920, -pkin(3) * t657 - t887, -pkin(3) * t604 + t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t842, t830, t811, t842, t824, -t813, -t697, -t718, 0, 0, -t813, -t811, -t824, -t842, t830, t842, -pkin(3) * t811 + t977, t670 + t1013, pkin(3) * t820 + t927 + (t893 - t930) * t865 + t877, -pkin(3) * t670 + t981, t701, t658, t693, -t699, t694, t725, t681 * t999 - t896, t688 * t999 - t911, t655 * t999 - t885, t1018 * t999 + t982; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t813, -t811, -t824, -t842, t830, t842, 0, t670, t677, 0, t701, t658, t693, -t699, t694, t725, -pkin(7) * t681 + t660, -pkin(7) * t688 + t947, -pkin(7) * t655 - t1018, -pkin(7) * t1018; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t843, -t814, -t832, -t812, -t823, -t843, -t670, 0, t678, 0, -t759, -t757, -t714, t759, -t712, -t796, t918, t884, -t654, -t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t834, -t821, -t925, t834, -t926, -qJDD(2), -t677, -t678, 0, 0, t700, t656, t691, -t698, t692, t724, -t919, t920, t887, -t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t748, -t710, t1016, -t788, t786, t788, 0, t665, t616, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t954, t1017, t924, t881, t737, t954, -t665, 0, t617, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t759, t757, t714, -t759, t712, t796, -t616, -t617, 0, 0;];
m_new_reg = t1;
