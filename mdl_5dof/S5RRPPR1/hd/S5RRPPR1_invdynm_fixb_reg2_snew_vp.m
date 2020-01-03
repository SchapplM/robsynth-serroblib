% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:13
% EndTime: 2020-01-03 11:56:24
% DurationCPUTime: 11.39s
% Computational Cost: add. (62787->490), mult. (88678->690), div. (0->0), fcn. (55375->10), ass. (0->336)
t936 = qJD(1) + qJD(2);
t933 = t936 ^ 2;
t934 = qJDD(1) + qJDD(2);
t939 = sin(pkin(8));
t941 = cos(pkin(8));
t900 = t933 * t939 - t934 * t941;
t937 = g(1) - qJDD(3);
t1026 = qJ(3) * t900 - t937 * t939;
t897 = t933 * t941 + t934 * t939;
t943 = sin(qJ(2));
t946 = cos(qJ(2));
t845 = t897 * t943 + t900 * t946;
t866 = qJ(3) * t897 - t937 * t941;
t1041 = pkin(6) * t845 + t1026 * t946 + t866 * t943;
t841 = t897 * t946 - t900 * t943;
t771 = pkin(6) * t841 - t1026 * t943 + t866 * t946;
t944 = sin(qJ(1));
t947 = cos(qJ(1));
t787 = t841 * t947 - t845 * t944;
t1055 = -pkin(5) * t787 + t1041 * t944 - t947 * t771;
t1025 = t841 * t944 + t845 * t947;
t1054 = pkin(5) * t1025 + t1041 * t947 + t944 * t771;
t917 = g(2) * t947 + g(3) * t944;
t956 = qJDD(1) * pkin(1) - t917;
t1014 = qJD(1) ^ 2;
t916 = g(2) * t944 - g(3) * t947;
t957 = t1014 * pkin(1) + t916;
t860 = t943 * t956 - t946 * t957;
t847 = -t933 * pkin(2) + t860;
t859 = -t943 * t957 - t946 * t956;
t953 = t934 * pkin(2) - t859;
t788 = t847 * t939 - t941 * t953;
t789 = t941 * t847 + t939 * t953;
t736 = t788 * t941 - t789 * t939;
t1006 = t736 * t943;
t979 = t788 * t939 + t941 * t789;
t1036 = t946 * t979 + t1006;
t1005 = t736 * t946;
t674 = -t943 * t979 + t1005;
t1053 = t1036 * t947 + t674 * t944;
t1052 = t1036 * t944 - t674 * t947;
t906 = t933 * t943 - t934 * t946;
t1027 = pkin(6) * t906 - g(1) * t943;
t903 = t933 * t946 + t934 * t943;
t849 = t903 * t947 - t906 * t944;
t876 = pkin(6) * t903 - g(1) * t946;
t1045 = -pkin(5) * t849 + t1027 * t944 - t947 * t876;
t964 = t903 * t944 + t906 * t947;
t1044 = pkin(5) * t964 + t1027 * t947 + t944 * t876;
t938 = sin(pkin(9));
t935 = t938 ^ 2;
t940 = cos(pkin(9));
t950 = t940 ^ 2;
t1015 = t933 * (t935 + t950);
t888 = t940 * t1015;
t990 = t940 * t941;
t855 = -t888 * t939 + t934 * t990;
t925 = t940 * t934;
t857 = t888 * t941 + t939 * t925;
t799 = t855 * t946 - t857 * t943;
t801 = t855 * t943 + t857 * t946;
t1043 = t799 * t947 - t801 * t944;
t1042 = t799 * t944 + t801 * t947;
t810 = t859 * t946 - t860 * t943;
t1001 = t810 * t947;
t978 = t859 * t943 + t946 * t860;
t1038 = t944 * t978 - t1001;
t1002 = t810 * t944;
t1037 = t947 * t978 + t1002;
t942 = sin(qJ(5));
t945 = cos(qJ(5));
t1031 = t938 * t942 - t940 * t945;
t881 = t1031 * t936;
t962 = t938 * t945 + t940 * t942;
t883 = t962 * t936;
t834 = t883 * t881;
t1016 = qJDD(5) - t834;
t1030 = t1016 * t942;
t1029 = t1016 * t945;
t975 = -t933 * pkin(3) + t934 * qJ(4) + 0.2e1 * qJD(4) * t936 + t789;
t1028 = pkin(7) * t934 + t975;
t993 = t938 * t940;
t861 = t897 * t993;
t912 = t933 * t993;
t983 = t938 * t925;
t862 = -t939 * t912 + t941 * t983;
t815 = t861 * t946 + t862 * t943;
t818 = t861 * t943 - t862 * t946;
t1024 = t815 * t947 - t818 * t944;
t1023 = t815 * t944 + t818 * t947;
t927 = t940 * t937;
t762 = t975 * t938 + t927;
t995 = t938 * t937;
t763 = t975 * t940 - t995;
t708 = t938 * t762 + t940 * t763;
t926 = t950 * t933;
t998 = t935 * t933;
t901 = t926 + t998;
t879 = t881 ^ 2;
t880 = t883 ^ 2;
t754 = -t927 + (pkin(4) * t933 * t940 - t1028) * t938;
t755 = -pkin(4) * t926 + t1028 * t940 - t995;
t684 = -t945 * t754 + t755 * t942;
t685 = t754 * t942 + t755 * t945;
t651 = -t684 * t945 + t685 * t942;
t1013 = pkin(4) * t651;
t877 = t1031 * t934;
t878 = t962 * t934;
t776 = -t877 * t942 - t878 * t945;
t1012 = pkin(4) * t776;
t1008 = t651 * t938;
t1007 = t651 * t940;
t775 = -t934 * pkin(3) - t933 * qJ(4) + qJDD(4) + t788;
t759 = -pkin(4) * t925 - t901 * pkin(7) + t775;
t1004 = t759 * t942;
t1003 = t759 * t945;
t827 = qJDD(5) + t834;
t1000 = t827 * t942;
t999 = t827 * t945;
t766 = t938 * t775;
t923 = t938 * t934;
t994 = t938 * t939;
t992 = t938 * t941;
t767 = t940 * t775;
t988 = -pkin(3) * t775 + qJ(4) * t708;
t872 = t881 * qJD(5);
t986 = t883 * qJD(5);
t985 = t939 * t834;
t984 = t941 * t834;
t982 = pkin(3) * t925 - qJ(4) * t888 - t767;
t677 = t708 * t939 - t775 * t941;
t981 = pkin(2) * t677 + t988;
t910 = -t944 * qJDD(1) - t947 * t1014;
t980 = pkin(5) * t910 + t947 * g(1);
t652 = t684 * t942 + t945 * t685;
t976 = -t944 * t916 - t917 * t947;
t778 = -t877 * t945 + t878 * t942;
t814 = -t879 - t880;
t635 = -pkin(4) * t814 + pkin(7) * t778 + t652;
t644 = -pkin(7) * t776 - t651;
t726 = -t776 * t938 + t778 * t940;
t974 = -pkin(3) * t814 + qJ(4) * t726 + t940 * t635 + t938 * t644;
t948 = qJD(5) ^ 2;
t825 = -t948 - t879;
t774 = t825 * t945 - t1030;
t829 = t877 + 0.2e1 * t986;
t689 = -pkin(4) * t829 + pkin(7) * t774 - t1003;
t773 = t825 * t942 + t1029;
t710 = -pkin(7) * t773 + t1004;
t723 = -t773 * t938 + t774 * t940;
t973 = -pkin(3) * t829 + qJ(4) * t723 + t940 * t689 + t938 * t710;
t869 = -t880 - t948;
t796 = -t869 * t942 - t999;
t831 = -0.2e1 * t872 + t878;
t702 = -pkin(4) * t831 + pkin(7) * t796 + t1004;
t793 = t869 * t945 - t1000;
t720 = -pkin(7) * t793 + t1003;
t744 = -t793 * t938 + t796 * t940;
t972 = -pkin(3) * t831 + qJ(4) * t744 + t940 * t702 + t938 * t720;
t922 = t935 * t934;
t924 = t950 * t934;
t895 = t924 + t922;
t971 = pkin(3) * t901 + qJ(4) * t895 + t708;
t970 = pkin(2) * t855 + t982;
t969 = -pkin(2) * t900 - t788;
t687 = t726 * t939 - t814 * t941;
t968 = pkin(2) * t687 + t974;
t695 = t723 * t939 - t829 * t941;
t967 = pkin(2) * t695 + t973;
t712 = t744 * t939 - t831 * t941;
t966 = pkin(2) * t712 + t972;
t837 = t895 * t939 + t901 * t941;
t965 = pkin(2) * t837 + t971;
t707 = t762 * t940 - t763 * t938;
t963 = t916 * t947 - t917 * t944;
t887 = t938 * t1015;
t961 = -pkin(3) * t923 + qJ(4) * t887 + t766;
t960 = pkin(4) * t773 - t684;
t853 = t887 * t939 - t934 * t992;
t959 = pkin(2) * t853 + t961;
t629 = t652 * t940 - t1008;
t640 = -pkin(4) * t759 + pkin(7) * t652;
t958 = -pkin(3) * t759 - pkin(7) * t1008 + qJ(4) * t629 + t940 * t640;
t625 = t629 * t939 - t759 * t941;
t955 = pkin(2) * t625 + t958;
t954 = pkin(4) * t793 - t685;
t952 = -pkin(2) * t897 - t789;
t911 = qJDD(1) * t947 - t1014 * t944;
t909 = 0.2e1 * t983;
t902 = -t926 + t998;
t896 = t924 - t922;
t893 = pkin(5) * t911 + g(1) * t944;
t868 = -t880 + t948;
t867 = t879 - t948;
t856 = t887 * t941 + t934 * t994;
t840 = t896 * t941 + t902 * t939;
t839 = t895 * t941 - t901 * t939;
t838 = t896 * t939 - t902 * t941;
t833 = t880 - t879;
t832 = -t872 + t878;
t830 = -t877 - t986;
t822 = -pkin(1) * t906 - t859;
t821 = -pkin(1) * t903 - t860;
t820 = (-t881 * t945 + t883 * t942) * qJD(5);
t819 = (-t881 * t942 - t883 * t945) * qJD(5);
t807 = pkin(1) * t810;
t806 = t832 * t945 - t942 * t986;
t805 = t832 * t942 + t945 * t986;
t804 = -t830 * t942 + t945 * t872;
t803 = t830 * t945 + t942 * t872;
t800 = -t853 * t943 + t856 * t946;
t797 = t853 * t946 + t856 * t943;
t795 = -t868 * t942 + t1029;
t794 = t867 * t945 - t1000;
t792 = t868 * t945 + t1030;
t791 = t867 * t942 + t999;
t790 = pkin(1) * g(1) + pkin(6) * t978;
t785 = -t838 * t943 + t840 * t946;
t784 = -t837 * t943 + t839 * t946;
t783 = t838 * t946 + t840 * t943;
t782 = t837 * t946 + t839 * t943;
t779 = -t829 * t945 - t831 * t942;
t777 = -t829 * t942 + t831 * t945;
t765 = -t819 * t938 + t820 * t940;
t764 = t819 * t940 + t820 * t938;
t758 = qJDD(5) * t939 + t765 * t941;
t757 = -qJDD(5) * t941 + t765 * t939;
t752 = -pkin(1) * t845 + t969;
t751 = -pkin(1) * t841 + t952;
t750 = -t805 * t938 + t806 * t940;
t749 = -t803 * t938 + t804 * t940;
t748 = t805 * t940 + t806 * t938;
t747 = t803 * t940 + t804 * t938;
t746 = t797 * t944 - t800 * t947;
t745 = t797 * t947 + t800 * t944;
t743 = -t792 * t938 + t795 * t940;
t742 = -t791 * t938 + t794 * t940;
t741 = t793 * t940 + t796 * t938;
t740 = t792 * t940 + t795 * t938;
t739 = t791 * t940 + t794 * t938;
t733 = pkin(2) * t736;
t732 = pkin(2) * t937 + qJ(3) * t979;
t731 = t743 * t941 + t878 * t939;
t730 = t742 * t941 - t877 * t939;
t729 = t743 * t939 - t878 * t941;
t728 = t742 * t939 + t877 * t941;
t727 = -t777 * t938 + t779 * t940;
t725 = t777 * t940 + t779 * t938;
t724 = t776 * t940 + t778 * t938;
t722 = t773 * t940 + t774 * t938;
t718 = t750 * t941 + t985;
t717 = t749 * t941 - t985;
t716 = t750 * t939 - t984;
t715 = t749 * t939 + t984;
t713 = t744 * t941 + t831 * t939;
t704 = pkin(1) * t799 + t970;
t703 = pkin(1) * t797 + t959;
t701 = -t757 * t943 + t758 * t946;
t700 = t757 * t946 + t758 * t943;
t699 = t727 * t941 + t833 * t939;
t698 = t727 * t939 - t833 * t941;
t696 = t723 * t941 + t829 * t939;
t693 = -qJ(3) * t853 - t763 * t939 + t775 * t990;
t692 = -qJ(3) * t855 - t762 * t939 + t775 * t992;
t691 = qJ(3) * t856 + t763 * t941 + t939 * t767;
t690 = -qJ(3) * t857 + t762 * t941 + t775 * t994;
t688 = t726 * t941 + t814 * t939;
t681 = -qJ(3) * t837 + t707 * t941;
t680 = qJ(3) * t839 + t707 * t939;
t679 = -pkin(3) * t724 - t1012;
t678 = t708 * t941 + t775 * t939;
t671 = -t729 * t943 + t731 * t946;
t670 = -t728 * t943 + t730 * t946;
t669 = t729 * t946 + t731 * t943;
t668 = t728 * t946 + t730 * t943;
t667 = pkin(1) * t782 + t965;
t666 = -t716 * t943 + t718 * t946;
t665 = -t715 * t943 + t717 * t946;
t664 = t716 * t946 + t718 * t943;
t663 = t715 * t946 + t717 * t943;
t662 = -t712 * t943 + t713 * t946;
t661 = t712 * t946 + t713 * t943;
t660 = -t698 * t943 + t699 * t946;
t659 = t698 * t946 + t699 * t943;
t658 = -pkin(3) * t741 - t954;
t657 = -t695 * t943 + t696 * t946;
t656 = t695 * t946 + t696 * t943;
t655 = -pkin(1) * t674 - t733;
t654 = -t687 * t943 + t688 * t946;
t653 = t687 * t946 + t688 * t943;
t650 = -pkin(3) * t722 - t960;
t649 = -pkin(6) * t797 - t691 * t943 + t693 * t946;
t648 = -pkin(6) * t799 - t690 * t943 + t692 * t946;
t647 = pkin(6) * t800 + t691 * t946 + t693 * t943;
t646 = -pkin(6) * t801 + t690 * t946 + t692 * t943;
t645 = -qJ(4) * t741 - t702 * t938 + t720 * t940;
t642 = -pkin(6) * t782 - t680 * t943 + t681 * t946;
t641 = pkin(6) * t784 + t680 * t946 + t681 * t943;
t638 = -t677 * t943 + t678 * t946;
t637 = t677 * t946 + t678 * t943;
t636 = -qJ(4) * t722 - t689 * t938 + t710 * t940;
t633 = pkin(6) * t674 + qJ(3) * t1005 - t732 * t943;
t632 = pkin(1) * t937 + pkin(6) * t1036 + qJ(3) * t1006 + t732 * t946;
t631 = -qJ(3) * t677 - (pkin(3) * t939 - qJ(4) * t941) * t707;
t630 = qJ(3) * t678 - (-pkin(3) * t941 - qJ(4) * t939 - pkin(2)) * t707;
t628 = t652 * t938 + t1007;
t626 = t629 * t941 + t759 * t939;
t623 = pkin(1) * t661 + t966;
t622 = -qJ(3) * t712 + t645 * t941 - t658 * t939;
t621 = pkin(1) * t637 + t981;
t620 = -pkin(2) * t741 + qJ(3) * t713 + t645 * t939 + t658 * t941;
t619 = pkin(1) * t656 + t967;
t618 = -qJ(3) * t695 + t636 * t941 - t650 * t939;
t617 = -qJ(4) * t724 - t635 * t938 + t644 * t940;
t616 = -pkin(3) * t628 - t1013;
t615 = -pkin(2) * t722 + qJ(3) * t696 + t636 * t939 + t650 * t941;
t614 = -qJ(3) * t687 + t617 * t941 - t679 * t939;
t613 = -pkin(7) * t1007 - qJ(4) * t628 - t640 * t938;
t612 = -pkin(2) * t724 + qJ(3) * t688 + t617 * t939 + t679 * t941;
t611 = pkin(1) * t653 + t968;
t610 = -t625 * t943 + t626 * t946;
t609 = t625 * t946 + t626 * t943;
t608 = -pkin(6) * t637 - t630 * t943 + t631 * t946;
t607 = pkin(1) * t707 + pkin(6) * t638 + t630 * t946 + t631 * t943;
t606 = -pkin(6) * t661 - t620 * t943 + t622 * t946;
t605 = -pkin(1) * t741 + pkin(6) * t662 + t620 * t946 + t622 * t943;
t604 = -pkin(6) * t656 - t615 * t943 + t618 * t946;
t603 = -pkin(1) * t722 + pkin(6) * t657 + t615 * t946 + t618 * t943;
t602 = -pkin(6) * t653 - t612 * t943 + t614 * t946;
t601 = -pkin(1) * t724 + pkin(6) * t654 + t612 * t946 + t614 * t943;
t600 = -qJ(3) * t625 + t613 * t941 - t616 * t939;
t599 = pkin(1) * t609 + t955;
t598 = -pkin(2) * t628 + qJ(3) * t626 + t613 * t939 + t616 * t941;
t597 = -pkin(6) * t609 - t598 * t943 + t600 * t946;
t596 = -pkin(1) * t628 + pkin(6) * t610 + t598 * t946 + t600 * t943;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t917, t916, 0, 0, 0, 0, 0, 0, 0, t934, t822, t821, 0, -t807, 0, 0, 0, 0, 0, t934, t752, t751, 0, t655, t922, t909, 0, t924, 0, 0, t704, t703, t667, t621, t748, t725, t740, t747, t739, t764, t619, t623, t611, t599; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t910, 0, t911, 0, t980, -t893, -t963, -pkin(5) * t963, 0, 0, t849, 0, -t964, 0, t1045, t1044, t1037, pkin(5) * t1037 + pkin(6) * t1002 + t947 * t790, 0, 0, t787, 0, -t1025, 0, t1055, t1054, t1053, pkin(5) * t1053 + t947 * t632 + t944 * t633, t1024, t783 * t947 + t785 * t944, t745, -t1024, -t1043, 0, -pkin(5) * t1042 + t947 * t646 + t944 * t648, -pkin(5) * t746 + t647 * t947 + t649 * t944, t944 * t642 + t947 * t641 - pkin(5) * (t782 * t944 - t784 * t947), t944 * t608 + t947 * t607 - pkin(5) * (t637 * t944 - t638 * t947), t664 * t947 + t666 * t944, t659 * t947 + t660 * t944, t669 * t947 + t671 * t944, t663 * t947 + t665 * t944, t668 * t947 + t670 * t944, t700 * t947 + t701 * t944, t944 * t604 + t947 * t603 - pkin(5) * (t656 * t944 - t657 * t947), t944 * t606 + t947 * t605 - pkin(5) * (t661 * t944 - t662 * t947), t944 * t602 + t947 * t601 - pkin(5) * (t653 * t944 - t654 * t947), t944 * t597 + t947 * t596 - pkin(5) * (t609 * t944 - t610 * t947); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t911, 0, -t910, 0, t893, t980, t976, pkin(5) * t976, 0, 0, t964, 0, t849, 0, -t1044, t1045, t1038, pkin(5) * t1038 - pkin(6) * t1001 + t944 * t790, 0, 0, t1025, 0, t787, 0, -t1054, t1055, t1052, pkin(5) * t1052 + t944 * t632 - t947 * t633, t1023, t783 * t944 - t785 * t947, t746, -t1023, -t1042, 0, pkin(5) * t1043 + t944 * t646 - t947 * t648, pkin(5) * t745 + t647 * t944 - t649 * t947, -t947 * t642 + t944 * t641 + pkin(5) * (t782 * t947 + t784 * t944), -t947 * t608 + t944 * t607 + pkin(5) * (t637 * t947 + t638 * t944), t664 * t944 - t666 * t947, t659 * t944 - t660 * t947, t669 * t944 - t671 * t947, t663 * t944 - t665 * t947, t668 * t944 - t670 * t947, t700 * t944 - t701 * t947, -t947 * t604 + t944 * t603 + pkin(5) * (t656 * t947 + t657 * t944), -t947 * t606 + t944 * t605 + pkin(5) * (t661 * t947 + t662 * t944), -t947 * t602 + t944 * t601 + pkin(5) * (t653 * t947 + t654 * t944), -t947 * t597 + t944 * t596 + pkin(5) * (t609 * t947 + t610 * t944); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1014, 0, 0, -g(1), t917, 0, 0, 0, -t906, 0, -t903, 0, t1027, t876, t810, pkin(6) * t810, 0, 0, -t845, 0, -t841, 0, t1041, t771, t674, t633, -t818, t785, t800, t818, t801, 0, t648, t649, t642, t608, t666, t660, t671, t665, t670, t701, t604, t606, t602, t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1014, 0, qJDD(1), 0, g(1), 0, -t916, 0, 0, 0, t903, 0, -t906, 0, -t876, t1027, t978, t790, 0, 0, t841, 0, -t845, 0, -t771, t1041, t1036, t632, t815, t783, t797, -t815, -t799, 0, t646, t647, t641, t607, t664, t659, t669, t663, t668, t700, t603, t605, t601, t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t917, t916, 0, 0, 0, 0, 0, 0, 0, t934, t822, t821, 0, -t807, 0, 0, 0, 0, 0, t934, t752, t751, 0, t655, t922, t909, 0, t924, 0, 0, t704, t703, t667, t621, t748, t725, t740, t747, t739, t764, t619, t623, t611, t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, 0, -t933, 0, 0, -g(1), t859, 0, 0, 0, -t900, 0, -t897, 0, t1026, t866, t736, qJ(3) * t736, t862, t840, t856, -t862, t857, 0, t692, t693, t681, t631, t718, t699, t731, t717, t730, t758, t618, t622, t614, t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, 0, t934, 0, g(1), 0, t860, 0, 0, 0, t897, 0, -t900, 0, -t866, t1026, t979, t732, t861, t838, t853, -t861, -t855, 0, t690, t691, t680, t630, t716, t698, t729, t715, t728, t757, t615, t620, t612, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, -t859, -t860, 0, 0, 0, 0, 0, 0, 0, t934, t969, t952, 0, -t733, t922, t909, 0, t924, 0, 0, t970, t959, t965, t981, t748, t725, t740, t747, t739, t764, t967, t966, t968, t955; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, 0, -t933, 0, 0, -t937, t788, 0, t983, t896, t887, -t983, t888, 0, t766, t767, t707, qJ(4) * t707, t750, t727, t743, t749, t742, t765, t636, t645, t617, t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, 0, t934, 0, t937, 0, t789, 0, t912, -t902, -t923, -t912, -t925, 0, t762, t763, 0, pkin(3) * t707, -t834, -t833, -t878, t834, t877, -qJDD(5), t650, t658, t679, t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, -t788, -t789, 0, 0, t922, t909, 0, t924, 0, 0, t982, t961, t971, t988, t748, t725, t740, t747, t739, t764, t973, t972, t974, t958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t923, t925, t912, 0, t926, 0, 0, t775, t762, 0, t806, t779, t795, t804, t794, t820, t710, t720, t644, -pkin(7) * t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t923, -t998, t925, -t912, 0, -t775, 0, t763, 0, t805, t777, t792, t803, t791, t819, t689, t702, t635, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t912, t902, t923, t912, t925, 0, -t762, -t763, 0, 0, t834, t833, t878, -t834, -t877, qJDD(5), t960, t954, t1012, t1013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t832, -t829, t1016, t872, t867, -t872, 0, t759, t684, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t986, t831, t868, t830, t827, -t986, -t759, 0, t685, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t834, t833, t878, -t834, -t877, qJDD(5), -t684, -t685, 0, 0;];
m_new_reg = t1;