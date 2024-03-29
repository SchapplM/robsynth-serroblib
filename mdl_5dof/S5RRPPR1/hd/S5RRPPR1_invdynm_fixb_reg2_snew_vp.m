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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:57
% EndTime: 2022-01-20 09:52:11
% DurationCPUTime: 13.69s
% Computational Cost: add. (62787->493), mult. (88678->687), div. (0->0), fcn. (55375->10), ass. (0->338)
t933 = qJDD(1) + qJDD(2);
t938 = sin(pkin(8));
t1003 = t938 * t933;
t935 = qJD(1) + qJD(2);
t932 = t935 ^ 2;
t940 = cos(pkin(8));
t898 = t932 * t940 + t1003;
t998 = t940 * t933;
t901 = t932 * t938 - t998;
t942 = sin(qJ(2));
t945 = cos(qJ(2));
t842 = t898 * t945 - t901 * t942;
t846 = t898 * t942 + t901 * t945;
t943 = sin(qJ(1));
t946 = cos(qJ(1));
t1027 = t842 * t946 - t846 * t943;
t936 = g(3) - qJDD(3);
t1028 = qJ(3) * t901 - t936 * t938;
t867 = qJ(3) * t898 - t936 * t940;
t1043 = pkin(6) * t846 + t1028 * t945 + t867 * t942;
t772 = pkin(6) * t842 - t1028 * t942 + t867 * t945;
t1056 = pkin(5) * t1027 - t1043 * t943 + t946 * t772;
t1044 = t842 * t943 + t846 * t946;
t1055 = pkin(5) * t1044 + t1043 * t946 + t943 * t772;
t1015 = qJD(1) ^ 2;
t918 = g(1) * t946 + g(2) * t943;
t955 = pkin(1) * t1015 + t918;
t917 = g(1) * t943 - g(2) * t946;
t959 = qJDD(1) * pkin(1) + t917;
t861 = t942 * t959 - t945 * t955;
t848 = -t932 * pkin(2) + t861;
t860 = -t942 * t955 - t945 * t959;
t952 = t933 * pkin(2) - t860;
t789 = t848 * t938 - t940 * t952;
t790 = t940 * t848 + t938 * t952;
t979 = t789 * t938 + t790 * t940;
t737 = t789 * t940 - t790 * t938;
t994 = t942 * t737;
t1036 = t945 * t979 + t994;
t990 = t945 * t737;
t674 = t942 * t979 - t990;
t1053 = t1036 * t946 - t674 * t943;
t1052 = t1036 * t943 + t674 * t946;
t907 = t932 * t942 - t933 * t945;
t1029 = pkin(6) * t907 - g(3) * t942;
t904 = t932 * t945 + t933 * t942;
t877 = pkin(6) * t904 - g(3) * t945;
t963 = t904 * t946 - t907 * t943;
t1047 = pkin(5) * t963 - t1029 * t943 + t946 * t877;
t937 = sin(pkin(9));
t934 = t937 ^ 2;
t939 = cos(pkin(9));
t949 = t939 ^ 2;
t1016 = t932 * (t934 + t949);
t889 = t939 * t1016;
t982 = t939 * t998;
t856 = -t889 * t938 + t982;
t926 = t939 * t933;
t858 = t889 * t940 + t926 * t938;
t800 = t856 * t945 - t858 * t942;
t802 = t856 * t942 + t858 * t945;
t1046 = t800 * t946 - t802 * t943;
t1045 = t800 * t943 + t802 * t946;
t1024 = t904 * t943 + t907 * t946;
t1042 = pkin(5) * t1024 + t1029 * t946 + t943 * t877;
t978 = t860 * t942 + t861 * t945;
t811 = t860 * t945 - t861 * t942;
t989 = t946 * t811;
t1038 = t943 * t978 - t989;
t993 = t943 * t811;
t1037 = t946 * t978 + t993;
t941 = sin(qJ(5));
t944 = cos(qJ(5));
t1030 = t937 * t941 - t939 * t944;
t882 = t1030 * t935;
t961 = t937 * t944 + t939 * t941;
t884 = t961 * t935;
t835 = t884 * t882;
t1017 = qJDD(5) - t835;
t1033 = t1017 * t941;
t1032 = t1017 * t944;
t975 = -t932 * pkin(3) + t933 * qJ(4) + 0.2e1 * qJD(4) * t935 + t790;
t1031 = pkin(7) * t933 + t975;
t1006 = t937 * t939;
t862 = t898 * t1006;
t913 = t932 * t1006;
t863 = -t913 * t938 + t937 * t982;
t816 = t862 * t945 + t863 * t942;
t819 = t862 * t942 - t863 * t945;
t1026 = t816 * t946 - t819 * t943;
t1025 = t816 * t943 + t819 * t946;
t928 = t939 * t936;
t763 = t937 * t975 + t928;
t1007 = t937 * t936;
t764 = t939 * t975 - t1007;
t709 = t763 * t937 + t764 * t939;
t1009 = t934 * t932;
t927 = t949 * t932;
t902 = t927 + t1009;
t880 = t882 ^ 2;
t881 = t884 ^ 2;
t755 = -t928 + (pkin(4) * t932 * t939 - t1031) * t937;
t756 = -pkin(4) * t927 + t1031 * t939 - t1007;
t685 = -t755 * t944 + t756 * t941;
t686 = t755 * t941 + t756 * t944;
t652 = -t685 * t944 + t686 * t941;
t1014 = pkin(4) * t652;
t878 = t1030 * t933;
t879 = t961 * t933;
t777 = -t878 * t941 - t879 * t944;
t1013 = pkin(4) * t777;
t1008 = t937 * t652;
t776 = -pkin(3) * t933 - qJ(4) * t932 + qJDD(4) + t789;
t767 = t937 * t776;
t924 = t937 * t933;
t1004 = t938 * t776;
t1001 = t939 * t652;
t768 = t939 * t776;
t999 = t940 * t776;
t760 = -pkin(4) * t926 - pkin(7) * t902 + t776;
t996 = t941 * t760;
t828 = qJDD(5) + t835;
t995 = t941 * t828;
t992 = t944 * t760;
t991 = t944 * t828;
t988 = -pkin(3) * t776 + qJ(4) * t709;
t873 = t882 * qJD(5);
t986 = t884 * qJD(5);
t985 = t937 * t926;
t984 = t938 * t835;
t983 = t940 * t835;
t981 = pkin(3) * t926 - qJ(4) * t889 - t768;
t678 = t709 * t938 - t999;
t980 = pkin(2) * t678 + t988;
t653 = t685 * t941 + t686 * t944;
t976 = -t917 * t943 - t918 * t946;
t779 = -t878 * t944 + t879 * t941;
t815 = -t880 - t881;
t636 = -pkin(4) * t815 + pkin(7) * t779 + t653;
t645 = -pkin(7) * t777 - t652;
t727 = -t777 * t937 + t779 * t939;
t974 = -pkin(3) * t815 + qJ(4) * t727 + t636 * t939 + t645 * t937;
t947 = qJD(5) ^ 2;
t826 = -t947 - t880;
t775 = t826 * t944 - t1033;
t830 = t878 + 0.2e1 * t986;
t690 = -pkin(4) * t830 + pkin(7) * t775 - t992;
t774 = t826 * t941 + t1032;
t711 = -pkin(7) * t774 + t996;
t724 = -t774 * t937 + t775 * t939;
t973 = -pkin(3) * t830 + qJ(4) * t724 + t690 * t939 + t711 * t937;
t870 = -t881 - t947;
t797 = -t870 * t941 - t991;
t832 = -0.2e1 * t873 + t879;
t703 = -pkin(4) * t832 + pkin(7) * t797 + t996;
t794 = t870 * t944 - t995;
t721 = -pkin(7) * t794 + t992;
t745 = -t794 * t937 + t797 * t939;
t972 = -pkin(3) * t832 + qJ(4) * t745 + t703 * t939 + t721 * t937;
t923 = t934 * t933;
t925 = t949 * t933;
t896 = t925 + t923;
t971 = pkin(3) * t902 + qJ(4) * t896 + t709;
t970 = pkin(2) * t856 + t981;
t912 = qJDD(1) * t946 - t1015 * t943;
t969 = -pkin(5) * t912 - g(3) * t943;
t968 = -pkin(2) * t901 - t789;
t688 = t727 * t938 - t815 * t940;
t967 = pkin(2) * t688 + t974;
t696 = t724 * t938 - t830 * t940;
t966 = pkin(2) * t696 + t973;
t713 = t745 * t938 - t832 * t940;
t965 = pkin(2) * t713 + t972;
t838 = t896 * t938 + t902 * t940;
t964 = pkin(2) * t838 + t971;
t708 = t763 * t939 - t764 * t937;
t962 = t917 * t946 - t918 * t943;
t888 = t937 * t1016;
t960 = -pkin(3) * t924 + qJ(4) * t888 + t767;
t958 = pkin(4) * t774 - t685;
t854 = t888 * t938 - t937 * t998;
t957 = pkin(2) * t854 + t960;
t630 = t653 * t939 - t1008;
t641 = -pkin(4) * t760 + pkin(7) * t653;
t956 = -pkin(3) * t760 - pkin(7) * t1008 + qJ(4) * t630 + t641 * t939;
t626 = t630 * t938 - t760 * t940;
t954 = pkin(2) * t626 + t956;
t953 = pkin(4) * t794 - t686;
t951 = -pkin(2) * t898 - t790;
t911 = qJDD(1) * t943 + t1015 * t946;
t910 = 0.2e1 * t985;
t903 = -t927 + t1009;
t897 = t925 - t923;
t894 = -pkin(5) * t911 + g(3) * t946;
t869 = -t881 + t947;
t868 = t880 - t947;
t857 = t1003 * t937 + t888 * t940;
t841 = t897 * t940 + t903 * t938;
t840 = t896 * t940 - t902 * t938;
t839 = t897 * t938 - t903 * t940;
t834 = t881 - t880;
t833 = -t873 + t879;
t831 = -t878 - t986;
t823 = -pkin(1) * t907 - t860;
t822 = -pkin(1) * t904 - t861;
t821 = (-t882 * t944 + t884 * t941) * qJD(5);
t820 = (-t882 * t941 - t884 * t944) * qJD(5);
t808 = pkin(1) * t811;
t807 = t833 * t944 - t941 * t986;
t806 = t833 * t941 + t944 * t986;
t805 = -t831 * t941 + t873 * t944;
t804 = t831 * t944 + t873 * t941;
t801 = -t854 * t942 + t857 * t945;
t798 = t854 * t945 + t857 * t942;
t796 = -t869 * t941 + t1032;
t795 = t868 * t944 - t995;
t793 = t869 * t944 + t1033;
t792 = t868 * t941 + t991;
t791 = pkin(1) * g(3) + pkin(6) * t978;
t786 = -t839 * t942 + t841 * t945;
t785 = -t838 * t942 + t840 * t945;
t784 = t839 * t945 + t841 * t942;
t783 = t838 * t945 + t840 * t942;
t780 = -t830 * t944 - t832 * t941;
t778 = -t830 * t941 + t832 * t944;
t766 = -t820 * t937 + t821 * t939;
t765 = t820 * t939 + t821 * t937;
t759 = qJDD(5) * t938 + t766 * t940;
t758 = -qJDD(5) * t940 + t766 * t938;
t753 = -pkin(1) * t846 + t968;
t752 = -pkin(1) * t842 + t951;
t751 = -t806 * t937 + t807 * t939;
t750 = -t804 * t937 + t805 * t939;
t749 = t806 * t939 + t807 * t937;
t748 = t804 * t939 + t805 * t937;
t747 = -t798 * t943 + t801 * t946;
t746 = t798 * t946 + t801 * t943;
t744 = -t793 * t937 + t796 * t939;
t743 = -t792 * t937 + t795 * t939;
t742 = t794 * t939 + t797 * t937;
t741 = t793 * t939 + t796 * t937;
t740 = t792 * t939 + t795 * t937;
t734 = pkin(2) * t737;
t733 = pkin(2) * t936 + qJ(3) * t979;
t732 = t744 * t940 + t879 * t938;
t731 = t743 * t940 - t878 * t938;
t730 = t744 * t938 - t879 * t940;
t729 = t743 * t938 + t878 * t940;
t728 = -t778 * t937 + t780 * t939;
t726 = t778 * t939 + t780 * t937;
t725 = t777 * t939 + t779 * t937;
t723 = t774 * t939 + t775 * t937;
t719 = t751 * t940 + t984;
t718 = t750 * t940 - t984;
t717 = t751 * t938 - t983;
t716 = t750 * t938 + t983;
t714 = t745 * t940 + t832 * t938;
t705 = pkin(1) * t800 + t970;
t704 = pkin(1) * t798 + t957;
t702 = -t758 * t942 + t759 * t945;
t701 = t758 * t945 + t759 * t942;
t700 = t728 * t940 + t834 * t938;
t699 = t728 * t938 - t834 * t940;
t697 = t724 * t940 + t830 * t938;
t694 = -qJ(3) * t854 - t764 * t938 + t939 * t999;
t693 = -qJ(3) * t856 - t763 * t938 + t937 * t999;
t692 = qJ(3) * t857 + t764 * t940 + t768 * t938;
t691 = -qJ(3) * t858 + t1004 * t937 + t763 * t940;
t689 = t727 * t940 + t815 * t938;
t682 = -qJ(3) * t838 + t708 * t940;
t681 = qJ(3) * t840 + t708 * t938;
t680 = -pkin(3) * t725 - t1013;
t679 = t709 * t940 + t1004;
t672 = -t730 * t942 + t732 * t945;
t671 = -t729 * t942 + t731 * t945;
t670 = t730 * t945 + t732 * t942;
t669 = t729 * t945 + t731 * t942;
t668 = pkin(1) * t783 + t964;
t667 = -t717 * t942 + t719 * t945;
t666 = -t716 * t942 + t718 * t945;
t665 = t717 * t945 + t719 * t942;
t664 = t716 * t945 + t718 * t942;
t663 = -t713 * t942 + t714 * t945;
t662 = t713 * t945 + t714 * t942;
t661 = -t699 * t942 + t700 * t945;
t660 = t699 * t945 + t700 * t942;
t659 = -pkin(3) * t742 - t953;
t658 = -t696 * t942 + t697 * t945;
t657 = t696 * t945 + t697 * t942;
t656 = pkin(1) * t674 - t734;
t655 = -t688 * t942 + t689 * t945;
t654 = t688 * t945 + t689 * t942;
t651 = -pkin(3) * t723 - t958;
t650 = -pkin(6) * t798 - t692 * t942 + t694 * t945;
t649 = -pkin(6) * t800 - t691 * t942 + t693 * t945;
t648 = pkin(6) * t801 + t692 * t945 + t694 * t942;
t647 = -pkin(6) * t802 + t691 * t945 + t693 * t942;
t646 = -qJ(4) * t742 - t703 * t937 + t721 * t939;
t643 = -pkin(6) * t783 - t681 * t942 + t682 * t945;
t642 = pkin(6) * t785 + t681 * t945 + t682 * t942;
t639 = -t678 * t942 + t679 * t945;
t638 = t678 * t945 + t679 * t942;
t637 = -qJ(4) * t723 - t690 * t937 + t711 * t939;
t634 = -pkin(6) * t674 + qJ(3) * t990 - t733 * t942;
t633 = pkin(1) * t936 + pkin(6) * t1036 + qJ(3) * t994 + t733 * t945;
t632 = -qJ(3) * t678 - (pkin(3) * t938 - qJ(4) * t940) * t708;
t631 = qJ(3) * t679 - (-pkin(3) * t940 - qJ(4) * t938 - pkin(2)) * t708;
t629 = t653 * t937 + t1001;
t627 = t630 * t940 + t760 * t938;
t624 = pkin(1) * t662 + t965;
t623 = -qJ(3) * t713 + t646 * t940 - t659 * t938;
t622 = pkin(1) * t638 + t980;
t621 = -pkin(2) * t742 + qJ(3) * t714 + t646 * t938 + t659 * t940;
t620 = pkin(1) * t657 + t966;
t619 = -qJ(3) * t696 + t637 * t940 - t651 * t938;
t618 = -qJ(4) * t725 - t636 * t937 + t645 * t939;
t617 = -pkin(3) * t629 - t1014;
t616 = -pkin(2) * t723 + qJ(3) * t697 + t637 * t938 + t651 * t940;
t615 = -qJ(3) * t688 + t618 * t940 - t680 * t938;
t614 = -pkin(7) * t1001 - qJ(4) * t629 - t641 * t937;
t613 = -pkin(2) * t725 + qJ(3) * t689 + t618 * t938 + t680 * t940;
t612 = pkin(1) * t654 + t967;
t611 = -t626 * t942 + t627 * t945;
t610 = t626 * t945 + t627 * t942;
t609 = -pkin(6) * t638 - t631 * t942 + t632 * t945;
t608 = pkin(1) * t708 + pkin(6) * t639 + t631 * t945 + t632 * t942;
t607 = -pkin(6) * t662 - t621 * t942 + t623 * t945;
t606 = -pkin(1) * t742 + pkin(6) * t663 + t621 * t945 + t623 * t942;
t605 = -pkin(6) * t657 - t616 * t942 + t619 * t945;
t604 = -pkin(1) * t723 + pkin(6) * t658 + t616 * t945 + t619 * t942;
t603 = -pkin(6) * t654 - t613 * t942 + t615 * t945;
t602 = -pkin(1) * t725 + pkin(6) * t655 + t613 * t945 + t615 * t942;
t601 = -qJ(3) * t626 + t614 * t940 - t617 * t938;
t600 = pkin(1) * t610 + t954;
t599 = -pkin(2) * t629 + qJ(3) * t627 + t614 * t938 + t617 * t940;
t598 = -pkin(6) * t610 - t599 * t942 + t601 * t945;
t597 = -pkin(1) * t629 + pkin(6) * t611 + t599 * t945 + t601 * t942;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t912, 0, -t911, 0, t969, -t894, -t962, -pkin(5) * t962, 0, 0, -t1024, 0, -t963, 0, t1042, t1047, -t1038, -pkin(5) * t1038 + pkin(6) * t989 - t943 * t791, 0, 0, -t1044, 0, -t1027, 0, t1055, t1056, -t1052, -pkin(5) * t1052 - t943 * t633 + t946 * t634, -t1025, -t784 * t943 + t786 * t946, t747, t1025, t1045, 0, -pkin(5) * t1046 - t943 * t647 + t946 * t649, -pkin(5) * t746 - t648 * t943 + t650 * t946, t946 * t643 - t943 * t642 - pkin(5) * (t783 * t946 + t785 * t943), t946 * t609 - t943 * t608 - pkin(5) * (t638 * t946 + t639 * t943), -t665 * t943 + t667 * t946, -t660 * t943 + t661 * t946, -t670 * t943 + t672 * t946, -t664 * t943 + t666 * t946, -t669 * t943 + t671 * t946, -t701 * t943 + t702 * t946, t946 * t605 - t943 * t604 - pkin(5) * (t657 * t946 + t658 * t943), t946 * t607 - t943 * t606 - pkin(5) * (t662 * t946 + t663 * t943), t946 * t603 - t943 * t602 - pkin(5) * (t654 * t946 + t655 * t943), t946 * t598 - t943 * t597 - pkin(5) * (t610 * t946 + t611 * t943); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t911, 0, t912, 0, t894, t969, t976, pkin(5) * t976, 0, 0, t963, 0, -t1024, 0, -t1047, t1042, t1037, pkin(5) * t1037 + pkin(6) * t993 + t946 * t791, 0, 0, t1027, 0, -t1044, 0, -t1056, t1055, t1053, pkin(5) * t1053 + t946 * t633 + t943 * t634, t1026, t784 * t946 + t786 * t943, t746, -t1026, -t1046, 0, -pkin(5) * t1045 + t946 * t647 + t943 * t649, pkin(5) * t747 + t648 * t946 + t650 * t943, t943 * t643 + t946 * t642 + pkin(5) * (-t783 * t943 + t785 * t946), t943 * t609 + t946 * t608 + pkin(5) * (-t638 * t943 + t639 * t946), t665 * t946 + t667 * t943, t660 * t946 + t661 * t943, t670 * t946 + t672 * t943, t664 * t946 + t666 * t943, t669 * t946 + t671 * t943, t701 * t946 + t702 * t943, t943 * t605 + t946 * t604 + pkin(5) * (-t657 * t943 + t658 * t946), t943 * t607 + t946 * t606 + pkin(5) * (-t662 * t943 + t663 * t946), t943 * t603 + t946 * t602 + pkin(5) * (-t654 * t943 + t655 * t946), t943 * t598 + t946 * t597 + pkin(5) * (-t610 * t943 + t611 * t946); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t917, t918, 0, 0, 0, 0, 0, 0, 0, t933, t823, t822, 0, -t808, 0, 0, 0, 0, 0, t933, t753, t752, 0, t656, t923, t910, 0, t925, 0, 0, t705, t704, t668, t622, t749, t726, t741, t748, t740, t765, t620, t624, t612, t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1015, 0, 0, -g(3), -t917, 0, 0, 0, -t907, 0, -t904, 0, t1029, t877, t811, pkin(6) * t811, 0, 0, -t846, 0, -t842, 0, t1043, t772, -t674, t634, -t819, t786, t801, t819, t802, 0, t649, t650, t643, t609, t667, t661, t672, t666, t671, t702, t605, t607, t603, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1015, 0, qJDD(1), 0, g(3), 0, -t918, 0, 0, 0, t904, 0, -t907, 0, -t877, t1029, t978, t791, 0, 0, t842, 0, -t846, 0, -t772, t1043, t1036, t633, t816, t784, t798, -t816, -t800, 0, t647, t648, t642, t608, t665, t660, t670, t664, t669, t701, t604, t606, t602, t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t917, t918, 0, 0, 0, 0, 0, 0, 0, t933, t823, t822, 0, -t808, 0, 0, 0, 0, 0, t933, t753, t752, 0, t656, t923, t910, 0, t925, 0, 0, t705, t704, t668, t622, t749, t726, t741, t748, t740, t765, t620, t624, t612, t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, 0, -t932, 0, 0, -g(3), t860, 0, 0, 0, -t901, 0, -t898, 0, t1028, t867, t737, qJ(3) * t737, t863, t841, t857, -t863, t858, 0, t693, t694, t682, t632, t719, t700, t732, t718, t731, t759, t619, t623, t615, t601; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t932, 0, t933, 0, g(3), 0, t861, 0, 0, 0, t898, 0, -t901, 0, -t867, t1028, t979, t733, t862, t839, t854, -t862, -t856, 0, t691, t692, t681, t631, t717, t699, t730, t716, t729, t758, t616, t621, t613, t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, -t860, -t861, 0, 0, 0, 0, 0, 0, 0, t933, t968, t951, 0, -t734, t923, t910, 0, t925, 0, 0, t970, t957, t964, t980, t749, t726, t741, t748, t740, t765, t966, t965, t967, t954; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, 0, -t932, 0, 0, -t936, t789, 0, t985, t897, t888, -t985, t889, 0, t767, t768, t708, qJ(4) * t708, t751, t728, t744, t750, t743, t766, t637, t646, t618, t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t932, 0, t933, 0, t936, 0, t790, 0, t913, -t903, -t924, -t913, -t926, 0, t763, t764, 0, pkin(3) * t708, -t835, -t834, -t879, t835, t878, -qJDD(5), t651, t659, t680, t617; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, -t789, -t790, 0, 0, t923, t910, 0, t925, 0, 0, t981, t960, t971, t988, t749, t726, t741, t748, t740, t765, t973, t972, t974, t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t924, t926, t913, 0, t927, 0, 0, t776, t763, 0, t807, t780, t796, t805, t795, t821, t711, t721, t645, -pkin(7) * t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t924, -t1009, t926, -t913, 0, -t776, 0, t764, 0, t806, t778, t793, t804, t792, t820, t690, t703, t636, t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t913, t903, t924, t913, t926, 0, -t763, -t764, 0, 0, t835, t834, t879, -t835, -t878, qJDD(5), t958, t953, t1013, t1014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t833, -t830, t1017, t873, t868, -t873, 0, t760, t685, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t986, t832, t869, t831, t828, -t986, -t760, 0, t686, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t835, t834, t879, -t835, -t878, qJDD(5), -t685, -t686, 0, 0;];
m_new_reg = t1;
