% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRRRR6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:31
% EndTime: 2019-12-05 17:10:44
% DurationCPUTime: 12.66s
% Computational Cost: add. (59596->546), mult. (78250->770), div. (0->0), fcn. (53006->10), ass. (0->367)
t947 = qJDD(2) + qJDD(3);
t959 = cos(qJ(3));
t1008 = t959 * t947;
t949 = qJD(2) + qJD(3);
t945 = t949 ^ 2;
t955 = sin(qJ(3));
t912 = t955 * t945 - t1008;
t951 = sin(pkin(9));
t952 = cos(pkin(9));
t926 = t951 * g(1) - t952 * g(2);
t1042 = pkin(6) * t912 - t955 * t926;
t1017 = t955 * t947;
t909 = t959 * t945 + t1017;
t956 = sin(qJ(2));
t960 = cos(qJ(2));
t851 = t956 * t909 + t960 * t912;
t862 = pkin(6) * t909 - t959 * t926;
t1048 = pkin(5) * t851 + t1042 * t960 + t956 * t862;
t846 = t960 * t909 - t956 * t912;
t764 = pkin(5) * t846 - t1042 * t956 + t960 * t862;
t1004 = t952 * g(1) + t951 * g(2);
t1006 = g(3) - qJDD(1);
t890 = -t960 * t1004 - t956 * t1006;
t962 = qJD(2) ^ 2;
t888 = -t962 * pkin(2) + t890;
t889 = -t956 * t1004 + t960 * t1006;
t963 = qJDD(2) * pkin(2) - t889;
t817 = t955 * t888 - t959 * t963;
t818 = t959 * t888 + t955 * t963;
t747 = t959 * t817 - t955 * t818;
t1016 = t956 * t747;
t989 = t955 * t817 + t959 * t818;
t1045 = t960 * t989 + t1016;
t1007 = t960 * t747;
t703 = -t956 * t989 + t1007;
t954 = sin(qJ(4));
t1027 = t949 * t954;
t953 = sin(qJ(5));
t957 = cos(qJ(5));
t958 = cos(qJ(4));
t893 = -t957 * t958 * t949 + t953 * t1027;
t895 = (t953 * t958 + t954 * t957) * t949;
t842 = t895 * t893;
t946 = qJDD(4) + qJDD(5);
t1039 = -t842 + t946;
t1044 = t1039 * t953;
t1043 = t1039 * t957;
t806 = -t945 * pkin(3) + t947 * pkin(7) + t818;
t791 = t954 * t806 + t958 * t926;
t792 = t958 * t806 - t954 * t926;
t729 = t954 * t791 + t958 * t792;
t1019 = t954 * t947;
t1003 = qJD(4) * t949;
t934 = t958 * t1003;
t905 = t934 + t1019;
t1010 = t958 * t947;
t996 = t954 * t1003;
t971 = t996 - t1010;
t810 = -t893 * qJD(5) + t957 * t905 - t953 * t971;
t948 = qJD(4) + qJD(5);
t886 = t948 * t893;
t1038 = -t886 + t810;
t915 = t952 * t926;
t1037 = t951 * t1004 - t915;
t1036 = t951 * t1006;
t1035 = t952 * t1006;
t987 = t953 * t905 + t957 * t971;
t784 = (qJD(5) - t948) * t895 + t987;
t891 = t893 ^ 2;
t892 = t895 ^ 2;
t944 = t948 ^ 2;
t1034 = t958 ^ 2;
t1033 = pkin(2) * t747;
t931 = t954 * t945 * t958;
t920 = qJDD(4) + t931;
t743 = (-t905 + t934) * pkin(8) + t920 * pkin(4) - t791;
t922 = qJD(4) * pkin(4) - pkin(8) * t1027;
t936 = t1034 * t945;
t744 = -pkin(4) * t936 - pkin(8) * t971 - qJD(4) * t922 + t792;
t697 = -t957 * t743 + t953 * t744;
t698 = t953 * t743 + t957 * t744;
t654 = -t957 * t697 + t953 * t698;
t1032 = pkin(4) * t654;
t788 = t886 + t810;
t723 = -t784 * t953 - t957 * t788;
t1031 = pkin(4) * t723;
t1030 = t895 * t948;
t1029 = t948 * t953;
t1028 = t948 * t957;
t950 = t954 ^ 2;
t1026 = t950 * t945;
t1025 = t951 * t926;
t805 = -t947 * pkin(3) - t945 * pkin(7) + t817;
t760 = pkin(4) * t971 - pkin(8) * t936 + t922 * t1027 + t805;
t1024 = t953 * t760;
t835 = t842 + t946;
t1023 = t953 * t835;
t1022 = t954 * t654;
t801 = t954 * t805;
t1021 = t954 * t920;
t921 = qJDD(4) - t931;
t1020 = t954 * t921;
t1015 = t957 * t760;
t1014 = t957 * t835;
t1013 = t958 * t654;
t802 = t958 * t805;
t1012 = t958 * t920;
t1011 = t958 * t921;
t1005 = -pkin(3) * t805 + pkin(7) * t729;
t1001 = t950 + t1034;
t1000 = t955 * t842;
t999 = t959 * t842;
t961 = qJD(4) ^ 2;
t928 = -t961 - t1026;
t870 = -t954 * t928 - t1011;
t904 = 0.2e1 * t934 + t1019;
t998 = -pkin(3) * t904 + pkin(7) * t870 + t801;
t930 = -t936 - t961;
t868 = t958 * t930 - t1021;
t906 = -0.2e1 * t996 + t1010;
t997 = pkin(3) * t906 + pkin(7) * t868 - t802;
t908 = t1001 * t947;
t913 = t936 + t1026;
t845 = t955 * t908 + t959 * t913;
t849 = t959 * t908 - t955 * t913;
t793 = t960 * t845 + t956 * t849;
t979 = pkin(3) * t913 + pkin(7) * t908 + t729;
t966 = pkin(2) * t845 + t979;
t675 = -pkin(1) * t793 - t966;
t794 = -t956 * t845 + t960 * t849;
t994 = qJ(1) * t794 + t675;
t977 = -pkin(2) * t909 - t818;
t758 = pkin(1) * t846 - t977;
t993 = qJ(1) * t851 + t758;
t973 = -pkin(2) * t912 - t817;
t759 = pkin(1) * t851 - t973;
t992 = -qJ(1) * t846 + t759;
t924 = t956 * qJDD(2) + t960 * t962;
t853 = pkin(1) * t924 + t890;
t925 = t960 * qJDD(2) - t956 * t962;
t991 = qJ(1) * t925 - t853;
t854 = -pkin(1) * t925 + t889;
t990 = qJ(1) * t924 - t854;
t655 = t953 * t697 + t957 * t698;
t988 = t956 * t889 + t960 * t890;
t985 = -t952 * t1004 - t1025;
t984 = t955 * t931;
t983 = t959 * t931;
t725 = -t784 * t957 + t953 * t788;
t813 = -t891 - t892;
t630 = -pkin(4) * t813 + pkin(8) * t725 + t655;
t637 = -pkin(8) * t723 - t654;
t673 = -t954 * t723 + t958 * t725;
t982 = -pkin(3) * t813 + pkin(7) * t673 + t958 * t630 + t954 * t637;
t829 = -t944 - t891;
t772 = t957 * t829 - t1044;
t783 = (qJD(5) + t948) * t895 + t987;
t679 = -pkin(4) * t783 + pkin(8) * t772 - t1015;
t771 = t953 * t829 + t1043;
t712 = -pkin(8) * t771 + t1024;
t715 = -t954 * t771 + t958 * t772;
t981 = -pkin(3) * t783 + pkin(7) * t715 + t958 * t679 + t954 * t712;
t871 = -t892 - t944;
t796 = -t953 * t871 - t1014;
t685 = -pkin(4) * t1038 + pkin(8) * t796 + t1024;
t795 = t957 * t871 - t1023;
t717 = -pkin(8) * t795 + t1015;
t732 = -t954 * t795 + t958 * t796;
t980 = -pkin(3) * t1038 + pkin(7) * t732 + t958 * t685 + t954 * t717;
t699 = t955 * t729 - t959 * t805;
t978 = pkin(2) * t699 + t1005;
t882 = pkin(5) * t924 - t960 * t926;
t976 = -pkin(5) * t925 - t956 * t926;
t728 = t958 * t791 - t954 * t792;
t827 = t960 * t889 - t956 * t890;
t822 = t955 * t870 - t959 * t904;
t975 = pkin(2) * t822 + t998;
t821 = t955 * t868 + t959 * t906;
t974 = pkin(2) * t821 + t997;
t972 = pkin(4) * t771 - t697;
t622 = t958 * t655 - t1022;
t640 = -pkin(4) * t760 + pkin(8) * t655;
t970 = -pkin(3) * t760 + pkin(7) * t622 - pkin(8) * t1022 + t958 * t640;
t665 = t955 * t673 - t959 * t813;
t969 = pkin(2) * t665 + t982;
t680 = t955 * t715 - t959 * t783;
t968 = pkin(2) * t680 + t981;
t688 = -t1038 * t959 + t955 * t732;
t967 = pkin(2) * t688 + t980;
t965 = pkin(4) * t795 - t698;
t618 = t955 * t622 - t959 * t760;
t964 = pkin(2) * t618 + t970;
t929 = t936 - t961;
t927 = t961 - t1026;
t923 = pkin(1) * t926;
t914 = -t936 + t1026;
t902 = t1001 * t1003;
t884 = -t892 + t944;
t883 = t891 - t944;
t878 = t955 * qJDD(4) + t959 * t902;
t877 = -t959 * qJDD(4) + t955 * t902;
t873 = -t950 * t1003 + t958 * t905;
t872 = -t1034 * t1003 + t954 * t971;
t869 = -t954 * t927 + t1012;
t867 = t958 * t929 - t1020;
t866 = t958 * t928 - t1020;
t865 = t958 * t927 + t1021;
t864 = t954 * t930 + t1012;
t863 = t954 * t929 + t1011;
t856 = (t905 + t934) * t954;
t855 = t954 * t934 + t958 * t971;
t844 = -t954 * t904 + t958 * t906;
t843 = t958 * t904 + t954 * t906;
t841 = t892 - t891;
t840 = t954 * t1017 + t959 * t869;
t839 = t955 * t1010 + t959 * t867;
t838 = -t954 * t1008 + t955 * t869;
t837 = -t958 * t1008 + t955 * t867;
t833 = t959 * t873 - t984;
t832 = t959 * t872 + t984;
t831 = t955 * t873 + t983;
t830 = t955 * t872 - t983;
t824 = t959 * t870 + t955 * t904;
t823 = t959 * t868 - t955 * t906;
t820 = (-t893 * t957 + t895 * t953) * t948;
t819 = (-t893 * t953 - t895 * t957) * t948;
t815 = -t956 * t877 + t960 * t878;
t814 = t960 * t877 + t956 * t878;
t812 = t959 * t844 + t955 * t914;
t811 = t955 * t844 - t959 * t914;
t809 = -t895 * qJD(5) - t987;
t808 = pkin(5) * t988 + t923;
t800 = t957 * t883 - t1023;
t799 = -t953 * t884 + t1043;
t798 = t953 * t883 + t1014;
t797 = t957 * t884 + t1044;
t780 = -t895 * t1029 + t957 * t810;
t779 = t895 * t1028 + t953 * t810;
t778 = t893 * t1028 - t953 * t809;
t777 = t893 * t1029 + t957 * t809;
t776 = -t956 * t838 + t960 * t840;
t775 = -t956 * t837 + t960 * t839;
t774 = t960 * t838 + t956 * t840;
t773 = t960 * t837 + t956 * t839;
t770 = -t956 * t831 + t960 * t833;
t769 = -t956 * t830 + t960 * t832;
t768 = t960 * t831 + t956 * t833;
t767 = t960 * t830 + t956 * t832;
t766 = -pkin(7) * t866 + t802;
t765 = -pkin(7) * t864 + t801;
t756 = -t956 * t822 + t960 * t824;
t755 = -t956 * t821 + t960 * t823;
t754 = t960 * t822 + t956 * t824;
t753 = t960 * t821 + t956 * t823;
t752 = -pkin(3) * t866 + t792;
t751 = -pkin(3) * t864 + t791;
t750 = -t954 * t819 + t958 * t820;
t749 = t958 * t819 + t954 * t820;
t741 = t959 * t750 + t955 * t946;
t740 = t955 * t750 - t959 * t946;
t739 = pkin(2) * t926 + pkin(6) * t989;
t738 = -t956 * t811 + t960 * t812;
t737 = t960 * t811 + t956 * t812;
t736 = -t954 * t798 + t958 * t800;
t735 = -t954 * t797 + t958 * t799;
t734 = t958 * t798 + t954 * t800;
t733 = t958 * t797 + t954 * t799;
t731 = t958 * t795 + t954 * t796;
t724 = -t1038 * t953 - t957 * t783;
t722 = t1038 * t957 - t953 * t783;
t721 = -t954 * t779 + t958 * t780;
t720 = -t954 * t777 + t958 * t778;
t719 = t958 * t779 + t954 * t780;
t718 = t958 * t777 + t954 * t778;
t714 = t958 * t771 + t954 * t772;
t710 = -pkin(6) * t845 + t959 * t728;
t709 = pkin(6) * t849 + t955 * t728;
t708 = t959 * t721 + t1000;
t707 = t959 * t720 - t1000;
t706 = t955 * t721 - t999;
t705 = t955 * t720 + t999;
t700 = t959 * t729 + t955 * t805;
t696 = t959 * t736 - t955 * t784;
t695 = t959 * t735 + t955 * t788;
t694 = t955 * t736 + t959 * t784;
t693 = t955 * t735 - t959 * t788;
t691 = -t956 * t740 + t960 * t741;
t690 = t960 * t740 + t956 * t741;
t689 = t1038 * t955 + t959 * t732;
t687 = -pkin(6) * t822 - t955 * t752 + t959 * t766;
t686 = -pkin(6) * t821 - t955 * t751 + t959 * t765;
t683 = -pkin(1) * t754 - t975;
t682 = -pkin(1) * t753 - t974;
t681 = t959 * t715 + t955 * t783;
t677 = -pkin(2) * t866 + pkin(6) * t824 + t959 * t752 + t955 * t766;
t676 = -pkin(2) * t864 + pkin(6) * t823 + t959 * t751 + t955 * t765;
t674 = pkin(1) * t703 + t1033;
t672 = -t954 * t722 + t958 * t724;
t671 = t958 * t723 + t954 * t725;
t670 = t958 * t722 + t954 * t724;
t668 = t959 * t672 + t955 * t841;
t667 = t955 * t672 - t959 * t841;
t666 = t959 * t673 + t955 * t813;
t664 = -t956 * t706 + t960 * t708;
t663 = -t956 * t705 + t960 * t707;
t662 = t960 * t706 + t956 * t708;
t661 = t960 * t705 + t956 * t707;
t660 = -t956 * t699 + t960 * t700;
t659 = t960 * t699 + t956 * t700;
t658 = pkin(5) * t703 + pkin(6) * t1007 - t956 * t739;
t657 = -pkin(3) * t731 - t965;
t656 = pkin(5) * t1045 + pkin(6) * t1016 + t960 * t739 + t923;
t653 = -t956 * t694 + t960 * t696;
t652 = -t956 * t693 + t960 * t695;
t651 = t960 * t694 + t956 * t696;
t650 = t960 * t693 + t956 * t695;
t649 = -pkin(5) * t793 - t956 * t709 + t960 * t710;
t648 = pkin(5) * t794 + t960 * t709 + t956 * t710;
t647 = -t956 * t688 + t960 * t689;
t646 = t960 * t688 + t956 * t689;
t645 = -pkin(3) * t671 - t1031;
t644 = -pkin(3) * t714 - t972;
t643 = -t956 * t680 + t960 * t681;
t642 = t960 * t680 + t956 * t681;
t641 = -pkin(6) * t699 - (pkin(3) * t955 - pkin(7) * t959) * t728;
t638 = -pkin(7) * t731 - t954 * t685 + t958 * t717;
t635 = -pkin(5) * t754 - t956 * t677 + t960 * t687;
t634 = -pkin(5) * t753 - t956 * t676 + t960 * t686;
t633 = -pkin(7) * t714 - t954 * t679 + t958 * t712;
t632 = -pkin(1) * t866 + pkin(5) * t756 + t960 * t677 + t956 * t687;
t631 = -pkin(1) * t864 + pkin(5) * t755 + t960 * t676 + t956 * t686;
t628 = pkin(6) * t700 - (-pkin(3) * t959 - pkin(7) * t955 - pkin(2)) * t728;
t627 = -t956 * t667 + t960 * t668;
t626 = t960 * t667 + t956 * t668;
t625 = -t956 * t665 + t960 * t666;
t624 = t960 * t665 + t956 * t666;
t623 = -pkin(1) * t659 - t978;
t621 = t954 * t655 + t1013;
t619 = t959 * t622 + t955 * t760;
t617 = -pkin(6) * t688 + t959 * t638 - t955 * t657;
t616 = -pkin(1) * t646 - t967;
t615 = -pkin(3) * t621 - t1032;
t614 = -pkin(6) * t680 + t959 * t633 - t955 * t644;
t613 = -pkin(2) * t731 + pkin(6) * t689 + t955 * t638 + t959 * t657;
t612 = -pkin(1) * t642 - t968;
t611 = -pkin(2) * t714 + pkin(6) * t681 + t955 * t633 + t959 * t644;
t610 = -pkin(7) * t671 - t954 * t630 + t958 * t637;
t609 = -pkin(5) * t659 - t956 * t628 + t960 * t641;
t608 = pkin(1) * t728 + pkin(5) * t660 + t960 * t628 + t956 * t641;
t607 = -pkin(7) * t621 - pkin(8) * t1013 - t954 * t640;
t606 = -t956 * t618 + t960 * t619;
t605 = t960 * t618 + t956 * t619;
t604 = -pkin(6) * t665 + t959 * t610 - t955 * t645;
t603 = -pkin(2) * t671 + pkin(6) * t666 + t955 * t610 + t959 * t645;
t602 = -pkin(1) * t624 - t969;
t601 = -pkin(5) * t646 - t956 * t613 + t960 * t617;
t600 = -pkin(1) * t731 + pkin(5) * t647 + t960 * t613 + t956 * t617;
t599 = -pkin(5) * t642 - t956 * t611 + t960 * t614;
t598 = -pkin(1) * t714 + pkin(5) * t643 + t960 * t611 + t956 * t614;
t597 = -pkin(6) * t618 + t959 * t607 - t955 * t615;
t596 = -pkin(1) * t605 - t964;
t595 = -pkin(2) * t621 + pkin(6) * t619 + t955 * t607 + t959 * t615;
t594 = -pkin(5) * t624 - t956 * t603 + t960 * t604;
t593 = -pkin(1) * t671 + pkin(5) * t625 + t960 * t603 + t956 * t604;
t592 = -pkin(5) * t605 - t956 * t595 + t960 * t597;
t591 = -pkin(1) * t621 + pkin(5) * t606 + t960 * t595 + t956 * t597;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1036, -t1035, t1037, qJ(1) * t1037, 0, 0, t952 * t925, 0, -t952 * t924, t951 * qJDD(2), t990 * t951 + t952 * t976, t952 * t882 + t991 * t951, t952 * t827, -qJ(1) * (t951 * t988 + t915) - (t951 * pkin(1) - t952 * pkin(5)) * t827, 0, 0, -t952 * t851, 0, -t952 * t846, t951 * t947, t1048 * t952 - t951 * t992, t952 * t764 - t951 * t993, t952 * t703, t952 * t658 - t951 * t674 - qJ(1) * (t1045 * t951 + t915), t952 * t770 + t951 * t856, t952 * t738 + t951 * t843, t952 * t776 + t951 * t865, t952 * t769 - t951 * t855, t952 * t775 + t951 * t863, t952 * t815, t952 * t634 - t951 * t682 - qJ(1) * (t951 * t755 - t952 * t864), t952 * t635 - t951 * t683 - qJ(1) * (t951 * t756 - t952 * t866), t952 * t649 - t951 * t994, t952 * t609 - t951 * t623 - qJ(1) * (t951 * t660 + t728 * t952), t952 * t664 + t951 * t719, t952 * t627 + t951 * t670, t952 * t652 + t951 * t733, t952 * t663 + t951 * t718, t952 * t653 + t951 * t734, t952 * t691 + t951 * t749, t952 * t599 - t951 * t612 - qJ(1) * (t951 * t643 - t952 * t714), t952 * t601 - t951 * t616 - qJ(1) * (t951 * t647 - t952 * t731), t952 * t594 - t951 * t602 - qJ(1) * (t951 * t625 - t952 * t671), t952 * t592 - t951 * t596 - qJ(1) * (t951 * t606 - t952 * t621); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1035, -t1036, t985, qJ(1) * t985, 0, 0, t951 * t925, 0, -t951 * t924, -t952 * qJDD(2), t951 * t976 - t990 * t952, t951 * t882 - t991 * t952, t951 * t827, qJ(1) * (t952 * t988 - t1025) - (-t952 * pkin(1) - t951 * pkin(5)) * t827, 0, 0, -t951 * t851, 0, -t951 * t846, -t952 * t947, t1048 * t951 + t952 * t992, t951 * t764 + t952 * t993, t951 * t703, t951 * t658 + t952 * t674 + qJ(1) * (t1045 * t952 - t1025), t951 * t770 - t952 * t856, t951 * t738 - t952 * t843, t951 * t776 - t952 * t865, t951 * t769 + t952 * t855, t951 * t775 - t952 * t863, t951 * t815, t951 * t634 + t952 * t682 + qJ(1) * (t952 * t755 + t951 * t864), t951 * t635 + t952 * t683 + qJ(1) * (t952 * t756 + t951 * t866), t951 * t649 + t952 * t994, t951 * t609 + t952 * t623 + qJ(1) * (t952 * t660 - t728 * t951), t951 * t664 - t952 * t719, t951 * t627 - t952 * t670, t951 * t652 - t952 * t733, t951 * t663 - t952 * t718, t951 * t653 - t952 * t734, t951 * t691 - t952 * t749, t951 * t599 + t952 * t612 + qJ(1) * (t952 * t643 + t951 * t714), t951 * t601 + t952 * t616 + qJ(1) * (t952 * t647 + t951 * t731), t951 * t594 + t952 * t602 + qJ(1) * (t952 * t625 + t951 * t671), t951 * t592 + t952 * t596 + qJ(1) * (t952 * t606 + t951 * t621); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t926, t1004, 0, 0, 0, 0, t924, 0, t925, 0, -t882, t976, t988, t808, 0, 0, t846, 0, -t851, 0, -t764, t1048, t1045, t656, t768, t737, t774, t767, t773, t814, t631, t632, t648, t608, t662, t626, t650, t661, t651, t690, t598, t600, t593, t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1006, -t926, 0, 0, 0, t925, 0, -t924, 0, t976, t882, t827, pkin(5) * t827, 0, 0, -t851, 0, -t846, 0, t1048, t764, t703, t658, t770, t738, t776, t769, t775, t815, t634, t635, t649, t609, t664, t627, t652, t663, t653, t691, t599, t601, t594, t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1006, 0, -t1004, 0, 0, 0, 0, 0, 0, -qJDD(2), t854, t853, 0, pkin(1) * t827, 0, 0, 0, 0, 0, -t947, t759, t758, 0, t674, -t856, -t843, -t865, t855, -t863, 0, t682, t683, t675, t623, -t719, -t670, -t733, -t718, -t734, -t749, t612, t616, t602, t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t926, t1004, 0, 0, 0, 0, t924, 0, t925, 0, -t882, t976, t988, t808, 0, 0, t846, 0, -t851, 0, -t764, t1048, t1045, t656, t768, t737, t774, t767, t773, t814, t631, t632, t648, t608, t662, t626, t650, t661, t651, t690, t598, t600, t593, t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t962, 0, 0, -t926, t889, 0, 0, 0, -t912, 0, -t909, 0, t1042, t862, t747, pkin(6) * t747, t833, t812, t840, t832, t839, t878, t686, t687, t710, t641, t708, t668, t695, t707, t696, t741, t614, t617, t604, t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t962, 0, qJDD(2), 0, t926, 0, t890, 0, 0, 0, t909, 0, -t912, 0, -t862, t1042, t989, t739, t831, t811, t838, t830, t837, t877, t676, t677, t709, t628, t706, t667, t693, t705, t694, t740, t611, t613, t603, t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t889, -t890, 0, 0, 0, 0, 0, 0, 0, t947, t973, t977, 0, -t1033, t856, t843, t865, -t855, t863, 0, t974, t975, t966, t978, t719, t670, t733, t718, t734, t749, t968, t967, t969, t964; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t947, 0, -t945, 0, 0, -t926, t817, 0, t873, t844, t869, t872, t867, t902, t765, t766, t728, pkin(7) * t728, t721, t672, t735, t720, t736, t750, t633, t638, t610, t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t945, 0, t947, 0, t926, 0, t818, 0, t931, -t914, -t1019, -t931, -t1010, -qJDD(4), t751, t752, 0, pkin(3) * t728, -t842, -t841, -t788, t842, t784, -t946, t644, t657, t645, t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t947, -t817, -t818, 0, 0, t856, t843, t865, -t855, t863, 0, t997, t998, t979, t1005, t719, t670, t733, t718, t734, t749, t981, t980, t982, t970; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t905, t906, t920, -t934, t929, t934, 0, t805, t791, 0, t780, t724, t799, t778, t800, t820, t712, t717, t637, -pkin(8) * t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t996, t904, t927, -t971, t921, -t996, -t805, 0, t792, 0, t779, t722, t797, t777, t798, t819, t679, t685, t630, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t931, t914, t1019, t931, t1010, qJDD(4), -t791, -t792, 0, 0, t842, t841, t788, -t842, -t784, t946, t972, t965, t1031, t1032; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t810, -t783, t1039, t886, t883, -t886, 0, t760, t697, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1030, t1038, t884, t809, t835, -t1030, -t760, 0, t698, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t842, t841, t788, -t842, -t784, t946, -t697, -t698, 0, 0;];
m_new_reg = t1;