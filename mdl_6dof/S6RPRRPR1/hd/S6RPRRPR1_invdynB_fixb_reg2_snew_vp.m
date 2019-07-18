% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRRPR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:57:00
% EndTime: 2019-05-05 21:57:31
% DurationCPUTime: 31.87s
% Computational Cost: add. (220113->879), mult. (476008->1358), div. (0->0), fcn. (342112->12), ass. (0->601)
t1041 = qJD(1) ^ 2;
t962 = sin(qJ(1));
t966 = cos(qJ(1));
t926 = g(1) * t966 + g(2) * t962;
t911 = -pkin(1) * t1041 - t926;
t955 = sin(pkin(10));
t957 = cos(pkin(10));
t925 = g(1) * t962 - t966 * g(2);
t973 = qJDD(1) * pkin(1) + t925;
t860 = t955 * t911 - t957 * t973;
t861 = t957 * t911 + t955 * t973;
t788 = t860 * t957 - t861 * t955;
t1026 = t788 * t962;
t982 = t860 * t955 + t957 * t861;
t717 = t966 * t982 + t1026;
t1025 = t788 * t966;
t1063 = -t962 * t982 + t1025;
t915 = qJDD(1) * t955 + t1041 * t957;
t916 = qJDD(1) * t957 - t1041 * t955;
t864 = -t915 * t962 + t966 * t916;
t952 = g(3) - qJDD(2);
t889 = qJ(2) * t915 - t952 * t957;
t975 = -qJ(2) * t916 - t952 * t955;
t1062 = -pkin(6) * t864 + t889 * t962 + t966 * t975;
t961 = sin(qJ(3));
t1004 = qJD(1) * t961;
t960 = sin(qJ(4));
t964 = cos(qJ(4));
t965 = cos(qJ(3));
t902 = qJD(1) * t964 * t965 - t1004 * t960;
t974 = t960 * t965 + t961 * t964;
t903 = t974 * qJD(1);
t862 = t902 * t903;
t949 = qJDD(3) + qJDD(4);
t1048 = t862 + t949;
t1061 = t1048 * t960;
t1060 = t1048 * t964;
t954 = sin(pkin(11));
t956 = cos(pkin(11));
t841 = -t956 * t902 + t903 * t954;
t843 = t954 * t902 + t956 * t903;
t779 = t843 * t841;
t1049 = -t779 + t949;
t1059 = t1049 * t954;
t1058 = t1049 * t956;
t950 = qJD(3) + qJD(4);
t959 = sin(qJ(6));
t963 = cos(qJ(6));
t813 = t843 * t959 - t963 * t950;
t815 = t843 * t963 + t950 * t959;
t757 = t815 * t813;
t999 = qJD(1) * qJD(3);
t987 = t965 * t999;
t998 = qJDD(1) * t961;
t913 = t987 + t998;
t988 = t961 * t999;
t997 = qJDD(1) * t965;
t972 = t988 - t997;
t981 = t960 * t913 + t964 * t972;
t820 = -qJD(4) * t903 - t981;
t821 = t902 * qJD(4) + t964 * t913 - t960 * t972;
t983 = -t956 * t820 + t821 * t954;
t976 = qJDD(6) + t983;
t1050 = -t757 + t976;
t1057 = t1050 * t959;
t1056 = t1050 * t963;
t1020 = t843 * t950;
t721 = t983 + t1020;
t1046 = t966 * t915 + t916 * t962;
t1055 = pkin(6) * t1046 + t889 * t966 - t962 * t975;
t755 = t820 * t954 + t821 * t956;
t833 = t950 * t841;
t725 = t755 - t833;
t704 = -t813 * qJD(6) + t963 * t755 + t959 * t949;
t838 = qJD(6) + t841;
t767 = t838 * t813;
t671 = -t767 + t704;
t891 = t950 * t902;
t798 = t891 - t821;
t1047 = t891 + t821;
t984 = t959 * t755 - t963 * t949;
t668 = (qJD(6) - t838) * t815 + t984;
t794 = (qJD(4) - t950) * t903 + t981;
t845 = -pkin(2) * t1041 + qJDD(1) * pkin(7) + t861;
t817 = t961 * t845 + t965 * t952;
t1045 = -t817 + (-t913 + t987) * pkin(8);
t810 = t813 ^ 2;
t811 = t815 ^ 2;
t837 = t838 ^ 2;
t839 = t841 ^ 2;
t840 = t843 ^ 2;
t897 = t902 ^ 2;
t898 = t903 ^ 2;
t1043 = t950 ^ 2;
t1042 = t965 ^ 2;
t1040 = 2 * qJD(5);
t1039 = pkin(5) * t954;
t1002 = t965 * t1041;
t934 = t961 * t1002;
t922 = qJDD(3) + t934;
t969 = t922 * pkin(3) + t1045;
t946 = t1042 * t1041;
t819 = t965 * t845 - t961 * t952;
t924 = qJD(3) * pkin(3) - pkin(8) * t1004;
t971 = -pkin(8) * t972 - qJD(3) * t924 + t819;
t970 = -pkin(3) * t946 + t971;
t706 = t960 * t969 + t964 * t970;
t879 = pkin(4) * t950 - qJ(5) * t903;
t667 = -pkin(4) * t897 + qJ(5) * t820 - t879 * t950 + t706;
t705 = t960 * t970 - t964 * t969;
t968 = pkin(4) * t1048 + qJ(5) * t798 - t705;
t986 = t954 * t667 - t956 * t968;
t581 = t1040 * t843 + t986;
t582 = -0.2e1 * qJD(5) * t841 + t956 * t667 + t954 * t968;
t526 = -t581 * t956 + t582 * t954;
t1038 = t526 * t960;
t1037 = t526 * t964;
t774 = pkin(5) * t841 - pkin(9) * t843;
t552 = -t949 * pkin(5) - t1043 * pkin(9) + (t1040 + t774) * t843 + t986;
t1036 = t552 * t959;
t1035 = t552 * t963;
t630 = -t705 * t964 + t706 * t960;
t1034 = t630 * t961;
t1033 = t630 * t965;
t687 = t757 + t976;
t1032 = t687 * t959;
t1031 = t687 * t963;
t844 = -qJDD(1) * pkin(2) - t1041 * pkin(7) + t860;
t790 = t972 * pkin(3) - pkin(8) * t946 + t1004 * t924 + t844;
t699 = -t820 * pkin(4) - t897 * qJ(5) + t903 * t879 + qJDD(5) + t790;
t1030 = t699 * t954;
t1029 = t699 * t956;
t772 = t779 + t949;
t1028 = t772 * t954;
t1027 = t772 * t956;
t1024 = t790 * t960;
t1023 = t790 * t964;
t1022 = t838 * t959;
t1021 = t838 * t963;
t1019 = t844 * t961;
t1018 = t844 * t965;
t857 = -t862 + t949;
t1017 = t857 * t960;
t1016 = t857 * t964;
t914 = -0.2e1 * t988 + t997;
t1015 = t914 * t965;
t1012 = t922 * t961;
t923 = qJDD(3) - t934;
t1011 = t923 * t961;
t1010 = t923 * t965;
t1009 = t949 * t955;
t1008 = t950 * t954;
t1007 = t950 * t956;
t1006 = t950 * t960;
t1005 = t950 * t964;
t553 = -pkin(5) * t1043 + pkin(9) * t949 - t774 * t841 + t582;
t605 = pkin(5) * t721 - pkin(9) * t725 + t699;
t525 = t963 * t553 + t959 * t605;
t951 = t961 ^ 2;
t1003 = t1041 * t951;
t996 = t951 + t1042;
t995 = t954 * t757;
t994 = t956 * t757;
t993 = t955 * t779;
t992 = t957 * t779;
t991 = t955 * t862;
t990 = t957 * t862;
t989 = -pkin(5) * t956 - pkin(4);
t524 = t553 * t959 - t963 * t605;
t527 = t581 * t954 + t956 * t582;
t631 = t705 * t960 + t964 * t706;
t753 = t817 * t961 + t965 * t819;
t871 = -t925 * t962 - t966 * t926;
t979 = t955 * t934;
t978 = t957 * t934;
t919 = qJDD(1) * t966 - t1041 * t962;
t977 = -pkin(6) * t919 - g(3) * t962;
t471 = -t524 * t963 + t525 * t959;
t472 = t524 * t959 + t525 * t963;
t752 = t817 * t965 - t819 * t961;
t870 = t925 * t966 - t926 * t962;
t722 = t983 - t1020;
t967 = qJD(3) ^ 2;
t932 = t957 * t949;
t931 = -t946 - t967;
t930 = t946 - t967;
t929 = -t967 - t1003;
t928 = t967 - t1003;
t921 = t946 - t1003;
t920 = t946 + t1003;
t918 = qJDD(1) * t962 + t1041 * t966;
t917 = t996 * qJDD(1);
t912 = 0.2e1 * t987 + t998;
t909 = t965 * t922;
t908 = t996 * t999;
t896 = -pkin(6) * t918 + g(3) * t966;
t885 = -t898 + t1043;
t884 = t897 - t1043;
t883 = t913 * t965 - t951 * t999;
t882 = -t1042 * t999 + t961 * t972;
t881 = qJDD(3) * t955 + t908 * t957;
t880 = -qJDD(3) * t957 + t908 * t955;
t878 = -t898 - t1043;
t877 = -t929 * t961 - t1010;
t876 = -t928 * t961 + t909;
t875 = t931 * t965 - t1012;
t874 = t930 * t965 - t1011;
t873 = t929 * t965 - t1011;
t872 = t931 * t961 + t909;
t869 = t917 * t957 - t920 * t955;
t868 = t917 * t955 + t920 * t957;
t863 = -t912 * t961 + t1015;
t859 = -t898 + t897;
t854 = -t1043 - t897;
t853 = t883 * t957 - t979;
t852 = t882 * t957 + t979;
t851 = t883 * t955 + t978;
t850 = t882 * t955 - t978;
t849 = t876 * t957 + t955 * t998;
t848 = t874 * t957 + t955 * t997;
t847 = t876 * t955 - t957 * t998;
t846 = t874 * t955 - t957 * t997;
t831 = t877 * t957 + t912 * t955;
t830 = t875 * t957 - t914 * t955;
t829 = t877 * t955 - t912 * t957;
t828 = t875 * t955 + t914 * t957;
t827 = -t840 + t1043;
t826 = t839 - t1043;
t825 = (t902 * t964 + t903 * t960) * t950;
t824 = (t902 * t960 - t903 * t964) * t950;
t823 = -t840 - t1043;
t822 = -t897 - t898;
t818 = t863 * t957 - t921 * t955;
t816 = t863 * t955 + t921 * t957;
t806 = t884 * t964 - t1017;
t805 = -t885 * t960 + t1060;
t804 = t884 * t960 + t1016;
t803 = t885 * t964 + t1061;
t802 = -t868 * t962 + t869 * t966;
t801 = t868 * t966 + t869 * t962;
t800 = -t878 * t960 - t1016;
t799 = t878 * t964 - t1017;
t793 = (qJD(4) + t950) * t903 + t981;
t792 = -pkin(7) * t873 + t1018;
t791 = -pkin(7) * t872 + t1019;
t785 = -t1006 * t903 + t821 * t964;
t784 = t1005 * t903 + t821 * t960;
t783 = -t1005 * t902 - t820 * t960;
t782 = -t1006 * t902 + t820 * t964;
t781 = -pkin(2) * t873 + t819;
t780 = -pkin(2) * t872 + t817;
t778 = t854 * t964 - t1061;
t777 = t854 * t960 + t1060;
t776 = -t840 + t839;
t775 = pkin(1) * t952 + qJ(2) * t982;
t769 = -t1043 - t839;
t766 = -t811 + t837;
t765 = t810 - t837;
t764 = (-t841 * t956 + t843 * t954) * t950;
t763 = (-t841 * t954 - t843 * t956) * t950;
t762 = -t829 * t962 + t831 * t966;
t761 = -t828 * t962 + t830 * t966;
t760 = t829 * t966 + t831 * t962;
t759 = t828 * t966 + t830 * t962;
t758 = -t824 * t961 + t825 * t965;
t756 = -t811 + t810;
t748 = t758 * t957 + t1009;
t747 = t758 * t955 - t932;
t746 = -t839 - t840;
t745 = -t811 - t837;
t744 = -t804 * t961 + t806 * t965;
t743 = -t803 * t961 + t805 * t965;
t742 = t826 * t956 - t1028;
t741 = -t827 * t954 + t1058;
t740 = t826 * t954 + t1027;
t739 = t827 * t956 + t1059;
t738 = -t823 * t954 - t1027;
t737 = t823 * t956 - t1028;
t736 = -qJ(2) * t868 + t752 * t957;
t735 = qJ(2) * t869 + t752 * t955;
t734 = -t837 - t810;
t733 = -t799 * t961 + t800 * t965;
t732 = t799 * t965 + t800 * t961;
t731 = -t794 * t964 - t798 * t960;
t730 = -t1047 * t960 - t793 * t964;
t729 = -t794 * t960 + t798 * t964;
t728 = t1047 * t964 - t793 * t960;
t727 = -pkin(8) * t799 + t1023;
t726 = -t755 - t833;
t720 = t810 + t811;
t719 = t753 * t957 + t844 * t955;
t718 = t753 * t955 - t844 * t957;
t715 = -t1008 * t843 + t755 * t956;
t714 = t1007 * t843 + t755 * t954;
t713 = t1007 * t841 + t954 * t983;
t712 = t1008 * t841 - t956 * t983;
t711 = -pkin(8) * t777 + t1024;
t710 = -t784 * t961 + t785 * t965;
t709 = -t782 * t961 + t783 * t965;
t708 = -t777 * t961 + t778 * t965;
t707 = t777 * t965 + t778 * t961;
t703 = -qJD(6) * t815 - t984;
t701 = t769 * t956 - t1059;
t700 = t769 * t954 + t1058;
t698 = (-t813 * t963 + t815 * t959) * t838;
t697 = (t813 * t959 + t815 * t963) * t838;
t696 = -t763 * t960 + t764 * t964;
t695 = t763 * t964 + t764 * t960;
t694 = t710 * t957 - t991;
t693 = t709 * t957 + t991;
t692 = t710 * t955 + t990;
t691 = t709 * t955 - t990;
t690 = -qJ(2) * t829 - t781 * t955 + t792 * t957;
t689 = -qJ(2) * t828 - t780 * t955 + t791 * t957;
t685 = t744 * t957 - t794 * t955;
t684 = t743 * t957 - t798 * t955;
t683 = t744 * t955 + t794 * t957;
t682 = t743 * t955 + t798 * t957;
t681 = -pkin(3) * t1047 + pkin(8) * t800 + t1024;
t680 = t1047 * t955 + t733 * t957;
t679 = -t1047 * t957 + t733 * t955;
t678 = -pkin(1) * t873 + qJ(2) * t831 + t781 * t957 + t792 * t955;
t677 = -pkin(1) * t872 + qJ(2) * t830 + t780 * t957 + t791 * t955;
t676 = -pkin(3) * t793 + pkin(8) * t778 - t1023;
t675 = t708 * t957 + t793 * t955;
t674 = t708 * t955 - t793 * t957;
t672 = -t767 - t704;
t669 = (-qJD(6) - t838) * t815 - t984;
t666 = -t1022 * t815 + t704 * t963;
t665 = -t1021 * t815 - t704 * t959;
t664 = t1021 * t813 - t703 * t959;
t663 = -t1022 * t813 - t703 * t963;
t661 = -t740 * t960 + t742 * t964;
t660 = -t739 * t960 + t741 * t964;
t659 = t740 * t964 + t742 * t960;
t658 = t739 * t964 + t741 * t960;
t655 = -t737 * t960 + t738 * t964;
t654 = t737 * t964 + t738 * t960;
t653 = -t729 * t961 + t731 * t965;
t652 = -t728 * t961 + t730 * t965;
t651 = t729 * t965 + t731 * t961;
t650 = t698 * t956 + t954 * t976;
t649 = t698 * t954 - t956 * t976;
t648 = -t722 * t956 - t726 * t954;
t647 = -t721 * t956 - t725 * t954;
t646 = -t722 * t954 + t726 * t956;
t645 = -t721 * t954 + t725 * t956;
t644 = t765 * t963 - t1032;
t643 = -t766 * t959 + t1056;
t642 = -t765 * t959 - t1031;
t641 = -t766 * t963 - t1057;
t640 = -t718 * t962 + t719 * t966;
t639 = t718 * t966 + t719 * t962;
t638 = -t714 * t960 + t715 * t964;
t637 = -t712 * t960 + t713 * t964;
t636 = t714 * t964 + t715 * t960;
t635 = t712 * t964 + t713 * t960;
t634 = -qJ(5) * t737 + t1029;
t633 = t652 * t957 - t859 * t955;
t632 = t652 * t955 + t859 * t957;
t629 = -t745 * t959 - t1031;
t628 = t745 * t963 - t1032;
t627 = t653 * t957 + t822 * t955;
t626 = t653 * t955 - t822 * t957;
t625 = -t700 * t960 + t701 * t964;
t624 = t700 * t964 + t701 * t960;
t623 = t734 * t963 - t1057;
t622 = t734 * t959 + t1056;
t621 = -qJ(5) * t700 + t1030;
t620 = t666 * t956 + t995;
t619 = t664 * t956 - t995;
t618 = t666 * t954 - t994;
t617 = t664 * t954 + t994;
t616 = -qJ(2) * t718 - (pkin(2) * t955 - pkin(7) * t957) * t752;
t615 = -pkin(2) * t732 - pkin(3) * t799 + t706;
t614 = -t695 * t961 + t696 * t965;
t613 = t614 * t957 + t1009;
t612 = t614 * t955 - t932;
t611 = -pkin(2) * t707 + t960 * t971 - t964 * t1045 + (-t964 * qJDD(3) - t1002 * t974 - t777) * pkin(3);
t610 = -pkin(3) * t790 + pkin(8) * t631;
t609 = -pkin(2) * t651 - pkin(3) * t729;
t608 = -t679 * t962 + t680 * t966;
t607 = t679 * t966 + t680 * t962;
t606 = -pkin(4) * t725 + qJ(5) * t738 + t1030;
t604 = qJ(2) * t719 - (-pkin(2) * t957 - pkin(7) * t955 - pkin(1)) * t752;
t601 = -pkin(8) * t729 - t630;
t600 = -t674 * t962 + t675 * t966;
t599 = t674 * t966 + t675 * t962;
t598 = -pkin(7) * t732 - t681 * t961 + t727 * t965;
t597 = -pkin(4) * t721 + qJ(5) * t701 - t1029;
t596 = -t668 * t963 - t672 * t959;
t595 = t669 * t963 - t671 * t959;
t594 = -t668 * t959 + t672 * t963;
t593 = -t669 * t959 - t671 * t963;
t592 = -pkin(3) * t822 + pkin(8) * t731 + t631;
t591 = -pkin(7) * t707 - t676 * t961 + t711 * t965;
t590 = -t659 * t961 + t661 * t965;
t589 = -t658 * t961 + t660 * t965;
t588 = t644 * t956 - t668 * t954;
t587 = t643 * t956 - t672 * t954;
t586 = t644 * t954 + t668 * t956;
t585 = t643 * t954 + t672 * t956;
t584 = -t654 * t961 + t655 * t965;
t583 = t654 * t965 + t655 * t961;
t579 = -t649 * t960 + t650 * t964;
t578 = t649 * t964 + t650 * t960;
t577 = t629 * t956 + t671 * t954;
t576 = t629 * t954 - t671 * t956;
t575 = -t646 * t960 + t648 * t964;
t574 = -t645 * t960 + t647 * t964;
t573 = t646 * t964 + t648 * t960;
t572 = t645 * t964 + t647 * t960;
t571 = t623 * t956 - t669 * t954;
t570 = t623 * t954 + t669 * t956;
t569 = t595 * t956 - t756 * t954;
t568 = t595 * t954 + t756 * t956;
t567 = -t636 * t961 + t638 * t965;
t566 = -t635 * t961 + t637 * t965;
t565 = t596 * t956 - t720 * t954;
t564 = t596 * t954 + t720 * t956;
t563 = t631 * t965 - t1034;
t562 = t631 * t961 + t1033;
t561 = -t626 * t962 + t627 * t966;
t560 = t626 * t966 + t627 * t962;
t559 = -t624 * t961 + t625 * t965;
t558 = t624 * t965 + t625 * t961;
t557 = t590 * t957 - t722 * t955;
t556 = t589 * t957 - t726 * t955;
t555 = t590 * t955 + t722 * t957;
t554 = t589 * t955 + t726 * t957;
t550 = t567 * t957 + t993;
t549 = t566 * t957 - t993;
t548 = t567 * t955 - t992;
t547 = t566 * t955 + t992;
t546 = t584 * t957 + t725 * t955;
t545 = t584 * t955 - t725 * t957;
t544 = -t618 * t960 + t620 * t964;
t543 = -t617 * t960 + t619 * t964;
t542 = t618 * t964 + t620 * t960;
t541 = t617 * t964 + t619 * t960;
t540 = t563 * t957 + t790 * t955;
t539 = t563 * t955 - t790 * t957;
t538 = t559 * t957 + t721 * t955;
t537 = t559 * t955 - t721 * t957;
t536 = -pkin(2) * t562 - pkin(3) * t630;
t535 = -pkin(9) * t628 + t1035;
t534 = -pkin(9) * t622 + t1036;
t533 = -pkin(8) * t654 - t606 * t960 + t634 * t964;
t532 = -qJ(2) * t679 + t598 * t957 - t615 * t955;
t531 = -t586 * t960 + t588 * t964;
t530 = -t585 * t960 + t587 * t964;
t529 = t586 * t964 + t588 * t960;
t528 = t585 * t964 + t587 * t960;
t523 = -qJ(2) * t674 + t591 * t957 - t611 * t955;
t522 = -t578 * t961 + t579 * t965;
t521 = -pkin(8) * t624 - t597 * t960 + t621 * t964;
t520 = -t576 * t960 + t577 * t964;
t519 = t576 * t964 + t577 * t960;
t518 = -t573 * t961 + t575 * t965;
t517 = -t572 * t961 + t574 * t965;
t516 = t573 * t965 + t575 * t961;
t515 = -pkin(3) * t725 + pkin(8) * t655 + t606 * t964 + t634 * t960;
t514 = -t570 * t960 + t571 * t964;
t513 = t570 * t964 + t571 * t960;
t512 = -pkin(1) * t732 + qJ(2) * t680 + t598 * t955 + t615 * t957;
t511 = -t568 * t960 + t569 * t964;
t510 = t568 * t964 + t569 * t960;
t509 = t517 * t957 - t776 * t955;
t508 = t517 * t955 + t776 * t957;
t507 = -pkin(7) * t651 - t592 * t961 + t601 * t965;
t506 = t518 * t957 + t746 * t955;
t505 = t518 * t955 - t746 * t957;
t504 = -pkin(4) * t699 + qJ(5) * t527;
t503 = -t564 * t960 + t565 * t964;
t502 = t564 * t964 + t565 * t960;
t501 = -pkin(3) * t721 + pkin(8) * t625 + t597 * t964 + t621 * t960;
t500 = -pkin(1) * t707 + qJ(2) * t675 + t591 * t955 + t611 * t957;
t499 = t522 * t957 - t697 * t955;
t498 = t522 * t955 + t697 * t957;
t497 = -pkin(7) * t562 - pkin(8) * t1033 - t610 * t961;
t496 = -t545 * t962 + t546 * t966;
t495 = t545 * t966 + t546 * t962;
t494 = -qJ(5) * t646 - t526;
t493 = -t542 * t961 + t544 * t965;
t492 = -t541 * t961 + t543 * t965;
t491 = -t539 * t962 + t540 * t966;
t490 = t539 * t966 + t540 * t962;
t489 = -pkin(5) * t628 + t525;
t488 = -pkin(5) * t622 + t524;
t487 = -pkin(2) * t583 - pkin(3) * t654 - pkin(4) * t737 + t582;
t486 = -pkin(4) * t746 + qJ(5) * t648 + t527;
t485 = -t537 * t962 + t538 * t966;
t484 = t537 * t966 + t538 * t962;
t483 = t493 * t957 - t665 * t955;
t482 = t492 * t957 - t663 * t955;
t481 = t493 * t955 + t665 * t957;
t480 = t492 * t955 + t663 * t957;
t479 = -pkin(2) * t558 - pkin(3) * t624 - pkin(4) * t700 + t581;
t478 = -qJ(2) * t626 + t507 * t957 - t609 * t955;
t477 = -pkin(2) * t516 - pkin(3) * t573 - pkin(4) * t646;
t476 = -t529 * t961 + t531 * t965;
t475 = -t528 * t961 + t530 * t965;
t474 = t527 * t964 - t1038;
t473 = t527 * t960 + t1037;
t470 = -pkin(1) * t651 + qJ(2) * t627 + t507 * t955 + t609 * t957;
t469 = -t519 * t961 + t520 * t965;
t468 = t519 * t965 + t520 * t961;
t467 = -t513 * t961 + t514 * t965;
t466 = t513 * t965 + t514 * t961;
t465 = -t510 * t961 + t511 * t965;
t464 = t476 * t957 - t642 * t955;
t463 = t475 * t957 - t641 * t955;
t462 = t476 * t955 + t642 * t957;
t461 = t475 * t955 + t641 * t957;
t460 = -t505 * t962 + t506 * t966;
t459 = t505 * t966 + t506 * t962;
t458 = -t502 * t961 + t503 * t965;
t457 = t502 * t965 + t503 * t961;
t456 = t469 * t957 + t628 * t955;
t455 = t469 * t955 - t628 * t957;
t454 = t467 * t957 + t622 * t955;
t453 = t467 * t955 - t622 * t957;
t452 = -pkin(7) * t583 - t515 * t961 + t533 * t965;
t451 = -pkin(9) * t594 - t471;
t450 = t465 * t957 - t593 * t955;
t449 = t465 * t955 + t593 * t957;
t448 = t472 * t956 + t552 * t954;
t447 = t472 * t954 - t552 * t956;
t446 = -qJ(5) * t576 - t489 * t954 + t535 * t956;
t445 = -qJ(5) * t570 - t488 * t954 + t534 * t956;
t444 = -qJ(2) * t539 + t497 * t957 - t536 * t955;
t443 = t458 * t957 + t594 * t955;
t442 = t458 * t955 - t594 * t957;
t441 = -pkin(7) * t558 - t501 * t961 + t521 * t965;
t440 = -pkin(4) * t628 + qJ(5) * t577 + t489 * t956 + t535 * t954;
t439 = -pkin(4) * t622 + qJ(5) * t571 + t488 * t956 + t534 * t954;
t438 = -pkin(8) * t573 - t486 * t960 + t494 * t964;
t437 = -pkin(3) * t746 + pkin(8) * t575 + t486 * t964 + t494 * t960;
t436 = -pkin(1) * t562 + qJ(2) * t540 + t497 * t955 + t536 * t957;
t435 = -qJ(5) * t564 + t1039 * t594 + t451 * t956;
t434 = qJ(5) * t565 + t954 * t451 + t594 * t989;
t433 = -t473 * t961 + t474 * t965;
t432 = t473 * t965 + t474 * t961;
t431 = -pkin(8) * t473 - qJ(5) * t1037 - t504 * t960;
t430 = t433 * t957 + t699 * t955;
t429 = t433 * t955 - t699 * t957;
t428 = -pkin(3) * t699 + pkin(8) * t474 - qJ(5) * t1038 + t504 * t964;
t427 = -qJ(2) * t545 + t452 * t957 - t487 * t955;
t426 = -pkin(2) * t468 - pkin(3) * t519 - pkin(4) * t576 + pkin(5) * t671 - pkin(9) * t629 - t1036;
t425 = -t455 * t962 + t456 * t966;
t424 = t455 * t966 + t456 * t962;
t423 = -pkin(2) * t466 - pkin(3) * t513 - pkin(4) * t570 - pkin(5) * t669 - pkin(9) * t623 + t1035;
t422 = -t453 * t962 + t454 * t966;
t421 = t453 * t966 + t454 * t962;
t420 = -pkin(1) * t583 + qJ(2) * t546 + t452 * t955 + t487 * t957;
t419 = -qJ(2) * t537 + t441 * t957 - t479 * t955;
t418 = -t447 * t960 + t448 * t964;
t417 = t447 * t964 + t448 * t960;
t416 = -t442 * t962 + t443 * t966;
t415 = t442 * t966 + t443 * t962;
t414 = -pkin(1) * t558 + qJ(2) * t538 + t441 * t955 + t479 * t957;
t413 = -pkin(2) * t457 - pkin(3) * t502 - pkin(4) * t564 - pkin(5) * t720 - pkin(9) * t596 - t472;
t412 = -qJ(5) * t447 + (-pkin(9) * t956 + t1039) * t471;
t411 = -pkin(8) * t519 - t440 * t960 + t446 * t964;
t410 = -pkin(8) * t513 - t439 * t960 + t445 * t964;
t409 = -pkin(2) * t432 - pkin(3) * t473 - pkin(4) * t526;
t408 = -pkin(3) * t628 + pkin(8) * t520 + t440 * t964 + t446 * t960;
t407 = -pkin(3) * t622 + pkin(8) * t514 + t439 * t964 + t445 * t960;
t406 = -pkin(7) * t516 - t437 * t961 + t438 * t965;
t405 = -t429 * t962 + t430 * t966;
t404 = t429 * t966 + t430 * t962;
t403 = -pkin(8) * t502 - t434 * t960 + t435 * t964;
t402 = qJ(5) * t448 + (-pkin(9) * t954 + t989) * t471;
t401 = -pkin(3) * t594 + pkin(8) * t503 + t434 * t964 + t435 * t960;
t400 = -qJ(2) * t505 + t406 * t957 - t477 * t955;
t399 = -t417 * t961 + t418 * t965;
t398 = t417 * t965 + t418 * t961;
t397 = -pkin(1) * t516 + qJ(2) * t506 + t406 * t955 + t477 * t957;
t396 = -pkin(7) * t432 - t428 * t961 + t431 * t965;
t395 = t399 * t957 + t471 * t955;
t394 = t399 * t955 - t471 * t957;
t393 = -pkin(7) * t468 - t408 * t961 + t411 * t965;
t392 = -pkin(7) * t466 - t407 * t961 + t410 * t965;
t391 = -pkin(7) * t457 - t401 * t961 + t403 * t965;
t390 = -pkin(8) * t417 - t402 * t960 + t412 * t964;
t389 = -pkin(2) * t398 - pkin(3) * t417 - pkin(4) * t447 + pkin(5) * t552 - pkin(9) * t472;
t388 = -pkin(3) * t471 + pkin(8) * t418 + t402 * t964 + t412 * t960;
t387 = -qJ(2) * t455 + t393 * t957 - t426 * t955;
t386 = -qJ(2) * t453 + t392 * t957 - t423 * t955;
t385 = -qJ(2) * t429 + t396 * t957 - t409 * t955;
t384 = -pkin(1) * t468 + qJ(2) * t456 + t393 * t955 + t426 * t957;
t383 = -t394 * t962 + t395 * t966;
t382 = t394 * t966 + t395 * t962;
t381 = -pkin(1) * t466 + qJ(2) * t454 + t392 * t955 + t423 * t957;
t380 = -pkin(1) * t432 + qJ(2) * t430 + t396 * t955 + t409 * t957;
t379 = -qJ(2) * t442 + t391 * t957 - t413 * t955;
t378 = -pkin(1) * t457 + qJ(2) * t443 + t391 * t955 + t413 * t957;
t377 = -pkin(7) * t398 - t388 * t961 + t390 * t965;
t376 = -qJ(2) * t394 + t377 * t957 - t389 * t955;
t375 = -pkin(1) * t398 + qJ(2) * t395 + t377 * t955 + t389 * t957;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t918, -t919, 0, t871, 0, 0, 0, 0, 0, 0, -t1046, -t864, 0, t717, 0, 0, 0, 0, 0, 0, t761, t762, t802, t640, 0, 0, 0, 0, 0, 0, t600, t608, t561, t491, 0, 0, 0, 0, 0, 0, t485, t496, t460, t405, 0, 0, 0, 0, 0, 0, t422, t425, t416, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t919, -t918, 0, t870, 0, 0, 0, 0, 0, 0, t864, -t1046, 0, -t1063, 0, 0, 0, 0, 0, 0, t759, t760, t801, t639, 0, 0, 0, 0, 0, 0, t599, t607, t560, t490, 0, 0, 0, 0, 0, 0, t484, t495, t459, t404, 0, 0, 0, 0, 0, 0, t421, t424, t415, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t952, 0, 0, 0, 0, 0, 0, t872, t873, 0, -t752, 0, 0, 0, 0, 0, 0, t707, t732, t651, t562, 0, 0, 0, 0, 0, 0, t558, t583, t516, t432, 0, 0, 0, 0, 0, 0, t466, t468, t457, t398; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t919, 0, -t918, 0, t977, -t896, -t870, -pkin(6) * t870, 0, 0, t864, 0, -t1046, 0, t1062, t1055, t1063, pkin(6) * t1063 + qJ(2) * t1025 - t775 * t962, -t851 * t962 + t853 * t966, -t816 * t962 + t818 * t966, -t847 * t962 + t849 * t966, -t850 * t962 + t852 * t966, -t846 * t962 + t848 * t966, -t880 * t962 + t881 * t966, -pkin(6) * t759 - t677 * t962 + t689 * t966, -pkin(6) * t760 - t678 * t962 + t690 * t966, -pkin(6) * t801 - t735 * t962 + t736 * t966, -pkin(6) * t639 - t604 * t962 + t616 * t966, -t692 * t962 + t694 * t966, -t632 * t962 + t633 * t966, -t682 * t962 + t684 * t966, -t691 * t962 + t693 * t966, -t683 * t962 + t685 * t966, -t747 * t962 + t748 * t966, -pkin(6) * t599 - t500 * t962 + t523 * t966, -pkin(6) * t607 - t512 * t962 + t532 * t966, -pkin(6) * t560 - t470 * t962 + t478 * t966, -pkin(6) * t490 - t436 * t962 + t444 * t966, -t548 * t962 + t550 * t966, -t508 * t962 + t509 * t966, -t554 * t962 + t556 * t966, -t547 * t962 + t549 * t966, -t555 * t962 + t557 * t966, -t612 * t962 + t613 * t966, -pkin(6) * t484 - t414 * t962 + t419 * t966, -pkin(6) * t495 - t420 * t962 + t427 * t966, -pkin(6) * t459 - t397 * t962 + t400 * t966, -pkin(6) * t404 - t380 * t962 + t385 * t966, -t481 * t962 + t483 * t966, -t449 * t962 + t450 * t966, -t461 * t962 + t463 * t966, -t480 * t962 + t482 * t966, -t462 * t962 + t464 * t966, -t498 * t962 + t499 * t966, -pkin(6) * t421 - t381 * t962 + t386 * t966, -pkin(6) * t424 - t384 * t962 + t387 * t966, -pkin(6) * t415 - t378 * t962 + t379 * t966, -pkin(6) * t382 - t375 * t962 + t376 * t966; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t918, 0, t919, 0, t896, t977, t871, pkin(6) * t871, 0, 0, t1046, 0, t864, 0, -t1055, t1062, t717, pkin(6) * t717 + qJ(2) * t1026 + t775 * t966, t851 * t966 + t853 * t962, t816 * t966 + t818 * t962, t847 * t966 + t849 * t962, t850 * t966 + t852 * t962, t846 * t966 + t848 * t962, t880 * t966 + t881 * t962, pkin(6) * t761 + t677 * t966 + t689 * t962, pkin(6) * t762 + t678 * t966 + t690 * t962, pkin(6) * t802 + t735 * t966 + t736 * t962, pkin(6) * t640 + t604 * t966 + t616 * t962, t692 * t966 + t694 * t962, t632 * t966 + t633 * t962, t682 * t966 + t684 * t962, t691 * t966 + t693 * t962, t683 * t966 + t685 * t962, t747 * t966 + t748 * t962, pkin(6) * t600 + t500 * t966 + t523 * t962, pkin(6) * t608 + t512 * t966 + t532 * t962, pkin(6) * t561 + t470 * t966 + t478 * t962, pkin(6) * t491 + t436 * t966 + t444 * t962, t548 * t966 + t550 * t962, t508 * t966 + t509 * t962, t554 * t966 + t556 * t962, t547 * t966 + t549 * t962, t555 * t966 + t557 * t962, t612 * t966 + t613 * t962, pkin(6) * t485 + t414 * t966 + t419 * t962, pkin(6) * t496 + t420 * t966 + t427 * t962, pkin(6) * t460 + t397 * t966 + t400 * t962, pkin(6) * t405 + t380 * t966 + t385 * t962, t481 * t966 + t483 * t962, t449 * t966 + t450 * t962, t461 * t966 + t463 * t962, t480 * t966 + t482 * t962, t462 * t966 + t464 * t962, t498 * t966 + t499 * t962, pkin(6) * t422 + t381 * t966 + t386 * t962, pkin(6) * t425 + t384 * t966 + t387 * t962, pkin(6) * t416 + t378 * t966 + t379 * t962, pkin(6) * t383 + t375 * t966 + t376 * t962; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t925, t926, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t916 - t860, -pkin(1) * t915 - t861, 0, -pkin(1) * t788, (t913 + t987) * t961, t912 * t965 + t914 * t961, t928 * t965 + t1012, t1015, t930 * t961 + t1010, 0, pkin(1) * t828 + pkin(2) * t914 + pkin(7) * t875 - t1018, pkin(1) * t829 - pkin(2) * t912 + pkin(7) * t877 + t1019, pkin(1) * t868 + pkin(2) * t920 + pkin(7) * t917 + t753, pkin(1) * t718 - pkin(2) * t844 + pkin(7) * t753, t784 * t965 + t785 * t961, t728 * t965 + t730 * t961, t803 * t965 + t805 * t961, t782 * t965 + t783 * t961, t804 * t965 + t806 * t961, t824 * t965 + t825 * t961, pkin(1) * t674 - pkin(2) * t793 + pkin(7) * t708 + t676 * t965 + t711 * t961, pkin(1) * t679 - pkin(2) * t1047 + pkin(7) * t733 + t681 * t965 + t727 * t961, pkin(1) * t626 - pkin(2) * t822 + pkin(7) * t653 + t592 * t965 + t601 * t961, pkin(1) * t539 - pkin(2) * t790 + pkin(7) * t563 - pkin(8) * t1034 + t610 * t965, t636 * t965 + t638 * t961, t572 * t965 + t574 * t961, t658 * t965 + t660 * t961, t635 * t965 + t637 * t961, t659 * t965 + t661 * t961, t695 * t965 + t696 * t961, pkin(1) * t537 - pkin(2) * t721 + pkin(7) * t559 + t501 * t965 + t521 * t961, pkin(1) * t545 - pkin(2) * t725 + pkin(7) * t584 + t515 * t965 + t533 * t961, pkin(1) * t505 - pkin(2) * t746 + pkin(7) * t518 + t437 * t965 + t438 * t961, pkin(1) * t429 - pkin(2) * t699 + pkin(7) * t433 + t428 * t965 + t431 * t961, t542 * t965 + t544 * t961, t510 * t965 + t511 * t961, t528 * t965 + t530 * t961, t541 * t965 + t543 * t961, t529 * t965 + t531 * t961, t578 * t965 + t579 * t961, pkin(1) * t453 - pkin(2) * t622 + pkin(7) * t467 + t407 * t965 + t410 * t961, pkin(1) * t455 - pkin(2) * t628 + pkin(7) * t469 + t408 * t965 + t411 * t961, pkin(1) * t442 - pkin(2) * t594 + pkin(7) * t458 + t401 * t965 + t403 * t961, pkin(1) * t394 - pkin(2) * t471 + pkin(7) * t399 + t388 * t965 + t390 * t961;];
tauB_reg  = t1;