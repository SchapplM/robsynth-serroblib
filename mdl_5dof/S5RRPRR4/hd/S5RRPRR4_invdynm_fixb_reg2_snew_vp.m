% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:22
% EndTime: 2020-01-03 12:02:35
% DurationCPUTime: 13.83s
% Computational Cost: add. (75545->550), mult. (100093->770), div. (0->0), fcn. (62071->10), ass. (0->377)
t1022 = sin(qJ(1));
t1026 = cos(qJ(1));
t1021 = sin(qJ(2));
t1025 = cos(qJ(2));
t1016 = g(1) - qJDD(3);
t1017 = sin(pkin(9));
t1014 = qJD(1) + qJD(2);
t1010 = t1014 ^ 2;
t1012 = qJDD(1) + qJDD(2);
t1018 = cos(pkin(9));
t974 = t1010 * t1017 - t1012 * t1018;
t1107 = qJ(3) * t974 - t1016 * t1017;
t971 = t1010 * t1018 + t1012 * t1017;
t912 = t1021 * t971 + t1025 * t974;
t944 = qJ(3) * t971 - t1016 * t1018;
t1120 = pkin(6) * t912 + t1021 * t944 + t1025 * t1107;
t909 = t1021 * t974 - t1025 * t971;
t1121 = t1022 * t909 - t1026 * t912;
t825 = pkin(6) * t909 + t1021 * t1107 - t1025 * t944;
t1131 = pkin(5) * t1121 + t1022 * t825 - t1026 * t1120;
t1122 = t1022 * t912 + t1026 * t909;
t1130 = pkin(5) * t1122 + t1022 * t1120 + t1026 * t825;
t996 = g(2) * t1026 + g(3) * t1022;
t1032 = qJDD(1) * pkin(1) - t996;
t1096 = qJD(1) ^ 2;
t995 = g(2) * t1022 - g(3) * t1026;
t1036 = pkin(1) * t1096 + t995;
t920 = -t1021 * t1036 - t1025 * t1032;
t1029 = t1012 * pkin(2) - t920;
t921 = t1021 * t1032 - t1025 * t1036;
t916 = -t1010 * pkin(2) + t921;
t858 = t1017 * t916 - t1018 * t1029;
t859 = t1017 * t1029 + t1018 * t916;
t1056 = t1017 * t858 + t1018 * t859;
t786 = t1017 * t859 - t1018 * t858;
t1083 = t1021 * t786;
t1114 = t1025 * t1056 - t1083;
t1077 = t1025 * t786;
t730 = -t1021 * t1056 - t1077;
t1127 = t1022 * t730 + t1026 * t1114;
t1126 = t1022 * t1114 - t1026 * t730;
t977 = t1010 * t1025 + t1012 * t1021;
t980 = t1010 * t1021 - t1012 * t1025;
t1038 = t1022 * t977 + t1026 * t980;
t1108 = pkin(6) * t980 - g(1) * t1021;
t950 = pkin(6) * t977 - g(1) * t1025;
t1123 = pkin(5) * t1038 + t1022 * t950 + t1026 * t1108;
t1105 = t1022 * t980 - t1026 * t977;
t1119 = pkin(5) * t1105 + t1022 * t1108 - t1026 * t950;
t1054 = t1021 * t920 + t1025 * t921;
t869 = t1021 * t921 - t1025 * t920;
t1076 = t1026 * t869;
t1115 = t1022 * t1054 + t1076;
t1082 = t1022 * t869;
t1113 = t1026 * t1054 - t1082;
t1019 = sin(qJ(5));
t1011 = qJDD(4) + qJDD(5);
t1023 = cos(qJ(5));
t1024 = cos(qJ(4));
t1071 = t1014 * t1024;
t1020 = sin(qJ(4));
t1072 = t1014 * t1020;
t953 = t1019 * t1072 - t1023 * t1071;
t955 = (t1019 * t1024 + t1020 * t1023) * t1014;
t904 = t955 * t953;
t1103 = -t904 + t1011;
t1110 = t1019 * t1103;
t1109 = t1023 * t1103;
t840 = -t1010 * pkin(3) + t1012 * pkin(7) + t859;
t819 = t1016 * t1024 + t1020 * t840;
t820 = -t1016 * t1020 + t1024 * t840;
t771 = t1020 * t819 + t1024 * t820;
t1058 = qJD(4) * t1072;
t1066 = t1024 * t1012;
t1033 = t1058 - t1066;
t1067 = t1020 * t1012;
t998 = qJD(4) * t1071;
t968 = t998 + t1067;
t876 = -qJD(5) * t953 - t1019 * t1033 + t1023 * t968;
t1013 = qJD(4) + qJD(5);
t946 = t1013 * t953;
t1104 = -t946 + t876;
t1055 = t1019 * t968 + t1023 * t1033;
t849 = (qJD(5) - t1013) * t955 + t1055;
t951 = t953 ^ 2;
t952 = t955 ^ 2;
t1009 = t1013 ^ 2;
t1095 = t1024 ^ 2;
t994 = t1024 * t1010 * t1020;
t984 = qJDD(4) + t994;
t797 = (-t968 + t998) * pkin(8) + t984 * pkin(4) - t819;
t1000 = t1095 * t1010;
t986 = qJD(4) * pkin(4) - pkin(8) * t1072;
t798 = -pkin(4) * t1000 - pkin(8) * t1033 - qJD(4) * t986 + t820;
t746 = t1019 * t798 - t1023 * t797;
t748 = t1019 * t797 + t1023 * t798;
t697 = t1019 * t748 - t1023 * t746;
t1094 = pkin(4) * t697;
t853 = t946 + t876;
t781 = -t1019 * t849 - t1023 * t853;
t1093 = pkin(4) * t781;
t839 = -pkin(3) * t1012 - pkin(7) * t1010 + t858;
t1090 = -pkin(3) * t839 + pkin(7) * t771;
t1089 = t1013 * t955;
t802 = pkin(4) * t1033 - pkin(8) * t1000 + t1072 * t986 + t839;
t1088 = t1019 * t802;
t896 = t904 + t1011;
t1087 = t1019 * t896;
t1086 = t1020 * t697;
t831 = t1020 * t839;
t1085 = t1020 * t984;
t985 = qJDD(4) - t994;
t1084 = t1020 * t985;
t1081 = t1023 * t802;
t1080 = t1023 * t896;
t1079 = t1024 * t697;
t832 = t1024 * t839;
t969 = -0.2e1 * t1058 + t1066;
t922 = t1024 * t969;
t1078 = t1024 * t985;
t1075 = qJD(4) * t1014;
t1074 = t1013 * t1019;
t1073 = t1013 * t1023;
t1015 = t1020 ^ 2;
t1070 = t1015 * t1010;
t1064 = t1015 + t1095;
t1027 = qJD(4) ^ 2;
t990 = -t1027 - t1070;
t933 = -t1020 * t990 - t1078;
t967 = 0.2e1 * t998 + t1067;
t1063 = -pkin(3) * t967 + pkin(7) * t933 + t831;
t992 = -t1000 - t1027;
t931 = t1024 * t992 - t1085;
t1062 = pkin(3) * t969 + pkin(7) * t931 - t832;
t737 = t1017 * t771 - t1018 * t839;
t1061 = pkin(2) * t737 + t1090;
t1060 = t1017 * t904;
t1059 = t1018 * t904;
t987 = -qJDD(1) * t1022 - t1026 * t1096;
t1057 = pkin(5) * t987 + g(1) * t1026;
t699 = t1019 * t746 + t1023 * t748;
t1052 = -t1022 * t995 - t1026 * t996;
t783 = t1019 * t853 - t1023 * t849;
t879 = -t951 - t952;
t677 = -pkin(4) * t879 + pkin(8) * t783 + t699;
t685 = -pkin(8) * t781 - t697;
t725 = -t1020 * t781 + t1024 * t783;
t1051 = -pkin(3) * t879 + pkin(7) * t725 + t1020 * t685 + t1024 * t677;
t894 = -t1009 - t951;
t830 = t1023 * t894 - t1110;
t848 = (qJD(5) + t1013) * t955 + t1055;
t733 = -pkin(4) * t848 + pkin(8) * t830 - t1081;
t829 = t1019 * t894 + t1109;
t765 = -pkin(8) * t829 + t1088;
t774 = -t1020 * t829 + t1024 * t830;
t1050 = -pkin(3) * t848 + pkin(7) * t774 + t1020 * t765 + t1024 * t733;
t934 = -t952 - t1009;
t861 = -t1019 * t934 - t1080;
t742 = -pkin(4) * t1104 + pkin(8) * t861 + t1088;
t860 = t1023 * t934 - t1087;
t767 = -pkin(8) * t860 + t1081;
t791 = -t1020 * t860 + t1024 * t861;
t1049 = -pkin(3) * t1104 + pkin(7) * t791 + t1020 * t767 + t1024 * t742;
t975 = t1064 * t1012;
t981 = t1000 + t1070;
t1048 = pkin(3) * t981 + pkin(7) * t975 + t771;
t885 = t1017 * t933 - t1018 * t967;
t1047 = pkin(2) * t885 + t1063;
t884 = t1017 * t931 + t1018 * t969;
t1046 = pkin(2) * t884 + t1062;
t1045 = t1017 * t994;
t1044 = t1018 * t994;
t1043 = -pkin(2) * t974 - t858;
t716 = t1017 * t725 - t1018 * t879;
t1042 = pkin(2) * t716 + t1051;
t745 = t1017 * t774 - t1018 * t848;
t1041 = pkin(2) * t745 + t1050;
t750 = t1017 * t791 - t1018 * t1104;
t1040 = pkin(2) * t750 + t1049;
t914 = t1017 * t975 + t1018 * t981;
t1039 = pkin(2) * t914 + t1048;
t769 = t1020 * t820 - t1024 * t819;
t1037 = t1022 * t996 - t1026 * t995;
t1035 = pkin(4) * t829 - t746;
t672 = t1024 * t699 - t1086;
t689 = -pkin(4) * t802 + pkin(8) * t699;
t1034 = -pkin(3) * t802 + pkin(7) * t672 - pkin(8) * t1086 + t1024 * t689;
t667 = t1017 * t672 - t1018 * t802;
t1031 = pkin(2) * t667 + t1034;
t1030 = pkin(4) * t860 - t748;
t1028 = -pkin(2) * t971 - t859;
t991 = t1000 - t1027;
t989 = t1027 - t1070;
t988 = qJDD(1) * t1026 - t1022 * t1096;
t982 = -t1000 + t1070;
t976 = t1024 * t984;
t965 = pkin(5) * t988 + g(1) * t1022;
t962 = t1064 * t1075;
t940 = -t952 + t1009;
t939 = t951 - t1009;
t938 = qJDD(4) * t1017 + t1018 * t962;
t937 = -qJDD(4) * t1018 + t1017 * t962;
t936 = -t1015 * t1075 + t1024 * t968;
t935 = t1020 * t1033 - t1075 * t1095;
t932 = -t1020 * t989 + t976;
t930 = t1024 * t991 - t1084;
t929 = t1024 * t990 - t1084;
t928 = t1024 * t989 + t1085;
t927 = t1020 * t992 + t976;
t926 = t1020 * t991 + t1078;
t923 = (t968 + t998) * t1020;
t915 = -t1017 * t981 + t1018 * t975;
t906 = -t1020 * t967 + t922;
t905 = t1020 * t969 + t1024 * t967;
t902 = t952 - t951;
t901 = t1017 * t1067 + t1018 * t932;
t900 = t1017 * t1066 + t1018 * t930;
t899 = t1017 * t932 - t1018 * t1067;
t898 = t1017 * t930 - t1018 * t1066;
t893 = t1018 * t936 - t1045;
t892 = t1018 * t935 + t1045;
t891 = t1017 * t936 + t1044;
t890 = t1017 * t935 - t1044;
t889 = -pkin(1) * t980 - t920;
t888 = -pkin(1) * t977 - t921;
t887 = t1017 * t967 + t1018 * t933;
t886 = -t1017 * t969 + t1018 * t931;
t883 = (t1019 * t955 - t1023 * t953) * t1013;
t882 = (-t1019 * t953 - t1023 * t955) * t1013;
t878 = -t1021 * t937 + t1025 * t938;
t877 = t1021 * t938 + t1025 * t937;
t875 = -qJD(5) * t955 - t1055;
t874 = t1017 * t982 + t1018 * t906;
t873 = t1017 * t906 - t1018 * t982;
t867 = pkin(1) * t869;
t866 = pkin(1) * g(1) + pkin(6) * t1054;
t865 = t1023 * t939 - t1087;
t864 = -t1019 * t940 + t1109;
t863 = t1019 * t939 + t1080;
t862 = t1023 * t940 + t1110;
t857 = -t1021 * t914 + t1025 * t915;
t856 = t1021 * t915 + t1025 * t914;
t844 = t1023 * t876 - t1074 * t955;
t843 = t1019 * t876 + t1073 * t955;
t842 = -t1019 * t875 + t1073 * t953;
t841 = t1023 * t875 + t1074 * t953;
t838 = -t1021 * t899 + t1025 * t901;
t837 = -t1021 * t898 + t1025 * t900;
t836 = t1021 * t901 + t1025 * t899;
t835 = t1021 * t900 + t1025 * t898;
t824 = -t1021 * t891 + t1025 * t893;
t823 = -t1021 * t890 + t1025 * t892;
t822 = t1021 * t893 + t1025 * t891;
t821 = t1021 * t892 + t1025 * t890;
t816 = -t1021 * t885 + t1025 * t887;
t815 = -t1021 * t884 + t1025 * t886;
t814 = t1021 * t887 + t1025 * t885;
t813 = t1021 * t886 + t1025 * t884;
t812 = -t1020 * t882 + t1024 * t883;
t811 = t1020 * t883 + t1024 * t882;
t810 = t1011 * t1017 + t1018 * t812;
t809 = -t1011 * t1018 + t1017 * t812;
t808 = -pkin(7) * t929 + t832;
t807 = -pkin(7) * t927 + t831;
t806 = -t1021 * t873 + t1025 * t874;
t805 = t1021 * t874 + t1025 * t873;
t804 = -pkin(3) * t929 + t820;
t803 = -pkin(3) * t927 + t819;
t800 = -pkin(1) * t912 + t1043;
t799 = pkin(1) * t909 + t1028;
t795 = -t1020 * t863 + t1024 * t865;
t794 = -t1020 * t862 + t1024 * t864;
t793 = t1020 * t865 + t1024 * t863;
t792 = t1020 * t864 + t1024 * t862;
t790 = t1020 * t861 + t1024 * t860;
t784 = pkin(2) * t786;
t782 = -t1019 * t1104 - t1023 * t848;
t780 = -t1019 * t848 + t1023 * t1104;
t779 = pkin(2) * t1016 + qJ(3) * t1056;
t778 = -t1020 * t843 + t1024 * t844;
t777 = -t1020 * t841 + t1024 * t842;
t776 = t1020 * t844 + t1024 * t843;
t775 = t1020 * t842 + t1024 * t841;
t773 = t1020 * t830 + t1024 * t829;
t764 = t1018 * t778 + t1060;
t763 = t1018 * t777 - t1060;
t762 = t1017 * t778 - t1059;
t761 = t1017 * t777 + t1059;
t759 = -t1017 * t849 + t1018 * t795;
t758 = t1017 * t853 + t1018 * t794;
t757 = t1017 * t795 + t1018 * t849;
t756 = t1017 * t794 - t1018 * t853;
t755 = -qJ(3) * t914 - t1018 * t769;
t754 = qJ(3) * t915 - t1017 * t769;
t753 = -t1021 * t809 + t1025 * t810;
t752 = t1017 * t1104 + t1018 * t791;
t751 = t1021 * t810 + t1025 * t809;
t747 = t1017 * t848 + t1018 * t774;
t740 = pkin(1) * t814 + t1047;
t739 = pkin(1) * t813 + t1046;
t738 = t1017 * t839 + t1018 * t771;
t735 = -qJ(3) * t885 - t1017 * t804 + t1018 * t808;
t734 = -qJ(3) * t884 - t1017 * t803 + t1018 * t807;
t727 = -pkin(2) * t929 + qJ(3) * t887 + t1017 * t808 + t1018 * t804;
t726 = -pkin(2) * t927 + qJ(3) * t886 + t1017 * t807 + t1018 * t803;
t724 = -t1020 * t780 + t1024 * t782;
t723 = t1020 * t783 + t1024 * t781;
t722 = t1020 * t782 + t1024 * t780;
t720 = pkin(1) * t856 + t1039;
t719 = t1017 * t902 + t1018 * t724;
t718 = t1017 * t724 - t1018 * t902;
t717 = t1017 * t879 + t1018 * t725;
t714 = -t1021 * t762 + t1025 * t764;
t713 = -t1021 * t761 + t1025 * t763;
t712 = t1021 * t764 + t1025 * t762;
t711 = t1021 * t763 + t1025 * t761;
t710 = -t1021 * t757 + t1025 * t759;
t709 = -t1021 * t756 + t1025 * t758;
t708 = t1021 * t759 + t1025 * t757;
t707 = t1021 * t758 + t1025 * t756;
t706 = -pkin(1) * t730 + t784;
t705 = -t1021 * t750 + t1025 * t752;
t704 = t1021 * t752 + t1025 * t750;
t703 = -pkin(3) * t790 - t1030;
t702 = -pkin(3) * t723 - t1093;
t701 = -pkin(3) * t773 - t1035;
t700 = -t1021 * t745 + t1025 * t747;
t698 = t1021 * t747 + t1025 * t745;
t696 = -pkin(6) * t856 - t1021 * t754 + t1025 * t755;
t695 = pkin(6) * t857 + t1021 * t755 + t1025 * t754;
t694 = -t1021 * t737 + t1025 * t738;
t693 = t1021 * t738 + t1025 * t737;
t692 = -pkin(7) * t790 - t1020 * t742 + t1024 * t767;
t691 = pkin(6) * t730 - qJ(3) * t1077 - t1021 * t779;
t690 = pkin(1) * t1016 + pkin(6) * t1114 - qJ(3) * t1083 + t1025 * t779;
t687 = -qJ(3) * t737 + (pkin(3) * t1017 - pkin(7) * t1018) * t769;
t686 = -pkin(7) * t773 - t1020 * t733 + t1024 * t765;
t683 = -pkin(6) * t814 - t1021 * t727 + t1025 * t735;
t682 = -pkin(6) * t813 - t1021 * t726 + t1025 * t734;
t681 = -pkin(1) * t929 + pkin(6) * t816 + t1021 * t735 + t1025 * t727;
t680 = -pkin(1) * t927 + pkin(6) * t815 + t1021 * t734 + t1025 * t726;
t679 = -t1021 * t718 + t1025 * t719;
t678 = t1021 * t719 + t1025 * t718;
t675 = -t1021 * t716 + t1025 * t717;
t674 = t1021 * t717 + t1025 * t716;
t673 = qJ(3) * t738 + (-pkin(3) * t1018 - pkin(7) * t1017 - pkin(2)) * t769;
t671 = t1020 * t699 + t1079;
t669 = pkin(1) * t693 + t1061;
t668 = t1017 * t802 + t1018 * t672;
t665 = -qJ(3) * t750 - t1017 * t703 + t1018 * t692;
t664 = pkin(1) * t704 + t1040;
t663 = -qJ(3) * t745 - t1017 * t701 + t1018 * t686;
t662 = -pkin(2) * t790 + qJ(3) * t752 + t1017 * t692 + t1018 * t703;
t661 = pkin(1) * t698 + t1041;
t660 = -pkin(3) * t671 - t1094;
t659 = -pkin(2) * t773 + qJ(3) * t747 + t1017 * t686 + t1018 * t701;
t658 = -pkin(7) * t723 - t1020 * t677 + t1024 * t685;
t657 = -pkin(6) * t693 - t1021 * t673 + t1025 * t687;
t656 = -pkin(7) * t671 - pkin(8) * t1079 - t1020 * t689;
t655 = -pkin(1) * t769 + pkin(6) * t694 + t1021 * t687 + t1025 * t673;
t654 = -t1021 * t667 + t1025 * t668;
t653 = t1021 * t668 + t1025 * t667;
t652 = -qJ(3) * t716 - t1017 * t702 + t1018 * t658;
t651 = -pkin(2) * t723 + qJ(3) * t717 + t1017 * t658 + t1018 * t702;
t650 = pkin(1) * t674 + t1042;
t649 = -pkin(6) * t704 - t1021 * t662 + t1025 * t665;
t648 = -pkin(1) * t790 + pkin(6) * t705 + t1021 * t665 + t1025 * t662;
t647 = -pkin(6) * t698 - t1021 * t659 + t1025 * t663;
t646 = -pkin(1) * t773 + pkin(6) * t700 + t1021 * t663 + t1025 * t659;
t645 = -qJ(3) * t667 - t1017 * t660 + t1018 * t656;
t644 = pkin(1) * t653 + t1031;
t643 = -pkin(2) * t671 + qJ(3) * t668 + t1017 * t656 + t1018 * t660;
t642 = -pkin(6) * t674 - t1021 * t651 + t1025 * t652;
t641 = -pkin(1) * t723 + pkin(6) * t675 + t1021 * t652 + t1025 * t651;
t640 = -pkin(6) * t653 - t1021 * t643 + t1025 * t645;
t639 = -pkin(1) * t671 + pkin(6) * t654 + t1021 * t645 + t1025 * t643;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t996, t995, 0, 0, 0, 0, 0, 0, 0, t1012, t889, t888, 0, t867, 0, 0, 0, 0, 0, t1012, t800, t799, 0, t706, t923, t905, t928, t922, t926, 0, t739, t740, t720, t669, t776, t722, t792, t775, t793, t811, t661, t664, t650, t644; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t987, 0, t988, 0, t1057, -t965, t1037, pkin(5) * t1037, 0, 0, -t1105, 0, -t1038, 0, t1119, t1123, t1113, pkin(5) * t1113 - pkin(6) * t1082 + t1026 * t866, 0, 0, -t1122, 0, t1121, 0, t1130, -t1131, t1127, pkin(5) * t1127 + t1022 * t691 + t1026 * t690, t1022 * t824 + t1026 * t822, t1022 * t806 + t1026 * t805, t1022 * t838 + t1026 * t836, t1022 * t823 + t1026 * t821, t1022 * t837 + t1026 * t835, t1022 * t878 + t1026 * t877, t1022 * t682 + t1026 * t680 - pkin(5) * (t1022 * t813 - t1026 * t815), t1022 * t683 + t1026 * t681 - pkin(5) * (t1022 * t814 - t1026 * t816), t1022 * t696 + t1026 * t695 - pkin(5) * (t1022 * t856 - t1026 * t857), t1022 * t657 + t1026 * t655 - pkin(5) * (t1022 * t693 - t1026 * t694), t1022 * t714 + t1026 * t712, t1022 * t679 + t1026 * t678, t1022 * t709 + t1026 * t707, t1022 * t713 + t1026 * t711, t1022 * t710 + t1026 * t708, t1022 * t753 + t1026 * t751, t1022 * t647 + t1026 * t646 - pkin(5) * (t1022 * t698 - t1026 * t700), t1022 * t649 + t1026 * t648 - pkin(5) * (t1022 * t704 - t1026 * t705), t1022 * t642 + t1026 * t641 - pkin(5) * (t1022 * t674 - t1026 * t675), t1022 * t640 + t1026 * t639 - pkin(5) * (t1022 * t653 - t1026 * t654); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t988, 0, -t987, 0, t965, t1057, t1052, pkin(5) * t1052, 0, 0, t1038, 0, -t1105, 0, -t1123, t1119, t1115, pkin(5) * t1115 + pkin(6) * t1076 + t1022 * t866, 0, 0, -t1121, 0, -t1122, 0, t1131, t1130, t1126, pkin(5) * t1126 + t1022 * t690 - t1026 * t691, t1022 * t822 - t1026 * t824, t1022 * t805 - t1026 * t806, t1022 * t836 - t1026 * t838, t1022 * t821 - t1026 * t823, t1022 * t835 - t1026 * t837, t1022 * t877 - t1026 * t878, -t1026 * t682 + t1022 * t680 + pkin(5) * (t1022 * t815 + t1026 * t813), -t1026 * t683 + t1022 * t681 + pkin(5) * (t1022 * t816 + t1026 * t814), -t1026 * t696 + t1022 * t695 + pkin(5) * (t1022 * t857 + t1026 * t856), -t1026 * t657 + t1022 * t655 + pkin(5) * (t1022 * t694 + t1026 * t693), t1022 * t712 - t1026 * t714, t1022 * t678 - t1026 * t679, t1022 * t707 - t1026 * t709, t1022 * t711 - t1026 * t713, t1022 * t708 - t1026 * t710, t1022 * t751 - t1026 * t753, -t1026 * t647 + t1022 * t646 + pkin(5) * (t1022 * t700 + t1026 * t698), -t1026 * t649 + t1022 * t648 + pkin(5) * (t1022 * t705 + t1026 * t704), -t1026 * t642 + t1022 * t641 + pkin(5) * (t1022 * t675 + t1026 * t674), -t1026 * t640 + t1022 * t639 + pkin(5) * (t1022 * t654 + t1026 * t653); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1096, 0, 0, -g(1), t996, 0, 0, 0, -t980, 0, -t977, 0, t1108, t950, -t869, -pkin(6) * t869, 0, 0, -t912, 0, t909, 0, t1120, -t825, t730, t691, t824, t806, t838, t823, t837, t878, t682, t683, t696, t657, t714, t679, t709, t713, t710, t753, t647, t649, t642, t640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1096, 0, qJDD(1), 0, g(1), 0, -t995, 0, 0, 0, t977, 0, -t980, 0, -t950, t1108, t1054, t866, 0, 0, -t909, 0, -t912, 0, t825, t1120, t1114, t690, t822, t805, t836, t821, t835, t877, t680, t681, t695, t655, t712, t678, t707, t711, t708, t751, t646, t648, t641, t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t996, t995, 0, 0, 0, 0, 0, 0, 0, t1012, t889, t888, 0, t867, 0, 0, 0, 0, 0, t1012, t800, t799, 0, t706, t923, t905, t928, t922, t926, 0, t739, t740, t720, t669, t776, t722, t792, t775, t793, t811, t661, t664, t650, t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1012, 0, -t1010, 0, 0, -g(1), t920, 0, 0, 0, -t974, 0, -t971, 0, t1107, t944, -t786, -qJ(3) * t786, t893, t874, t901, t892, t900, t938, t734, t735, t755, t687, t764, t719, t758, t763, t759, t810, t663, t665, t652, t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1010, 0, t1012, 0, g(1), 0, t921, 0, 0, 0, t971, 0, -t974, 0, -t944, t1107, t1056, t779, t891, t873, t899, t890, t898, t937, t726, t727, t754, t673, t762, t718, t756, t761, t757, t809, t659, t662, t651, t643; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1012, -t920, -t921, 0, 0, 0, 0, 0, 0, 0, t1012, t1043, t1028, 0, t784, t923, t905, t928, t922, t926, 0, t1046, t1047, t1039, t1061, t776, t722, t792, t775, t793, t811, t1041, t1040, t1042, t1031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1012, 0, -t1010, 0, 0, -t1016, t858, 0, t936, t906, t932, t935, t930, t962, t807, t808, -t769, -pkin(7) * t769, t778, t724, t794, t777, t795, t812, t686, t692, t658, t656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1010, 0, t1012, 0, t1016, 0, t859, 0, t994, -t982, -t1067, -t994, -t1066, -qJDD(4), t803, t804, 0, -pkin(3) * t769, -t904, -t902, -t853, t904, t849, -t1011, t701, t703, t702, t660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1012, -t858, -t859, 0, 0, t923, t905, t928, t922, t926, 0, t1062, t1063, t1048, t1090, t776, t722, t792, t775, t793, t811, t1050, t1049, t1051, t1034; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t968, t969, t984, -t998, t991, t998, 0, t839, t819, 0, t844, t782, t864, t842, t865, t883, t765, t767, t685, -pkin(8) * t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, t967, t989, -t1033, t985, -t1058, -t839, 0, t820, 0, t843, t780, t862, t841, t863, t882, t733, t742, t677, t689; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t994, t982, t1067, t994, t1066, qJDD(4), -t819, -t820, 0, 0, t904, t902, t853, -t904, -t849, t1011, t1035, t1030, t1093, t1094; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876, -t848, t1103, t946, t939, -t946, 0, t802, t746, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1089, t1104, t940, t875, t896, -t1089, -t802, 0, t748, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t904, t902, t853, -t904, -t849, t1011, -t746, -t748, 0, 0;];
m_new_reg = t1;
