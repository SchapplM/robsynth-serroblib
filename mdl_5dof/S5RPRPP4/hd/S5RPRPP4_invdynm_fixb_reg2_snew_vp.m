% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPP4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:15:05
% EndTime: 2019-12-31 18:15:17
% DurationCPUTime: 12.82s
% Computational Cost: add. (20955->517), mult. (46568->536), div. (0->0), fcn. (28502->6), ass. (0->345)
t940 = sin(pkin(7));
t941 = cos(pkin(7));
t943 = sin(qJ(3));
t945 = cos(qJ(3));
t898 = (-t940 * t945 - t941 * t943) * qJD(1);
t1028 = qJD(3) * t898;
t1023 = qJD(1) * qJD(3);
t1011 = t943 * t1023;
t927 = t945 * qJDD(1);
t907 = t927 - t1011;
t1010 = t945 * t1023;
t1021 = t943 * qJDD(1);
t971 = -t1010 - t1021;
t851 = t941 * t907 + t940 * t971;
t1085 = t851 + t1028;
t1029 = qJD(1) * t945;
t900 = -t940 * t943 * qJD(1) + t941 * t1029;
t1059 = t900 * t898;
t1084 = qJDD(3) - t1059;
t1052 = t940 * t1084;
t896 = t900 ^ 2;
t947 = qJD(3) ^ 2;
t837 = t947 + t896;
t777 = t941 * t837 + t1052;
t1047 = t941 * t1084;
t779 = -t940 * t837 + t1047;
t721 = t945 * t777 + t943 * t779;
t944 = sin(qJ(1));
t946 = cos(qJ(1));
t1153 = pkin(5) * (t1085 * t944 + t946 * t721);
t1152 = pkin(5) * (-t1085 * t946 + t944 * t721);
t1151 = pkin(2) * t721;
t1150 = pkin(6) * t721;
t742 = t943 * t777 - t945 * t779;
t1149 = pkin(6) * t742;
t1148 = qJ(2) * t742;
t1067 = pkin(6) + pkin(1);
t1147 = t1067 * t742;
t1068 = t898 ^ 2;
t875 = t1068 - t947;
t788 = t940 * t875 + t1047;
t795 = -t941 * t875 + t1052;
t731 = t945 * t788 - t943 * t795;
t850 = t940 * t907 - t941 * t971;
t888 = qJD(3) * t900;
t820 = t850 - t888;
t1146 = t944 * t731 - t946 * t820;
t1145 = t946 * t731 + t944 * t820;
t1073 = t896 - t1068;
t819 = t850 + t888;
t757 = t1085 * t941 - t940 * t819;
t1048 = t941 * t819;
t1054 = t940 * t1085;
t759 = t1048 + t1054;
t705 = -t945 * t757 + t943 * t759;
t1144 = -t1073 * t946 + t944 * t705;
t1143 = t1073 * t944 + t946 * t705;
t1140 = qJ(2) * t1085 + t1067 * t721;
t775 = pkin(3) * t777;
t1139 = qJ(4) * t777;
t1138 = qJ(4) * t779;
t709 = t943 * t757 + t945 * t759;
t1086 = -t1028 + t851;
t1083 = qJDD(3) + t1059;
t1051 = t940 * t1083;
t876 = t896 - t947;
t1118 = -t941 * t876 + t1051;
t831 = t941 * t1083;
t1119 = t940 * t876 + t831;
t1125 = t1118 * t945 + t1119 * t943;
t1137 = t1086 * t946 + t1125 * t944;
t1136 = t1086 * t944 - t1125 * t946;
t738 = t943 * t788 + t945 * t795;
t1091 = t940 * t1086 - t941 * t820;
t1092 = -t941 * t1086 - t940 * t820;
t1111 = t1091 * t943 + t1092 * t945;
t1135 = pkin(2) * t1111;
t1134 = pkin(6) * t1111;
t1112 = t1091 * t945 - t1092 * t943;
t1133 = pkin(6) * t1112;
t1130 = qJ(2) * t1112;
t1129 = t1067 * t1112;
t818 = t896 + t1068;
t1128 = -qJ(2) * t818 - t1067 * t1111;
t1127 = pkin(5) * (t1111 * t944 - t946 * t818);
t1126 = pkin(5) * (-t1111 * t946 - t944 * t818);
t1124 = -t1118 * t943 + t1119 * t945;
t1074 = -t1068 - t947;
t1077 = t1074 * t941 - t1051;
t1080 = t1074 * t940 + t831;
t1089 = t1077 * t943 + t1080 * t945;
t1066 = pkin(2) * t1089;
t754 = pkin(3) * t1092;
t1123 = pkin(6) * t1089;
t1090 = t1077 * t945 - t1080 * t943;
t1122 = pkin(6) * t1090;
t1121 = qJ(2) * t1090;
t1120 = qJ(4) * t1092;
t1117 = t1067 * t1090;
t1116 = pkin(3) * t818 + qJ(4) * t1091;
t1115 = qJ(2) * t819 - t1067 * t1089;
t1114 = pkin(5) * (-t1089 * t946 + t944 * t819);
t1113 = pkin(5) * (t1089 * t944 + t946 * t819);
t1110 = pkin(2) * t818;
t1108 = pkin(2) * t1085;
t785 = pkin(3) * t1080;
t1099 = qJ(4) * t1077;
t1098 = qJ(4) * t1080;
t1097 = t1085 * qJ(5);
t970 = (t898 * t940 - t900 * t941) * qJD(3);
t1026 = qJD(3) * t941;
t1012 = t898 * t1026;
t1027 = qJD(3) * t940;
t874 = t900 * t1027;
t986 = t874 + t1012;
t1072 = t943 * t986 + t945 * t970;
t926 = t944 * qJDD(3);
t1082 = -t1072 * t946 + t926;
t928 = t946 * qJDD(3);
t1081 = t1072 * t944 + t928;
t1016 = t946 * t1059;
t972 = t940 * t850 - t1012;
t987 = -t898 * t1027 - t941 * t850;
t1069 = t943 * t972 + t945 * t987;
t1079 = t1069 * t944 + t1016;
t1017 = t944 * t1059;
t1078 = -t1069 * t946 + t1017;
t1076 = -pkin(4) * t1083 - qJ(5) * t1074;
t1022 = qJD(5) * qJD(3);
t1025 = qJD(4) * t898;
t883 = 0.2e1 * t1025;
t1075 = t883 + 0.2e1 * t1022;
t948 = qJD(1) ^ 2;
t916 = t944 * g(1) - t946 * g(2);
t989 = qJDD(2) - t916;
t966 = -t948 * qJ(2) + t989;
t955 = -t1067 * qJDD(1) + t966;
t953 = t945 * t955;
t854 = t943 * g(3) + t953;
t855 = t945 * g(3) - t943 * t955;
t783 = t945 * t854 - t943 * t855;
t937 = t943 ^ 2;
t1058 = t937 * t948;
t920 = -t947 - t1058;
t1071 = -t943 * t970 + t945 * t986;
t1070 = -t943 * t987 + t945 * t972;
t1065 = pkin(2) * t783;
t935 = qJDD(1) * qJ(2);
t917 = t946 * g(1) + t944 * g(2);
t960 = -0.2e1 * qJD(2) * qJD(1) + t917;
t958 = -t935 + t960;
t870 = t1067 * t948 + t958;
t1064 = pkin(2) * t870;
t938 = t945 ^ 2;
t1030 = t937 + t938;
t909 = t1030 * qJDD(1);
t1063 = pkin(2) * t909;
t1062 = pkin(4) * t941;
t1061 = t850 * pkin(4);
t1060 = qJDD(1) * pkin(1);
t1057 = t938 * t948;
t802 = -(qJD(3) * pkin(3) - qJ(4) * t1029) * t1029 + t971 * pkin(3) - qJDD(4) + (t937 * qJ(4) + t1067) * t948 + t958;
t1056 = t940 * t802;
t1049 = t941 * t802;
t1024 = qJD(4) * t900;
t1018 = 0.2e1 * t1024;
t1036 = t945 * t948;
t800 = t953 - t907 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t1036 - qJ(4) * t1023 + g(3)) * t943;
t801 = pkin(3) * t920 - qJ(4) * t1021 - t855;
t999 = -t941 * t800 + t940 * t801;
t727 = t999 + t1018;
t1032 = t940 * t800 + t941 * t801;
t728 = t883 + t1032;
t699 = -t941 * t727 + t940 * t728;
t1046 = t943 * t699;
t1044 = t943 * t870;
t906 = 0.2e1 * t1010 + t1021;
t857 = t943 * t906;
t923 = t943 * t1036;
t914 = qJDD(3) + t923;
t1043 = t943 * t914;
t915 = qJDD(3) - t923;
t1042 = t943 * t915;
t1040 = t944 * t909;
t1039 = t945 * t699;
t856 = t945 * t870;
t1038 = t945 * t914;
t1037 = t945 * t915;
t1034 = t946 * t909;
t835 = -t898 * pkin(4) - t900 * qJ(5);
t983 = -t947 * pkin(4) + qJDD(3) * qJ(5) + t898 * t835 + t1032;
t714 = t983 + t1075;
t968 = -qJDD(3) * pkin(4) - t947 * qJ(5) + qJDD(5) + t999;
t716 = (0.2e1 * qJD(4) + t835) * t900 + t968;
t1033 = -pkin(4) * t716 + qJ(5) * t714;
t1031 = -pkin(4) * t1086 - qJ(5) * t820;
t1020 = t944 * qJDD(1);
t1019 = t946 * qJDD(1);
t687 = t940 * t714 - t941 * t716;
t1015 = pkin(3) * t687 + t1033;
t1014 = -t1032 - t775;
t1013 = t1031 + t754;
t700 = t940 * t727 + t941 * t728;
t679 = t943 * t700 + t1039;
t698 = pkin(3) * t699;
t1009 = -pkin(2) * t679 - t698;
t1008 = -t754 - t1135;
t1007 = -qJ(5) * t940 - pkin(3);
t688 = t941 * t714 + t940 * t716;
t950 = -pkin(4) * t888 + 0.2e1 * qJD(5) * t900 + t802;
t949 = t950 + t1097;
t718 = t949 - t1061;
t665 = qJ(4) * t688 + (-t1007 + t1062) * t718;
t669 = -qJ(4) * t687 + (-pkin(4) * t940 + qJ(5) * t941) * t718;
t1006 = -t943 * t665 + t945 * t669;
t701 = pkin(4) * t818 + t714;
t702 = qJ(5) * t818 + t716;
t673 = t941 * t701 + t940 * t702 + t1116;
t676 = -t940 * t701 + t941 * t702 - t1120;
t1005 = -t943 * t673 + t945 * t676;
t685 = t1116 + t700;
t691 = -t699 - t1120;
t1004 = -t943 * t685 + t945 * t691;
t703 = -t1061 + t950 + 0.2e1 * t1097;
t689 = t1138 + t940 * t703 + (pkin(3) + t1062) * t1085;
t694 = -pkin(4) * t1054 + t941 * t703 - t1139;
t1003 = -t943 * t689 + t945 * t694;
t704 = (-t850 - t819) * pkin(4) + t949;
t692 = t1007 * t819 + t941 * t704 + t1099;
t697 = -qJ(5) * t1048 - t940 * t704 - t1098;
t1002 = -t943 * t692 + t945 * t697;
t717 = -pkin(3) * t819 + t1049 + t1099;
t730 = -t1056 - t1098;
t1001 = -t943 * t717 + t945 * t730;
t719 = -pkin(3) * t1085 - t1056 - t1138;
t745 = -t1049 + t1139;
t1000 = -t943 * t719 + t945 * t745;
t880 = t948 * pkin(1) + t958;
t889 = -t966 + t1060;
t997 = -t946 * t880 - t944 * t889;
t996 = -t944 * t916 - t946 * t917;
t994 = t944 * t923;
t993 = t946 * t923;
t910 = -t944 * t948 + t1019;
t992 = pkin(5) * t910 + t944 * g(3);
t911 = t946 * t948 + t1020;
t991 = -pkin(5) * t911 + t946 * g(3);
t807 = t900 * t1026 + t940 * t851;
t808 = t941 * t851 - t874;
t749 = t945 * t807 + t943 * t808;
t990 = -t946 * t749 - t1017;
t988 = -t999 + t785;
t985 = pkin(2) * t906 - t856;
t908 = t927 - 0.2e1 * t1011;
t984 = pkin(2) * t908 + t1044;
t784 = -t943 * t854 - t945 * t855;
t982 = t944 * t880 - t946 * t889;
t981 = t946 * t916 - t944 * t917;
t666 = t945 * t687 + t943 * t688;
t980 = -pkin(2) * t666 - t1015;
t979 = -t1013 - t1135;
t978 = -t1014 + t1151;
t922 = -t947 - t1057;
t861 = t945 * t922 - t1043;
t977 = -pkin(2) * t861 - t855;
t695 = pkin(3) * t802 + qJ(4) * t700;
t976 = -qJ(4) * t1039 - t943 * t695;
t975 = t944 * t749 - t1016;
t816 = pkin(2) * t819;
t974 = -t945 * t692 - t943 * t697 + t816;
t973 = -t945 * t717 - t943 * t730 + t816;
t969 = -t988 - t1066;
t967 = pkin(4) * t837 + qJ(5) * t1084 + t983;
t965 = -pkin(2) * t718 - t945 * t665 - t943 * t669;
t964 = -t945 * t673 - t943 * t676 - t1110;
t963 = -t945 * t685 - t943 * t691 - t1110;
t962 = -t945 * t689 - t943 * t694 - t1108;
t961 = -t945 * t719 - t943 * t745 + t1108;
t959 = -pkin(2) * t802 + qJ(4) * t1046 - t945 * t695;
t957 = t967 + t1075;
t956 = -t775 - t967 - t1151;
t885 = -0.2e1 * t1024;
t954 = -t900 * t835 - t1076 + t885 - t968;
t952 = t785 + t954;
t859 = t943 * t920 + t1037;
t951 = -pkin(2) * t859 - t854;
t921 = t947 - t1057;
t919 = -t947 + t1058;
t913 = (-t937 + t938) * t948;
t912 = t1030 * t948;
t904 = t1030 * t1023;
t903 = t989 - 0.2e1 * t1060;
t897 = 0.2e1 * t935 - t960;
t884 = -0.2e1 * t1025;
t872 = t938 * t1023 + t943 * t907;
t871 = t937 * t1023 + t945 * t971;
t866 = -t943 * t922 - t1038;
t865 = -t943 * t921 + t1037;
t864 = (t907 - t1011) * t945;
t863 = t945 * t920 - t1042;
t862 = t945 * t919 - t1043;
t860 = t945 * t921 + t1042;
t858 = t943 * t919 + t1038;
t853 = -t945 * t906 - t943 * t908;
t852 = t945 * t908 - t857;
t830 = pkin(1) * t889 - qJ(2) * t880;
t780 = pkin(2) * t912 + t784;
t773 = -qJ(2) * t866 - t977;
t772 = -qJ(2) * t863 - t951;
t771 = -t1067 * t863 + t985;
t770 = -t1067 * t866 + t984;
t769 = qJ(2) * t908 - t1067 * t861 - t856;
t768 = qJ(2) * t906 - t1067 * t859 - t1044;
t763 = -qJ(2) * t912 + t1067 * t909 - t783;
t752 = -t943 * t807 + t945 * t808;
t743 = -qJ(2) * t784 + t1065;
t725 = -t1067 * t784 - t1064;
t724 = -qJ(2) * t870 - t1067 * t783;
t684 = t884 - t978 - t1148;
t683 = -t1008 - t1130;
t682 = t885 - t969 - t1121;
t681 = t1066 + t952 - t1121;
t680 = t945 * t700 - t1046;
t678 = t961 - t1147;
t677 = t1000 + t1140;
t674 = -t979 - t1130;
t672 = t1075 - t956 + t1148;
t671 = t973 - t1117;
t670 = t1001 + t1115;
t667 = -t943 * t687 + t945 * t688;
t664 = t974 - t1117;
t663 = t1002 + t1115;
t662 = t962 + t1147;
t661 = t1003 - t1140;
t660 = t963 - t1129;
t659 = t1004 + t1128;
t658 = -qJ(2) * t680 - t1009;
t657 = t964 - t1129;
t656 = t1005 + t1128;
t655 = -t1067 * t680 + t959;
t654 = -qJ(2) * t802 - t1067 * t679 + t976;
t653 = -qJ(2) * t667 - t980;
t652 = -t1067 * t667 + t965;
t651 = -qJ(2) * t718 - t1067 * t666 + t1006;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t910, 0, -t911, 0, -t992, -t991, -t981, -pkin(5) * t981, 0, -t910, t911, 0, 0, 0, t982, t992, t991, pkin(5) * t982 + (-t944 * pkin(1) + t946 * qJ(2)) * g(3), t944 * t872 + t993, t944 * t852 + t946 * t913, t1019 * t945 + t944 * t860, t944 * t871 - t993, -t1019 * t943 + t944 * t858, -t944 * t904 + t928, t946 * t772 - t944 * t771 - pkin(5) * (-t946 * t859 + t944 * t906), t946 * t773 - t944 * t770 - pkin(5) * (-t946 * t861 + t944 * t908), -pkin(2) * t1034 + t944 * t780 - pkin(5) * (-t944 * t912 + t1034), t946 * t743 - t944 * t725 - pkin(5) * (-t946 * t783 - t944 * t870), t975, -t1144, t1137, t1079, t1146, t1081, -t944 * t671 + t946 * t682 - t1114, -t944 * t678 + t946 * t684 - t1153, -t944 * t660 + t946 * t683 - t1126, t946 * t658 - t944 * t655 - pkin(5) * (-t946 * t679 - t944 * t802), t975, t1137, t1144, t1081, -t1146, t1079, -t944 * t664 + t946 * t681 - t1114, -t944 * t657 + t946 * t674 - t1126, -t944 * t662 + t946 * t672 + t1153, t946 * t653 - t944 * t652 - pkin(5) * (-t946 * t666 - t944 * t718); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t911, 0, t910, 0, t991, -t992, t996, pkin(5) * t996, 0, -t911, -t910, 0, 0, 0, t997, -t991, t992, pkin(5) * t997 + (t946 * pkin(1) + t944 * qJ(2)) * g(3), -t946 * t872 + t994, -t946 * t852 + t944 * t913, -t946 * t860 + t927 * t944, -t946 * t871 - t994, -t1020 * t943 - t946 * t858, t946 * t904 + t926, t944 * t772 + t946 * t771 + pkin(5) * (t944 * t859 + t946 * t906), t944 * t773 + t946 * t770 + pkin(5) * (t944 * t861 + t946 * t908), -pkin(2) * t1040 - t946 * t780 + pkin(5) * (-t946 * t912 - t1040), t944 * t743 + t946 * t725 + pkin(5) * (t944 * t783 - t946 * t870), t990, t1143, t1136, t1078, -t1145, t1082, t946 * t671 + t944 * t682 + t1113, t946 * t678 + t944 * t684 - t1152, t946 * t660 + t944 * t683 + t1127, t944 * t658 + t946 * t655 + pkin(5) * (t944 * t679 - t946 * t802), t990, t1136, -t1143, t1082, t1145, t1078, t946 * t664 + t944 * t681 + t1113, t946 * t657 + t944 * t674 + t1127, t946 * t662 + t944 * t672 + t1152, t944 * t653 + t946 * t652 + pkin(5) * (t944 * t666 - t946 * t718); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t916, t917, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t903, t897, t830, t864, t853, t865, t857, t862, 0, t768, t769, t763, t724, t752, -t709, t1124, t1070, -t738, t1071, t670, t677, t659, t654, t752, t1124, t709, t1071, t738, t1070, t663, t656, t661, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t948, 0, 0, -g(3), -t916, 0, 0, -qJDD(1), t948, 0, 0, 0, -t889, 0, g(3), qJ(2) * g(3), t923, t913, t927, -t923, -t1021, qJDD(3), t772, t773, -t1063, t743, -t1059, t1073, t1086, t1059, -t820, qJDD(3), t682, t684, t683, t658, -t1059, t1086, -t1073, qJDD(3), t820, t1059, t681, t674, t672, t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t948, 0, qJDD(1), 0, g(3), 0, -t917, 0, 0, -t948, -qJDD(1), 0, 0, 0, -t880, -g(3), 0, pkin(1) * g(3), -t872, -t852, -t860, -t871, -t858, t904, t771, t770, -t780, t725, -t749, t705, -t1125, -t1069, -t731, -t1072, t671, t678, t660, t655, -t749, -t1125, -t705, -t1072, t731, -t1069, t664, t657, t662, t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t916, t917, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t903, t897, t830, t864, t853, t865, t857, t862, 0, t768, t769, t763, t724, t752, -t709, t1124, t1070, -t738, t1071, t670, t677, t659, t654, t752, t1124, t709, t1071, t738, t1070, t663, t656, t661, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t889, -t880, 0, t864, t853, t865, t857, t862, 0, -pkin(6) * t859 - t1044, -pkin(6) * t861 - t856, pkin(6) * t909 - t783, -pkin(6) * t783, t752, -t709, t1124, t1070, -t738, t1071, t1001 - t1123, t1000 + t1150, t1004 - t1134, -pkin(6) * t679 + t976, t752, t1124, t709, t1071, t738, t1070, t1002 - t1123, t1005 - t1134, t1003 - t1150, -pkin(6) * t666 + t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t948, 0, 0, 0, t889, 0, -g(3), 0, -t923, -t913, -t927, t923, t1021, -qJDD(3), t951, t977, t1063, -t1065, t1059, -t1073, -t1086, -t1059, t820, -qJDD(3), t969 + t1018, t883 + t978, t1008, t1009, t1059, -t1086, t1073, -qJDD(3), -t820, -t1059, t716 - t785 - t1066 + t1076, t979, t884 + t956 - 0.2e1 * t1022, t980; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t948, qJDD(1), 0, 0, 0, t880, g(3), 0, 0, t872, t852, t860, t871, t858, -t904, pkin(6) * t863 - t985, pkin(6) * t866 - t984, t780, pkin(6) * t784 + t1064, t749, -t705, t1125, t1069, t731, t1072, -t973 + t1122, -t961 + t1149, -t963 + t1133, pkin(6) * t680 - t959, t749, t1125, t705, t1072, -t731, t1069, -t974 + t1122, -t964 + t1133, -t962 - t1149, pkin(6) * t667 - t965; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t907, -t906, t915, t1011, t919, -t1011, 0, -t870, -t854, 0, t808, -t759, t1119, t972, -t795, t986, t730, t745, t691, -qJ(4) * t699, t808, t1119, t759, t986, t795, t972, t697, t676, t694, t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1010, t908, t921, t971, t914, -t1010, t870, 0, -t855, 0, t807, t757, t1118, t987, t788, t970, t717, t719, t685, t695, t807, t1118, -t757, t970, -t788, t987, t692, t673, t689, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t923, t913, t927, -t923, -t1021, qJDD(3), t854, t855, 0, 0, -t1059, t1073, t1086, t1059, -t820, qJDD(3), t885 + t988, t884 + t1014, t754, t698, -t1059, t1086, -t1073, qJDD(3), t820, t1059, t952, t1013, t775 + t957, t1015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, -t819, t1083, -t1028, t875, t1028, 0, -t802, t727, 0, t851, t1083, t819, t1028, -t875, -t1028, -qJ(5) * t819, t702, t703, qJ(5) * t718; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t888, t1085, -t876, -t850, t1084, -t888, t802, 0, t728, 0, t888, -t876, -t1085, -t888, -t1084, -t850, t704, t701, pkin(4) * t1085, pkin(4) * t718; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1059, t1073, t1086, t1059, -t820, qJDD(3), -t727, -t728, 0, 0, -t1059, t1086, -t1073, qJDD(3), t820, t1059, t954, t1031, t957, t1033; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, t1083, t819, t1028, -t875, -t1028, 0, t716, t718, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1059, t1086, -t1073, qJDD(3), t820, t1059, -t716, 0, t714, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t888, t876, t1085, t888, t1084, t850, -t718, -t714, 0, 0;];
m_new_reg = t1;