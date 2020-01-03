% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRRP13_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:53
% EndTime: 2019-12-31 19:00:04
% DurationCPUTime: 11.46s
% Computational Cost: add. (21035->536), mult. (40837->566), div. (0->0), fcn. (24530->6), ass. (0->341)
t980 = cos(qJ(3));
t1060 = qJD(1) * t980;
t976 = sin(qJ(4));
t979 = cos(qJ(4));
t935 = -t979 * qJD(3) + t976 * t1060;
t937 = t976 * qJD(3) + t979 * t1060;
t890 = t937 * t935;
t1058 = qJD(1) * qJD(3);
t1046 = t980 * t1058;
t977 = sin(qJ(3));
t1056 = t977 * qJDD(1);
t940 = -t1046 - t1056;
t933 = qJDD(4) - t940;
t1123 = t890 + t933;
t1136 = t976 * t1123;
t963 = t977 * qJD(1) + qJD(4);
t1106 = t963 ^ 2;
t1107 = t935 ^ 2;
t908 = t1107 - t1106;
t820 = t979 * t908 - t1136;
t1047 = t977 * t1058;
t967 = t980 * qJDD(1);
t941 = t967 - t1047;
t1035 = t979 * qJDD(3) - t976 * t941;
t875 = t937 * qJD(4) - t1035;
t916 = t963 * t937;
t846 = t875 - t916;
t765 = t977 * t820 + t980 * t846;
t1135 = t979 * t1123;
t816 = t976 * t908 + t1135;
t978 = sin(qJ(1));
t981 = cos(qJ(1));
t1189 = t978 * t765 + t981 * t816;
t932 = t937 ^ 2;
t1117 = t932 - t1107;
t1094 = t935 * t963;
t1010 = -t976 * qJDD(3) - t979 * t941;
t876 = -t935 * qJD(4) - t1010;
t1119 = t876 - t1094;
t1088 = t976 * t1119;
t1125 = t875 + t916;
t782 = t979 * t1125 + t1088;
t746 = t1117 * t980 + t977 * t782;
t1074 = t979 * t1119;
t780 = -t976 * t1125 + t1074;
t1188 = t978 * t746 - t981 * t780;
t1187 = t981 * t765 - t978 * t816;
t1186 = t981 * t746 + t978 * t780;
t886 = -t932 - t1106;
t810 = t979 * t886 - t1136;
t1185 = pkin(2) * t810;
t1184 = pkin(3) * t810;
t803 = t976 * t886 + t1135;
t1183 = pkin(7) * t803;
t1182 = pkin(7) * t810;
t1181 = qJ(2) * t810;
t1180 = t977 * t803;
t1178 = t978 * t810;
t1176 = t980 * t803;
t1174 = t981 * t810;
t769 = t980 * t820 - t977 * t846;
t748 = -t1117 * t977 + t980 * t782;
t1124 = -t890 + t933;
t1085 = t976 * t1124;
t909 = t932 - t1106;
t1150 = -t979 * t909 + t1085;
t1118 = t876 + t1094;
t1071 = t979 * t1124;
t1149 = t976 * t909 + t1071;
t1163 = -t1118 * t980 + t1149 * t977;
t1172 = t1150 * t981 + t1163 * t978;
t1171 = t1150 * t978 - t1163 * t981;
t1115 = -t1106 - t1107;
t1132 = t1115 * t979 - t1085;
t1148 = -t1125 * t980 + t1132 * t977;
t1170 = pkin(2) * t1148;
t1147 = t1125 * t977 + t1132 * t980;
t1169 = pkin(6) * t1147;
t1168 = pkin(6) * t1148;
t1167 = qJ(2) * t1147;
t1105 = pkin(6) + pkin(1);
t1164 = t1105 * t1147;
t1162 = t1118 * t977 + t1149 * t980;
t1133 = t1115 * t976 + t1071;
t1161 = qJ(2) * t1133 - t1105 * t1148;
t1160 = pkin(5) * (t1133 * t981 + t1148 * t978);
t1159 = pkin(5) * (t1133 * t978 - t1148 * t981);
t1158 = pkin(2) * t1133;
t1157 = pkin(3) * t1133;
t1156 = pkin(7) * t1132;
t1155 = pkin(7) * t1133;
t1154 = -qJ(5) * t976 - pkin(3);
t1116 = t932 + t1107;
t1146 = pkin(3) * t1116;
t1145 = t1116 * t977;
t1144 = t1116 * t980;
t1139 = t1119 * qJ(5);
t953 = t981 * g(1) + t978 * g(2);
t973 = qJDD(1) * qJ(2);
t1008 = t953 - t973;
t1018 = -t941 + t1047;
t1019 = -t940 + t1046;
t983 = qJD(1) ^ 2;
t1120 = t1105 * t983;
t1057 = qJD(2) * qJD(1);
t971 = 0.2e1 * t1057;
t833 = t1019 * pkin(3) + t1018 * pkin(7) - t1008 - t1120 + t971;
t1099 = pkin(7) * t980;
t1100 = pkin(3) * t977;
t1005 = t983 * (-t1099 + t1100);
t952 = t978 * g(1) - t981 * g(2);
t1025 = qJDD(2) - t952;
t999 = -t983 * qJ(2) + t1025;
t912 = -t1105 * qJDD(1) + t999;
t884 = t980 * g(3) - t977 * t912;
t982 = qJD(3) ^ 2;
t858 = -t982 * pkin(3) + qJDD(3) * pkin(7) - t1005 * t977 - t884;
t775 = t976 * t833 + t979 * t858;
t882 = t935 * pkin(4) - t937 * qJ(5);
t1017 = -pkin(4) * t1106 + t933 * qJ(5) - t935 * t882 + t775;
t1059 = qJD(5) * t963;
t950 = 0.2e1 * t1059;
t736 = t950 + t1017;
t774 = -t979 * t833 + t976 * t858;
t737 = -t933 * pkin(4) - qJ(5) * t1106 + t937 * t882 + qJDD(5) + t774;
t710 = t979 * t736 + t976 * t737;
t1098 = t875 * pkin(4);
t883 = t977 * g(3) + t980 * t912;
t857 = qJDD(3) * pkin(3) + t982 * pkin(7) - t980 * t1005 + t883;
t985 = -pkin(4) * t916 + 0.2e1 * qJD(5) * t937 + t857;
t984 = t985 + t1139;
t739 = t984 - t1098;
t704 = t977 * t710 + t980 * t739;
t989 = pkin(7) * t710 + (pkin(4) * t979 - t1154) * t739;
t1134 = -pkin(2) * t704 - t989;
t1092 = t963 * t979;
t1052 = t935 * t1092;
t1006 = t976 * t875 + t1052;
t1051 = t980 * t890;
t1112 = t1006 * t977 + t1051;
t1093 = t963 * t976;
t834 = t935 * t1093 - t979 * t875;
t1131 = -t1112 * t981 + t978 * t834;
t1130 = t1112 * t978 + t981 * t834;
t903 = t937 * t1093;
t1026 = t903 - t1052;
t1110 = t1026 * t977 - t980 * t933;
t1121 = (t935 * t976 + t937 * t979) * t963;
t1129 = -t1110 * t981 - t1121 * t978;
t1128 = t1110 * t978 - t1121 * t981;
t724 = t976 * t774 + t979 * t775;
t1048 = pkin(3) * t980 + pkin(2);
t723 = -t979 * t774 + t976 * t775;
t1122 = (pkin(7) * t977 + t1048) * t723;
t831 = t980 * t883 - t977 * t884;
t839 = t979 * t876 - t903;
t1007 = t977 * t839 - t1051;
t838 = t937 * t1092 + t976 * t876;
t1114 = -t1007 * t981 + t978 * t838;
t1113 = t1007 * t978 + t981 * t838;
t1053 = t977 * t890;
t1111 = t1006 * t980 - t1053;
t1109 = t1026 * t980 + t977 * t933;
t1073 = t979 * t1118;
t847 = (-qJD(4) + t963) * t937 + t1035;
t781 = t976 * t847 - t1073;
t716 = -pkin(7) * t781 - t723;
t1108 = -t1048 * t781 + t977 * t716;
t1103 = pkin(2) * t831;
t992 = t1008 - 0.2e1 * t1057;
t905 = t992 + t1120;
t1102 = pkin(2) * t905;
t974 = t977 ^ 2;
t975 = t980 ^ 2;
t1061 = t974 + t975;
t943 = t1061 * qJDD(1);
t1101 = pkin(2) * t943;
t1096 = qJ(5) * t979;
t1095 = qJDD(1) * pkin(1);
t1091 = t974 * t983;
t1090 = t975 * t983;
t1087 = t976 * t1118;
t853 = t976 * t857;
t1080 = t977 * t905;
t961 = t977 * t983 * t980;
t948 = qJDD(3) + t961;
t1078 = t977 * t948;
t949 = qJDD(3) - t961;
t1077 = t977 * t949;
t1076 = t978 * t943;
t854 = t979 * t857;
t891 = t980 * t905;
t1067 = t980 * t948;
t1066 = t980 * t949;
t1065 = t981 * t943;
t1063 = t781 * t1100 + t980 * t716;
t1062 = pkin(3) * t857 + pkin(7) * t724;
t1055 = t978 * qJDD(1);
t1054 = t981 * qJDD(1);
t852 = (qJD(4) + t963) * t935 + t1010;
t1050 = pkin(3) * t852 - t1183 - t853;
t1049 = -pkin(3) * t1125 + t1156 + t854;
t1022 = -pkin(4) * t737 + qJ(5) * t736;
t709 = t976 * t736 - t979 * t737;
t692 = -pkin(3) * t709 - t1022;
t694 = -pkin(7) * t709 + (-pkin(4) * t976 + t1096) * t739;
t1045 = -t977 * t692 + t980 * t694;
t730 = pkin(4) * t1116 + t736;
t733 = qJ(5) * t1116 + t737;
t779 = -t976 * t846 - t1073;
t703 = -pkin(7) * t779 - t976 * t730 + t979 * t733;
t1021 = -pkin(4) * t1118 - qJ(5) * t846;
t731 = -pkin(3) * t779 - t1021;
t1044 = t980 * t703 - t977 * t731;
t728 = -t1098 + t985 + 0.2e1 * t1139;
t712 = -pkin(4) * t1088 + t979 * t728 + t1182;
t988 = -pkin(4) * t886 + qJ(5) * t1123 + t1017;
t717 = -0.2e1 * t1059 - t988 + t1184;
t1043 = t980 * t712 - t977 * t717;
t729 = (-t1125 - t875) * pkin(4) + t984;
t714 = -t1096 * t1125 - t976 * t729 - t1155;
t986 = pkin(4) * t1124 + qJ(5) * t1115 - t737;
t718 = -t986 - t1157;
t1042 = t980 * t714 - t977 * t718;
t734 = t774 - t1157;
t759 = -t853 - t1155;
t1041 = -t977 * t734 + t980 * t759;
t735 = t775 - t1184;
t762 = -t854 - t1182;
t1040 = -t977 * t735 + t980 * t762;
t917 = t983 * pkin(1) + t992;
t918 = -t999 + t1095;
t1037 = -t981 * t917 - t978 * t918;
t1036 = -t978 * t952 - t981 * t953;
t1034 = t978 * t961;
t1033 = t981 * t961;
t783 = -t979 * t846 + t1087;
t1032 = pkin(7) * t783 + t979 * t730 + t976 * t733 + t1146;
t785 = t979 * t847 + t1087;
t1031 = pkin(7) * t785 + t1146 + t724;
t719 = t977 * t724 + t980 * t857;
t1030 = -pkin(2) * t719 - t1062;
t944 = -t978 * t983 + t1054;
t1029 = pkin(5) * t944 + t978 * g(3);
t945 = t981 * t983 + t1055;
t1028 = -pkin(5) * t945 + t981 * g(3);
t1027 = t980 * t839 + t1053;
t939 = 0.2e1 * t1046 + t1056;
t1024 = pkin(2) * t939 - t891;
t942 = t967 - 0.2e1 * t1047;
t1023 = pkin(2) * t942 + t1080;
t832 = -t977 * t883 - t980 * t884;
t1016 = t978 * t917 - t981 * t918;
t1014 = t981 * t952 - t978 * t953;
t1012 = -t1049 - t1170;
t755 = t980 * t852 - t1180;
t1011 = -pkin(2) * t755 - t1050;
t959 = -t982 - t1090;
t896 = t980 * t959 - t1078;
t1009 = -pkin(2) * t896 - t884;
t1004 = pkin(3) * t1119 + pkin(4) * t1074 + t976 * t728 + t1183;
t740 = t977 * t783 + t1144;
t1003 = -pkin(2) * t740 - t1032;
t741 = t977 * t785 + t1144;
t1002 = -pkin(2) * t741 - t1031;
t1001 = t1125 * t1154 + t979 * t729 + t1156;
t957 = -t982 - t1091;
t894 = t977 * t957 + t1066;
t1000 = -pkin(2) * t894 - t883;
t998 = pkin(2) * t709 - t980 * t692 - t977 * t694;
t997 = pkin(2) * t779 - t977 * t703 - t980 * t731;
t996 = -t977 * t712 - t980 * t717 - t1185;
t995 = -t977 * t714 - t980 * t718 + t1158;
t994 = -t980 * t734 - t977 * t759 + t1158;
t993 = -t980 * t735 - t977 * t762 + t1185;
t750 = t1119 * t980 + t1180;
t991 = -pkin(2) * t750 - t1004;
t990 = -t1001 - t1170;
t958 = t982 - t1090;
t956 = -t982 + t1091;
t947 = (-t974 + t975) * t983;
t946 = t1061 * t983;
t934 = t1061 * t1058;
t929 = t1025 - 0.2e1 * t1095;
t924 = -t953 + t971 + 0.2e1 * t973;
t907 = t975 * t1058 + t977 * t941;
t906 = t974 * t1058 + t980 * t940;
t901 = -t977 * t959 - t1067;
t900 = -t977 * t958 + t1066;
t899 = t1018 * t980;
t898 = t980 * t957 - t1077;
t897 = t980 * t956 - t1078;
t895 = t980 * t958 + t1077;
t893 = t977 * t956 + t1067;
t892 = t1019 * t977;
t881 = -t980 * t939 - t977 * t942;
t880 = -t977 * t939 + t980 * t942;
t867 = pkin(1) * t918 - qJ(2) * t917;
t809 = pkin(2) * t946 + t832;
t806 = -qJ(2) * t901 - t1009;
t805 = -qJ(2) * t898 - t1000;
t798 = -t1105 * t898 + t1024;
t797 = -t1105 * t901 + t1023;
t796 = qJ(2) * t942 - t1105 * t896 - t891;
t795 = qJ(2) * t939 - t1105 * t894 - t1080;
t786 = -qJ(2) * t946 + t1105 * t943 - t831;
t761 = -qJ(2) * t832 + t1103;
t757 = -t977 * t852 - t1176;
t752 = -t1119 * t977 + t1176;
t745 = -t1105 * t832 - t1102;
t744 = -qJ(2) * t905 - t1105 * t831;
t743 = t980 * t785 - t1145;
t742 = t980 * t783 - t1145;
t721 = t723 * t1100;
t720 = t980 * t724 - t977 * t857;
t707 = -qJ(2) * t757 - t1011;
t706 = -t1012 - t1167;
t705 = t980 * t710 - t977 * t739;
t701 = -t1105 * t757 + t993;
t700 = -t1105 * t755 + t1040 + t1181;
t699 = -t990 - t1167;
t698 = t994 - t1164;
t697 = t1041 + t1161;
t696 = -qJ(2) * t752 - t991;
t695 = -qJ(2) * t743 - t1002;
t691 = -qJ(2) * t720 - t1030;
t690 = -t1105 * t743 - t1108;
t689 = qJ(2) * t781 - t1105 * t741 + t1063;
t688 = -qJ(2) * t742 - t1003;
t687 = t995 - t1164;
t686 = t1042 + t1161;
t685 = -t1105 * t752 + t996;
t684 = -t1105 * t750 + t1043 - t1181;
t683 = -t1105 * t742 + t997;
t682 = qJ(2) * t779 - t1105 * t740 + t1044;
t681 = -t1105 * t720 + t1122;
t680 = t721 + (qJ(2) - t1099) * t723 - t1105 * t719;
t679 = -qJ(2) * t705 - t1134;
t678 = -t1105 * t705 + t998;
t677 = qJ(2) * t709 - t1105 * t704 + t1045;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t944, 0, -t945, 0, -t1029, -t1028, -t1014, -pkin(5) * t1014, 0, -t944, t945, 0, 0, 0, t1016, t1029, t1028, pkin(5) * t1016 + (-t978 * pkin(1) + t981 * qJ(2)) * g(3), t978 * t907 + t1033, t978 * t880 + t981 * t947, t1054 * t980 + t978 * t895, t978 * t906 - t1033, -t1054 * t977 + t978 * t893, t981 * qJDD(3) - t978 * t934, t981 * t805 - t978 * t798 - pkin(5) * (-t981 * t894 + t978 * t939), t981 * t806 - t978 * t797 - pkin(5) * (-t981 * t896 + t978 * t942), -pkin(2) * t1065 + t978 * t809 - pkin(5) * (-t978 * t946 + t1065), t981 * t761 - t978 * t745 - pkin(5) * (-t981 * t831 - t978 * t905), t1113, -t1188, t1172, t1130, t1189, t1128, -t978 * t698 + t981 * t706 - t1159, t981 * t707 - t978 * t701 - pkin(5) * (-t981 * t755 + t1178), t981 * t695 - t978 * t690 - pkin(5) * (-t981 * t741 + t978 * t781), t981 * t691 - t978 * t681 - pkin(5) * (-t981 * t719 + t978 * t723), t1113, t1172, t1188, t1128, -t1189, t1130, -t978 * t687 + t981 * t699 - t1159, t981 * t688 - t978 * t683 - pkin(5) * (-t981 * t740 + t978 * t779), t981 * t696 - t978 * t685 - pkin(5) * (-t981 * t750 - t1178), t981 * t679 - t978 * t678 - pkin(5) * (-t981 * t704 + t978 * t709); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t945, 0, t944, 0, t1028, -t1029, t1036, pkin(5) * t1036, 0, -t945, -t944, 0, 0, 0, t1037, -t1028, t1029, pkin(5) * t1037 + (t981 * pkin(1) + t978 * qJ(2)) * g(3), -t981 * t907 + t1034, -t981 * t880 + t978 * t947, -t981 * t895 + t967 * t978, -t981 * t906 - t1034, -t1055 * t977 - t981 * t893, t978 * qJDD(3) + t981 * t934, t978 * t805 + t981 * t798 + pkin(5) * (t978 * t894 + t981 * t939), t978 * t806 + t981 * t797 + pkin(5) * (t978 * t896 + t981 * t942), -pkin(2) * t1076 - t981 * t809 + pkin(5) * (-t981 * t946 - t1076), t978 * t761 + t981 * t745 + pkin(5) * (t978 * t831 - t981 * t905), t1114, t1186, t1171, t1131, -t1187, t1129, t981 * t698 + t978 * t706 + t1160, t978 * t707 + t981 * t701 + pkin(5) * (t978 * t755 + t1174), t978 * t695 + t981 * t690 + pkin(5) * (t978 * t741 + t981 * t781), t978 * t691 + t981 * t681 + pkin(5) * (t978 * t719 + t981 * t723), t1114, t1171, -t1186, t1129, t1187, t1131, t981 * t687 + t978 * t699 + t1160, t978 * t688 + t981 * t683 + pkin(5) * (t978 * t740 + t981 * t779), t978 * t696 + t981 * t685 + pkin(5) * (t978 * t750 - t1174), t978 * t679 + t981 * t678 + pkin(5) * (t978 * t704 + t981 * t709); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t952, t953, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t929, t924, t867, -t899, t881, t900, t892, t897, 0, t795, t796, t786, t744, t1027, -t748, t1162, t1111, t769, t1109, t697, t700, t689, t680, t1027, t1162, t748, t1109, -t769, t1111, t686, t682, t684, t677; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t983, 0, 0, -g(3), -t952, 0, 0, -qJDD(1), t983, 0, 0, 0, -t918, 0, g(3), qJ(2) * g(3), t961, t947, t967, -t961, -t1056, qJDD(3), t805, t806, -t1101, t761, t838, t780, t1150, t834, t816, -t1121, t706, t707, t695, t691, t838, t1150, -t780, -t1121, -t816, t834, t699, t688, t696, t679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t983, 0, qJDD(1), 0, g(3), 0, -t953, 0, 0, -t983, -qJDD(1), 0, 0, 0, -t917, -g(3), 0, pkin(1) * g(3), -t907, -t880, -t895, -t906, -t893, t934, t798, t797, -t809, t745, -t1007, t746, -t1163, -t1112, -t765, -t1110, t698, t701, t690, t681, -t1007, -t1163, -t746, -t1110, t765, -t1112, t687, t683, t685, t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t952, t953, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t929, t924, t867, -t899, t881, t900, t892, t897, 0, t795, t796, t786, t744, t1027, -t748, t1162, t1111, t769, t1109, t697, t700, t689, t680, t1027, t1162, t748, t1109, -t769, t1111, t686, t682, t684, t677; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t918, -t917, 0, -t899, t881, t900, t892, t897, 0, -pkin(6) * t894 - t1080, -pkin(6) * t896 - t891, pkin(6) * t943 - t831, -pkin(6) * t831, t1027, -t748, t1162, t1111, t769, t1109, t1041 - t1168, -pkin(6) * t755 + t1040, -pkin(6) * t741 + t1063, -pkin(6) * t719 - t1099 * t723 + t721, t1027, t1162, t748, t1109, -t769, t1111, t1042 - t1168, -pkin(6) * t740 + t1044, -pkin(6) * t750 + t1043, -pkin(6) * t704 + t1045; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t983, 0, 0, 0, t918, 0, -g(3), 0, -t961, -t947, -t967, t961, t1056, -qJDD(3), t1000, t1009, t1101, -t1103, -t838, -t780, -t1150, -t834, -t816, t1121, t1012, t1011, t1002, t1030, -t838, -t1150, t780, t1121, t816, -t834, t990, t1003, t991, t1134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t983, qJDD(1), 0, 0, 0, t917, g(3), 0, 0, t907, t880, t895, t906, t893, -t934, pkin(6) * t898 - t1024, pkin(6) * t901 - t1023, t809, pkin(6) * t832 + t1102, t1007, -t746, t1163, t1112, t765, t1110, -t994 + t1169, pkin(6) * t757 - t993, pkin(6) * t743 + t1108, pkin(6) * t720 - t1122, t1007, t1163, t746, t1110, -t765, t1112, -t995 + t1169, pkin(6) * t742 - t997, pkin(6) * t752 - t996, pkin(6) * t705 - t998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t941, -t939, t949, t1047, t956, -t1047, 0, -t905, -t883, 0, t839, -t782, t1149, t1006, t820, t1026, t759, t762, t716, -pkin(7) * t723, t839, t1149, t782, t1026, -t820, t1006, t714, t703, t712, t694; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1046, t942, t958, t940, t948, -t1046, t905, 0, -t884, 0, -t890, -t1117, -t1118, t890, t846, -t933, t734, t735, -pkin(3) * t781, -pkin(3) * t723, -t890, -t1118, t1117, -t933, -t846, t890, t718, t731, t717, t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t961, t947, t967, -t961, -t1056, qJDD(3), t883, t884, 0, 0, t838, t780, t1150, t834, t816, -t1121, t1049, t1050, t1031, t1062, t838, t1150, -t780, -t1121, -t816, t834, t1001, t1032, t1004, t989; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876, -t1125, t1124, t1094, t908, -t1094, 0, -t857, t774, 0, t876, t1124, t1125, -t1094, -t908, t1094, -qJ(5) * t1125, t733, t728, qJ(5) * t739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t916, t1119, -t909, -t875, t1123, -t916, t857, 0, t775, 0, t916, -t909, -t1119, -t916, -t1123, -t875, t729, t730, pkin(4) * t1119, pkin(4) * t739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t890, t1117, t1118, -t890, -t846, t933, -t774, -t775, 0, 0, t890, t1118, -t1117, t933, t846, -t890, t986, t1021, t950 + t988, t1022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876, t1124, t1125, -t1094, -t908, t1094, 0, t737, t739, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t890, t1118, -t1117, t933, t846, -t890, -t737, 0, t736, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t916, t909, t1119, t916, t1123, t875, -t739, -t736, 0, 0;];
m_new_reg = t1;
