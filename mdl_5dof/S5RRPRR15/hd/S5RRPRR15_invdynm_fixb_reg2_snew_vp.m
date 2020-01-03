% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRR15_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:31
% EndTime: 2019-12-31 20:43:45
% DurationCPUTime: 14.85s
% Computational Cost: add. (58084->648), mult. (120906->777), div. (0->0), fcn. (73884->8), ass. (0->420)
t1076 = sin(qJ(2));
t1071 = t1076 ^ 2;
t1082 = qJD(1) ^ 2;
t1065 = t1071 * t1082;
t1211 = qJD(2) ^ 2;
t1051 = -t1065 - t1211;
t1080 = cos(qJ(2));
t1162 = t1080 * t1082;
t1056 = t1076 * t1162;
t1044 = qJDD(2) - t1056;
t1163 = t1080 * t1044;
t1000 = t1076 * t1051 + t1163;
t1184 = qJD(1) * t1080;
t1062 = qJD(2) * t1184;
t1155 = t1076 * qJDD(1);
t1034 = 0.2e1 * t1062 + t1155;
t1077 = sin(qJ(1));
t1081 = cos(qJ(1));
t1253 = pkin(5) * (t1081 * t1000 - t1077 * t1034);
t1252 = pkin(5) * (t1077 * t1000 + t1081 * t1034);
t1209 = t1080 ^ 2;
t1066 = t1209 * t1082;
t1053 = -t1066 - t1211;
t1043 = t1056 + qJDD(2);
t1164 = t1080 * t1043;
t990 = t1076 * t1053 + t1164;
t1251 = pkin(1) * t990;
t1250 = pkin(6) * t990;
t1064 = t1080 * qJDD(1);
t1176 = t1076 * qJD(1);
t1144 = qJD(2) * t1176;
t1036 = t1064 - 0.2e1 * t1144;
t1173 = t1076 * t1043;
t999 = -t1080 * t1053 + t1173;
t1249 = pkin(5) * (t1077 * t1036 + t1081 * t999);
t1248 = pkin(5) * (-t1081 * t1036 + t1077 * t999);
t1247 = pkin(6) * t1000;
t1153 = t1081 * qJDD(1);
t1050 = -t1065 + t1211;
t994 = -t1076 * t1050 + t1164;
t1246 = -t1076 * t1153 + t1077 * t994;
t1154 = t1077 * qJDD(1);
t1245 = t1076 * t1154 + t1081 * t994;
t1242 = 2 * qJD(3);
t1172 = t1076 * t1044;
t992 = -t1080 * t1051 + t1172;
t1241 = pkin(1) * t992;
t1240 = pkin(6) * t992;
t1239 = pkin(6) * t999;
t1074 = sin(qJ(5));
t1035 = t1062 + t1155;
t1020 = qJDD(4) + t1035;
t1017 = qJDD(5) + t1020;
t1075 = sin(qJ(4));
t1079 = cos(qJ(4));
t1027 = t1075 * qJD(2) + t1079 * t1184;
t1029 = t1079 * qJD(2) - t1075 * t1184;
t1078 = cos(qJ(5));
t972 = t1078 * t1027 + t1074 * t1029;
t974 = -t1074 * t1027 + t1078 * t1029;
t921 = t974 * t972;
t1226 = t1017 - t921;
t1238 = t1074 * t1226;
t982 = t1029 * t1027;
t1225 = t1020 - t982;
t1237 = t1075 * t1225;
t1236 = t1078 * t1226;
t1235 = t1079 * t1225;
t1068 = t1076 * g(3);
t1047 = t1081 * g(1) + t1077 * g(2);
t1015 = -t1082 * pkin(1) + qJDD(1) * pkin(6) - t1047;
t1196 = qJ(3) * t1076;
t1201 = pkin(2) * t1080;
t1120 = -t1196 - t1201;
t1131 = t1082 * t1120 + t1015;
t1089 = -t1211 * pkin(2) + t1131 * t1080 - t1068;
t1086 = qJD(2) * t1242 + t1089;
t1234 = t1086 + (qJDD(2) + t1044) * qJ(3) - pkin(2) * t1051;
t1052 = t1066 - t1211;
t997 = -t1080 * t1052 + t1172;
t1233 = t1077 * t997 + t1080 * t1153;
t1231 = -t1077 * t1064 + t1081 * t997;
t1117 = -t1064 + t1144;
t1129 = t1075 * qJDD(2) - t1079 * t1117;
t1098 = t1029 * qJD(4) + t1129;
t968 = -t1027 * qJD(4) + t1079 * qJDD(2) + t1075 * t1117;
t883 = -t972 * qJD(5) - t1074 * t1098 + t1078 * t968;
t1058 = qJD(4) + t1176;
t1049 = qJD(5) + t1058;
t956 = t1049 * t972;
t1229 = -t956 + t883;
t1118 = t1035 + t1062;
t1228 = qJ(3) * t1118;
t1011 = t1058 * t1027;
t1227 = -t1011 + t968;
t1045 = pkin(3) * t1176 - qJD(2) * pkin(7);
t1046 = t1077 * g(1) - t1081 * g(2);
t1108 = -qJDD(1) * pkin(1) - t1046;
t1222 = -pkin(2) * t1144 + t1176 * t1242;
t1084 = t1108 - t1222 - t1228;
t1103 = t1117 * pkin(2);
t889 = t1103 + t1117 * pkin(7) - t1045 * t1176 + (-t1209 * pkin(3) - pkin(6)) * t1082 + t1084;
t1200 = t1080 * g(3);
t1110 = -qJDD(2) * pkin(2) - t1211 * qJ(3) + qJDD(3) + t1200;
t906 = -qJDD(2) * pkin(7) + (t1035 - t1062) * pkin(3) + (-pkin(7) * t1162 + t1131) * t1076 + t1110;
t835 = t1075 * t889 - t1079 * t906;
t836 = t1075 * t906 + t1079 * t889;
t1224 = t1075 * t836 - t1079 * t835;
t1223 = -t968 - t1011;
t988 = t1076 * t1052 + t1163;
t945 = t1076 * t1131 + t1110;
t1221 = -pkin(2) * t1043 - qJ(3) * t1053 + t945;
t820 = pkin(4) * t1225 + t1223 * pkin(8) - t835;
t1004 = t1058 * pkin(4) - t1029 * pkin(8);
t1210 = t1027 ^ 2;
t825 = -t1210 * pkin(4) - pkin(8) * t1098 - t1058 * t1004 + t836;
t774 = t1074 * t825 - t1078 * t820;
t775 = t1074 * t820 + t1078 * t825;
t754 = t1074 * t775 - t1078 * t774;
t1187 = t1079 * t754;
t755 = t1074 * t774 + t1078 * t775;
t1182 = qJDD(2) * qJ(3);
t905 = t1182 - t1117 * pkin(3) - pkin(7) * t1066 + (t1242 + t1045) * qJD(2) + t1089;
t837 = pkin(4) * t1098 - t1210 * pkin(8) + t1029 * t1004 + t905;
t749 = -pkin(4) * t837 + pkin(8) * t755;
t1109 = pkin(8) * t1187 + t1075 * t749;
t1208 = pkin(2) + pkin(7);
t736 = t1075 * t755 + t1187;
t1219 = qJ(3) * t837 - t1208 * t736 - t1109;
t934 = (qJD(4) - t1058) * t1029 + t1129;
t879 = -t1075 * t934 + t1079 * t1223;
t1019 = t1029 ^ 2;
t959 = -t1019 - t1210;
t1218 = qJ(3) * t959 - t1208 * t879 - t1224;
t1195 = t1074 * t837;
t901 = t1017 + t921;
t1188 = t1078 * t901;
t1048 = t1049 ^ 2;
t971 = t974 ^ 2;
t941 = -t971 - t1048;
t862 = -t1074 * t941 - t1188;
t779 = -pkin(4) * t1229 + pkin(8) * t862 + t1195;
t1189 = t1078 * t837;
t1194 = t1074 * t901;
t861 = t1078 * t941 - t1194;
t813 = -pkin(8) * t861 + t1189;
t1133 = t1075 * t779 - t1079 * t813;
t817 = t1075 * t862 + t1079 * t861;
t1217 = qJ(3) * t1229 - t1208 * t817 - t1133;
t970 = t972 ^ 2;
t909 = -t1048 - t970;
t845 = t1078 * t909 - t1238;
t1136 = -t1074 * t968 - t1078 * t1098;
t1102 = t974 * qJD(5) - t1136;
t957 = t1049 * t974;
t848 = t1102 + t957;
t778 = -pkin(4) * t848 + pkin(8) * t845 - t1189;
t844 = t1074 * t909 + t1236;
t798 = -pkin(8) * t844 + t1195;
t1134 = t1075 * t778 - t1079 * t798;
t800 = t1075 * t845 + t1079 * t844;
t1216 = qJ(3) * t848 - t1208 * t800 - t1134;
t850 = (-qJD(5) + t1049) * t974 + t1136;
t852 = t956 + t883;
t812 = t1074 * t852 + t1078 * t850;
t886 = -t970 - t971;
t743 = -pkin(4) * t886 + pkin(8) * t812 + t755;
t810 = t1074 * t850 - t1078 * t852;
t746 = -pkin(8) * t810 - t754;
t1135 = t1075 * t743 - t1079 * t746;
t770 = t1075 * t812 + t1079 * t810;
t1215 = qJ(3) * t886 - t1208 * t770 - t1135;
t1191 = t1075 * t905;
t1054 = t1058 ^ 2;
t969 = -t1054 - t1210;
t903 = t1075 * t969 + t1235;
t933 = (qJD(4) + t1058) * t1029 + t1129;
t1214 = qJ(3) * t933 - t1208 * t903 + t1191;
t897 = t1079 * t905;
t1151 = -t1019 - t1054;
t961 = t1020 + t982;
t1190 = t1075 * t961;
t911 = t1079 * t1151 - t1190;
t1213 = qJ(3) * t1227 - t1208 * t911 + t897;
t1212 = qJ(3) * t905 - t1208 * t1224;
t1207 = pkin(3) * t1224;
t1206 = pkin(3) * t879;
t1205 = pkin(3) * t905;
t1150 = t1071 + t1209;
t1038 = t1150 * qJDD(1);
t1041 = t1065 + t1066;
t1204 = pkin(5) * (t1077 * t1038 + t1081 * t1041);
t1199 = t1082 * pkin(6);
t1193 = t1075 * t754;
t1186 = t1079 * t961;
t1185 = qJD(1) * qJD(2);
t1181 = t1049 * t1074;
t1180 = t1049 * t1078;
t1179 = t1058 * t1029;
t1178 = t1058 * t1075;
t1177 = t1058 * t1079;
t1014 = -t1108 + t1199;
t1175 = t1076 * t1014;
t1174 = t1076 * t1034;
t1165 = t1080 * t1014;
t984 = t1080 * t1036;
t1158 = -t1019 + t1054;
t1157 = pkin(1) * t1041 + pkin(6) * t1038;
t1149 = t1076 * t921;
t1148 = t1080 * t921;
t753 = pkin(4) * t754;
t1147 = -pkin(3) * t736 - t753;
t808 = pkin(4) * t810;
t1146 = -pkin(3) * t770 - t808;
t1143 = t1076 * t982;
t1142 = t1080 * t982;
t1002 = t1076 * t1015 + t1200;
t1003 = t1080 * t1015 - t1068;
t940 = t1076 * t1002 + t1080 * t1003;
t1130 = -t1077 * t1046 - t1081 * t1047;
t1128 = -pkin(3) * t911 + t836;
t1127 = t1077 * t1056;
t1126 = t1081 * t1056;
t1125 = pkin(4) * t844 - t774;
t942 = t1086 + t1182;
t1124 = -pkin(2) * t945 + qJ(3) * t942;
t1040 = -t1077 * t1082 + t1153;
t1123 = -pkin(5) * t1040 - t1077 * g(3);
t1122 = pkin(3) * t1227 - t1191;
t1121 = pkin(3) * t933 + t897;
t1119 = pkin(2) * t1076 - qJ(3) * t1080;
t1115 = -t1074 * t848 + t1078 * t1229;
t953 = -t971 + t1048;
t1114 = t1078 * t953 + t1238;
t952 = t970 - t1048;
t1113 = t1074 * t952 + t1188;
t795 = t1075 * t835 + t1079 * t836;
t939 = t1080 * t1002 - t1076 * t1003;
t985 = t1080 * t1050 + t1173;
t1112 = t1081 * t1046 - t1077 * t1047;
t1111 = -pkin(3) * t903 + t835;
t1101 = pkin(4) * t861 - t775;
t1100 = -t1078 * t1102 + t972 * t1181;
t1099 = t1074 * t883 + t974 * t1180;
t1097 = -pkin(3) * t800 - t1125;
t1096 = (-t1074 * t972 - t1078 * t974) * t1049;
t1095 = pkin(3) * t886 - t1075 * t746 - t1079 * t743;
t1094 = pkin(3) * t848 - t1075 * t798 - t1079 * t778;
t1093 = pkin(3) * t1229 - t1075 * t813 - t1079 * t779;
t1092 = pkin(3) * t959 - t795;
t1090 = pkin(3) * t837 + pkin(8) * t1193 - t1079 * t749;
t1087 = -pkin(3) * t817 - t1101;
t1083 = t1014 - t1103 + t1222;
t1042 = t1065 - t1066;
t1039 = t1081 * t1082 + t1154;
t1030 = t1119 * qJDD(1);
t1022 = t1150 * t1185;
t1012 = -pkin(5) * t1039 + t1081 * g(3);
t1009 = -t1054 + t1210;
t1008 = t1077 * qJDD(2) + t1081 * t1022;
t1007 = t1080 * t1035 - t1071 * t1185;
t1006 = -t1081 * qJDD(2) + t1077 * t1022;
t1005 = t1076 * t1117 - t1209 * t1185;
t987 = t1118 * t1076;
t980 = t1019 - t1210;
t978 = pkin(5) * (t1081 * t1038 - t1077 * t1041);
t977 = t984 - t1174;
t976 = t1080 * t1034 + t1076 * t1036;
t967 = t1081 * t1007 - t1127;
t966 = t1081 * t1005 + t1127;
t965 = t1077 * t1007 + t1126;
t964 = t1077 * t1005 - t1126;
t951 = -t1165 + t1240;
t950 = -t1175 - t1250;
t949 = (-t1027 * t1079 + t1029 * t1075) * t1058;
t948 = (-t1027 * t1075 - t1029 * t1079) * t1058;
t947 = t1077 * t1042 + t1081 * t977;
t946 = -t1081 * t1042 + t1077 * t977;
t944 = t1003 + t1241;
t943 = t1002 - t1251;
t932 = pkin(1) * t1036 + t1165 - t1239;
t931 = -pkin(1) * t1034 - t1175 - t1247;
t928 = qJ(3) * t1041 + t945;
t927 = pkin(2) * t1041 + t942;
t926 = -t1029 * t1178 + t1079 * t968;
t925 = t1029 * t1177 + t1075 * t968;
t924 = -t1027 * t1177 - t1075 * t1098;
t923 = -t1027 * t1178 + t1079 * t1098;
t922 = t1083 + t1228;
t920 = t1080 * t1020 + t1076 * t948;
t919 = t1076 * t1020 - t1080 * t948;
t918 = t971 - t970;
t917 = t1079 * t1009 - t1190;
t916 = -t1075 * t1158 + t1235;
t915 = t1075 * t1009 + t1186;
t914 = t1079 * t1158 + t1237;
t913 = pkin(1) * t1014 + pkin(6) * t940;
t912 = -t1075 * t1151 - t1186;
t910 = -t1199 + (-t1036 + t1117) * pkin(2) + t1084;
t908 = (t1034 + t1118) * qJ(3) + t1083;
t907 = t940 + t1157;
t904 = t1079 * t969 - t1237;
t896 = (t1074 * t974 - t1078 * t972) * t1049;
t894 = -t1221 + t1251;
t893 = t1076 * t925 + t1142;
t892 = -t1076 * t923 - t1142;
t891 = -t1080 * t925 + t1143;
t890 = t1080 * t923 - t1143;
t888 = -t1234 - t1241;
t885 = t1076 * t945 + t1080 * t942;
t884 = t1076 * t942 - t1080 * t945;
t881 = -t1075 * t1223 - t1079 * t934;
t880 = -t1075 * t1227 - t1079 * t933;
t878 = -t1075 * t933 + t1079 * t1227;
t877 = -pkin(2) * t1174 + t1080 * t908 - t1240;
t876 = -qJ(3) * t984 - t1076 * t910 + t1250;
t875 = -t1076 * t927 + t1080 * t928;
t874 = t1076 * t914 - t1080 * t1223;
t873 = t1076 * t915 - t1080 * t934;
t872 = -t1076 * t1223 - t1080 * t914;
t871 = -t1076 * t934 - t1080 * t915;
t870 = t1078 * t952 - t1194;
t869 = -t1074 * t953 + t1236;
t866 = t1247 + t1076 * t908 + (pkin(1) + t1201) * t1034;
t865 = t1239 + t1080 * t910 + (-pkin(1) - t1196) * t1036;
t864 = t1076 * t911 + t1080 * t1227;
t863 = t1076 * t1227 - t1080 * t911;
t858 = t1076 * t903 + t1080 * t933;
t857 = t1076 * t933 - t1080 * t903;
t856 = t1076 * t928 + t1080 * t927 + t1157;
t855 = t1076 * t878 + t1080 * t980;
t854 = t1076 * t980 - t1080 * t878;
t849 = t1102 - t957;
t847 = t1078 * t883 - t974 * t1181;
t846 = t1074 * t1102 + t972 * t1180;
t840 = t1076 * t879 + t1080 * t959;
t839 = t1076 * t959 - t1080 * t879;
t832 = -t1075 * t1096 + t1079 * t896;
t831 = t1075 * t896 + t1079 * t1096;
t830 = t1080 * t1017 + t1076 * t831;
t829 = t1076 * t1017 - t1080 * t831;
t828 = -pkin(1) * t884 - t1124;
t827 = -qJ(3) * t881 + t1206;
t826 = -pkin(6) * t884 - t1119 * t922;
t824 = -t1075 * t1113 + t1079 * t870;
t823 = -t1075 * t1114 + t1079 * t869;
t822 = t1075 * t870 + t1079 * t1113;
t821 = t1075 * t869 + t1079 * t1114;
t818 = -t1075 * t861 + t1079 * t862;
t816 = -t1208 * t912 + t1122;
t815 = pkin(6) * t885 + (pkin(1) - t1120) * t922;
t814 = -t1208 * t904 + t1121;
t811 = -t1074 * t1229 - t1078 * t848;
t805 = -t1075 * t1099 + t1079 * t847;
t804 = -t1075 * t1100 + t1079 * t846;
t803 = t1075 * t847 + t1079 * t1099;
t802 = t1075 * t846 + t1079 * t1100;
t801 = -t1075 * t844 + t1079 * t845;
t799 = -qJ(3) * t912 - t1128;
t796 = -qJ(3) * t904 - t1111;
t793 = t1076 * t803 + t1148;
t792 = t1076 * t802 - t1148;
t791 = -t1080 * t803 + t1149;
t790 = -t1080 * t802 - t1149;
t789 = t1076 * t821 + t1080 * t852;
t788 = t1076 * t822 - t1080 * t849;
t787 = t1076 * t852 - t1080 * t821;
t786 = -t1076 * t849 - t1080 * t822;
t785 = -pkin(1) * t863 - t1213;
t784 = t1076 * t817 + t1080 * t1229;
t783 = t1076 * t1229 - t1080 * t817;
t782 = -pkin(1) * t857 - t1214;
t781 = t1076 * t1224 + t1080 * t905;
t780 = t1076 * t905 - t1080 * t1224;
t777 = t1076 * t800 + t1080 * t848;
t776 = t1076 * t848 - t1080 * t800;
t772 = -t1075 * t810 + t1079 * t812;
t771 = -t1075 * t1115 + t1079 * t811;
t769 = t1075 * t811 + t1079 * t1115;
t768 = -t1208 * t881 + t1092;
t767 = t1076 * t769 + t1080 * t918;
t766 = t1076 * t918 - t1080 * t769;
t765 = -qJ(3) * t795 + t1207;
t764 = t1076 * t770 + t1080 * t886;
t763 = t1076 * t886 - t1080 * t770;
t762 = -pkin(6) * t863 - t1076 * t816 + t1080 * t799;
t761 = -pkin(6) * t857 - t1076 * t814 + t1080 * t796;
t760 = -pkin(1) * t839 - t1218;
t759 = -t1208 * t795 + t1205;
t758 = -pkin(1) * t912 + pkin(6) * t864 + t1076 * t799 + t1080 * t816;
t757 = -pkin(1) * t904 + pkin(6) * t858 + t1076 * t796 + t1080 * t814;
t756 = -pkin(6) * t839 - t1076 * t768 + t1080 * t827;
t751 = -pkin(1) * t881 + pkin(6) * t840 + t1076 * t827 + t1080 * t768;
t750 = -qJ(3) * t818 - t1087;
t748 = -pkin(1) * t780 - t1212;
t747 = -qJ(3) * t801 - t1097;
t744 = -t1208 * t818 + t1093;
t742 = -qJ(3) * t772 - t1146;
t741 = -t1208 * t801 + t1094;
t740 = -pkin(1) * t783 - t1217;
t739 = -pkin(6) * t780 - t1076 * t759 + t1080 * t765;
t738 = -pkin(1) * t776 - t1216;
t737 = t1079 * t755 - t1193;
t735 = t1076 * t736 + t1080 * t837;
t734 = t1076 * t837 - t1080 * t736;
t733 = -pkin(1) * t795 + pkin(6) * t781 + t1076 * t765 + t1080 * t759;
t732 = -pkin(6) * t783 - t1076 * t744 + t1080 * t750;
t731 = -pkin(1) * t818 + pkin(6) * t784 + t1076 * t750 + t1080 * t744;
t730 = -pkin(6) * t776 - t1076 * t741 + t1080 * t747;
t729 = -pkin(1) * t801 + pkin(6) * t777 + t1076 * t747 + t1080 * t741;
t728 = -t1208 * t772 + t1095;
t727 = -pkin(1) * t763 - t1215;
t726 = -qJ(3) * t737 - t1147;
t725 = -pkin(6) * t763 - t1076 * t728 + t1080 * t742;
t724 = -t1208 * t737 + t1090;
t723 = -pkin(1) * t772 + pkin(6) * t764 + t1076 * t742 + t1080 * t728;
t722 = -pkin(1) * t734 - t1219;
t721 = -pkin(6) * t734 - t1076 * t724 + t1080 * t726;
t720 = -pkin(1) * t737 + pkin(6) * t735 + t1076 * t726 + t1080 * t724;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1040, 0, -t1039, 0, t1123, -t1012, -t1112, -pkin(5) * t1112, t967, t947, t1245, t966, -t1231, t1008, -t1077 * t943 + t1081 * t950 + t1248, -t1077 * t944 + t1081 * t951 + t1252, t1081 * t939 - t1204, -pkin(5) * (t1081 * t1014 + t1077 * t940) - (t1077 * pkin(1) - t1081 * pkin(6)) * t939, t1008, -t1245, t1231, t967, t947, t966, -t1077 * t1030 + t1081 * t875 - t1204, -t1077 * t894 + t1081 * t876 - t1248, -t1077 * t888 + t1081 * t877 - t1252, t1081 * t826 - t1077 * t828 - pkin(5) * (t1077 * t885 + t1081 * t922), t1077 * t926 + t1081 * t893, t1077 * t880 + t1081 * t855, t1077 * t916 + t1081 * t874, -t1077 * t924 + t1081 * t892, t1077 * t917 + t1081 * t873, t1077 * t949 + t1081 * t920, t1081 * t761 - t1077 * t782 - pkin(5) * (t1077 * t858 - t1081 * t904), t1081 * t762 - t1077 * t785 - pkin(5) * (t1077 * t864 - t1081 * t912), t1081 * t756 - t1077 * t760 - pkin(5) * (t1077 * t840 - t1081 * t881), t1081 * t739 - t1077 * t748 - pkin(5) * (t1077 * t781 - t1081 * t795), t1077 * t805 + t1081 * t793, t1077 * t771 + t1081 * t767, t1077 * t823 + t1081 * t789, t1077 * t804 + t1081 * t792, t1077 * t824 + t1081 * t788, t1077 * t832 + t1081 * t830, t1081 * t730 - t1077 * t738 - pkin(5) * (t1077 * t777 - t1081 * t801), t1081 * t732 - t1077 * t740 - pkin(5) * (t1077 * t784 - t1081 * t818), t1081 * t725 - t1077 * t727 - pkin(5) * (t1077 * t764 - t1081 * t772), t1081 * t721 - t1077 * t722 - pkin(5) * (t1077 * t735 - t1081 * t737); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1039, 0, t1040, 0, t1012, t1123, t1130, pkin(5) * t1130, t965, t946, t1246, t964, -t1233, t1006, t1077 * t950 + t1081 * t943 - t1249, t1077 * t951 + t1081 * t944 - t1253, t1077 * t939 + t978, pkin(5) * (-t1077 * t1014 + t1081 * t940) - (-t1081 * pkin(1) - t1077 * pkin(6)) * t939, t1006, -t1246, t1233, t965, t946, t964, t1081 * t1030 + t1077 * t875 + t978, t1077 * t876 + t1081 * t894 + t1249, t1077 * t877 + t1081 * t888 + t1253, t1077 * t826 + t1081 * t828 + pkin(5) * (-t1077 * t922 + t1081 * t885), t1077 * t893 - t1081 * t926, t1077 * t855 - t1081 * t880, t1077 * t874 - t1081 * t916, t1077 * t892 + t1081 * t924, t1077 * t873 - t1081 * t917, t1077 * t920 - t1081 * t949, t1077 * t761 + t1081 * t782 + pkin(5) * (t1077 * t904 + t1081 * t858), t1077 * t762 + t1081 * t785 + pkin(5) * (t1077 * t912 + t1081 * t864), t1077 * t756 + t1081 * t760 + pkin(5) * (t1077 * t881 + t1081 * t840), t1077 * t739 + t1081 * t748 + pkin(5) * (t1077 * t795 + t1081 * t781), t1077 * t793 - t1081 * t805, t1077 * t767 - t1081 * t771, t1077 * t789 - t1081 * t823, t1077 * t792 - t1081 * t804, t1077 * t788 - t1081 * t824, t1077 * t830 - t1081 * t832, t1077 * t730 + t1081 * t738 + pkin(5) * (t1077 * t801 + t1081 * t777), t1077 * t732 + t1081 * t740 + pkin(5) * (t1077 * t818 + t1081 * t784), t1077 * t725 + t1081 * t727 + pkin(5) * (t1077 * t772 + t1081 * t764), t1077 * t721 + t1081 * t722 + pkin(5) * (t1077 * t737 + t1081 * t735); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1046, t1047, 0, 0, t987, t976, t985, t984, t988, 0, t932, t931, t907, t913, 0, -t985, -t988, t987, t976, t984, t856, t865, t866, t815, t891, t854, t872, t890, t871, t919, t757, t758, t751, t733, t791, t766, t787, t790, t786, t829, t729, t731, t723, t720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1082, 0, 0, -g(3), -t1046, 0, t1007, t977, t994, t1005, -t997, t1022, t950, t951, t939, pkin(6) * t939, t1022, -t994, t997, t1007, t977, t1005, t875, t876, t877, t826, t893, t855, t874, t892, t873, t920, t761, t762, t756, t739, t793, t767, t789, t792, t788, t830, t730, t732, t725, t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1082, 0, qJDD(1), 0, g(3), 0, -t1047, 0, t1056, -t1042, -t1155, -t1056, -t1064, -qJDD(2), t943, t944, 0, pkin(1) * t939, -qJDD(2), t1155, t1064, t1056, -t1042, -t1056, t1030, t894, t888, t828, -t926, -t880, -t916, t924, -t917, -t949, t782, t785, t760, t748, -t805, -t771, -t823, -t804, -t824, -t832, t738, t740, t727, t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1046, t1047, 0, 0, t987, t976, t985, t984, t988, 0, t932, t931, t907, t913, 0, -t985, -t988, t987, t976, t984, t856, t865, t866, t815, t891, t854, t872, t890, t871, t919, t757, t758, t751, t733, t791, t766, t787, t790, t786, t829, t729, t731, t723, t720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1035, t1036, t1043, -t1062, t1052, t1062, 0, -t1014, t1002, 0, t1062, -t1043, -t1052, t1035, t1036, -t1062, t928, -qJ(3) * t1036, t908, qJ(3) * t922, t982, t980, -t1223, -t982, -t934, t1020, t796, t799, t827, t765, t921, t918, t852, -t921, -t849, t1017, t747, t750, t742, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1144, t1034, t1050, -t1117, t1044, -t1144, t1014, 0, t1003, 0, -t1144, -t1050, -t1044, t1144, t1034, -t1117, t927, t910, pkin(2) * t1034, pkin(2) * t922, -t925, -t878, -t914, t923, -t915, -t948, t814, t816, t768, t759, -t803, -t769, -t821, -t802, -t822, -t831, t741, t744, t728, t724; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1056, t1042, t1155, t1056, t1064, qJDD(2), -t1002, -t1003, 0, 0, qJDD(2), -t1155, -t1064, -t1056, t1042, t1056, -t1030, t1221, t1234, t1124, t926, t880, t916, -t924, t917, t949, t1214, t1213, t1218, t1212, t805, t771, t823, t804, t824, t832, t1216, t1217, t1215, t1219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1155, -t1064, -t1056, t1042, t1056, 0, t945, t942, 0, t926, t880, t916, -t924, t917, t949, -pkin(7) * t903 + t1191, -pkin(7) * t911 + t897, -pkin(7) * t879 - t1224, -pkin(7) * t1224, t805, t771, t823, t804, t824, t832, -pkin(7) * t800 - t1134, -pkin(7) * t817 - t1133, -pkin(7) * t770 - t1135, -pkin(7) * t736 - t1109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1062, t1043, t1052, -t1035, -t1036, t1062, -t945, 0, -t922, 0, -t982, -t980, t1223, t982, t934, -t1020, t1111, t1128, -t1206, -t1207, -t921, -t918, -t852, t921, t849, -t1017, t1097, t1087, t1146, t1147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1144, t1050, t1044, -t1144, -t1034, t1117, -t942, t922, 0, 0, t925, t878, t914, -t923, t915, t948, pkin(7) * t904 - t1121, pkin(7) * t912 - t1122, pkin(7) * t881 - t1092, pkin(7) * t795 - t1205, t803, t769, t821, t802, t822, t831, pkin(7) * t801 - t1094, pkin(7) * t818 - t1093, pkin(7) * t772 - t1095, pkin(7) * t737 - t1090; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t968, -t933, t1225, t1011, t1009, -t1011, 0, t905, t835, 0, t847, t811, t869, t846, t870, t896, t798, t813, t746, -pkin(8) * t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1179, t1227, t1158, -t1098, t961, -t1179, -t905, 0, t836, 0, t1099, t1115, t1114, t1100, t1113, t1096, t778, t779, t743, t749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t982, t980, -t1223, -t982, -t934, t1020, -t835, -t836, 0, 0, t921, t918, t852, -t921, -t849, t1017, t1125, t1101, t808, t753; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t883, -t848, t1226, t956, t952, -t956, 0, t837, t774, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t957, t1229, t953, -t1102, t901, -t957, -t837, 0, t775, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t921, t918, t852, -t921, -t849, t1017, -t774, -t775, 0, 0;];
m_new_reg = t1;
