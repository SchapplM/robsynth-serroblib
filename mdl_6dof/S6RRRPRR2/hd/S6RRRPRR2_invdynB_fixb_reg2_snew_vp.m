% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:59:42
% EndTime: 2019-05-07 10:00:28
% DurationCPUTime: 47.53s
% Computational Cost: add. (440386->958), mult. (993289->1478), div. (0->0), fcn. (757904->12), ass. (0->656)
t1155 = cos(qJ(2));
t1138 = t1155 * qJDD(1);
t1150 = sin(qJ(2));
t1228 = qJD(1) * t1150;
t1192 = qJD(2) * t1228;
t1113 = t1138 - t1192;
t1144 = t1155 ^ 2;
t1158 = qJD(1) ^ 2;
t1151 = sin(qJ(1));
t1156 = cos(qJ(1));
t1122 = g(1) * t1151 - t1156 * g(2);
t1171 = qJDD(1) * pkin(1) + t1122;
t1172 = qJD(2) * pkin(2) - pkin(8) * t1228;
t1044 = pkin(2) * t1113 - t1172 * t1228 + t1171 + (pkin(8) * t1144 + pkin(7)) * t1158;
t1145 = sin(pkin(11));
t1149 = sin(qJ(3));
t1154 = cos(qJ(3));
t1227 = qJD(1) * t1155;
t1103 = -t1149 * t1228 + t1154 * t1227;
t1104 = (t1149 * t1155 + t1150 * t1154) * qJD(1);
t1146 = cos(pkin(11));
t1067 = -t1146 * t1103 + t1104 * t1145;
t1069 = t1145 * t1103 + t1146 * t1104;
t1007 = t1069 * t1067;
t1197 = qJDD(2) + qJDD(3);
t1252 = -t1007 + t1197;
t1263 = t1145 * t1252;
t1262 = t1146 * t1252;
t1147 = sin(qJ(6));
t1185 = qJD(2) * t1227;
t1200 = qJDD(1) * t1150;
t1112 = t1185 + t1200;
t1179 = t1112 * t1149 - t1154 * t1113;
t1041 = -qJD(3) * t1104 - t1179;
t1042 = qJD(3) * t1103 + t1154 * t1112 + t1149 * t1113;
t1180 = -t1146 * t1041 + t1042 * t1145;
t987 = qJDD(5) + t1180;
t1170 = qJDD(6) + t987;
t1142 = qJD(2) + qJD(3);
t1148 = sin(qJ(5));
t1153 = cos(qJ(5));
t1036 = t1069 * t1148 - t1153 * t1142;
t1038 = t1069 * t1153 + t1142 * t1148;
t1152 = cos(qJ(6));
t984 = t1152 * t1036 + t1038 * t1147;
t986 = -t1036 * t1147 + t1038 * t1152;
t920 = t986 * t984;
t1253 = t1170 - t920;
t1261 = t1147 * t1253;
t991 = t1038 * t1036;
t1254 = t987 - t991;
t1260 = t1148 * t1254;
t1075 = t1103 * t1104;
t1250 = t1075 + t1197;
t1259 = t1149 * t1250;
t1258 = t1152 * t1253;
t1257 = t1153 * t1254;
t1256 = t1154 * t1250;
t1101 = t1103 ^ 2;
t1173 = pkin(3) * t1142 - qJ(4) * t1104;
t942 = pkin(3) * t1041 + qJ(4) * t1101 - t1104 * t1173 - qJDD(4) + t1044;
t1218 = t1069 * t1142;
t952 = t1180 + t1218;
t989 = t1041 * t1145 + t1042 * t1146;
t1181 = t1148 * t989 - t1153 * t1197;
t936 = -qJD(5) * t1038 - t1181;
t937 = -t1036 * qJD(5) + t1148 * t1197 + t1153 * t989;
t847 = -qJD(6) * t984 + t1147 * t936 + t1152 * t937;
t1063 = qJD(5) + t1067;
t1059 = qJD(6) + t1063;
t950 = t1059 * t984;
t1255 = t847 - t950;
t998 = t1063 * t1036;
t907 = -t998 - t937;
t906 = -t998 + t937;
t1057 = t1142 * t1067;
t1182 = -t989 + t1057;
t1096 = t1142 * t1103;
t1019 = -t1042 + t1096;
t1251 = t1042 + t1096;
t1140 = t1144 * t1158;
t1157 = qJD(2) ^ 2;
t1128 = -t1140 - t1157;
t1183 = t1147 * t937 - t1152 * t936;
t808 = (qJD(6) - t1059) * t986 + t1183;
t903 = (qJD(5) - t1063) * t1038 + t1181;
t1015 = (qJD(3) - t1142) * t1104 + t1179;
t1123 = g(1) * t1156 + g(2) * t1151;
t1167 = pkin(1) * t1158 - qJDD(1) * pkin(7) + t1123;
t1166 = t1150 * t1167;
t1163 = qJDD(2) * pkin(2) - t1112 * pkin(8) + t1166;
t1246 = t1150 * g(3);
t1164 = -pkin(2) * t1140 + t1113 * pkin(8) - qJD(2) * t1172 - t1246;
t1205 = t1150 * t1158;
t1229 = qJD(1) * qJD(2);
t1169 = pkin(2) * t1205 + pkin(8) * t1229 - g(3);
t978 = t1149 * t1164 - t1154 * t1163 + (-t1149 * t1167 - t1154 * t1169) * t1155;
t982 = t984 ^ 2;
t983 = t986 ^ 2;
t1249 = t1036 ^ 2;
t1035 = t1038 ^ 2;
t1058 = t1059 ^ 2;
t1062 = t1063 ^ 2;
t1064 = t1067 ^ 2;
t1065 = t1069 ^ 2;
t1102 = t1104 ^ 2;
t1124 = t1142 ^ 2;
t1248 = 2 * qJD(4);
t1247 = pkin(4) * t1145;
t1003 = pkin(4) * t1067 - pkin(9) * t1069;
t1159 = pkin(3) * t1250 + qJ(4) * t1019 - t978;
t1165 = t1155 * t1167;
t979 = t1154 * (-t1165 + t1164) + t1149 * (t1155 * t1169 + t1163);
t921 = -t1101 * pkin(3) + t1041 * qJ(4) - t1142 * t1173 + t979;
t837 = -0.2e1 * qJD(4) * t1067 + t1145 * t1159 + t1146 * t921;
t792 = -pkin(4) * t1124 + pkin(9) * t1197 - t1067 * t1003 + t837;
t845 = pkin(4) * t952 + t1182 * pkin(9) - t942;
t754 = t1148 * t792 - t1153 * t845;
t704 = pkin(5) * t1254 + pkin(10) * t907 - t754;
t755 = t1148 * t845 + t1153 * t792;
t993 = pkin(5) * t1063 - pkin(10) * t1038;
t717 = -pkin(5) * t1249 + pkin(10) * t936 - t1063 * t993 + t755;
t654 = t1147 * t704 + t1152 * t717;
t1245 = t1145 * t942;
t1244 = t1146 * t942;
t1184 = t1145 * t921 - t1146 * t1159;
t791 = -t1197 * pkin(4) - t1124 * pkin(9) + (t1248 + t1003) * t1069 + t1184;
t760 = -t936 * pkin(5) - pkin(10) * t1249 + t1038 * t993 + t791;
t1243 = t1147 * t760;
t870 = t1170 + t920;
t1242 = t1147 * t870;
t653 = t1147 * t717 - t1152 * t704;
t610 = t1147 * t654 - t1152 * t653;
t1241 = t1148 * t610;
t1240 = t1148 * t791;
t923 = t987 + t991;
t1239 = t1148 * t923;
t836 = t1069 * t1248 + t1184;
t761 = t1145 * t837 - t1146 * t836;
t1238 = t1149 * t761;
t912 = t1149 * t979 - t1154 * t978;
t1237 = t1150 * t912;
t1236 = t1152 * t760;
t1235 = t1152 * t870;
t1234 = t1153 * t610;
t1233 = t1153 * t791;
t1232 = t1153 * t923;
t1231 = t1154 * t761;
t1230 = t1155 * t912;
t1001 = t1007 + t1197;
t1226 = t1001 * t1145;
t1225 = t1001 * t1146;
t1224 = t1044 * t1149;
t1223 = t1044 * t1154;
t1222 = t1059 * t1147;
t1221 = t1059 * t1152;
t1220 = t1063 * t1148;
t1219 = t1063 * t1153;
t1072 = -t1075 + t1197;
t1217 = t1072 * t1149;
t1216 = t1072 * t1154;
t1105 = pkin(7) * t1158 + t1171;
t1215 = t1105 * t1150;
t1214 = t1105 * t1155;
t1130 = t1155 * t1205;
t1120 = qJDD(2) + t1130;
t1213 = t1120 * t1150;
t1121 = qJDD(2) - t1130;
t1212 = t1121 * t1150;
t1211 = t1121 * t1155;
t1210 = t1142 * t1145;
t1209 = t1142 * t1146;
t1208 = t1142 * t1149;
t1207 = t1142 * t1154;
t1143 = t1150 ^ 2;
t1206 = t1143 * t1158;
t1201 = t1143 + t1144;
t1199 = qJDD(1) * t1151;
t1198 = qJDD(1) * t1156;
t1196 = t1145 * t920;
t1195 = t1146 * t920;
t1194 = -pkin(4) * t1146 - pkin(3);
t1191 = t1145 * t991;
t1190 = t1146 * t991;
t1189 = t1151 * t1007;
t1188 = t1156 * t1007;
t1187 = t1151 * t1075;
t1186 = t1156 * t1075;
t762 = t1145 * t836 + t1146 * t837;
t611 = t1147 * t653 + t1152 * t654;
t913 = t1149 * t978 + t1154 * t979;
t1089 = t1155 * g(3) - t1166;
t1090 = -t1165 - t1246;
t1030 = t1089 * t1150 + t1155 * t1090;
t1081 = -t1122 * t1151 - t1156 * t1123;
t1178 = t1151 * t1197;
t1177 = t1151 * t1130;
t1176 = t1156 * t1130;
t1117 = -t1151 * t1158 + t1198;
t1174 = -pkin(6) * t1117 - g(3) * t1151;
t686 = t1148 * t755 - t1153 * t754;
t687 = t1148 * t754 + t1153 * t755;
t1029 = t1089 * t1155 - t1090 * t1150;
t1080 = t1122 * t1156 - t1123 * t1151;
t953 = t1180 - t1218;
t1133 = t1156 * t1197;
t1127 = t1140 - t1157;
t1126 = -t1157 - t1206;
t1125 = t1157 - t1206;
t1119 = t1140 - t1206;
t1118 = t1140 + t1206;
t1116 = t1156 * t1158 + t1199;
t1115 = t1201 * qJDD(1);
t1114 = t1138 - 0.2e1 * t1192;
t1111 = 0.2e1 * t1185 + t1200;
t1109 = t1155 * t1120;
t1108 = t1201 * t1229;
t1100 = -pkin(6) * t1116 + g(3) * t1156;
t1094 = -t1102 + t1124;
t1093 = t1101 - t1124;
t1092 = t1112 * t1155 - t1143 * t1229;
t1091 = -t1113 * t1150 - t1144 * t1229;
t1088 = -t1102 - t1124;
t1087 = -t1126 * t1150 - t1211;
t1086 = -t1125 * t1150 + t1109;
t1085 = t1128 * t1155 - t1213;
t1084 = t1127 * t1155 - t1212;
t1083 = t1126 * t1155 - t1212;
t1082 = t1128 * t1150 + t1109;
t1078 = t1115 * t1156 - t1118 * t1151;
t1077 = t1115 * t1151 + t1118 * t1156;
t1076 = -t1111 * t1150 + t1114 * t1155;
t1074 = -t1102 + t1101;
t1070 = -t1124 - t1101;
t1055 = t1087 * t1156 + t1111 * t1151;
t1054 = t1085 * t1156 - t1114 * t1151;
t1053 = t1087 * t1151 - t1111 * t1156;
t1052 = t1085 * t1151 + t1114 * t1156;
t1051 = -t1065 + t1124;
t1050 = t1064 - t1124;
t1049 = (t1103 * t1154 + t1104 * t1149) * t1142;
t1048 = (t1103 * t1149 - t1104 * t1154) * t1142;
t1047 = -pkin(7) * t1083 - t1214;
t1046 = -pkin(7) * t1082 - t1215;
t1045 = -t1065 - t1124;
t1043 = -t1101 - t1102;
t1040 = -pkin(1) * t1083 + t1090;
t1039 = -pkin(1) * t1082 + t1089;
t1025 = t1093 * t1154 - t1217;
t1024 = -t1094 * t1149 + t1256;
t1023 = t1093 * t1149 + t1216;
t1022 = t1094 * t1154 + t1259;
t1021 = -t1088 * t1149 - t1216;
t1020 = t1088 * t1154 - t1217;
t1014 = (qJD(3) + t1142) * t1104 + t1179;
t1013 = t1042 * t1154 - t1104 * t1208;
t1012 = t1042 * t1149 + t1104 * t1207;
t1011 = -t1041 * t1149 - t1103 * t1207;
t1010 = t1041 * t1154 - t1103 * t1208;
t1009 = t1030 * t1156 - t1105 * t1151;
t1008 = t1030 * t1151 + t1105 * t1156;
t1006 = t1070 * t1154 - t1259;
t1005 = t1070 * t1149 + t1256;
t1004 = -t1065 + t1064;
t999 = -t1124 - t1064;
t997 = -t1035 + t1062;
t996 = -t1062 + t1249;
t995 = (-t1067 * t1146 + t1069 * t1145) * t1142;
t994 = (-t1067 * t1145 - t1069 * t1146) * t1142;
t992 = -t1048 * t1150 + t1049 * t1155;
t990 = -t1035 + t1249;
t977 = -t1064 - t1065;
t975 = -pkin(8) * t1020 - t1223;
t974 = -t1035 - t1062;
t973 = -pkin(8) * t1005 - t1224;
t972 = -t1023 * t1150 + t1025 * t1155;
t971 = -t1022 * t1150 + t1024 * t1155;
t970 = t1050 * t1146 - t1226;
t969 = -t1051 * t1145 + t1262;
t968 = t1050 * t1145 + t1225;
t967 = t1051 * t1146 + t1263;
t966 = -t1045 * t1145 - t1225;
t965 = t1045 * t1146 - t1226;
t964 = -t1062 - t1249;
t963 = -t1020 * t1150 + t1021 * t1155;
t962 = t1020 * t1155 + t1021 * t1150;
t961 = -t1015 * t1154 - t1019 * t1149;
t960 = -t1014 * t1154 - t1149 * t1251;
t959 = -t1015 * t1149 + t1019 * t1154;
t958 = -t1014 * t1149 + t1154 * t1251;
t957 = -t1057 - t989;
t949 = -t983 + t1058;
t948 = t982 - t1058;
t947 = t1035 + t1249;
t946 = -t1069 * t1210 + t1146 * t989;
t945 = t1069 * t1209 + t1145 * t989;
t944 = t1067 * t1209 + t1145 * t1180;
t943 = t1067 * t1210 - t1146 * t1180;
t941 = -t1012 * t1150 + t1013 * t1155;
t940 = -t1010 * t1150 + t1011 * t1155;
t939 = -t1005 * t1150 + t1006 * t1155;
t938 = t1005 * t1155 + t1006 * t1150;
t933 = t1146 * t999 - t1263;
t932 = t1145 * t999 + t1262;
t931 = (-t1036 * t1153 + t1038 * t1148) * t1063;
t930 = (t1036 * t1148 + t1038 * t1153) * t1063;
t929 = -pkin(2) * t1251 + pkin(8) * t1021 - t1224;
t928 = -t983 - t1058;
t927 = -t1149 * t994 + t1154 * t995;
t926 = t1149 * t995 + t1154 * t994;
t925 = -pkin(2) * t1014 + pkin(8) * t1006 + t1223;
t918 = -t983 + t982;
t917 = t1151 * t1251 + t1156 * t963;
t916 = t1151 * t963 - t1156 * t1251;
t911 = t1014 * t1151 + t1156 * t939;
t910 = -t1014 * t1156 + t1151 * t939;
t909 = -t1058 - t982;
t904 = (-qJD(5) - t1063) * t1038 - t1181;
t902 = -t1038 * t1220 + t1153 * t937;
t901 = -t1038 * t1219 - t1148 * t937;
t900 = t1036 * t1219 - t1148 * t936;
t899 = -t1036 * t1220 - t1153 * t936;
t898 = -t1149 * t968 + t1154 * t970;
t897 = -t1149 * t967 + t1154 * t969;
t896 = t1149 * t970 + t1154 * t968;
t895 = t1149 * t969 + t1154 * t967;
t894 = -t1149 * t965 + t1154 * t966;
t893 = t1149 * t966 + t1154 * t965;
t892 = pkin(2) * t1044 + pkin(8) * t913;
t891 = -t1150 * t959 + t1155 * t961;
t890 = -t1150 * t958 + t1155 * t960;
t889 = t1150 * t961 + t1155 * t959;
t888 = -qJ(4) * t965 - t1244;
t887 = t1145 * t987 + t1146 * t931;
t886 = t1145 * t931 - t1146 * t987;
t885 = -t1145 * t957 - t1146 * t953;
t884 = t1145 * t1182 - t1146 * t952;
t883 = -t1145 * t953 + t1146 * t957;
t882 = -t1145 * t952 - t1146 * t1182;
t881 = (t1147 * t986 - t1152 * t984) * t1059;
t880 = (-t1147 * t984 - t1152 * t986) * t1059;
t879 = t1153 * t996 - t1239;
t878 = -t1148 * t997 + t1257;
t877 = -t1148 * t996 - t1232;
t876 = -t1153 * t997 - t1260;
t875 = -t1149 * t945 + t1154 * t946;
t874 = -t1149 * t943 + t1154 * t944;
t873 = t1149 * t946 + t1154 * t945;
t872 = t1149 * t944 + t1154 * t943;
t868 = -t982 - t983;
t867 = -pkin(1) * t962 - pkin(2) * t1020 + t979;
t866 = -qJ(4) * t932 - t1245;
t865 = -t1148 * t974 - t1232;
t864 = t1153 * t974 - t1239;
t863 = t1043 * t1151 + t1156 * t891;
t862 = -t1043 * t1156 + t1151 * t891;
t861 = -t1149 * t932 + t1154 * t933;
t860 = t1149 * t933 + t1154 * t932;
t859 = t1153 * t964 - t1260;
t858 = t1148 * t964 + t1257;
t857 = -pkin(1) * t938 + t1149 * t1090 + t1154 * t1089 + (t1149 * (t1113 + t1192) - t1154 * (-t1112 + t1185)) * pkin(8) + (-t1154 * t1120 + t1128 * t1149 - t1005) * pkin(2);
t856 = t1146 * t902 + t1191;
t855 = t1146 * t900 - t1191;
t854 = t1145 * t902 - t1190;
t853 = t1145 * t900 + t1190;
t852 = -pkin(8) * t959 - t912;
t851 = -t1150 * t926 + t1155 * t927;
t850 = -pkin(2) * t1043 + pkin(8) * t961 + t913;
t849 = -pkin(1) * t889 - pkin(2) * t959;
t848 = pkin(3) * t1182 + qJ(4) * t966 - t1245;
t846 = -qJD(6) * t986 - t1183;
t842 = -pkin(7) * t962 - t1150 * t929 + t1155 * t975;
t841 = t1152 * t948 - t1242;
t840 = -t1147 * t949 + t1258;
t839 = t1147 * t948 + t1235;
t838 = t1152 * t949 + t1261;
t835 = -pkin(3) * t952 + qJ(4) * t933 + t1244;
t834 = t1155 * t913 - t1237;
t833 = t1150 * t913 + t1230;
t831 = -pkin(7) * t938 - t1150 * t925 + t1155 * t973;
t830 = -t1147 * t928 - t1235;
t829 = t1152 * t928 - t1242;
t828 = -t1148 * t907 - t1153 * t903;
t827 = -t1148 * t906 + t1153 * t904;
t826 = -t1148 * t903 + t1153 * t907;
t825 = -t1148 * t904 - t1153 * t906;
t824 = -t1044 * t1151 + t1156 * t834;
t823 = t1044 * t1156 + t1151 * t834;
t822 = -t1150 * t896 + t1155 * t898;
t821 = -t1150 * t895 + t1155 * t897;
t820 = -t1145 * t903 + t1146 * t879;
t819 = -t1145 * t907 + t1146 * t878;
t818 = t1145 * t879 + t1146 * t903;
t817 = t1145 * t878 + t1146 * t907;
t816 = -t1150 * t893 + t1155 * t894;
t815 = t1150 * t894 + t1155 * t893;
t814 = t1152 * t909 - t1261;
t813 = t1147 * t909 + t1258;
t812 = -t847 - t950;
t807 = (qJD(6) + t1059) * t986 + t1183;
t806 = t1152 * t847 - t1222 * t986;
t805 = t1147 * t847 + t1221 * t986;
t804 = -t1147 * t846 + t1221 * t984;
t803 = t1152 * t846 + t1222 * t984;
t802 = -t1149 * t886 + t1154 * t887;
t801 = t1149 * t887 + t1154 * t886;
t800 = t1145 * t906 + t1146 * t865;
t799 = t1145 * t865 - t1146 * t906;
t798 = -t1149 * t883 + t1154 * t885;
t797 = -t1149 * t882 + t1154 * t884;
t796 = t1149 * t885 + t1154 * t883;
t795 = t1149 * t884 + t1154 * t882;
t794 = -t1148 * t880 + t1153 * t881;
t793 = -t1148 * t881 - t1153 * t880;
t789 = -t1145 * t904 + t1146 * t859;
t788 = t1145 * t859 + t1146 * t904;
t787 = -t1145 * t990 + t1146 * t827;
t786 = t1145 * t827 + t1146 * t990;
t785 = -t1150 * t873 + t1155 * t875;
t784 = -t1150 * t872 + t1155 * t874;
t783 = -t1145 * t947 + t1146 * t828;
t782 = t1145 * t828 + t1146 * t947;
t781 = -t1150 * t860 + t1155 * t861;
t780 = t1150 * t861 + t1155 * t860;
t779 = -t1151 * t1182 + t1156 * t816;
t778 = t1151 * t816 + t1156 * t1182;
t777 = -pkin(1) * t833 - pkin(2) * t912;
t776 = t1145 * t1170 + t1146 * t794;
t775 = t1145 * t794 - t1146 * t1170;
t774 = -t1149 * t854 + t1154 * t856;
t773 = -t1149 * t853 + t1154 * t855;
t772 = t1149 * t856 + t1154 * t854;
t771 = t1149 * t855 + t1154 * t853;
t770 = t1151 * t952 + t1156 * t781;
t769 = t1151 * t781 - t1156 * t952;
t768 = -t1148 * t839 + t1153 * t841;
t767 = -t1148 * t838 + t1153 * t840;
t766 = -t1148 * t841 - t1153 * t839;
t765 = -t1148 * t840 - t1153 * t838;
t764 = -pkin(9) * t864 + t1233;
t763 = -pkin(9) * t858 + t1240;
t759 = -pkin(7) * t833 - pkin(8) * t1230 - t1150 * t892;
t758 = -pkin(8) * t893 - t1149 * t848 + t1154 * t888;
t757 = -t1148 * t829 + t1153 * t830;
t756 = t1148 * t830 + t1153 * t829;
t753 = -pkin(7) * t889 - t1150 * t850 + t1155 * t852;
t752 = pkin(3) * t942 + qJ(4) * t762;
t751 = -t1149 * t818 + t1154 * t820;
t750 = -t1149 * t817 + t1154 * t819;
t749 = t1149 * t820 + t1154 * t818;
t748 = t1149 * t819 + t1154 * t817;
t747 = -t1148 * t813 + t1153 * t814;
t746 = t1148 * t814 + t1153 * t813;
t745 = -pkin(8) * t860 - t1149 * t835 + t1154 * t866;
t744 = pkin(2) * t1182 + pkin(8) * t894 + t1149 * t888 + t1154 * t848;
t743 = -t1147 * t812 - t1152 * t808;
t742 = -t1147 * t1255 - t1152 * t807;
t741 = -t1147 * t808 + t1152 * t812;
t740 = -t1147 * t807 + t1152 * t1255;
t739 = -t1148 * t805 + t1153 * t806;
t738 = -t1148 * t803 + t1153 * t804;
t737 = -t1148 * t806 - t1153 * t805;
t736 = -t1148 * t804 - t1153 * t803;
t735 = -t1150 * t801 + t1155 * t802;
t734 = -t1149 * t799 + t1154 * t800;
t733 = t1149 * t800 + t1154 * t799;
t732 = -t1150 * t796 + t1155 * t798;
t731 = -t1150 * t795 + t1155 * t797;
t730 = t1150 * t798 + t1155 * t796;
t729 = -t1149 * t788 + t1154 * t789;
t728 = t1149 * t789 + t1154 * t788;
t727 = -t1149 * t786 + t1154 * t787;
t726 = t1149 * t787 + t1154 * t786;
t725 = -qJ(4) * t883 - t761;
t724 = -pkin(2) * t952 + pkin(8) * t861 + t1149 * t866 + t1154 * t835;
t723 = t1151 * t977 + t1156 * t732;
t722 = t1151 * t732 - t1156 * t977;
t721 = -t1149 * t782 + t1154 * t783;
t720 = t1149 * t783 + t1154 * t782;
t719 = -pkin(3) * t977 + qJ(4) * t885 + t762;
t718 = -pkin(4) * t864 + t755;
t716 = -pkin(4) * t858 + t754;
t714 = t1146 * t739 + t1196;
t713 = t1146 * t738 - t1196;
t712 = t1145 * t739 - t1195;
t711 = t1145 * t738 + t1195;
t710 = -t1149 * t775 + t1154 * t776;
t709 = t1149 * t776 + t1154 * t775;
t708 = -pkin(1) * t815 - pkin(2) * t893 - pkin(3) * t965 + t837;
t707 = -pkin(10) * t829 + t1236;
t706 = -t1150 * t772 + t1155 * t774;
t705 = -t1150 * t771 + t1155 * t773;
t701 = -pkin(10) * t813 + t1243;
t700 = -t1145 * t808 + t1146 * t768;
t699 = -t1145 * t812 + t1146 * t767;
t698 = t1145 * t768 + t1146 * t808;
t697 = t1145 * t767 + t1146 * t812;
t696 = t1145 * t1255 + t1146 * t757;
t695 = t1145 * t757 - t1146 * t1255;
t694 = -pkin(1) * t780 - pkin(2) * t860 - pkin(3) * t932 + t836;
t693 = t1145 * t807 + t1146 * t747;
t692 = t1145 * t747 - t1146 * t807;
t691 = t1154 * t762 - t1238;
t690 = t1149 * t762 + t1231;
t689 = -pkin(5) * t1255 + pkin(10) * t830 + t1243;
t688 = -pkin(1) * t730 - pkin(2) * t796 - pkin(3) * t883;
t685 = -pkin(5) * t807 + pkin(10) * t814 - t1236;
t684 = -t1150 * t749 + t1155 * t751;
t683 = -t1150 * t748 + t1155 * t750;
t682 = -t1148 * t741 + t1153 * t743;
t681 = -t1148 * t740 + t1153 * t742;
t680 = t1148 * t743 + t1153 * t741;
t679 = -t1148 * t742 - t1153 * t740;
t678 = -t1150 * t733 + t1155 * t734;
t677 = t1150 * t734 + t1155 * t733;
t676 = -t1150 * t728 + t1155 * t729;
t675 = t1150 * t729 + t1155 * t728;
t674 = -t1150 * t726 + t1155 * t727;
t673 = -t1145 * t918 + t1146 * t681;
t672 = t1145 * t681 + t1146 * t918;
t671 = -pkin(9) * t826 - t686;
t670 = -t1150 * t720 + t1155 * t721;
t669 = t1145 * t868 + t1146 * t682;
t668 = t1150 * t721 + t1155 * t720;
t667 = t1145 * t682 - t1146 * t868;
t666 = t1145 * t791 + t1146 * t687;
t665 = t1145 * t687 - t1146 * t791;
t664 = -pkin(7) * t815 - t1150 * t744 + t1155 * t758;
t663 = t1151 * t864 + t1156 * t678;
t662 = t1151 * t678 - t1156 * t864;
t661 = t1151 * t858 + t1156 * t676;
t660 = t1151 * t676 - t1156 * t858;
t659 = -t1149 * t712 + t1154 * t714;
t658 = -t1149 * t711 + t1154 * t713;
t657 = t1149 * t714 + t1154 * t712;
t656 = t1149 * t713 + t1154 * t711;
t655 = -t1150 * t709 + t1155 * t710;
t651 = -qJ(4) * t799 - t1145 * t718 + t1146 * t764;
t650 = -qJ(4) * t788 - t1145 * t716 + t1146 * t763;
t649 = -t1149 * t698 + t1154 * t700;
t648 = -t1149 * t697 + t1154 * t699;
t647 = t1149 * t700 + t1154 * t698;
t646 = t1149 * t699 + t1154 * t697;
t645 = -pkin(7) * t780 - t1150 * t724 + t1155 * t745;
t644 = t1151 * t826 + t1156 * t670;
t643 = t1151 * t670 - t1156 * t826;
t642 = -pkin(8) * t796 - t1149 * t719 + t1154 * t725;
t641 = -t1149 * t695 + t1154 * t696;
t640 = t1149 * t696 + t1154 * t695;
t639 = -pkin(3) * t864 + qJ(4) * t800 + t1145 * t764 + t1146 * t718;
t638 = -pkin(3) * t858 + qJ(4) * t789 + t1145 * t763 + t1146 * t716;
t637 = -pkin(2) * t977 + pkin(8) * t798 + t1149 * t725 + t1154 * t719;
t636 = -t1149 * t692 + t1154 * t693;
t635 = t1149 * t693 + t1154 * t692;
t634 = -pkin(4) * t680 - pkin(5) * t741;
t633 = -qJ(4) * t782 + t1146 * t671 + t1247 * t826;
t632 = -t1150 * t690 + t1155 * t691;
t631 = t1150 * t691 + t1155 * t690;
t630 = -pkin(8) * t690 - qJ(4) * t1231 - t1149 * t752;
t629 = -t1151 * t942 + t1156 * t632;
t628 = t1151 * t632 + t1156 * t942;
t627 = pkin(2) * t942 + pkin(8) * t691 - qJ(4) * t1238 + t1154 * t752;
t626 = qJ(4) * t783 + t1145 * t671 + t1194 * t826;
t625 = -pkin(4) * t756 - pkin(5) * t829 + t654;
t624 = -pkin(9) * t756 - t1148 * t689 + t1153 * t707;
t623 = -pkin(4) * t746 - pkin(5) * t813 + t653;
t622 = -pkin(9) * t746 - t1148 * t685 + t1153 * t701;
t621 = -t1149 * t672 + t1154 * t673;
t620 = t1149 * t673 + t1154 * t672;
t619 = -t1149 * t667 + t1154 * t669;
t618 = t1149 * t669 + t1154 * t667;
t617 = -pkin(1) * t677 - pkin(2) * t733 - pkin(3) * t799 + pkin(4) * t906 - pkin(9) * t865 - t1240;
t616 = -t1149 * t665 + t1154 * t666;
t615 = t1149 * t666 + t1154 * t665;
t614 = -pkin(1) * t675 - pkin(2) * t728 - pkin(3) * t788 - pkin(4) * t904 - pkin(9) * t859 + t1233;
t613 = -t1150 * t657 + t1155 * t659;
t612 = -t1150 * t656 + t1155 * t658;
t609 = -t1150 * t647 + t1155 * t649;
t608 = -t1150 * t646 + t1155 * t648;
t607 = -t1150 * t640 + t1155 * t641;
t606 = t1150 * t641 + t1155 * t640;
t605 = -pkin(1) * t631 - pkin(2) * t690 - pkin(3) * t761;
t604 = -pkin(5) * t760 + pkin(10) * t611;
t603 = -qJ(4) * t665 + (-pkin(9) * t1146 + t1247) * t686;
t602 = -t1150 * t635 + t1155 * t636;
t601 = t1150 * t636 + t1155 * t635;
t600 = -pkin(10) * t741 - t610;
t599 = -pkin(5) * t868 + pkin(10) * t743 + t611;
t598 = -pkin(8) * t733 - t1149 * t639 + t1154 * t651;
t597 = -pkin(1) * t668 - pkin(2) * t720 - pkin(3) * t782 - pkin(4) * t947 - pkin(9) * t828 - t687;
t596 = -pkin(8) * t728 - t1149 * t638 + t1154 * t650;
t595 = t1151 * t756 + t1156 * t607;
t594 = t1151 * t607 - t1156 * t756;
t593 = -pkin(2) * t864 + pkin(8) * t734 + t1149 * t651 + t1154 * t639;
t592 = -pkin(7) * t730 - t1150 * t637 + t1155 * t642;
t591 = -pkin(2) * t858 + pkin(8) * t729 + t1149 * t650 + t1154 * t638;
t590 = t1151 * t746 + t1156 * t602;
t589 = t1151 * t602 - t1156 * t746;
t588 = qJ(4) * t666 + (-pkin(9) * t1145 + t1194) * t686;
t587 = -pkin(8) * t720 - t1149 * t626 + t1154 * t633;
t586 = -pkin(2) * t826 + pkin(8) * t721 + t1149 * t633 + t1154 * t626;
t585 = -t1150 * t620 + t1155 * t621;
t584 = -qJ(4) * t695 - t1145 * t625 + t1146 * t624;
t583 = -t1150 * t618 + t1155 * t619;
t582 = t1150 * t619 + t1155 * t618;
t581 = -t1150 * t615 + t1155 * t616;
t580 = t1150 * t616 + t1155 * t615;
t579 = -qJ(4) * t692 - t1145 * t623 + t1146 * t622;
t578 = -pkin(3) * t756 + qJ(4) * t696 + t1145 * t624 + t1146 * t625;
t577 = t1153 * t611 - t1241;
t576 = t1148 * t611 + t1234;
t575 = -pkin(7) * t631 - t1150 * t627 + t1155 * t630;
t574 = -pkin(3) * t746 + qJ(4) * t693 + t1145 * t622 + t1146 * t623;
t573 = t1145 * t760 + t1146 * t577;
t572 = t1145 * t577 - t1146 * t760;
t571 = t1151 * t686 + t1156 * t581;
t570 = t1151 * t581 - t1156 * t686;
t569 = t1151 * t680 + t1156 * t583;
t568 = t1151 * t583 - t1156 * t680;
t567 = -pkin(1) * t606 - pkin(2) * t640 - pkin(3) * t695 + pkin(4) * t1255 - pkin(9) * t757 - t1148 * t707 - t1153 * t689;
t566 = -pkin(9) * t680 - t1148 * t599 + t1153 * t600;
t565 = -pkin(4) * t576 - pkin(5) * t610;
t564 = -pkin(7) * t677 - t1150 * t593 + t1155 * t598;
t563 = -pkin(1) * t601 - pkin(2) * t635 - pkin(3) * t692 + pkin(4) * t807 - pkin(9) * t747 - t1148 * t701 - t1153 * t685;
t562 = -pkin(7) * t675 - t1150 * t591 + t1155 * t596;
t561 = -pkin(8) * t615 - t1149 * t588 + t1154 * t603;
t560 = -pkin(1) * t580 - pkin(2) * t615 - pkin(3) * t665 + pkin(4) * t791 - pkin(9) * t687;
t559 = -pkin(7) * t668 - t1150 * t586 + t1155 * t587;
t558 = -qJ(4) * t667 - t1145 * t634 + t1146 * t566;
t557 = -pkin(2) * t686 + pkin(8) * t616 + t1149 * t603 + t1154 * t588;
t556 = -pkin(9) * t576 - pkin(10) * t1234 - t1148 * t604;
t555 = -t1149 * t572 + t1154 * t573;
t554 = t1149 * t573 + t1154 * t572;
t553 = -pkin(3) * t680 + qJ(4) * t669 + t1145 * t566 + t1146 * t634;
t552 = -pkin(8) * t640 - t1149 * t578 + t1154 * t584;
t551 = -pkin(2) * t756 + pkin(8) * t641 + t1149 * t584 + t1154 * t578;
t550 = -pkin(8) * t635 - t1149 * t574 + t1154 * t579;
t549 = -pkin(2) * t746 + pkin(8) * t636 + t1149 * t579 + t1154 * t574;
t548 = -pkin(1) * t582 - pkin(2) * t618 - pkin(3) * t667 + pkin(4) * t868 - pkin(9) * t682 - t1148 * t600 - t1153 * t599;
t547 = -t1150 * t554 + t1155 * t555;
t546 = t1150 * t555 + t1155 * t554;
t545 = -pkin(8) * t618 - t1149 * t553 + t1154 * t558;
t544 = -qJ(4) * t572 - t1145 * t565 + t1146 * t556;
t543 = -pkin(7) * t580 - t1150 * t557 + t1155 * t561;
t542 = -pkin(2) * t680 + pkin(8) * t619 + t1149 * t558 + t1154 * t553;
t541 = -pkin(7) * t606 - t1150 * t551 + t1155 * t552;
t540 = -pkin(7) * t601 - t1150 * t549 + t1155 * t550;
t539 = t1151 * t576 + t1156 * t547;
t538 = t1151 * t547 - t1156 * t576;
t537 = -pkin(3) * t576 + qJ(4) * t573 + t1145 * t556 + t1146 * t565;
t536 = -pkin(1) * t546 - pkin(2) * t554 - pkin(3) * t572 + pkin(4) * t760 - pkin(9) * t577 + pkin(10) * t1241 - t1153 * t604;
t535 = -pkin(7) * t582 - t1150 * t542 + t1155 * t545;
t534 = -pkin(8) * t554 - t1149 * t537 + t1154 * t544;
t533 = -pkin(2) * t576 + pkin(8) * t555 + t1149 * t544 + t1154 * t537;
t532 = -pkin(7) * t546 - t1150 * t533 + t1155 * t534;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1116, -t1117, 0, t1081, 0, 0, 0, 0, 0, 0, t1054, t1055, t1078, t1009, 0, 0, 0, 0, 0, 0, t911, t917, t863, t824, 0, 0, 0, 0, 0, 0, t770, t779, t723, t629, 0, 0, 0, 0, 0, 0, t661, t663, t644, t571, 0, 0, 0, 0, 0, 0, t590, t595, t569, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1117, -t1116, 0, t1080, 0, 0, 0, 0, 0, 0, t1052, t1053, t1077, t1008, 0, 0, 0, 0, 0, 0, t910, t916, t862, t823, 0, 0, 0, 0, 0, 0, t769, t778, t722, t628, 0, 0, 0, 0, 0, 0, t660, t662, t643, t570, 0, 0, 0, 0, 0, 0, t589, t594, t568, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1082, t1083, 0, -t1029, 0, 0, 0, 0, 0, 0, t938, t962, t889, t833, 0, 0, 0, 0, 0, 0, t780, t815, t730, t631, 0, 0, 0, 0, 0, 0, t675, t677, t668, t580, 0, 0, 0, 0, 0, 0, t601, t606, t582, t546; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1117, 0, -t1116, 0, t1174, -t1100, -t1080, -pkin(6) * t1080, t1092 * t1156 - t1177, t1076 * t1156 - t1119 * t1151, t1086 * t1156 + t1150 * t1199, t1091 * t1156 + t1177, t1084 * t1156 + t1138 * t1151, qJDD(2) * t1151 + t1108 * t1156, -pkin(6) * t1052 - t1039 * t1151 + t1046 * t1156, -pkin(6) * t1053 - t1040 * t1151 + t1047 * t1156, -pkin(6) * t1077 + t1029 * t1156, -pkin(6) * t1008 - (pkin(1) * t1151 - pkin(7) * t1156) * t1029, t1156 * t941 - t1187, -t1074 * t1151 + t1156 * t890, -t1019 * t1151 + t1156 * t971, t1156 * t940 + t1187, -t1015 * t1151 + t1156 * t972, t1156 * t992 + t1178, -pkin(6) * t910 - t1151 * t857 + t1156 * t831, -pkin(6) * t916 - t1151 * t867 + t1156 * t842, -pkin(6) * t862 - t1151 * t849 + t1156 * t753, -pkin(6) * t823 - t1151 * t777 + t1156 * t759, t1156 * t785 + t1189, -t1004 * t1151 + t1156 * t731, -t1151 * t957 + t1156 * t821, t1156 * t784 - t1189, -t1151 * t953 + t1156 * t822, t1156 * t851 + t1178, -pkin(6) * t769 - t1151 * t694 + t1156 * t645, -pkin(6) * t778 - t1151 * t708 + t1156 * t664, -pkin(6) * t722 - t1151 * t688 + t1156 * t592, -pkin(6) * t628 - t1151 * t605 + t1156 * t575, -t1151 * t901 + t1156 * t706, -t1151 * t825 + t1156 * t674, -t1151 * t876 + t1156 * t683, -t1151 * t899 + t1156 * t705, -t1151 * t877 + t1156 * t684, -t1151 * t930 + t1156 * t735, -pkin(6) * t660 - t1151 * t614 + t1156 * t562, -pkin(6) * t662 - t1151 * t617 + t1156 * t564, -pkin(6) * t643 - t1151 * t597 + t1156 * t559, -pkin(6) * t570 - t1151 * t560 + t1156 * t543, -t1151 * t737 + t1156 * t613, -t1151 * t679 + t1156 * t585, -t1151 * t765 + t1156 * t608, -t1151 * t736 + t1156 * t612, -t1151 * t766 + t1156 * t609, -t1151 * t793 + t1156 * t655, -pkin(6) * t589 - t1151 * t563 + t1156 * t540, -pkin(6) * t594 - t1151 * t567 + t1156 * t541, -pkin(6) * t568 - t1151 * t548 + t1156 * t535, -pkin(6) * t538 - t1151 * t536 + t1156 * t532; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1116, 0, t1117, 0, t1100, t1174, t1081, pkin(6) * t1081, t1092 * t1151 + t1176, t1076 * t1151 + t1119 * t1156, t1086 * t1151 - t1150 * t1198, t1091 * t1151 - t1176, t1084 * t1151 - t1138 * t1156, -qJDD(2) * t1156 + t1108 * t1151, pkin(6) * t1054 + t1039 * t1156 + t1046 * t1151, pkin(6) * t1055 + t1040 * t1156 + t1047 * t1151, pkin(6) * t1078 + t1029 * t1151, pkin(6) * t1009 - (-pkin(1) * t1156 - pkin(7) * t1151) * t1029, t1151 * t941 + t1186, t1074 * t1156 + t1151 * t890, t1019 * t1156 + t1151 * t971, t1151 * t940 - t1186, t1015 * t1156 + t1151 * t972, t1151 * t992 - t1133, pkin(6) * t911 + t1151 * t831 + t1156 * t857, pkin(6) * t917 + t1151 * t842 + t1156 * t867, pkin(6) * t863 + t1151 * t753 + t1156 * t849, pkin(6) * t824 + t1151 * t759 + t1156 * t777, t1151 * t785 - t1188, t1004 * t1156 + t1151 * t731, t1151 * t821 + t1156 * t957, t1151 * t784 + t1188, t1151 * t822 + t1156 * t953, t1151 * t851 - t1133, pkin(6) * t770 + t1151 * t645 + t1156 * t694, pkin(6) * t779 + t1151 * t664 + t1156 * t708, pkin(6) * t723 + t1151 * t592 + t1156 * t688, pkin(6) * t629 + t1151 * t575 + t1156 * t605, t1151 * t706 + t1156 * t901, t1151 * t674 + t1156 * t825, t1151 * t683 + t1156 * t876, t1151 * t705 + t1156 * t899, t1151 * t684 + t1156 * t877, t1151 * t735 + t1156 * t930, pkin(6) * t661 + t1151 * t562 + t1156 * t614, pkin(6) * t663 + t1151 * t564 + t1156 * t617, pkin(6) * t644 + t1151 * t559 + t1156 * t597, pkin(6) * t571 + t1151 * t543 + t1156 * t560, t1151 * t613 + t1156 * t737, t1151 * t585 + t1156 * t679, t1151 * t608 + t1156 * t765, t1151 * t612 + t1156 * t736, t1151 * t609 + t1156 * t766, t1151 * t655 + t1156 * t793, pkin(6) * t590 + t1151 * t540 + t1156 * t563, pkin(6) * t595 + t1151 * t541 + t1156 * t567, pkin(6) * t569 + t1151 * t535 + t1156 * t548, pkin(6) * t539 + t1151 * t532 + t1156 * t536; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1122, t1123, 0, 0, (t1112 + t1185) * t1150, t1111 * t1155 + t1114 * t1150, t1125 * t1155 + t1213, (t1113 - t1192) * t1155, t1127 * t1150 + t1211, 0, pkin(1) * t1114 + pkin(7) * t1085 + t1214, -pkin(1) * t1111 + pkin(7) * t1087 - t1215, pkin(1) * t1118 + pkin(7) * t1115 + t1030, pkin(1) * t1105 + pkin(7) * t1030, t1012 * t1155 + t1013 * t1150, t1150 * t960 + t1155 * t958, t1022 * t1155 + t1024 * t1150, t1010 * t1155 + t1011 * t1150, t1023 * t1155 + t1025 * t1150, t1048 * t1155 + t1049 * t1150, -pkin(1) * t1014 + pkin(7) * t939 + t1150 * t973 + t1155 * t925, -pkin(1) * t1251 + pkin(7) * t963 + t1150 * t975 + t1155 * t929, -pkin(1) * t1043 + pkin(7) * t891 + t1150 * t852 + t1155 * t850, pkin(1) * t1044 + pkin(7) * t834 - pkin(8) * t1237 + t1155 * t892, t1150 * t875 + t1155 * t873, t1150 * t797 + t1155 * t795, t1150 * t897 + t1155 * t895, t1150 * t874 + t1155 * t872, t1150 * t898 + t1155 * t896, t1150 * t927 + t1155 * t926, -pkin(1) * t952 + pkin(7) * t781 + t1150 * t745 + t1155 * t724, pkin(1) * t1182 + pkin(7) * t816 + t1150 * t758 + t1155 * t744, -pkin(1) * t977 + pkin(7) * t732 + t1150 * t642 + t1155 * t637, pkin(1) * t942 + pkin(7) * t632 + t1150 * t630 + t1155 * t627, t1150 * t774 + t1155 * t772, t1150 * t727 + t1155 * t726, t1150 * t750 + t1155 * t748, t1150 * t773 + t1155 * t771, t1150 * t751 + t1155 * t749, t1150 * t802 + t1155 * t801, -pkin(1) * t858 + pkin(7) * t676 + t1150 * t596 + t1155 * t591, -pkin(1) * t864 + pkin(7) * t678 + t1150 * t598 + t1155 * t593, -pkin(1) * t826 + pkin(7) * t670 + t1150 * t587 + t1155 * t586, -pkin(1) * t686 + pkin(7) * t581 + t1150 * t561 + t1155 * t557, t1150 * t659 + t1155 * t657, t1150 * t621 + t1155 * t620, t1150 * t648 + t1155 * t646, t1150 * t658 + t1155 * t656, t1150 * t649 + t1155 * t647, t1150 * t710 + t1155 * t709, -pkin(1) * t746 + pkin(7) * t602 + t1150 * t550 + t1155 * t549, -pkin(1) * t756 + pkin(7) * t607 + t1150 * t552 + t1155 * t551, -pkin(1) * t680 + pkin(7) * t583 + t1150 * t545 + t1155 * t542, -pkin(1) * t576 + pkin(7) * t547 + t1150 * t534 + t1155 * t533;];
tauB_reg  = t1;