% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:38
% EndTime: 2020-01-03 11:44:00
% DurationCPUTime: 23.12s
% Computational Cost: add. (115607->713), mult. (302492->995), div. (0->0), fcn. (213627->10), ass. (0->488)
t1200 = cos(pkin(8));
t1181 = qJD(1) * t1200 - qJD(3);
t1310 = t1181 ^ 2;
t1197 = sin(pkin(9));
t1199 = cos(pkin(9));
t1202 = sin(qJ(3));
t1205 = cos(qJ(3));
t1198 = sin(pkin(8));
t1294 = qJD(1) * t1198;
t1136 = (t1197 * t1205 + t1199 * t1202) * t1294;
t1311 = t1136 ^ 2;
t1076 = -t1310 - t1311;
t1293 = qJD(1) * t1205;
t1254 = t1198 * t1293;
t1269 = t1198 * t1202;
t1255 = qJD(1) * t1269;
t1137 = -t1197 * t1255 + t1199 * t1254;
t1093 = t1137 * t1136;
t1260 = t1200 * qJDD(1);
t1180 = -qJDD(3) + t1260;
t1313 = -t1093 - t1180;
t1318 = t1199 * t1313;
t1013 = t1076 * t1197 + t1318;
t1203 = sin(qJ(1));
t1206 = cos(qJ(1));
t1170 = g(2) * t1203 - t1206 * g(3);
t1207 = qJD(1) ^ 2;
t1151 = -pkin(1) * t1207 + qJDD(1) * qJ(2) - t1170;
t1304 = pkin(2) * t1200;
t1234 = -pkin(6) * t1198 - t1304;
t1308 = 2 * qJD(2);
t1264 = t1234 * qJD(1) + t1308;
t1231 = t1264 * qJD(1) + t1151;
t1302 = t1198 * g(1);
t1087 = t1200 * t1231 - t1302;
t1171 = g(2) * t1206 + g(3) * t1203;
t1196 = qJDD(1) * pkin(1);
t1144 = -t1207 * qJ(2) + qJDD(2) + t1171 - t1196;
t1124 = t1234 * qJDD(1) + t1144;
t1111 = t1205 * t1124;
t1185 = t1198 * qJDD(1);
t1177 = t1205 * t1185;
t1143 = -qJD(3) * t1255 + t1177;
t1193 = t1198 ^ 2;
t1256 = t1181 * t1294;
t1265 = t1205 * t1207;
t1005 = -pkin(3) * t1180 - qJ(4) * t1143 + t1111 + (-pkin(3) * t1193 * t1265 + qJ(4) * t1256 - t1087) * t1202;
t1040 = t1205 * t1087 + t1202 * t1124;
t1140 = -pkin(3) * t1181 - qJ(4) * t1254;
t1272 = t1193 * t1207;
t1309 = t1202 ^ 2;
t1179 = t1309 * t1272;
t1261 = qJDD(1) * t1202;
t1227 = qJD(3) * t1293 + t1261;
t1221 = t1227 * t1198;
t1008 = -pkin(3) * t1179 - qJ(4) * t1221 + t1181 * t1140 + t1040;
t1233 = -0.2e1 * qJD(4) * t1137 + t1199 * t1005 - t1008 * t1197;
t1320 = pkin(3) * t1013 + t1233;
t1319 = t1197 * t1313;
t1201 = sin(qJ(5));
t1204 = cos(qJ(5));
t1083 = t1204 * t1136 + t1137 * t1201;
t1085 = -t1136 * t1201 + t1137 * t1204;
t1021 = t1085 * t1083;
t1169 = -qJDD(5) + t1180;
t1315 = -t1021 - t1169;
t1317 = t1201 * t1315;
t1316 = t1204 * t1315;
t1153 = t1181 * t1255;
t1115 = t1143 + t1153;
t1104 = t1115 * t1205;
t1096 = t1199 * t1143 - t1197 * t1221;
t1240 = t1143 * t1197 + t1199 * t1221;
t1001 = -t1083 * qJD(5) + t1204 * t1096 - t1201 * t1240;
t1174 = -qJD(5) + t1181;
t1070 = t1083 * t1174;
t1314 = t1070 + t1001;
t1123 = t1136 * t1181;
t1056 = -t1123 + t1096;
t1312 = t1123 + t1096;
t1241 = t1201 * t1096 + t1204 * t1240;
t969 = (qJD(5) + t1174) * t1085 + t1241;
t1080 = t1083 ^ 2;
t1081 = t1085 ^ 2;
t1133 = t1137 ^ 2;
t1168 = t1174 ^ 2;
t1292 = qJD(4) * t1136;
t1128 = -0.2e1 * t1292;
t1262 = t1197 * t1005 + t1199 * t1008;
t936 = t1128 + t1262;
t885 = t1197 * t936 + t1199 * t1233;
t1307 = pkin(3) * t885;
t1275 = t1181 * t1137;
t1052 = t1240 + t1275;
t997 = -t1052 * t1197 - t1056 * t1199;
t1306 = pkin(3) * t997;
t1305 = pkin(2) * t1198;
t1301 = t1200 * g(1);
t900 = pkin(4) * t1313 - pkin(7) * t1056 + t1233;
t1109 = -pkin(4) * t1181 - pkin(7) * t1137;
t912 = -t1311 * pkin(4) - pkin(7) * t1240 + t1181 * t1109 + t936;
t864 = t1201 * t912 - t1204 * t900;
t865 = t1201 * t900 + t1204 * t912;
t829 = t1201 * t865 - t1204 * t864;
t1300 = t1197 * t829;
t1299 = t1199 * t829;
t1020 = t1301 + qJDD(4) + pkin(3) * t1221 - qJ(4) * t1179 + (t1151 + (t1140 * t1205 + t1264) * qJD(1)) * t1198;
t955 = pkin(4) * t1240 - t1311 * pkin(7) + t1109 * t1137 + t1020;
t1298 = t1201 * t955;
t1297 = t1202 * t885;
t1296 = t1204 * t955;
t1295 = t1205 * t885;
t1016 = -t1021 + t1169;
t1290 = t1016 * t1201;
t1289 = t1016 * t1204;
t1288 = t1020 * t1197;
t1287 = t1020 * t1199;
t1077 = -t1093 + t1180;
t1286 = t1077 * t1197;
t1285 = t1077 * t1199;
t1086 = t1198 * t1231 + t1301;
t1284 = t1086 * t1202;
t1283 = t1086 * t1205;
t1248 = t1202 * t1265;
t1167 = t1193 * t1248;
t1141 = -t1167 + t1180;
t1282 = t1141 * t1202;
t1281 = t1141 * t1205;
t1142 = -t1167 - t1180;
t1280 = t1142 * t1202;
t1279 = t1142 * t1205;
t1278 = t1174 * t1085;
t1277 = t1174 * t1201;
t1276 = t1174 * t1204;
t1274 = t1181 * t1197;
t1273 = t1181 * t1199;
t1271 = t1198 * t1180;
t1270 = t1198 * t1200;
t1268 = t1198 * t1203;
t1267 = t1198 * t1206;
t1266 = t1200 * t1144;
t1259 = t1203 * qJDD(1);
t830 = t1201 * t864 + t1204 * t865;
t806 = t1197 * t830 + t1299;
t828 = pkin(4) * t829;
t1258 = pkin(3) * t806 + t828;
t972 = -t1070 + t1001;
t909 = -t1201 * t969 - t1204 * t972;
t911 = t1201 * t972 - t1204 * t969;
t868 = t1197 * t911 + t1199 * t909;
t906 = pkin(4) * t909;
t1257 = pkin(3) * t868 + t906;
t1253 = t1198 * t1021;
t1252 = t1200 * t1021;
t1251 = t1198 * t1093;
t1250 = t1200 * t1093;
t1195 = t1205 ^ 2;
t1249 = t1195 * t1272;
t1247 = t1198 * t1260;
t1229 = t1206 * t1207 + t1259;
t1246 = -pkin(5) * t1229 + t1206 * g(1);
t1245 = -t1144 + t1196;
t886 = -t1197 * t1233 + t1199 * t936;
t1244 = qJD(1) * (qJD(3) - t1181);
t1243 = qJD(1) * t1308 + t1151;
t1039 = t1087 * t1202 - t1111;
t988 = t1039 * t1202 + t1205 * t1040;
t1119 = t1243 * t1198 + t1301;
t1120 = t1243 * t1200 - t1302;
t1064 = t1119 * t1198 + t1200 * t1120;
t1239 = -t1203 * t1170 - t1171 * t1206;
t1192 = t1198 * t1193;
t1238 = t1192 * t1248;
t1237 = -pkin(2) * t1086 + pkin(6) * t988;
t1105 = -t1133 - t1310;
t1023 = t1105 * t1199 + t1286;
t1236 = pkin(3) * t1023 - t1262;
t1015 = -t1168 - t1080;
t953 = t1015 * t1201 + t1316;
t1235 = pkin(4) * t953 - t864;
t1232 = t1200 * t1167;
t987 = -t1039 * t1205 + t1040 * t1202;
t1063 = t1119 * t1200 - t1120 * t1198;
t1230 = t1170 * t1206 - t1171 * t1203;
t1166 = qJDD(1) * t1206 - t1203 * t1207;
t1061 = -t1081 - t1168;
t977 = t1061 * t1204 + t1290;
t1228 = pkin(4) * t977 - t865;
t1209 = t1200 ^ 2;
t1156 = (t1193 + t1209) * t1200 * t1207;
t1226 = -t1156 * t1203 + t1206 * t1260;
t1225 = t1156 * t1206 + t1200 * t1259;
t954 = t1015 * t1204 - t1317;
t894 = t1197 * t954 + t1199 * t953;
t1224 = pkin(3) * t894 + t1235;
t1132 = -t1310 - t1249;
t1089 = -t1132 * t1202 + t1281;
t1116 = t1244 * t1269 - t1177;
t1223 = pkin(2) * t1116 + pkin(6) * t1089 + t1284;
t1147 = -t1179 - t1310;
t1102 = t1147 * t1205 - t1280;
t1154 = t1181 * t1254;
t1114 = t1154 - t1221;
t1222 = pkin(2) * t1114 + pkin(6) * t1102 - t1283;
t980 = -t1061 * t1201 + t1289;
t915 = t1197 * t980 + t1199 * t977;
t1220 = pkin(3) * t915 + t1228;
t1112 = -t1143 + t1153;
t1113 = t1154 + t1221;
t1059 = -t1112 * t1202 - t1113 * t1205;
t1149 = t1179 + t1249;
t1219 = pkin(2) * t1149 + pkin(6) * t1059 + t988;
t807 = t1199 * t830 - t1300;
t818 = -pkin(4) * t955 + pkin(7) * t830;
t790 = -pkin(3) * t955 - pkin(7) * t1300 + qJ(4) * t807 + t1199 * t818;
t793 = -pkin(7) * t1299 - qJ(4) * t806 - t1197 * t818;
t797 = -t1202 * t806 + t1205 * t807;
t1218 = -pkin(2) * t955 + pkin(6) * t797 + t1202 * t793 + t1205 * t790;
t993 = -t1080 - t1081;
t814 = -pkin(4) * t993 + pkin(7) * t911 + t830;
t816 = -pkin(7) * t909 - t829;
t870 = -t1197 * t909 + t1199 * t911;
t800 = -pkin(3) * t993 + qJ(4) * t870 + t1197 * t816 + t1199 * t814;
t801 = -qJ(4) * t868 - t1197 * t814 + t1199 * t816;
t834 = -t1202 * t868 + t1205 * t870;
t1217 = -pkin(2) * t993 + pkin(6) * t834 + t1202 * t801 + t1205 * t800;
t968 = (qJD(5) - t1174) * t1085 + t1241;
t882 = -pkin(4) * t968 + pkin(7) * t954 - t1296;
t895 = -t1197 * t953 + t1199 * t954;
t897 = -pkin(7) * t953 + t1298;
t823 = -pkin(3) * t968 + qJ(4) * t895 + t1197 * t897 + t1199 * t882;
t836 = -qJ(4) * t894 - t1197 * t882 + t1199 * t897;
t858 = -t1202 * t894 + t1205 * t895;
t1216 = -pkin(2) * t968 + pkin(6) * t858 + t1202 * t836 + t1205 * t823;
t884 = -pkin(4) * t1314 + pkin(7) * t980 + t1298;
t907 = -pkin(7) * t977 + t1296;
t916 = -t1197 * t977 + t1199 * t980;
t835 = -pkin(3) * t1314 + qJ(4) * t916 + t1197 * t907 + t1199 * t884;
t838 = -qJ(4) * t915 - t1197 * t884 + t1199 * t907;
t873 = -t1202 * t915 + t1205 * t916;
t1215 = -pkin(2) * t1314 + pkin(6) * t873 + t1202 * t838 + t1205 * t835;
t1051 = t1240 - t1275;
t1014 = t1076 * t1199 - t1319;
t937 = -pkin(3) * t1051 + qJ(4) * t1014 - t1287;
t952 = -t1013 * t1202 + t1014 * t1205;
t957 = -qJ(4) * t1013 + t1288;
t1214 = -pkin(2) * t1051 + pkin(6) * t952 + t1202 * t957 + t1205 * t937;
t1024 = -t1105 * t1197 + t1285;
t940 = -pkin(3) * t1312 + qJ(4) * t1024 + t1288;
t967 = -qJ(4) * t1023 + t1287;
t975 = -t1023 * t1202 + t1024 * t1205;
t1213 = -pkin(2) * t1312 + pkin(6) * t975 + t1202 * t967 + t1205 * t940;
t1065 = -t1133 - t1311;
t999 = -t1052 * t1199 + t1056 * t1197;
t866 = -pkin(3) * t1065 + qJ(4) * t999 + t886;
t871 = -qJ(4) * t997 - t885;
t931 = -t1202 * t997 + t1205 * t999;
t1212 = -pkin(2) * t1065 + pkin(6) * t931 + t1202 * t871 + t1205 * t866;
t842 = t1205 * t886 - t1297;
t880 = -pkin(3) * t1020 + qJ(4) * t886;
t1211 = -pkin(2) * t1020 + pkin(6) * t842 - qJ(4) * t1297 + t1205 * t880;
t1187 = t1209 * t1207;
t1186 = t1209 * qJDD(1);
t1184 = t1193 * qJDD(1);
t1176 = t1207 * t1270;
t1172 = 0.2e1 * t1247;
t1163 = -t1187 + t1272;
t1162 = t1187 + t1272;
t1161 = t1200 * t1180;
t1160 = t1186 - t1184;
t1159 = t1186 + t1184;
t1155 = (t1198 * t1209 + t1192) * t1207;
t1150 = -t1179 + t1249;
t1148 = t1310 - t1249;
t1146 = t1179 - t1310;
t1145 = pkin(5) * t1166 + g(1) * t1203;
t1135 = t1166 * t1270;
t1134 = t1229 * t1270;
t1126 = -t1155 * t1206 - t1198 * t1259;
t1125 = t1155 * t1203 - t1206 * t1185;
t1121 = (-t1195 - t1309) * t1256;
t1118 = -t1133 + t1310;
t1117 = -t1310 + t1311;
t1107 = t1143 * t1202 - t1195 * t1256;
t1106 = (-t1309 * t1181 * qJD(1) - t1205 * t1227) * t1198;
t1103 = (t1205 * t1244 + t1261) * t1269;
t1101 = t1146 * t1205 + t1282;
t1100 = -t1148 * t1202 + t1279;
t1099 = t1147 * t1202 + t1279;
t1098 = t1146 * t1202 - t1281;
t1097 = t1148 * t1205 + t1280;
t1095 = -qJ(2) * t1156 + t1245 * t1200;
t1094 = qJ(2) * t1155 - t1245 * t1198;
t1091 = t1133 - t1311;
t1088 = t1132 * t1205 + t1282;
t1075 = t1104 * t1200 + t1238;
t1074 = t1103 * t1200 - t1238;
t1073 = t1104 * t1198 - t1232;
t1072 = t1103 * t1198 + t1232;
t1069 = -t1081 + t1168;
t1068 = t1080 - t1168;
t1067 = (t1136 * t1199 - t1137 * t1197) * t1181;
t1066 = (t1136 * t1197 + t1137 * t1199) * t1181;
t1060 = t1114 * t1205 - t1115 * t1202;
t1058 = t1114 * t1202 + t1104;
t1057 = t1112 * t1205 - t1113 * t1202;
t1050 = t1096 * t1199 + t1137 * t1274;
t1049 = t1096 * t1197 - t1137 * t1273;
t1048 = -t1136 * t1273 + t1197 * t1240;
t1047 = -t1136 * t1274 - t1199 * t1240;
t1046 = t1102 * t1200 - t1114 * t1198;
t1045 = t1101 * t1200 - t1113 * t1198;
t1044 = t1100 * t1200 - t1112 * t1198;
t1043 = t1102 * t1198 + t1114 * t1200;
t1042 = t1101 * t1198 + t1113 * t1200;
t1041 = t1100 * t1198 + t1112 * t1200;
t1037 = t1089 * t1200 - t1116 * t1198;
t1036 = t1089 * t1198 + t1116 * t1200;
t1035 = -pkin(1) * t1144 + qJ(2) * t1064;
t1034 = t1117 * t1199 + t1286;
t1033 = -t1118 * t1197 + t1318;
t1032 = t1117 * t1197 - t1285;
t1031 = t1118 * t1199 + t1319;
t1030 = pkin(1) * t1162 + qJ(2) * t1159 + t1064;
t1029 = -pkin(6) * t1099 + t1284;
t1028 = t1060 * t1200 + t1150 * t1198;
t1027 = t1059 * t1200 - t1149 * t1198;
t1026 = t1060 * t1198 - t1150 * t1200;
t1025 = t1059 * t1198 + t1149 * t1200;
t1022 = -pkin(6) * t1088 + t1283;
t1019 = t1081 - t1080;
t1012 = (t1083 * t1204 - t1085 * t1201) * t1174;
t1011 = (t1083 * t1201 + t1085 * t1204) * t1174;
t1010 = -pkin(2) * t1099 + t1039;
t1009 = -pkin(2) * t1088 + t1040;
t1007 = -t1066 * t1202 + t1067 * t1205;
t1006 = t1066 * t1205 + t1067 * t1202;
t1000 = -qJD(5) * t1085 - t1241;
t998 = -t1051 * t1199 - t1197 * t1312;
t996 = -t1051 * t1197 + t1199 * t1312;
t995 = t1007 * t1200 - t1271;
t994 = t1007 * t1198 + t1161;
t992 = -t1049 * t1202 + t1050 * t1205;
t991 = -t1047 * t1202 + t1048 * t1205;
t990 = t1049 * t1205 + t1050 * t1202;
t989 = t1047 * t1205 + t1048 * t1202;
t986 = t1068 * t1204 + t1290;
t985 = -t1069 * t1201 + t1316;
t984 = t1068 * t1201 - t1289;
t983 = t1069 * t1204 + t1317;
t982 = -t1032 * t1202 + t1034 * t1205;
t981 = -t1031 * t1202 + t1033 * t1205;
t979 = t1032 * t1205 + t1034 * t1202;
t978 = t1031 * t1205 + t1033 * t1202;
t974 = t1023 * t1205 + t1024 * t1202;
t966 = -pkin(1) * t1043 - t1222;
t965 = t1001 * t1204 + t1085 * t1277;
t964 = t1001 * t1201 - t1085 * t1276;
t963 = -t1000 * t1201 - t1083 * t1276;
t962 = t1000 * t1204 - t1083 * t1277;
t961 = t1200 * t992 + t1251;
t960 = t1200 * t991 - t1251;
t959 = t1198 * t992 - t1250;
t958 = t1198 * t991 + t1250;
t956 = -pkin(1) * t1036 - t1223;
t951 = t1086 * t1198 + t1200 * t988;
t950 = t1013 * t1205 + t1014 * t1202;
t949 = -t1086 * t1200 + t1198 * t988;
t947 = -pkin(6) * t1057 - t987;
t946 = -t1011 * t1197 + t1012 * t1199;
t945 = t1011 * t1199 + t1012 * t1197;
t944 = -t1052 * t1198 + t1200 * t982;
t943 = t1056 * t1198 + t1200 * t981;
t942 = t1052 * t1200 + t1198 * t982;
t941 = -t1056 * t1200 + t1198 * t981;
t939 = t1198 * t1312 + t1200 * t975;
t938 = t1198 * t975 - t1200 * t1312;
t934 = t1051 * t1198 + t1200 * t952;
t933 = -t1051 * t1200 + t1198 * t952;
t932 = -qJ(2) * t1043 - t1010 * t1198 + t1029 * t1200;
t930 = -t1202 * t996 + t1205 * t998;
t929 = t1202 * t999 + t1205 * t997;
t928 = t1202 * t998 + t1205 * t996;
t926 = -qJ(2) * t1036 - t1009 * t1198 + t1022 * t1200;
t925 = -pkin(1) * t1099 + qJ(2) * t1046 + t1010 * t1200 + t1029 * t1198;
t924 = -pkin(1) * t1025 - t1219;
t923 = t1091 * t1198 + t1200 * t930;
t922 = -t1091 * t1200 + t1198 * t930;
t921 = -t1197 * t984 + t1199 * t986;
t920 = -t1197 * t983 + t1199 * t985;
t919 = t1197 * t986 + t1199 * t984;
t918 = t1197 * t985 + t1199 * t983;
t917 = -pkin(1) * t1088 + qJ(2) * t1037 + t1009 * t1200 + t1022 * t1198;
t914 = t1065 * t1198 + t1200 * t931;
t913 = -t1065 * t1200 + t1198 * t931;
t910 = -t1201 * t1314 - t1204 * t968;
t908 = -t1201 * t968 + t1204 * t1314;
t905 = -qJ(2) * t1025 + t1057 * t1305 + t1200 * t947;
t904 = -t1197 * t964 + t1199 * t965;
t903 = -t1197 * t962 + t1199 * t963;
t902 = t1197 * t965 + t1199 * t964;
t901 = t1197 * t963 + t1199 * t962;
t898 = -pkin(1) * t949 - t1237;
t896 = -pkin(2) * t929 - t1306;
t893 = qJ(2) * t1027 + t1198 * t947 + (-pkin(1) - t1304) * t1057;
t892 = -t1202 * t945 + t1205 * t946;
t891 = t1202 * t946 + t1205 * t945;
t890 = -t1169 * t1198 + t1200 * t892;
t889 = t1169 * t1200 + t1198 * t892;
t888 = -qJ(2) * t949 + (-pkin(6) * t1200 + t1305) * t987;
t887 = -pkin(2) * t974 + t1128 - t1236;
t883 = -pkin(2) * t950 - t1320;
t881 = -pkin(6) * t974 - t1202 * t940 + t1205 * t967;
t879 = -t1202 * t919 + t1205 * t921;
t878 = -t1202 * t918 + t1205 * t920;
t877 = t1202 * t921 + t1205 * t919;
t876 = t1202 * t920 + t1205 * t918;
t875 = -pkin(6) * t950 - t1202 * t937 + t1205 * t957;
t874 = qJ(2) * t951 + (-pkin(1) + t1234) * t987;
t872 = t1202 * t916 + t1205 * t915;
t869 = -t1197 * t908 + t1199 * t910;
t867 = t1197 * t910 + t1199 * t908;
t862 = -t1202 * t902 + t1205 * t904;
t861 = -t1202 * t901 + t1205 * t903;
t860 = t1202 * t904 + t1205 * t902;
t859 = t1202 * t903 + t1205 * t901;
t857 = t1202 * t895 + t1205 * t894;
t856 = -t1198 * t969 + t1200 * t879;
t855 = t1198 * t972 + t1200 * t878;
t854 = t1198 * t879 + t1200 * t969;
t853 = t1198 * t878 - t1200 * t972;
t852 = t1200 * t862 + t1253;
t851 = t1200 * t861 - t1253;
t850 = t1198 * t862 - t1252;
t849 = t1198 * t861 + t1252;
t848 = t1198 * t1314 + t1200 * t873;
t847 = t1198 * t873 - t1200 * t1314;
t846 = -pkin(1) * t938 - t1213;
t845 = -pkin(1) * t933 - t1214;
t844 = t1198 * t968 + t1200 * t858;
t843 = t1198 * t858 - t1200 * t968;
t841 = t1202 * t886 + t1295;
t840 = t1020 * t1198 + t1200 * t842;
t839 = -t1020 * t1200 + t1198 * t842;
t837 = -qJ(2) * t938 - t1198 * t887 + t1200 * t881;
t833 = -t1202 * t867 + t1205 * t869;
t832 = t1202 * t870 + t1205 * t868;
t831 = t1202 * t869 + t1205 * t867;
t827 = -qJ(2) * t933 - t1198 * t883 + t1200 * t875;
t826 = -pkin(1) * t974 + qJ(2) * t939 + t1198 * t881 + t1200 * t887;
t825 = t1019 * t1198 + t1200 * t833;
t824 = -t1019 * t1200 + t1198 * t833;
t822 = t1198 * t993 + t1200 * t834;
t821 = t1198 * t834 - t1200 * t993;
t820 = -pkin(2) * t841 - t1307;
t819 = -pkin(6) * t929 - t1202 * t866 + t1205 * t871;
t817 = -pkin(1) * t950 + qJ(2) * t934 + t1198 * t875 + t1200 * t883;
t815 = -pkin(2) * t872 - t1220;
t813 = -pkin(1) * t913 - t1212;
t812 = -pkin(2) * t857 - t1224;
t811 = -pkin(6) * t841 - qJ(4) * t1295 - t1202 * t880;
t810 = -qJ(2) * t913 - t1198 * t896 + t1200 * t819;
t809 = -pkin(2) * t832 - t1257;
t808 = -pkin(1) * t929 + qJ(2) * t914 + t1198 * t819 + t1200 * t896;
t805 = -pkin(6) * t872 - t1202 * t835 + t1205 * t838;
t804 = -pkin(1) * t839 - t1211;
t803 = -pkin(6) * t857 - t1202 * t823 + t1205 * t836;
t802 = -pkin(1) * t847 - t1215;
t799 = -pkin(1) * t843 - t1216;
t798 = -qJ(2) * t839 - t1198 * t820 + t1200 * t811;
t796 = t1202 * t807 + t1205 * t806;
t795 = -qJ(2) * t847 - t1198 * t815 + t1200 * t805;
t794 = -pkin(1) * t841 + qJ(2) * t840 + t1198 * t811 + t1200 * t820;
t792 = t1198 * t955 + t1200 * t797;
t791 = t1198 * t797 - t1200 * t955;
t789 = -pkin(1) * t872 + qJ(2) * t848 + t1198 * t805 + t1200 * t815;
t788 = -qJ(2) * t843 - t1198 * t812 + t1200 * t803;
t787 = -pkin(1) * t857 + qJ(2) * t844 + t1198 * t803 + t1200 * t812;
t786 = -pkin(2) * t796 - t1258;
t785 = -pkin(6) * t832 - t1202 * t800 + t1205 * t801;
t784 = -pkin(1) * t821 - t1217;
t783 = -qJ(2) * t821 - t1198 * t809 + t1200 * t785;
t782 = -pkin(1) * t832 + qJ(2) * t822 + t1198 * t785 + t1200 * t809;
t781 = -pkin(6) * t796 - t1202 * t790 + t1205 * t793;
t780 = -pkin(1) * t791 - t1218;
t779 = -qJ(2) * t791 - t1198 * t786 + t1200 * t781;
t778 = -pkin(1) * t796 + qJ(2) * t792 + t1198 * t781 + t1200 * t786;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t1171, t1170, 0, 0, t1184, t1172, 0, t1186, 0, 0, t1095, t1094, t1030, t1035, t1073, t1026, t1041, t1072, t1042, t1161, t925, t917, t893, t874, t959, t922, t941, t958, t942, t994, t817, t826, t808, t794, t850, t824, t853, t849, t854, t889, t787, t789, t782, t778; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1229, 0, t1166, 0, t1246, -t1145, -t1230, -pkin(5) * t1230, t1134, t1160 * t1203 - t1163 * t1206, t1125, -t1134, -t1226, 0, -pkin(5) * t1225 + t1206 * t1119 + t1144 * t1268, -pkin(5) * t1126 + t1120 * t1206 + t1203 * t1266, t1203 * t1063 - pkin(5) * (-t1159 * t1206 + t1162 * t1203), -pkin(5) * (-t1064 * t1206 - t1144 * t1203) - (-pkin(1) * t1206 - qJ(2) * t1203) * t1063, t1075 * t1203 - t1107 * t1206, t1028 * t1203 - t1058 * t1206, t1044 * t1203 - t1097 * t1206, t1074 * t1203 - t1106 * t1206, t1045 * t1203 - t1098 * t1206, t1121 * t1206 - t1180 * t1268, t1203 * t932 + t1206 * t966 - pkin(5) * (-t1046 * t1206 - t1099 * t1203), t1203 * t926 + t1206 * t956 - pkin(5) * (-t1037 * t1206 - t1088 * t1203), t1203 * t905 + t1206 * t924 - pkin(5) * (-t1027 * t1206 - t1057 * t1203), t1203 * t888 + t1206 * t898 - pkin(5) * (-t1203 * t987 - t1206 * t951), t1203 * t961 - t1206 * t990, t1203 * t923 - t1206 * t928, t1203 * t943 - t1206 * t978, t1203 * t960 - t1206 * t989, t1203 * t944 - t1206 * t979, -t1006 * t1206 + t1203 * t995, t1203 * t827 + t1206 * t845 - pkin(5) * (-t1203 * t950 - t1206 * t934), t1203 * t837 + t1206 * t846 - pkin(5) * (-t1203 * t974 - t1206 * t939), t1203 * t810 + t1206 * t813 - pkin(5) * (-t1203 * t929 - t1206 * t914), t1203 * t798 + t1206 * t804 - pkin(5) * (-t1203 * t841 - t1206 * t840), t1203 * t852 - t1206 * t860, t1203 * t825 - t1206 * t831, t1203 * t855 - t1206 * t876, t1203 * t851 - t1206 * t859, t1203 * t856 - t1206 * t877, t1203 * t890 - t1206 * t891, t1203 * t788 + t1206 * t799 - pkin(5) * (-t1203 * t857 - t1206 * t844), t1203 * t795 + t1206 * t802 - pkin(5) * (-t1203 * t872 - t1206 * t848), t1203 * t783 + t1206 * t784 - pkin(5) * (-t1203 * t832 - t1206 * t822), t1203 * t779 + t1206 * t780 - pkin(5) * (-t1203 * t796 - t1206 * t792); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t1166, 0, t1229, 0, t1145, t1246, t1239, pkin(5) * t1239, -t1135, -t1160 * t1206 - t1163 * t1203, t1126, t1135, -t1225, 0, pkin(5) * t1226 + t1203 * t1119 - t1144 * t1267, pkin(5) * t1125 + t1120 * t1203 - t1206 * t1266, -t1206 * t1063 + pkin(5) * (t1159 * t1203 + t1162 * t1206), pkin(5) * (t1064 * t1203 - t1144 * t1206) - (-pkin(1) * t1203 + qJ(2) * t1206) * t1063, -t1075 * t1206 - t1107 * t1203, -t1028 * t1206 - t1058 * t1203, -t1044 * t1206 - t1097 * t1203, -t1074 * t1206 - t1106 * t1203, -t1045 * t1206 - t1098 * t1203, t1121 * t1203 + t1180 * t1267, -t1206 * t932 + t1203 * t966 + pkin(5) * (t1046 * t1203 - t1099 * t1206), -t1206 * t926 + t1203 * t956 + pkin(5) * (t1037 * t1203 - t1088 * t1206), -t1206 * t905 + t1203 * t924 + pkin(5) * (t1027 * t1203 - t1057 * t1206), -t1206 * t888 + t1203 * t898 + pkin(5) * (t1203 * t951 - t1206 * t987), -t1203 * t990 - t1206 * t961, -t1203 * t928 - t1206 * t923, -t1203 * t978 - t1206 * t943, -t1203 * t989 - t1206 * t960, -t1203 * t979 - t1206 * t944, -t1006 * t1203 - t1206 * t995, -t1206 * t827 + t1203 * t845 + pkin(5) * (t1203 * t934 - t1206 * t950), -t1206 * t837 + t1203 * t846 + pkin(5) * (t1203 * t939 - t1206 * t974), -t1206 * t810 + t1203 * t813 + pkin(5) * (t1203 * t914 - t1206 * t929), -t1206 * t798 + t1203 * t804 + pkin(5) * (t1203 * t840 - t1206 * t841), -t1203 * t860 - t1206 * t852, -t1203 * t831 - t1206 * t825, -t1203 * t876 - t1206 * t855, -t1203 * t859 - t1206 * t851, -t1203 * t877 - t1206 * t856, -t1203 * t891 - t1206 * t890, -t1206 * t788 + t1203 * t799 + pkin(5) * (t1203 * t844 - t1206 * t857), -t1206 * t795 + t1203 * t802 + pkin(5) * (t1203 * t848 - t1206 * t872), -t1206 * t783 + t1203 * t784 + pkin(5) * (t1203 * t822 - t1206 * t832), -t1206 * t779 + t1203 * t780 + pkin(5) * (t1203 * t792 - t1206 * t796); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1207, 0, 0, -g(1), t1171, 0, t1247, t1160, t1155, -t1247, t1156, 0, t1198 * t1144, t1266, t1063, qJ(2) * t1063, t1075, t1028, t1044, t1074, t1045, -t1271, t932, t926, t905, t888, t961, t923, t943, t960, t944, t995, t827, t837, t810, t798, t852, t825, t855, t851, t856, t890, t788, t795, t783, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1207, 0, qJDD(1), 0, g(1), 0, -t1170, 0, t1176, -t1163, -t1185, -t1176, -t1260, 0, t1119, t1120, 0, pkin(1) * t1063, -t1107, -t1058, -t1097, -t1106, -t1098, t1121, t966, t956, t924, t898, -t990, -t928, -t978, -t989, -t979, -t1006, t845, t846, t813, t804, -t860, -t831, -t876, -t859, -t877, -t891, t799, t802, t784, t780; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1171, t1170, 0, 0, t1184, t1172, 0, t1186, 0, 0, t1095, t1094, t1030, t1035, t1073, t1026, t1041, t1072, t1042, t1161, t925, t917, t893, t874, t959, t922, t941, t958, t942, t994, t817, t826, t808, t794, t850, t824, t853, t849, t854, t889, t787, t789, t782, t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, t1260, t1176, 0, t1187, 0, 0, t1144, t1119, 0, t1104, t1060, t1100, t1103, t1101, 0, t1029, t1022, t947, -pkin(6) * t987, t992, t930, t981, t991, t982, t1007, t875, t881, t819, t811, t862, t833, t878, t861, t879, t892, t803, t805, t785, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, -t1272, t1260, -t1176, 0, -t1144, 0, t1120, 0, -t1167, -t1150, t1112, t1167, t1113, t1180, t1010, t1009, -pkin(2) * t1057, -pkin(2) * t987, -t1093, -t1091, -t1056, t1093, t1052, t1180, t883, t887, t896, t820, -t1021, -t1019, -t972, t1021, t969, t1169, t812, t815, t809, t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1176, t1163, t1185, t1176, t1260, 0, -t1119, -t1120, 0, 0, t1107, t1058, t1097, t1106, t1098, -t1121, t1222, t1223, t1219, t1237, t990, t928, t978, t989, t979, t1006, t1214, t1213, t1212, t1211, t860, t831, t876, t859, t877, t891, t1216, t1215, t1217, t1218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1143, t1114, t1142, -t1153, t1146, t1153, 0, t1086, t1039, 0, t1050, t998, t1033, t1048, t1034, t1067, t957, t967, t871, -qJ(4) * t885, t904, t869, t920, t903, t921, t946, t836, t838, t801, t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1154, t1115, t1148, -t1221, -t1141, t1154, -t1086, 0, t1040, 0, t1049, t996, t1031, t1047, t1032, t1066, t937, t940, t866, t880, t902, t867, t918, t901, t919, t945, t823, t835, t800, t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1167, t1150, -t1112, -t1167, -t1113, -t1180, -t1039, -t1040, 0, 0, t1093, t1091, t1056, -t1093, -t1052, -t1180, t1320, t1236 + 0.2e1 * t1292, t1306, t1307, t1021, t1019, t972, -t1021, -t969, -t1169, t1224, t1220, t1257, t1258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1096, -t1051, t1313, -t1123, t1117, t1123, 0, t1020, -t1233, 0, t965, t910, t985, t963, t986, t1012, t897, t907, t816, -pkin(7) * t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1275, t1312, t1118, -t1240, -t1077, t1275, -t1020, 0, t936, 0, t964, t908, t983, t962, t984, t1011, t882, t884, t814, t818; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1093, t1091, t1056, -t1093, -t1052, -t1180, t1233, -t936, 0, 0, t1021, t1019, t972, -t1021, -t969, -t1169, t1235, t1228, t906, t828; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1001, -t968, t1315, -t1070, t1068, t1070, 0, t955, t864, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1278, t1314, t1069, t1000, -t1016, t1278, -t955, 0, t865, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1021, t1019, t972, -t1021, -t969, -t1169, -t864, -t865, 0, 0;];
m_new_reg = t1;
