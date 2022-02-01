% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRRR6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:52
% EndTime: 2022-01-20 12:09:09
% DurationCPUTime: 18.02s
% Computational Cost: add. (163384->651), mult. (211504->931), div. (0->0), fcn. (143913->10), ass. (0->452)
t1224 = cos(qJ(2));
t1212 = qJDD(1) + qJDD(2);
t1219 = sin(qJ(2));
t1276 = t1219 * t1212;
t1214 = qJD(1) + qJD(2);
t1305 = t1214 ^ 2;
t1179 = t1224 * t1305 + t1276;
t1159 = pkin(6) * t1179 - g(3) * t1224;
t1220 = sin(qJ(1));
t1225 = cos(qJ(1));
t1269 = t1224 * t1212;
t1182 = t1219 * t1305 - t1269;
t1237 = t1179 * t1225 - t1182 * t1220;
t1314 = pkin(6) * t1182 - g(3) * t1219;
t1325 = pkin(5) * t1237 + t1225 * t1159 - t1220 * t1314;
t1313 = t1179 * t1220 + t1182 * t1225;
t1324 = pkin(5) * t1313 + t1220 * t1159 + t1225 * t1314;
t1199 = g(1) * t1225 + g(2) * t1220;
t1227 = qJD(1) ^ 2;
t1186 = -pkin(1) * t1227 - t1199;
t1198 = g(1) * t1220 - g(2) * t1225;
t1235 = qJDD(1) * pkin(1) + t1198;
t1131 = t1186 * t1219 - t1224 * t1235;
t1132 = t1224 * t1186 + t1219 * t1235;
t1250 = t1131 * t1219 + t1132 * t1224;
t1073 = t1131 * t1224 - t1132 * t1219;
t1268 = t1225 * t1073;
t1322 = -t1220 * t1250 + t1268;
t1275 = t1220 * t1073;
t1321 = t1225 * t1250 + t1275;
t1216 = sin(qJ(5));
t1217 = sin(qJ(4));
t1222 = cos(qJ(4));
t1223 = cos(qJ(3));
t1285 = t1214 * t1223;
t1218 = sin(qJ(3));
t1286 = t1214 * t1218;
t1162 = t1217 * t1286 - t1222 * t1285;
t1163 = (t1217 * t1223 + t1218 * t1222) * t1214;
t1221 = cos(qJ(5));
t1101 = t1162 * t1221 + t1163 * t1216;
t1103 = -t1162 * t1216 + t1163 * t1221;
t1041 = t1103 * t1101;
t1211 = qJDD(3) + qJDD(4);
t1205 = qJDD(5) + t1211;
t1309 = -t1041 + t1205;
t1318 = t1216 * t1309;
t1123 = t1163 * t1162;
t1308 = -t1123 + t1211;
t1317 = t1217 * t1308;
t1316 = t1221 * t1309;
t1315 = t1222 * t1308;
t1213 = qJD(3) + qJD(4);
t1208 = qJD(5) + t1213;
t1093 = t1208 * t1101;
t1201 = qJD(3) * t1285;
t1277 = t1218 * t1212;
t1173 = t1201 + t1277;
t1261 = qJD(3) * t1286;
t1270 = t1223 * t1212;
t1231 = -t1261 + t1270;
t1078 = -qJD(4) * t1162 + t1173 * t1222 + t1217 * t1231;
t1249 = t1217 * t1173 - t1222 * t1231;
t1230 = qJD(4) * t1163 + t1249;
t993 = -qJD(5) * t1101 + t1078 * t1221 - t1216 * t1230;
t1312 = -t1093 + t993;
t1121 = -pkin(2) * t1305 + pkin(7) * t1212 + t1132;
t1280 = t1218 * t1121;
t1095 = t1223 * g(3) + t1280;
t1096 = -t1218 * g(3) + t1121 * t1223;
t1031 = t1095 * t1218 + t1096 * t1223;
t1154 = t1213 * t1162;
t1307 = -t1154 + t1078;
t1058 = t1154 + t1078;
t1251 = t1216 * t1078 + t1221 * t1230;
t965 = (qJD(5) - t1208) * t1103 + t1251;
t1054 = (qJD(4) - t1213) * t1163 + t1249;
t1099 = t1101 ^ 2;
t1100 = t1103 ^ 2;
t1306 = t1162 ^ 2;
t1161 = t1163 ^ 2;
t1204 = t1208 ^ 2;
t1210 = t1213 ^ 2;
t1304 = t1223 ^ 2;
t1047 = qJDD(3) * pkin(3) - t1173 * pkin(8) - t1280 + (-g(3) + (pkin(3) * t1286 + pkin(8) * qJD(3)) * t1214) * t1223;
t1189 = qJD(3) * pkin(3) - pkin(8) * t1286;
t1203 = t1304 * t1305;
t1048 = -pkin(3) * t1203 + pkin(8) * t1231 - qJD(3) * t1189 + t1096;
t989 = -t1047 * t1222 + t1217 * t1048;
t990 = t1047 * t1217 + t1048 * t1222;
t923 = t1217 * t990 - t1222 * t989;
t1303 = pkin(3) * t923;
t995 = -t1054 * t1217 - t1058 * t1222;
t1302 = pkin(3) * t995;
t1120 = -pkin(2) * t1212 - pkin(7) * t1305 + t1131;
t1052 = -pkin(3) * t1231 - pkin(8) * t1203 + t1189 * t1286 + t1120;
t1148 = pkin(4) * t1213 - pkin(9) * t1163;
t975 = pkin(4) * t1230 - pkin(9) * t1306 + t1148 * t1163 + t1052;
t1299 = t1216 * t975;
t937 = pkin(4) * t1308 - pkin(9) * t1058 - t989;
t939 = -pkin(4) * t1306 - pkin(9) * t1230 - t1148 * t1213 + t990;
t887 = t1216 * t939 - t1221 * t937;
t888 = t1216 * t937 + t1221 * t939;
t844 = t1216 * t888 - t1221 * t887;
t1298 = t1217 * t844;
t1297 = t1218 * t923;
t1296 = t1221 * t975;
t1295 = t1222 * t844;
t1294 = t1223 * t923;
t1293 = qJD(3) * t1214;
t1292 = t1208 * t1103;
t1291 = t1208 * t1216;
t1290 = t1208 * t1221;
t1289 = t1213 * t1163;
t1288 = t1213 * t1217;
t1287 = t1213 * t1222;
t1215 = t1218 ^ 2;
t1284 = t1215 * t1305;
t1034 = t1041 + t1205;
t1283 = t1216 * t1034;
t1282 = t1217 * t1052;
t1113 = t1123 + t1211;
t1281 = t1217 * t1113;
t1104 = t1218 * t1120;
t1197 = t1223 * t1305 * t1218;
t1187 = qJDD(3) + t1197;
t1279 = t1218 * t1187;
t1188 = qJDD(3) - t1197;
t1278 = t1218 * t1188;
t1274 = t1221 * t1034;
t1273 = t1222 * t1052;
t1272 = t1222 * t1113;
t1105 = t1223 * t1120;
t1174 = -0.2e1 * t1261 + t1270;
t1133 = t1223 * t1174;
t1271 = t1223 * t1188;
t1265 = -pkin(2) * t1120 + pkin(7) * t1031;
t1264 = t1215 + t1304;
t845 = t1216 * t887 + t1221 * t888;
t811 = t1217 * t845 + t1295;
t843 = pkin(4) * t844;
t1263 = pkin(3) * t811 + t843;
t968 = t1093 + t993;
t908 = -t1216 * t965 - t1221 * t968;
t910 = t1216 * t968 - t1221 * t965;
t860 = t1217 * t910 + t1222 * t908;
t906 = pkin(4) * t908;
t1262 = pkin(3) * t860 + t906;
t1260 = t1219 * t1041;
t1259 = t1219 * t1123;
t1258 = t1224 * t1041;
t1257 = t1224 * t1123;
t1226 = qJD(3) ^ 2;
t1193 = -t1226 - t1284;
t1144 = -t1193 * t1218 - t1271;
t1172 = 0.2e1 * t1201 + t1277;
t1256 = -pkin(2) * t1172 + pkin(7) * t1144 + t1104;
t1195 = -t1203 - t1226;
t1142 = t1195 * t1223 - t1279;
t1255 = pkin(2) * t1174 + pkin(7) * t1142 - t1105;
t924 = t1217 * t989 + t1222 * t990;
t812 = t1222 * t845 - t1298;
t837 = -pkin(4) * t975 + pkin(9) * t845;
t790 = -pkin(3) * t975 + pkin(8) * t812 - pkin(9) * t1298 + t1222 * t837;
t795 = -pkin(8) * t811 - pkin(9) * t1295 - t1217 * t837;
t798 = -t1218 * t811 + t1223 * t812;
t1254 = -pkin(2) * t975 + pkin(7) * t798 + t1218 * t795 + t1223 * t790;
t964 = (qJD(5) + t1208) * t1103 + t1251;
t1032 = -t1204 - t1099;
t978 = t1032 * t1221 - t1318;
t885 = -pkin(4) * t964 + pkin(9) * t978 - t1296;
t977 = t1032 * t1216 + t1316;
t916 = -pkin(9) * t977 + t1299;
t918 = -t1217 * t977 + t1222 * t978;
t828 = -pkin(3) * t964 + pkin(8) * t918 + t1217 * t916 + t1222 * t885;
t917 = t1217 * t978 + t1222 * t977;
t836 = -pkin(8) * t917 - t1217 * t885 + t1222 * t916;
t873 = -t1218 * t917 + t1223 * t918;
t1253 = -pkin(2) * t964 + pkin(7) * t873 + t1218 * t836 + t1223 * t828;
t1082 = -t1100 - t1204;
t1005 = -t1082 * t1216 - t1274;
t891 = -pkin(4) * t1312 + pkin(9) * t1005 + t1299;
t1004 = t1082 * t1221 - t1283;
t921 = -pkin(9) * t1004 + t1296;
t931 = -t1004 * t1217 + t1005 * t1222;
t833 = -pkin(3) * t1312 + pkin(8) * t931 + t1217 * t921 + t1222 * t891;
t930 = t1004 * t1222 + t1005 * t1217;
t841 = -pkin(8) * t930 - t1217 * t891 + t1222 * t921;
t879 = -t1218 * t930 + t1223 * t931;
t1252 = -pkin(2) * t1312 + pkin(7) * t879 + t1218 * t841 + t1223 * t833;
t1248 = -t1198 * t1220 - t1199 * t1225;
t1021 = -t1099 - t1100;
t825 = -pkin(4) * t1021 + pkin(9) * t910 + t845;
t829 = -pkin(9) * t908 - t844;
t862 = -t1217 * t908 + t1222 * t910;
t801 = -pkin(3) * t1021 + pkin(8) * t862 + t1217 * t829 + t1222 * t825;
t803 = -pkin(8) * t860 - t1217 * t825 + t1222 * t829;
t824 = -t1218 * t860 + t1223 * t862;
t1247 = -pkin(2) * t1021 + pkin(7) * t824 + t1218 * t803 + t1223 * t801;
t1081 = -t1161 - t1306;
t997 = -t1054 * t1222 + t1058 * t1217;
t890 = -pkin(3) * t1081 + pkin(8) * t997 + t924;
t894 = -pkin(8) * t995 - t923;
t929 = -t1218 * t995 + t1223 * t997;
t1246 = -pkin(2) * t1081 + pkin(7) * t929 + t1218 * t894 + t1223 * t890;
t1053 = (qJD(4) + t1213) * t1163 + t1249;
t1107 = -t1210 - t1306;
t1037 = t1107 * t1222 - t1317;
t943 = -pkin(3) * t1053 + pkin(8) * t1037 - t1273;
t1036 = t1107 * t1217 + t1315;
t981 = -t1036 * t1218 + t1037 * t1223;
t988 = -pkin(8) * t1036 + t1282;
t1245 = -pkin(2) * t1053 + pkin(7) * t981 + t1218 * t988 + t1223 * t943;
t1145 = -t1161 - t1210;
t1061 = t1145 * t1222 - t1281;
t1062 = -t1145 * t1217 - t1272;
t1002 = -t1061 * t1218 + t1062 * t1223;
t947 = -pkin(3) * t1307 + pkin(8) * t1062 + t1282;
t999 = -pkin(8) * t1061 + t1273;
t1244 = -pkin(2) * t1307 + pkin(7) * t1002 + t1218 * t999 + t1223 * t947;
t1243 = t1219 * t1197;
t1242 = t1224 * t1197;
t1177 = t1264 * t1212;
t1183 = t1203 + t1284;
t1241 = pkin(2) * t1183 + pkin(7) * t1177 + t1031;
t1240 = pkin(3) * t1061 - t990;
t1239 = pkin(4) * t977 - t887;
t1191 = qJDD(1) * t1225 - t1220 * t1227;
t1238 = -pkin(5) * t1191 - g(3) * t1220;
t1030 = t1095 * t1223 - t1096 * t1218;
t1236 = t1198 * t1225 - t1199 * t1220;
t1234 = pkin(3) * t1036 - t989;
t1233 = pkin(4) * t1004 - t888;
t876 = t1223 * t924 - t1297;
t905 = -pkin(3) * t1052 + pkin(8) * t924;
t1232 = -pkin(2) * t1052 + pkin(7) * t876 - pkin(8) * t1297 + t1223 * t905;
t1229 = pkin(3) * t917 + t1239;
t1228 = pkin(3) * t930 + t1233;
t1194 = t1203 - t1226;
t1192 = t1226 - t1284;
t1190 = qJDD(1) * t1220 + t1225 * t1227;
t1184 = -t1203 + t1284;
t1178 = t1223 * t1187;
t1170 = -pkin(5) * t1190 + g(3) * t1225;
t1169 = t1264 * t1293;
t1152 = -t1161 + t1210;
t1151 = -t1210 + t1306;
t1150 = qJDD(3) * t1219 + t1169 * t1224;
t1149 = -qJDD(3) * t1224 + t1169 * t1219;
t1147 = t1173 * t1223 - t1215 * t1293;
t1146 = -t1218 * t1231 - t1293 * t1304;
t1143 = -t1192 * t1218 + t1178;
t1141 = t1194 * t1223 - t1278;
t1140 = t1193 * t1223 - t1278;
t1139 = t1192 * t1223 + t1279;
t1138 = t1195 * t1218 + t1178;
t1137 = t1194 * t1218 + t1271;
t1134 = (t1173 + t1201) * t1218;
t1128 = t1177 * t1224 - t1183 * t1219;
t1126 = t1177 * t1219 + t1183 * t1224;
t1125 = -t1172 * t1218 + t1133;
t1124 = t1172 * t1223 + t1174 * t1218;
t1122 = t1161 - t1306;
t1119 = t1143 * t1224 + t1218 * t1276;
t1118 = t1141 * t1224 + t1219 * t1270;
t1117 = t1143 * t1219 - t1218 * t1269;
t1116 = t1141 * t1219 - t1223 * t1269;
t1111 = t1147 * t1224 - t1243;
t1110 = t1146 * t1224 + t1243;
t1109 = t1147 * t1219 + t1242;
t1108 = t1146 * t1219 - t1242;
t1098 = -pkin(1) * t1179 - t1132;
t1097 = -pkin(1) * t1182 - t1131;
t1092 = t1144 * t1224 + t1172 * t1219;
t1091 = t1142 * t1224 - t1174 * t1219;
t1090 = t1144 * t1219 - t1172 * t1224;
t1089 = t1142 * t1219 + t1174 * t1224;
t1088 = -t1100 + t1204;
t1087 = t1099 - t1204;
t1086 = (-t1162 * t1222 + t1163 * t1217) * t1213;
t1085 = (-t1162 * t1217 - t1163 * t1222) * t1213;
t1080 = t1125 * t1224 + t1184 * t1219;
t1079 = t1125 * t1219 - t1184 * t1224;
t1070 = pkin(1) * t1073;
t1069 = pkin(1) * g(3) + pkin(6) * t1250;
t1068 = t1151 * t1222 - t1281;
t1067 = -t1152 * t1217 + t1315;
t1066 = t1151 * t1217 + t1272;
t1065 = t1152 * t1222 + t1317;
t1064 = -pkin(7) * t1140 + t1105;
t1063 = -pkin(7) * t1138 + t1104;
t1060 = -pkin(2) * t1140 + t1096;
t1059 = -pkin(2) * t1138 + t1095;
t1046 = t1078 * t1222 - t1163 * t1288;
t1045 = t1078 * t1217 + t1163 * t1287;
t1044 = t1162 * t1287 + t1217 * t1230;
t1043 = t1162 * t1288 - t1222 * t1230;
t1038 = t1100 - t1099;
t1027 = (-t1101 * t1221 + t1103 * t1216) * t1208;
t1026 = (-t1101 * t1216 - t1103 * t1221) * t1208;
t1025 = -t1085 * t1218 + t1086 * t1223;
t1024 = t1085 * t1223 + t1086 * t1218;
t1023 = t1025 * t1224 + t1211 * t1219;
t1022 = t1025 * t1219 - t1211 * t1224;
t1019 = pkin(1) * t1089 + t1255;
t1018 = pkin(1) * t1090 + t1256;
t1017 = -pkin(6) * t1126 + t1030 * t1224;
t1016 = pkin(6) * t1128 + t1030 * t1219;
t1015 = t1031 * t1224 + t1120 * t1219;
t1014 = t1031 * t1219 - t1120 * t1224;
t1013 = -t1066 * t1218 + t1068 * t1223;
t1012 = -t1065 * t1218 + t1067 * t1223;
t1011 = t1066 * t1223 + t1068 * t1218;
t1010 = t1065 * t1223 + t1067 * t1218;
t1009 = t1087 * t1221 - t1283;
t1008 = -t1088 * t1216 + t1316;
t1007 = t1087 * t1216 + t1274;
t1006 = t1088 * t1221 + t1318;
t1001 = t1061 * t1223 + t1062 * t1218;
t996 = -t1053 * t1222 - t1217 * t1307;
t994 = -t1053 * t1217 + t1222 * t1307;
t992 = -qJD(5) * t1103 - t1251;
t991 = pkin(1) * t1126 + t1241;
t985 = -t1045 * t1218 + t1046 * t1223;
t984 = -t1043 * t1218 + t1044 * t1223;
t983 = t1045 * t1223 + t1046 * t1218;
t982 = t1043 * t1223 + t1044 * t1218;
t980 = t1036 * t1223 + t1037 * t1218;
t974 = -pkin(6) * t1090 - t1060 * t1219 + t1064 * t1224;
t973 = -pkin(6) * t1089 - t1059 * t1219 + t1063 * t1224;
t971 = -t1026 * t1217 + t1027 * t1222;
t970 = t1026 * t1222 + t1027 * t1217;
t961 = t1224 * t985 + t1259;
t960 = t1224 * t984 - t1259;
t959 = -t1103 * t1291 + t1221 * t993;
t958 = t1219 * t985 - t1257;
t957 = t1219 * t984 + t1257;
t956 = t1103 * t1290 + t1216 * t993;
t955 = t1101 * t1290 - t1216 * t992;
t954 = t1101 * t1291 + t1221 * t992;
t953 = -pkin(1) * t1140 + pkin(6) * t1092 + t1060 * t1224 + t1064 * t1219;
t952 = -pkin(1) * t1138 + pkin(6) * t1091 + t1059 * t1224 + t1063 * t1219;
t951 = t1013 * t1224 - t1054 * t1219;
t950 = t1012 * t1224 + t1058 * t1219;
t949 = t1013 * t1219 + t1054 * t1224;
t948 = t1012 * t1219 - t1058 * t1224;
t946 = t1002 * t1224 + t1219 * t1307;
t945 = t1002 * t1219 - t1224 * t1307;
t941 = t1053 * t1219 + t1224 * t981;
t940 = -t1053 * t1224 + t1219 * t981;
t938 = pkin(1) * t1014 + t1265;
t935 = -t1007 * t1217 + t1009 * t1222;
t934 = -t1006 * t1217 + t1008 * t1222;
t933 = t1007 * t1222 + t1009 * t1217;
t932 = t1006 * t1222 + t1008 * t1217;
t928 = -t1218 * t994 + t1223 * t996;
t927 = t1218 * t997 + t1223 * t995;
t926 = t1218 * t996 + t1223 * t994;
t922 = -pkin(6) * t1014 - (pkin(2) * t1219 - pkin(7) * t1224) * t1030;
t920 = t1122 * t1219 + t1224 * t928;
t919 = -t1122 * t1224 + t1219 * t928;
t915 = t1081 * t1219 + t1224 * t929;
t914 = -t1081 * t1224 + t1219 * t929;
t913 = -pkin(2) * t1001 - t1240;
t912 = -t1218 * t970 + t1223 * t971;
t911 = t1218 * t971 + t1223 * t970;
t909 = -t1216 * t1312 - t1221 * t964;
t907 = -t1216 * t964 + t1221 * t1312;
t903 = t1205 * t1219 + t1224 * t912;
t902 = -t1205 * t1224 + t1219 * t912;
t901 = -t1217 * t956 + t1222 * t959;
t900 = -t1217 * t954 + t1222 * t955;
t899 = t1217 * t959 + t1222 * t956;
t898 = t1217 * t955 + t1222 * t954;
t897 = -pkin(2) * t980 - t1234;
t896 = pkin(6) * t1015 - (-pkin(2) * t1224 - pkin(7) * t1219 - pkin(1)) * t1030;
t895 = -pkin(2) * t927 - t1302;
t892 = -pkin(7) * t1001 - t1218 * t947 + t1223 * t999;
t884 = -t1218 * t933 + t1223 * t935;
t883 = -t1218 * t932 + t1223 * t934;
t882 = t1218 * t935 + t1223 * t933;
t881 = t1218 * t934 + t1223 * t932;
t880 = -pkin(7) * t980 - t1218 * t943 + t1223 * t988;
t878 = t1218 * t931 + t1223 * t930;
t875 = t1218 * t924 + t1294;
t872 = t1218 * t918 + t1223 * t917;
t870 = t1052 * t1219 + t1224 * t876;
t869 = -t1052 * t1224 + t1219 * t876;
t868 = -t1219 * t965 + t1224 * t884;
t867 = t1219 * t968 + t1224 * t883;
t866 = t1219 * t884 + t1224 * t965;
t865 = t1219 * t883 - t1224 * t968;
t864 = t1219 * t1312 + t1224 * t879;
t863 = t1219 * t879 - t1224 * t1312;
t861 = -t1217 * t907 + t1222 * t909;
t859 = t1217 * t909 + t1222 * t907;
t858 = pkin(1) * t945 + t1244;
t857 = -t1218 * t899 + t1223 * t901;
t856 = -t1218 * t898 + t1223 * t900;
t855 = t1218 * t901 + t1223 * t899;
t854 = t1218 * t900 + t1223 * t898;
t853 = pkin(1) * t940 + t1245;
t852 = t1224 * t857 + t1260;
t851 = t1224 * t856 - t1260;
t850 = t1219 * t857 - t1258;
t849 = t1219 * t856 + t1258;
t848 = t1219 * t964 + t1224 * t873;
t847 = t1219 * t873 - t1224 * t964;
t846 = -pkin(2) * t875 - t1303;
t842 = -pkin(6) * t945 - t1219 * t913 + t1224 * t892;
t839 = -pkin(6) * t940 - t1219 * t897 + t1224 * t880;
t838 = -pkin(1) * t1001 + pkin(6) * t946 + t1219 * t892 + t1224 * t913;
t834 = -pkin(7) * t927 - t1218 * t890 + t1223 * t894;
t831 = -pkin(7) * t875 - pkin(8) * t1294 - t1218 * t905;
t830 = -pkin(1) * t980 + pkin(6) * t941 + t1219 * t880 + t1224 * t897;
t826 = -pkin(2) * t878 - t1228;
t823 = -t1218 * t859 + t1223 * t861;
t822 = t1218 * t862 + t1223 * t860;
t821 = t1218 * t861 + t1223 * t859;
t819 = t1038 * t1219 + t1224 * t823;
t818 = -t1038 * t1224 + t1219 * t823;
t817 = t1021 * t1219 + t1224 * t824;
t816 = -t1021 * t1224 + t1219 * t824;
t815 = pkin(1) * t914 + t1246;
t814 = -pkin(2) * t872 - t1229;
t813 = -pkin(6) * t914 - t1219 * t895 + t1224 * t834;
t810 = pkin(1) * t869 + t1232;
t809 = -pkin(1) * t927 + pkin(6) * t915 + t1219 * t834 + t1224 * t895;
t808 = -pkin(2) * t822 - t1262;
t807 = -pkin(7) * t878 - t1218 * t833 + t1223 * t841;
t806 = -pkin(6) * t869 - t1219 * t846 + t1224 * t831;
t805 = -pkin(7) * t872 - t1218 * t828 + t1223 * t836;
t804 = -pkin(1) * t875 + pkin(6) * t870 + t1219 * t831 + t1224 * t846;
t799 = pkin(1) * t863 + t1252;
t797 = t1218 * t812 + t1223 * t811;
t793 = t1219 * t975 + t1224 * t798;
t792 = t1219 * t798 - t1224 * t975;
t791 = pkin(1) * t847 + t1253;
t788 = -pkin(6) * t863 - t1219 * t826 + t1224 * t807;
t787 = -pkin(1) * t878 + pkin(6) * t864 + t1219 * t807 + t1224 * t826;
t786 = -pkin(6) * t847 - t1219 * t814 + t1224 * t805;
t785 = -pkin(1) * t872 + pkin(6) * t848 + t1219 * t805 + t1224 * t814;
t784 = -pkin(2) * t797 - t1263;
t783 = -pkin(7) * t822 - t1218 * t801 + t1223 * t803;
t782 = pkin(1) * t816 + t1247;
t781 = -pkin(7) * t797 - t1218 * t790 + t1223 * t795;
t780 = -pkin(6) * t816 - t1219 * t808 + t1224 * t783;
t779 = -pkin(1) * t822 + pkin(6) * t817 + t1219 * t783 + t1224 * t808;
t778 = pkin(1) * t792 + t1254;
t777 = -pkin(6) * t792 - t1219 * t784 + t1224 * t781;
t776 = -pkin(1) * t797 + pkin(6) * t793 + t1219 * t781 + t1224 * t784;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1191, 0, -t1190, 0, t1238, -t1170, -t1236, -pkin(5) * t1236, 0, 0, -t1313, 0, -t1237, 0, t1324, t1325, t1322, pkin(5) * t1322 + pkin(6) * t1268 - t1220 * t1069, -t1109 * t1220 + t1111 * t1225, -t1079 * t1220 + t1080 * t1225, -t1117 * t1220 + t1119 * t1225, -t1108 * t1220 + t1110 * t1225, -t1116 * t1220 + t1118 * t1225, -t1149 * t1220 + t1150 * t1225, t1225 * t973 - t1220 * t952 - pkin(5) * (t1089 * t1225 + t1091 * t1220), t1225 * t974 - t1220 * t953 - pkin(5) * (t1090 * t1225 + t1092 * t1220), t1225 * t1017 - t1220 * t1016 - pkin(5) * (t1126 * t1225 + t1128 * t1220), t1225 * t922 - t1220 * t896 - pkin(5) * (t1014 * t1225 + t1015 * t1220), -t1220 * t958 + t1225 * t961, -t1220 * t919 + t1225 * t920, -t1220 * t948 + t1225 * t950, -t1220 * t957 + t1225 * t960, -t1220 * t949 + t1225 * t951, -t1022 * t1220 + t1023 * t1225, t1225 * t839 - t1220 * t830 - pkin(5) * (t1220 * t941 + t1225 * t940), t1225 * t842 - t1220 * t838 - pkin(5) * (t1220 * t946 + t1225 * t945), t1225 * t813 - t1220 * t809 - pkin(5) * (t1220 * t915 + t1225 * t914), t1225 * t806 - t1220 * t804 - pkin(5) * (t1220 * t870 + t1225 * t869), -t1220 * t850 + t1225 * t852, -t1220 * t818 + t1225 * t819, -t1220 * t865 + t1225 * t867, -t1220 * t849 + t1225 * t851, -t1220 * t866 + t1225 * t868, -t1220 * t902 + t1225 * t903, t1225 * t786 - t1220 * t785 - pkin(5) * (t1220 * t848 + t1225 * t847), t1225 * t788 - t1220 * t787 - pkin(5) * (t1220 * t864 + t1225 * t863), t1225 * t780 - t1220 * t779 - pkin(5) * (t1220 * t817 + t1225 * t816), t1225 * t777 - t1220 * t776 - pkin(5) * (t1220 * t793 + t1225 * t792); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1190, 0, t1191, 0, t1170, t1238, t1248, pkin(5) * t1248, 0, 0, t1237, 0, -t1313, 0, -t1325, t1324, t1321, pkin(5) * t1321 + pkin(6) * t1275 + t1225 * t1069, t1109 * t1225 + t1111 * t1220, t1079 * t1225 + t1080 * t1220, t1117 * t1225 + t1119 * t1220, t1108 * t1225 + t1110 * t1220, t1116 * t1225 + t1118 * t1220, t1149 * t1225 + t1150 * t1220, t1220 * t973 + t1225 * t952 + pkin(5) * (-t1089 * t1220 + t1091 * t1225), t1220 * t974 + t1225 * t953 + pkin(5) * (-t1090 * t1220 + t1092 * t1225), t1220 * t1017 + t1225 * t1016 + pkin(5) * (-t1126 * t1220 + t1128 * t1225), t1220 * t922 + t1225 * t896 + pkin(5) * (-t1014 * t1220 + t1015 * t1225), t1220 * t961 + t1225 * t958, t1220 * t920 + t1225 * t919, t1220 * t950 + t1225 * t948, t1220 * t960 + t1225 * t957, t1220 * t951 + t1225 * t949, t1022 * t1225 + t1023 * t1220, t1220 * t839 + t1225 * t830 + pkin(5) * (-t1220 * t940 + t1225 * t941), t1220 * t842 + t1225 * t838 + pkin(5) * (-t1220 * t945 + t1225 * t946), t1220 * t813 + t1225 * t809 + pkin(5) * (-t1220 * t914 + t1225 * t915), t1220 * t806 + t1225 * t804 + pkin(5) * (-t1220 * t869 + t1225 * t870), t1220 * t852 + t1225 * t850, t1220 * t819 + t1225 * t818, t1220 * t867 + t1225 * t865, t1220 * t851 + t1225 * t849, t1220 * t868 + t1225 * t866, t1220 * t903 + t1225 * t902, t1220 * t786 + t1225 * t785 + pkin(5) * (-t1220 * t847 + t1225 * t848), t1220 * t788 + t1225 * t787 + pkin(5) * (-t1220 * t863 + t1225 * t864), t1220 * t780 + t1225 * t779 + pkin(5) * (-t1220 * t816 + t1225 * t817), t1220 * t777 + t1225 * t776 + pkin(5) * (-t1220 * t792 + t1225 * t793); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1198, t1199, 0, 0, 0, 0, 0, 0, 0, t1212, t1097, t1098, 0, -t1070, t1134, t1124, t1139, t1133, t1137, 0, t1019, t1018, t991, t938, t983, t926, t1010, t982, t1011, t1024, t853, t858, t815, t810, t855, t821, t881, t854, t882, t911, t791, t799, t782, t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1227, 0, 0, -g(3), -t1198, 0, 0, 0, -t1182, 0, -t1179, 0, t1314, t1159, t1073, pkin(6) * t1073, t1111, t1080, t1119, t1110, t1118, t1150, t973, t974, t1017, t922, t961, t920, t950, t960, t951, t1023, t839, t842, t813, t806, t852, t819, t867, t851, t868, t903, t786, t788, t780, t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1227, 0, qJDD(1), 0, g(3), 0, -t1199, 0, 0, 0, t1179, 0, -t1182, 0, -t1159, t1314, t1250, t1069, t1109, t1079, t1117, t1108, t1116, t1149, t952, t953, t1016, t896, t958, t919, t948, t957, t949, t1022, t830, t838, t809, t804, t850, t818, t865, t849, t866, t902, t785, t787, t779, t776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1198, t1199, 0, 0, 0, 0, 0, 0, 0, t1212, t1097, t1098, 0, -t1070, t1134, t1124, t1139, t1133, t1137, 0, t1019, t1018, t991, t938, t983, t926, t1010, t982, t1011, t1024, t853, t858, t815, t810, t855, t821, t881, t854, t882, t911, t791, t799, t782, t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1212, 0, -t1305, 0, 0, -g(3), t1131, 0, t1147, t1125, t1143, t1146, t1141, t1169, t1063, t1064, t1030, pkin(7) * t1030, t985, t928, t1012, t984, t1013, t1025, t880, t892, t834, t831, t857, t823, t883, t856, t884, t912, t805, t807, t783, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1305, 0, t1212, 0, g(3), 0, t1132, 0, t1197, -t1184, -t1277, -t1197, -t1270, -qJDD(3), t1059, t1060, 0, pkin(2) * t1030, -t1123, -t1122, -t1058, t1123, t1054, -t1211, t897, t913, t895, t846, -t1041, -t1038, -t968, t1041, t965, -t1205, t814, t826, t808, t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1212, -t1131, -t1132, 0, 0, t1134, t1124, t1139, t1133, t1137, 0, t1255, t1256, t1241, t1265, t983, t926, t1010, t982, t1011, t1024, t1245, t1244, t1246, t1232, t855, t821, t881, t854, t882, t911, t1253, t1252, t1247, t1254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1173, t1174, t1187, -t1201, t1194, t1201, 0, t1120, t1095, 0, t1046, t996, t1067, t1044, t1068, t1086, t988, t999, t894, -pkin(8) * t923, t901, t861, t934, t900, t935, t971, t836, t841, t803, t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1261, t1172, t1192, t1231, t1188, -t1261, -t1120, 0, t1096, 0, t1045, t994, t1065, t1043, t1066, t1085, t943, t947, t890, t905, t899, t859, t932, t898, t933, t970, t828, t833, t801, t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1197, t1184, t1277, t1197, t1270, qJDD(3), -t1095, -t1096, 0, 0, t1123, t1122, t1058, -t1123, -t1054, t1211, t1234, t1240, t1302, t1303, t1041, t1038, t968, -t1041, -t965, t1205, t1229, t1228, t1262, t1263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1078, -t1053, t1308, t1154, t1151, -t1154, 0, t1052, t989, 0, t959, t909, t1008, t955, t1009, t1027, t916, t921, t829, -pkin(9) * t844; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1289, t1307, t1152, -t1230, t1113, -t1289, -t1052, 0, t990, 0, t956, t907, t1006, t954, t1007, t1026, t885, t891, t825, t837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1123, t1122, t1058, -t1123, -t1054, t1211, -t989, -t990, 0, 0, t1041, t1038, t968, -t1041, -t965, t1205, t1239, t1233, t906, t843; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t993, -t964, t1309, t1093, t1087, -t1093, 0, t975, t887, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1292, t1312, t1088, t992, t1034, -t1292, -t975, 0, t888, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1041, t1038, t968, -t1041, -t965, t1205, -t887, -t888, 0, 0;];
m_new_reg = t1;
