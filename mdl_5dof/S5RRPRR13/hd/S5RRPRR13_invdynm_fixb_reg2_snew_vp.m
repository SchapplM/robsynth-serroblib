% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRR13
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRR13_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:34:27
% EndTime: 2019-12-31 20:34:47
% DurationCPUTime: 20.82s
% Computational Cost: add. (181506->746), mult. (397862->1029), div. (0->0), fcn. (284412->10), ass. (0->499)
t1277 = sin(qJ(5));
t1275 = sin(pkin(9));
t1276 = cos(pkin(9));
t1279 = sin(qJ(2));
t1362 = qJD(1) * t1279;
t1235 = -t1276 * qJD(2) + t1275 * t1362;
t1236 = t1275 * qJD(2) + t1276 * t1362;
t1278 = sin(qJ(4));
t1282 = cos(qJ(4));
t1194 = t1282 * t1235 + t1278 * t1236;
t1195 = -t1278 * t1235 + t1282 * t1236;
t1281 = cos(qJ(5));
t1141 = t1281 * t1194 + t1277 * t1195;
t1143 = -t1277 * t1194 + t1281 * t1195;
t1078 = t1143 * t1141;
t1266 = qJD(2) * t1362;
t1283 = cos(qJ(2));
t1329 = t1283 * qJDD(1);
t1243 = -t1266 + t1329;
t1237 = -qJDD(4) + t1243;
t1233 = -qJDD(5) + t1237;
t1379 = -t1078 - t1233;
t1385 = t1277 * t1379;
t1146 = t1195 * t1194;
t1377 = -t1146 - t1237;
t1384 = t1278 * t1377;
t1383 = t1281 * t1379;
t1382 = t1282 * t1377;
t1361 = t1236 * t1235;
t1301 = -t1243 - t1361;
t1381 = t1275 * t1301;
t1380 = t1276 * t1301;
t1337 = t1283 * qJD(1);
t1267 = qJD(2) * t1337;
t1331 = t1279 * qJDD(1);
t1315 = t1267 + t1331;
t1217 = t1275 * qJDD(2) + t1276 * t1315;
t1294 = t1276 * qJDD(2) - t1275 * t1315;
t1123 = -t1194 * qJD(4) + t1282 * t1217 + t1278 * t1294;
t1317 = t1278 * t1217 - t1282 * t1294;
t1302 = t1195 * qJD(4) + t1317;
t1056 = -t1141 * qJD(5) + t1281 * t1123 - t1277 * t1302;
t1263 = -qJD(4) + t1337;
t1255 = -qJD(5) + t1263;
t1126 = t1141 * t1255;
t1378 = t1126 + t1056;
t1182 = t1194 * t1263;
t1098 = -t1182 + t1123;
t1376 = t1182 + t1123;
t1224 = t1236 * t1337;
t1185 = t1294 - t1224;
t1318 = t1277 * t1123 + t1281 * t1302;
t1015 = (qJD(5) + t1255) * t1143 + t1318;
t1094 = (qJD(4) + t1263) * t1195 + t1317;
t1375 = qJD(2) ^ 2;
t1139 = t1141 ^ 2;
t1140 = t1143 ^ 2;
t1374 = t1194 ^ 2;
t1193 = t1195 ^ 2;
t1373 = t1235 ^ 2;
t1234 = t1236 ^ 2;
t1254 = t1255 ^ 2;
t1261 = t1263 ^ 2;
t1280 = sin(qJ(1));
t1284 = cos(qJ(1));
t1252 = t1280 * g(1) - t1284 * g(2);
t1285 = qJD(1) ^ 2;
t1229 = qJDD(1) * pkin(1) + t1285 * pkin(6) + t1252;
t1297 = t1315 + t1267;
t1159 = -t1297 * qJ(3) + (-t1243 + t1266) * pkin(2) - t1229;
t1253 = t1284 * g(1) + t1280 * g(2);
t1230 = -t1285 * pkin(1) + qJDD(1) * pkin(6) - t1253;
t1213 = -t1279 * g(3) + t1283 * t1230;
t1370 = pkin(2) * t1283;
t1307 = -t1279 * qJ(3) - t1370;
t1241 = t1307 * qJD(1);
t1168 = -t1375 * pkin(2) + qJDD(2) * qJ(3) + t1241 * t1337 + t1213;
t1101 = 0.2e1 * qJD(3) * t1236 - t1276 * t1159 + t1275 * t1168;
t1322 = t1235 * t1337;
t1306 = -t1217 + t1322;
t1060 = pkin(3) * t1301 + pkin(7) * t1306 - t1101;
t1102 = -0.2e1 * qJD(3) * t1235 + t1275 * t1159 + t1276 * t1168;
t1218 = -pkin(3) * t1337 - t1236 * pkin(7);
t1066 = -t1373 * pkin(3) + pkin(7) * t1294 + t1218 * t1337 + t1102;
t1000 = -t1282 * t1060 + t1278 * t1066;
t1001 = t1278 * t1060 + t1282 * t1066;
t945 = -t1282 * t1000 + t1278 * t1001;
t1372 = pkin(3) * t945;
t1371 = pkin(2) * t1279;
t1046 = -t1094 * t1278 - t1282 * t1098;
t1369 = pkin(3) * t1046;
t1368 = t1283 * g(3);
t1367 = t1275 * t945;
t1366 = t1276 * t945;
t969 = pkin(4) * t1377 - t1098 * pkin(8) - t1000;
t1175 = -t1263 * pkin(4) - t1195 * pkin(8);
t973 = -t1374 * pkin(4) - pkin(8) * t1302 + t1263 * t1175 + t1001;
t926 = t1277 * t973 - t1281 * t969;
t927 = t1277 * t969 + t1281 * t973;
t892 = t1277 * t927 - t1281 * t926;
t1365 = t1278 * t892;
t1364 = t1282 * t892;
t1363 = qJD(1) * qJD(2);
t1360 = t1255 * t1143;
t1359 = t1255 * t1277;
t1358 = t1255 * t1281;
t1357 = t1263 * t1195;
t1356 = t1263 * t1278;
t1355 = t1263 * t1282;
t1272 = t1279 ^ 2;
t1354 = t1272 * t1285;
t1167 = t1368 + qJDD(3) - t1375 * qJ(3) - qJDD(2) * pkin(2) + (qJD(1) * t1241 + t1230) * t1279;
t1353 = t1275 * t1167;
t1188 = t1243 - t1361;
t1352 = t1275 * t1188;
t1351 = t1276 * t1167;
t1350 = t1276 * t1188;
t1105 = -t1294 * pkin(3) - t1373 * pkin(7) + t1236 * t1218 + t1167;
t1025 = pkin(4) * t1302 - t1374 * pkin(8) + t1195 * t1175 + t1105;
t1349 = t1277 * t1025;
t1071 = -t1078 + t1233;
t1348 = t1277 * t1071;
t1347 = t1278 * t1105;
t1128 = -t1146 + t1237;
t1346 = t1278 * t1128;
t1345 = t1279 * t1229;
t1344 = t1279 * t1243;
t1262 = t1283 * t1285 * t1279;
t1250 = -t1262 + qJDD(2);
t1343 = t1279 * t1250;
t1251 = qJDD(2) + t1262;
t1342 = t1279 * t1251;
t1341 = t1281 * t1025;
t1340 = t1281 * t1071;
t1339 = t1282 * t1105;
t1338 = t1282 * t1128;
t1336 = t1283 * t1229;
t1335 = t1283 * t1250;
t1273 = t1283 ^ 2;
t1332 = t1272 + t1273;
t1330 = t1280 * qJDD(1);
t1328 = t1284 * qJDD(1);
t893 = t1277 * t926 + t1281 * t927;
t864 = t1278 * t893 + t1364;
t891 = pkin(4) * t892;
t1327 = pkin(3) * t864 + t891;
t1018 = -t1126 + t1056;
t958 = -t1015 * t1277 - t1281 * t1018;
t960 = -t1015 * t1281 + t1277 * t1018;
t920 = t1278 * t960 + t1282 * t958;
t956 = pkin(4) * t958;
t1326 = pkin(3) * t920 + t956;
t1325 = t1279 * t1078;
t1324 = t1279 * t1146;
t1323 = t1279 * t1361;
t1321 = t1283 * t1078;
t1320 = t1283 * t1146;
t1319 = t1283 * t1361;
t946 = t1278 * t1000 + t1282 * t1001;
t1050 = t1275 * t1101 + t1276 * t1102;
t1212 = t1279 * t1230 + t1368;
t1162 = t1279 * t1212 + t1283 * t1213;
t1316 = -t1280 * t1252 - t1284 * t1253;
t1314 = t1275 * t1322;
t1313 = t1280 * t1262;
t1312 = t1284 * t1262;
t1164 = -t1193 - t1261;
t1079 = t1282 * t1164 + t1346;
t1311 = pkin(3) * t1079 - t1001;
t1247 = -t1280 * t1285 + t1328;
t1310 = -pkin(5) * t1247 - t1280 * g(3);
t1076 = -t1254 - t1139;
t1023 = t1277 * t1076 + t1383;
t1309 = pkin(4) * t1023 - t926;
t1308 = -pkin(2) * t1167 + qJ(3) * t1050;
t1049 = -t1276 * t1101 + t1275 * t1102;
t1161 = t1283 * t1212 - t1279 * t1213;
t1305 = t1284 * t1252 - t1280 * t1253;
t1138 = -t1261 - t1374;
t1074 = t1278 * t1138 + t1382;
t1304 = pkin(3) * t1074 - t1000;
t1114 = -t1140 - t1254;
t1042 = t1281 * t1114 + t1348;
t1303 = pkin(4) * t1042 - t927;
t1024 = t1281 * t1076 - t1385;
t961 = t1282 * t1023 + t1278 * t1024;
t1300 = pkin(3) * t961 + t1309;
t1269 = t1273 * t1285;
t1196 = -t1269 - t1373;
t1136 = t1276 * t1196 - t1381;
t1184 = t1224 + t1294;
t1299 = pkin(2) * t1184 + qJ(3) * t1136 - t1351;
t1222 = -t1234 - t1269;
t1150 = -t1275 * t1222 + t1350;
t1187 = t1217 + t1322;
t1298 = -pkin(2) * t1187 + qJ(3) * t1150 + t1353;
t1043 = -t1277 * t1114 + t1340;
t976 = t1282 * t1042 + t1278 * t1043;
t1296 = pkin(3) * t976 + t1303;
t1122 = t1276 * t1185 - t1275 * t1306;
t1180 = t1234 + t1373;
t1295 = pkin(2) * t1180 + qJ(3) * t1122 + t1050;
t1014 = (qJD(5) - t1255) * t1143 + t1318;
t938 = -pkin(4) * t1014 + pkin(8) * t1024 - t1341;
t962 = -t1278 * t1023 + t1282 * t1024;
t964 = -pkin(8) * t1023 + t1349;
t886 = -pkin(3) * t1014 + pkin(7) * t962 + t1278 * t964 + t1282 * t938;
t894 = -pkin(7) * t961 - t1278 * t938 + t1282 * t964;
t924 = -t1275 * t961 + t1276 * t962;
t1293 = -pkin(2) * t1014 + qJ(3) * t924 + t1275 * t894 + t1276 * t886;
t941 = -pkin(4) * t1378 + pkin(8) * t1043 + t1349;
t972 = -pkin(8) * t1042 + t1341;
t977 = -t1278 * t1042 + t1282 * t1043;
t890 = -pkin(3) * t1378 + pkin(7) * t977 + t1278 * t972 + t1282 * t941;
t896 = -pkin(7) * t976 - t1278 * t941 + t1282 * t972;
t931 = -t1275 * t976 + t1276 * t977;
t1292 = -pkin(2) * t1378 + qJ(3) * t931 + t1275 * t896 + t1276 * t890;
t865 = t1282 * t893 - t1365;
t885 = -pkin(4) * t1025 + pkin(8) * t893;
t848 = -pkin(3) * t1025 + pkin(7) * t865 - pkin(8) * t1365 + t1282 * t885;
t852 = -pkin(7) * t864 - pkin(8) * t1364 - t1278 * t885;
t854 = -t1275 * t864 + t1276 * t865;
t1291 = -pkin(2) * t1025 + qJ(3) * t854 + t1275 * t852 + t1276 * t848;
t1057 = -t1139 - t1140;
t872 = -pkin(4) * t1057 + pkin(8) * t960 + t893;
t873 = -pkin(8) * t958 - t892;
t922 = -t1278 * t958 + t1282 * t960;
t857 = -pkin(3) * t1057 + pkin(7) * t922 + t1278 * t873 + t1282 * t872;
t858 = -pkin(7) * t920 - t1278 * t872 + t1282 * t873;
t884 = -t1275 * t920 + t1276 * t922;
t1290 = -pkin(2) * t1057 + qJ(3) * t884 + t1275 * t858 + t1276 * t857;
t1113 = -t1193 - t1374;
t1048 = -t1094 * t1282 + t1278 * t1098;
t928 = -pkin(3) * t1113 + pkin(7) * t1048 + t946;
t932 = -pkin(7) * t1046 - t945;
t981 = -t1275 * t1046 + t1276 * t1048;
t1289 = -pkin(2) * t1113 + qJ(3) * t981 + t1275 * t932 + t1276 * t928;
t1075 = t1282 * t1138 - t1384;
t1021 = -t1275 * t1074 + t1276 * t1075;
t1029 = -pkin(7) * t1074 + t1347;
t1093 = (qJD(4) - t1263) * t1195 + t1317;
t992 = -pkin(3) * t1093 + pkin(7) * t1075 - t1339;
t1288 = -pkin(2) * t1093 + qJ(3) * t1021 + t1275 * t1029 + t1276 * t992;
t1080 = -t1278 * t1164 + t1338;
t1002 = -pkin(3) * t1376 + pkin(7) * t1080 + t1347;
t1031 = -t1275 * t1079 + t1276 * t1080;
t1044 = -pkin(7) * t1079 + t1339;
t1287 = -pkin(2) * t1376 + qJ(3) * t1031 + t1276 * t1002 + t1275 * t1044;
t903 = t1276 * t946 - t1367;
t940 = -pkin(3) * t1105 + pkin(7) * t946;
t1286 = -pkin(2) * t1105 - pkin(7) * t1367 + qJ(3) * t903 + t1276 * t940;
t1260 = -t1269 - t1375;
t1259 = t1269 - t1375;
t1258 = -t1354 - t1375;
t1257 = -t1354 + t1375;
t1249 = -t1269 + t1354;
t1248 = t1269 + t1354;
t1246 = t1284 * t1285 + t1330;
t1245 = t1332 * qJDD(1);
t1244 = -0.2e1 * t1266 + t1329;
t1242 = 0.2e1 * t1267 + t1331;
t1239 = t1283 * t1251;
t1238 = t1332 * t1363;
t1231 = t1283 * t1243;
t1225 = -pkin(5) * t1246 + t1284 * g(3);
t1221 = -t1234 + t1269;
t1220 = -t1269 + t1373;
t1219 = t1276 * t1224;
t1216 = -t1272 * t1363 + t1283 * t1315;
t1215 = -t1273 * t1363 - t1344;
t1211 = -t1279 * t1258 - t1335;
t1210 = -t1279 * t1257 + t1239;
t1209 = t1283 * t1260 - t1342;
t1208 = t1283 * t1259 - t1343;
t1207 = t1283 * t1258 - t1343;
t1206 = t1283 * t1257 + t1342;
t1205 = t1279 * t1260 + t1239;
t1204 = t1279 * t1259 + t1335;
t1203 = t1297 * t1279;
t1202 = -t1279 * t1267 + t1231;
t1199 = t1234 - t1373;
t1198 = -t1279 * t1242 + t1283 * t1244;
t1197 = t1283 * t1242 + t1279 * t1244;
t1179 = (t1235 * t1276 - t1236 * t1275) * t1337;
t1178 = -t1219 - t1314;
t1177 = -t1193 + t1261;
t1176 = -t1261 + t1374;
t1174 = -pkin(6) * t1207 - t1336;
t1173 = -pkin(6) * t1205 - t1345;
t1172 = t1276 * t1217 + t1275 * t1224;
t1171 = t1275 * t1217 - t1219;
t1170 = -t1275 * t1294 - t1276 * t1322;
t1169 = t1276 * t1294 - t1314;
t1166 = -pkin(1) * t1207 + t1213;
t1165 = -pkin(1) * t1205 + t1212;
t1158 = pkin(1) * t1244 + pkin(6) * t1209 + t1336;
t1157 = -pkin(1) * t1242 + pkin(6) * t1211 - t1345;
t1154 = t1283 * t1179 - t1344;
t1153 = t1279 * t1179 + t1231;
t1152 = t1276 * t1220 + t1352;
t1151 = -t1275 * t1221 + t1380;
t1149 = t1275 * t1220 - t1350;
t1148 = t1276 * t1221 + t1381;
t1147 = t1276 * t1222 + t1352;
t1145 = t1193 - t1374;
t1144 = pkin(1) * t1229 + pkin(6) * t1162;
t1137 = pkin(1) * t1248 + pkin(6) * t1245 + t1162;
t1135 = t1275 * t1196 + t1380;
t1134 = t1283 * t1172 + t1323;
t1133 = t1283 * t1170 - t1323;
t1132 = t1279 * t1172 - t1319;
t1131 = t1279 * t1170 + t1319;
t1125 = -t1140 + t1254;
t1124 = t1139 - t1254;
t1121 = t1276 * t1184 - t1275 * t1187;
t1120 = t1275 * t1185 + t1276 * t1306;
t1119 = t1275 * t1184 + t1276 * t1187;
t1117 = (t1194 * t1282 - t1195 * t1278) * t1263;
t1116 = (t1194 * t1278 + t1195 * t1282) * t1263;
t1112 = t1283 * t1152 + t1185 * t1279;
t1111 = t1283 * t1151 - t1279 * t1306;
t1110 = t1283 * t1150 + t1279 * t1187;
t1109 = t1279 * t1152 - t1185 * t1283;
t1108 = t1279 * t1151 + t1283 * t1306;
t1107 = t1279 * t1150 - t1283 * t1187;
t1106 = -qJ(3) * t1147 + t1351;
t1104 = t1283 * t1121 + t1279 * t1199;
t1103 = t1279 * t1121 - t1283 * t1199;
t1100 = t1283 * t1136 - t1279 * t1184;
t1099 = t1279 * t1136 + t1283 * t1184;
t1091 = t1282 * t1176 + t1346;
t1090 = -t1278 * t1177 + t1382;
t1089 = t1278 * t1176 - t1338;
t1088 = t1282 * t1177 + t1384;
t1087 = t1282 * t1123 + t1195 * t1356;
t1086 = t1278 * t1123 - t1195 * t1355;
t1085 = -t1194 * t1355 + t1278 * t1302;
t1084 = -t1194 * t1356 - t1282 * t1302;
t1083 = -qJ(3) * t1135 + t1353;
t1082 = t1283 * t1122 - t1279 * t1180;
t1081 = t1279 * t1122 + t1283 * t1180;
t1077 = t1140 - t1139;
t1070 = (t1141 * t1281 - t1143 * t1277) * t1255;
t1069 = (t1141 * t1277 + t1143 * t1281) * t1255;
t1068 = -t1275 * t1116 + t1276 * t1117;
t1067 = t1276 * t1116 + t1275 * t1117;
t1064 = -pkin(2) * t1147 + t1102;
t1063 = t1283 * t1068 - t1279 * t1237;
t1062 = t1279 * t1068 + t1283 * t1237;
t1061 = -pkin(2) * t1135 + t1101;
t1055 = -t1143 * qJD(5) - t1318;
t1054 = t1281 * t1124 + t1348;
t1053 = -t1277 * t1125 + t1383;
t1052 = t1277 * t1124 - t1340;
t1051 = t1281 * t1125 + t1385;
t1047 = -t1282 * t1093 - t1278 * t1376;
t1045 = -t1278 * t1093 + t1282 * t1376;
t1040 = -t1275 * t1089 + t1276 * t1091;
t1039 = -t1275 * t1088 + t1276 * t1090;
t1038 = t1276 * t1089 + t1275 * t1091;
t1037 = t1276 * t1088 + t1275 * t1090;
t1036 = -t1275 * t1086 + t1276 * t1087;
t1035 = -t1275 * t1084 + t1276 * t1085;
t1034 = t1276 * t1086 + t1275 * t1087;
t1033 = t1276 * t1084 + t1275 * t1085;
t1032 = -pkin(1) * t1107 - t1298;
t1030 = t1276 * t1079 + t1275 * t1080;
t1028 = -pkin(1) * t1099 - t1299;
t1027 = t1283 * t1050 + t1279 * t1167;
t1026 = t1279 * t1050 - t1283 * t1167;
t1020 = t1276 * t1074 + t1275 * t1075;
t1013 = t1281 * t1056 + t1143 * t1359;
t1012 = t1277 * t1056 - t1143 * t1358;
t1011 = -t1277 * t1055 - t1141 * t1358;
t1010 = t1281 * t1055 - t1141 * t1359;
t1009 = -qJ(3) * t1120 - t1049;
t1008 = t1283 * t1036 + t1324;
t1007 = t1283 * t1035 - t1324;
t1006 = t1279 * t1036 - t1320;
t1005 = t1279 * t1035 + t1320;
t1004 = -t1278 * t1069 + t1282 * t1070;
t1003 = t1282 * t1069 + t1278 * t1070;
t999 = t1283 * t1040 - t1279 * t1094;
t998 = t1283 * t1039 + t1279 * t1098;
t997 = t1279 * t1040 + t1283 * t1094;
t996 = t1279 * t1039 - t1283 * t1098;
t994 = t1283 * t1031 + t1279 * t1376;
t993 = t1279 * t1031 - t1283 * t1376;
t991 = -pkin(6) * t1107 - t1279 * t1064 + t1283 * t1106;
t990 = t1283 * t1021 + t1279 * t1093;
t989 = t1279 * t1021 - t1283 * t1093;
t988 = -pkin(6) * t1099 - t1279 * t1061 + t1283 * t1083;
t987 = -pkin(1) * t1147 + pkin(6) * t1110 + t1283 * t1064 + t1279 * t1106;
t986 = -t1278 * t1052 + t1282 * t1054;
t985 = -t1278 * t1051 + t1282 * t1053;
t984 = t1282 * t1052 + t1278 * t1054;
t983 = t1282 * t1051 + t1278 * t1053;
t982 = -pkin(1) * t1081 - t1295;
t980 = -t1275 * t1045 + t1276 * t1047;
t979 = t1276 * t1046 + t1275 * t1048;
t978 = t1276 * t1045 + t1275 * t1047;
t975 = -pkin(1) * t1135 + pkin(6) * t1100 + t1283 * t1061 + t1279 * t1083;
t974 = -pkin(6) * t1081 + t1283 * t1009 + t1120 * t1371;
t971 = t1279 * t1145 + t1283 * t980;
t970 = -t1283 * t1145 + t1279 * t980;
t967 = -pkin(1) * t1026 - t1308;
t966 = t1279 * t1113 + t1283 * t981;
t965 = -t1283 * t1113 + t1279 * t981;
t963 = pkin(6) * t1082 + t1279 * t1009 + (-pkin(1) - t1370) * t1120;
t959 = -t1281 * t1014 - t1277 * t1378;
t957 = -t1277 * t1014 + t1281 * t1378;
t955 = -t1278 * t1012 + t1282 * t1013;
t954 = -t1278 * t1010 + t1282 * t1011;
t953 = t1282 * t1012 + t1278 * t1013;
t952 = t1282 * t1010 + t1278 * t1011;
t951 = -t1275 * t1003 + t1276 * t1004;
t950 = t1276 * t1003 + t1275 * t1004;
t949 = -t1279 * t1233 + t1283 * t951;
t948 = t1283 * t1233 + t1279 * t951;
t947 = -pkin(2) * t979 - t1369;
t944 = -pkin(2) * t1030 - t1311;
t943 = -pkin(6) * t1026 + (-qJ(3) * t1283 + t1371) * t1049;
t942 = -pkin(2) * t1020 - t1304;
t939 = -qJ(3) * t1030 - t1275 * t1002 + t1276 * t1044;
t937 = -t1275 * t984 + t1276 * t986;
t936 = -t1275 * t983 + t1276 * t985;
t935 = t1275 * t986 + t1276 * t984;
t934 = t1275 * t985 + t1276 * t983;
t933 = -qJ(3) * t1020 + t1276 * t1029 - t1275 * t992;
t930 = t1275 * t977 + t1276 * t976;
t929 = pkin(6) * t1027 + (-pkin(1) + t1307) * t1049;
t923 = t1275 * t962 + t1276 * t961;
t921 = -t1278 * t957 + t1282 * t959;
t919 = t1278 * t959 + t1282 * t957;
t918 = -t1279 * t1015 + t1283 * t937;
t917 = t1279 * t1018 + t1283 * t936;
t916 = t1283 * t1015 + t1279 * t937;
t915 = -t1283 * t1018 + t1279 * t936;
t914 = -t1275 * t953 + t1276 * t955;
t913 = -t1275 * t952 + t1276 * t954;
t912 = t1275 * t955 + t1276 * t953;
t911 = t1275 * t954 + t1276 * t952;
t910 = t1279 * t1378 + t1283 * t931;
t909 = t1279 * t931 - t1283 * t1378;
t908 = -pkin(1) * t993 - t1287;
t907 = t1283 * t914 + t1325;
t906 = t1283 * t913 - t1325;
t905 = t1279 * t914 - t1321;
t904 = t1279 * t913 + t1321;
t902 = t1275 * t946 + t1366;
t901 = -pkin(1) * t989 - t1288;
t900 = t1279 * t1014 + t1283 * t924;
t899 = -t1283 * t1014 + t1279 * t924;
t898 = t1279 * t1105 + t1283 * t903;
t897 = -t1283 * t1105 + t1279 * t903;
t895 = -pkin(6) * t993 - t1279 * t944 + t1283 * t939;
t889 = -pkin(6) * t989 - t1279 * t942 + t1283 * t933;
t888 = -pkin(2) * t902 - t1372;
t887 = -pkin(1) * t1030 + pkin(6) * t994 + t1279 * t939 + t1283 * t944;
t883 = -t1275 * t919 + t1276 * t921;
t882 = t1275 * t922 + t1276 * t920;
t881 = t1275 * t921 + t1276 * t919;
t880 = -pkin(1) * t1020 + pkin(6) * t990 + t1279 * t933 + t1283 * t942;
t879 = -qJ(3) * t979 - t1275 * t928 + t1276 * t932;
t878 = t1279 * t1077 + t1283 * t883;
t877 = -t1283 * t1077 + t1279 * t883;
t876 = t1279 * t1057 + t1283 * t884;
t875 = -t1283 * t1057 + t1279 * t884;
t874 = -pkin(2) * t930 - t1296;
t871 = -pkin(7) * t1366 - qJ(3) * t902 - t1275 * t940;
t870 = -pkin(1) * t965 - t1289;
t869 = -pkin(2) * t923 - t1300;
t868 = -pkin(6) * t965 - t1279 * t947 + t1283 * t879;
t867 = -pkin(1) * t979 + pkin(6) * t966 + t1279 * t879 + t1283 * t947;
t866 = -pkin(2) * t882 - t1326;
t863 = -pkin(1) * t897 - t1286;
t862 = -qJ(3) * t930 - t1275 * t890 + t1276 * t896;
t861 = -qJ(3) * t923 - t1275 * t886 + t1276 * t894;
t860 = -pkin(1) * t909 - t1292;
t859 = -pkin(6) * t897 - t1279 * t888 + t1283 * t871;
t856 = -pkin(1) * t899 - t1293;
t855 = -pkin(1) * t902 + pkin(6) * t898 + t1279 * t871 + t1283 * t888;
t853 = t1275 * t865 + t1276 * t864;
t851 = t1279 * t1025 + t1283 * t854;
t850 = -t1283 * t1025 + t1279 * t854;
t849 = -pkin(6) * t909 - t1279 * t874 + t1283 * t862;
t847 = -pkin(1) * t930 + pkin(6) * t910 + t1279 * t862 + t1283 * t874;
t846 = -pkin(6) * t899 - t1279 * t869 + t1283 * t861;
t845 = -pkin(1) * t923 + pkin(6) * t900 + t1279 * t861 + t1283 * t869;
t844 = -pkin(2) * t853 - t1327;
t843 = -qJ(3) * t882 - t1275 * t857 + t1276 * t858;
t842 = -pkin(1) * t875 - t1290;
t841 = -pkin(6) * t875 - t1279 * t866 + t1283 * t843;
t840 = -pkin(1) * t882 + pkin(6) * t876 + t1279 * t843 + t1283 * t866;
t839 = -qJ(3) * t853 - t1275 * t848 + t1276 * t852;
t838 = -pkin(1) * t850 - t1291;
t837 = -pkin(6) * t850 - t1279 * t844 + t1283 * t839;
t836 = -pkin(1) * t853 + pkin(6) * t851 + t1279 * t839 + t1283 * t844;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1247, 0, -t1246, 0, t1310, -t1225, -t1305, -pkin(5) * t1305, t1284 * t1216 - t1313, t1284 * t1198 + t1280 * t1249, t1284 * t1210 + t1279 * t1330, t1284 * t1215 + t1313, t1284 * t1208 + t1280 * t1329, t1280 * qJDD(2) + t1284 * t1238, t1284 * t1173 - t1280 * t1165 - pkin(5) * (t1280 * t1209 + t1284 * t1244), t1284 * t1174 - t1280 * t1166 - pkin(5) * (t1280 * t1211 - t1284 * t1242), t1284 * t1161 - pkin(5) * (t1280 * t1245 + t1284 * t1248), -pkin(5) * (t1280 * t1162 + t1284 * t1229) - (t1280 * pkin(1) - t1284 * pkin(6)) * t1161, t1284 * t1134 + t1280 * t1171, t1284 * t1104 + t1280 * t1119, t1284 * t1111 + t1280 * t1148, t1284 * t1133 + t1280 * t1169, t1284 * t1112 + t1280 * t1149, t1284 * t1154 - t1280 * t1178, t1284 * t988 - t1280 * t1028 - pkin(5) * (t1280 * t1100 - t1284 * t1135), t1284 * t991 - t1280 * t1032 - pkin(5) * (t1280 * t1110 - t1284 * t1147), t1284 * t974 - t1280 * t982 - pkin(5) * (t1280 * t1082 - t1284 * t1120), t1284 * t943 - t1280 * t967 - pkin(5) * (t1280 * t1027 - t1284 * t1049), t1284 * t1008 + t1280 * t1034, t1280 * t978 + t1284 * t971, t1280 * t1037 + t1284 * t998, t1284 * t1007 + t1280 * t1033, t1280 * t1038 + t1284 * t999, t1284 * t1063 + t1280 * t1067, t1284 * t889 - t1280 * t901 - pkin(5) * (-t1284 * t1020 + t1280 * t990), t1284 * t895 - t1280 * t908 - pkin(5) * (-t1284 * t1030 + t1280 * t994), t1284 * t868 - t1280 * t870 - pkin(5) * (t1280 * t966 - t1284 * t979), t1284 * t859 - t1280 * t863 - pkin(5) * (t1280 * t898 - t1284 * t902), t1280 * t912 + t1284 * t907, t1280 * t881 + t1284 * t878, t1280 * t934 + t1284 * t917, t1280 * t911 + t1284 * t906, t1280 * t935 + t1284 * t918, t1280 * t950 + t1284 * t949, t1284 * t846 - t1280 * t856 - pkin(5) * (t1280 * t900 - t1284 * t923), t1284 * t849 - t1280 * t860 - pkin(5) * (t1280 * t910 - t1284 * t930), t1284 * t841 - t1280 * t842 - pkin(5) * (t1280 * t876 - t1284 * t882), t1284 * t837 - t1280 * t838 - pkin(5) * (t1280 * t851 - t1284 * t853); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1246, 0, t1247, 0, t1225, t1310, t1316, pkin(5) * t1316, t1280 * t1216 + t1312, t1280 * t1198 - t1284 * t1249, t1280 * t1210 - t1279 * t1328, t1280 * t1215 - t1312, t1280 * t1208 - t1283 * t1328, -t1284 * qJDD(2) + t1280 * t1238, t1280 * t1173 + t1284 * t1165 + pkin(5) * (t1284 * t1209 - t1280 * t1244), t1280 * t1174 + t1284 * t1166 + pkin(5) * (t1284 * t1211 + t1280 * t1242), t1280 * t1161 + pkin(5) * (t1284 * t1245 - t1280 * t1248), pkin(5) * (t1284 * t1162 - t1280 * t1229) - (-t1284 * pkin(1) - t1280 * pkin(6)) * t1161, t1280 * t1134 - t1284 * t1171, t1280 * t1104 - t1284 * t1119, t1280 * t1111 - t1284 * t1148, t1280 * t1133 - t1284 * t1169, t1280 * t1112 - t1284 * t1149, t1280 * t1154 + t1284 * t1178, t1280 * t988 + t1284 * t1028 + pkin(5) * (t1284 * t1100 + t1280 * t1135), t1280 * t991 + t1284 * t1032 + pkin(5) * (t1284 * t1110 + t1280 * t1147), t1280 * t974 + t1284 * t982 + pkin(5) * (t1284 * t1082 + t1280 * t1120), t1280 * t943 + t1284 * t967 + pkin(5) * (t1284 * t1027 + t1280 * t1049), t1280 * t1008 - t1284 * t1034, t1280 * t971 - t1284 * t978, -t1284 * t1037 + t1280 * t998, t1280 * t1007 - t1284 * t1033, -t1284 * t1038 + t1280 * t999, t1280 * t1063 - t1284 * t1067, t1280 * t889 + t1284 * t901 + pkin(5) * (t1280 * t1020 + t1284 * t990), t1280 * t895 + t1284 * t908 + pkin(5) * (t1280 * t1030 + t1284 * t994), t1280 * t868 + t1284 * t870 + pkin(5) * (t1280 * t979 + t1284 * t966), t1280 * t859 + t1284 * t863 + pkin(5) * (t1280 * t902 + t1284 * t898), t1280 * t907 - t1284 * t912, t1280 * t878 - t1284 * t881, t1280 * t917 - t1284 * t934, t1280 * t906 - t1284 * t911, t1280 * t918 - t1284 * t935, t1280 * t949 - t1284 * t950, t1280 * t846 + t1284 * t856 + pkin(5) * (t1280 * t923 + t1284 * t900), t1280 * t849 + t1284 * t860 + pkin(5) * (t1280 * t930 + t1284 * t910), t1280 * t841 + t1284 * t842 + pkin(5) * (t1280 * t882 + t1284 * t876), t1280 * t837 + t1284 * t838 + pkin(5) * (t1280 * t853 + t1284 * t851); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1252, t1253, 0, 0, t1203, t1197, t1206, t1202, t1204, 0, t1158, t1157, t1137, t1144, t1132, t1103, t1108, t1131, t1109, t1153, t975, t987, t963, t929, t1006, t970, t996, t1005, t997, t1062, t880, t887, t867, t855, t905, t877, t915, t904, t916, t948, t845, t847, t840, t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1285, 0, 0, -g(3), -t1252, 0, t1216, t1198, t1210, t1215, t1208, t1238, t1173, t1174, t1161, pkin(6) * t1161, t1134, t1104, t1111, t1133, t1112, t1154, t988, t991, t974, t943, t1008, t971, t998, t1007, t999, t1063, t889, t895, t868, t859, t907, t878, t917, t906, t918, t949, t846, t849, t841, t837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1285, 0, qJDD(1), 0, g(3), 0, -t1253, 0, t1262, -t1249, -t1331, -t1262, -t1329, -qJDD(2), t1165, t1166, 0, pkin(1) * t1161, -t1171, -t1119, -t1148, -t1169, -t1149, t1178, t1028, t1032, t982, t967, -t1034, -t978, -t1037, -t1033, -t1038, -t1067, t901, t908, t870, t863, -t912, -t881, -t934, -t911, -t935, -t950, t856, t860, t842, t838; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1252, t1253, 0, 0, t1203, t1197, t1206, t1202, t1204, 0, t1158, t1157, t1137, t1144, t1132, t1103, t1108, t1131, t1109, t1153, t975, t987, t963, t929, t1006, t970, t996, t1005, t997, t1062, t880, t887, t867, t855, t905, t877, t915, t904, t916, t948, t845, t847, t840, t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1315, t1244, t1251, -t1267, t1259, t1267, 0, -t1229, t1212, 0, t1172, t1121, t1151, t1170, t1152, t1179, t1083, t1106, t1009, -qJ(3) * t1049, t1036, t980, t1039, t1035, t1040, t1068, t933, t939, t879, t871, t914, t883, t936, t913, t937, t951, t861, t862, t843, t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1266, t1242, t1257, t1243, t1250, -t1266, t1229, 0, t1213, 0, -t1361, -t1199, t1306, t1361, -t1185, t1243, t1061, t1064, -pkin(2) * t1120, -pkin(2) * t1049, -t1146, -t1145, -t1098, t1146, t1094, t1237, t942, t944, t947, t888, -t1078, -t1077, -t1018, t1078, t1015, t1233, t869, t874, t866, t844; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1262, t1249, t1331, t1262, t1329, qJDD(2), -t1212, -t1213, 0, 0, t1171, t1119, t1148, t1169, t1149, -t1178, t1299, t1298, t1295, t1308, t1034, t978, t1037, t1033, t1038, t1067, t1288, t1287, t1289, t1286, t912, t881, t934, t911, t935, t950, t1293, t1292, t1290, t1291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1217, t1184, t1301, -t1322, t1220, t1322, 0, t1167, t1101, 0, t1087, t1047, t1090, t1085, t1091, t1117, t1029, t1044, t932, -pkin(7) * t945, t955, t921, t985, t954, t986, t1004, t894, t896, t858, t852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1224, t1187, t1221, t1294, -t1188, t1224, -t1167, 0, t1102, 0, t1086, t1045, t1088, t1084, t1089, t1116, t992, t1002, t928, t940, t953, t919, t983, t952, t984, t1003, t886, t890, t857, t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1361, t1199, -t1306, -t1361, t1185, -t1243, -t1101, -t1102, 0, 0, t1146, t1145, t1098, -t1146, -t1094, -t1237, t1304, t1311, t1369, t1372, t1078, t1077, t1018, -t1078, -t1015, -t1233, t1300, t1296, t1326, t1327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1123, -t1093, t1377, -t1182, t1176, t1182, 0, t1105, t1000, 0, t1013, t959, t1053, t1011, t1054, t1070, t964, t972, t873, -pkin(8) * t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1357, t1376, t1177, -t1302, -t1128, t1357, -t1105, 0, t1001, 0, t1012, t957, t1051, t1010, t1052, t1069, t938, t941, t872, t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1146, t1145, t1098, -t1146, -t1094, -t1237, -t1000, -t1001, 0, 0, t1078, t1077, t1018, -t1078, -t1015, -t1233, t1309, t1303, t956, t891; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1056, -t1014, t1379, -t1126, t1124, t1126, 0, t1025, t926, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1360, t1378, t1125, t1055, -t1071, t1360, -t1025, 0, t927, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1078, t1077, t1018, -t1078, -t1015, -t1233, -t926, -t927, 0, 0;];
m_new_reg = t1;