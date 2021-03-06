% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPPR6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:39
% EndTime: 2019-12-31 19:33:57
% DurationCPUTime: 19.00s
% Computational Cost: add. (137763->734), mult. (328238->1029), div. (0->0), fcn. (231563->10), ass. (0->505)
t1303 = sin(pkin(9));
t1304 = sin(pkin(8));
t1306 = cos(pkin(8));
t1311 = cos(qJ(2));
t1396 = qJD(1) * t1311;
t1308 = sin(qJ(2));
t1397 = qJD(1) * t1308;
t1266 = t1304 * t1396 + t1306 * t1397;
t1305 = cos(pkin(9));
t1234 = -t1305 * qJD(2) + t1303 * t1266;
t1235 = t1303 * qJD(2) + t1305 * t1266;
t1190 = t1235 * t1234;
t1294 = qJD(2) * t1396;
t1361 = t1308 * qJDD(1);
t1274 = t1294 + t1361;
t1350 = qJD(2) * t1397;
t1359 = t1311 * qJDD(1);
t1322 = t1350 - t1359;
t1224 = t1304 * t1274 + t1306 * t1322;
t1413 = -t1190 + t1224;
t1421 = t1303 * t1413;
t1264 = t1304 * t1397 - t1306 * t1396;
t1223 = t1266 * t1264;
t1411 = qJDD(2) - t1223;
t1420 = t1304 * t1411;
t1419 = t1305 * t1413;
t1418 = t1306 * t1411;
t1307 = sin(qJ(5));
t1310 = cos(qJ(5));
t1177 = t1310 * t1234 + t1307 * t1235;
t1179 = -t1307 * t1234 + t1310 * t1235;
t1128 = t1179 * t1177;
t1221 = qJDD(5) + t1224;
t1415 = -t1128 + t1221;
t1417 = t1307 * t1415;
t1416 = t1310 * t1415;
t1394 = qJD(3) * t1264;
t1253 = -0.2e1 * t1394;
t1309 = sin(qJ(1));
t1312 = cos(qJ(1));
t1284 = t1312 * g(1) + t1309 * g(2);
t1313 = qJD(1) ^ 2;
t1317 = -t1313 * pkin(1) + qJDD(1) * pkin(6) - t1284;
t1247 = -t1308 * g(3) + t1311 * t1317;
t1301 = t1311 ^ 2;
t1297 = t1301 * t1313;
t1410 = qJD(2) ^ 2;
t1291 = -t1297 - t1410;
t1185 = t1291 * pkin(2) + qJ(3) * t1359 + t1247;
t1316 = t1308 * t1317;
t1370 = t1308 * t1313;
t1398 = qJD(1) * qJD(2);
t1314 = -t1316 - t1274 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t1370 + qJ(3) * t1398 - g(3)) * t1311;
t1363 = t1306 * t1185 + t1304 * t1314;
t1107 = t1253 + t1363;
t1214 = t1264 * pkin(3) - t1266 * qJ(4);
t1073 = -t1410 * pkin(3) + qJDD(2) * qJ(4) - t1264 * t1214 + t1107;
t1283 = t1309 * g(1) - t1312 * g(2);
t1328 = qJDD(1) * pkin(1) + t1283;
t1191 = (t1301 * qJ(3) + pkin(6)) * t1313 - pkin(2) * t1322 - qJDD(3) - (qJD(2) * pkin(2) - qJ(3) * t1397) * t1397 + t1328;
t1395 = qJD(2) * t1266;
t1193 = t1224 + t1395;
t1225 = t1306 * t1274 - t1304 * t1322;
t1257 = qJD(2) * t1264;
t1196 = -t1257 + t1225;
t1092 = t1193 * pkin(3) - t1196 * qJ(4) - t1191;
t1010 = 0.2e1 * qJD(4) * t1235 + t1303 * t1073 - t1305 * t1092;
t1011 = -0.2e1 * qJD(4) * t1234 + t1305 * t1073 + t1303 * t1092;
t948 = t1303 * t1010 + t1305 * t1011;
t1213 = t1303 * qJDD(2) + t1305 * t1225;
t1331 = t1305 * qJDD(2) - t1303 * t1225;
t1100 = -t1177 * qJD(5) + t1310 * t1213 + t1307 * t1331;
t1259 = qJD(5) + t1264;
t1160 = t1259 * t1177;
t1414 = -t1160 + t1100;
t1211 = t1264 * t1234;
t1412 = -t1213 - t1211;
t1149 = -t1213 + t1211;
t1212 = t1264 * t1235;
t1146 = t1331 + t1212;
t1343 = t1307 * t1213 - t1310 * t1331;
t1055 = (qJD(5) - t1259) * t1179 + t1343;
t1173 = t1177 ^ 2;
t1174 = t1179 ^ 2;
t1409 = t1234 ^ 2;
t1233 = t1235 ^ 2;
t1258 = t1259 ^ 2;
t1408 = t1264 ^ 2;
t1263 = t1266 ^ 2;
t1407 = 0.2e1 * qJD(3);
t971 = pkin(4) * t1413 + t1412 * pkin(7) - t1010;
t1200 = t1264 * pkin(4) - t1235 * pkin(7);
t977 = -t1409 * pkin(4) + pkin(7) * t1331 - t1264 * t1200 + t1011;
t926 = t1307 * t977 - t1310 * t971;
t927 = t1307 * t971 + t1310 * t977;
t888 = t1307 * t927 - t1310 * t926;
t1406 = pkin(4) * t888;
t1058 = t1160 + t1100;
t998 = -t1055 * t1307 - t1310 * t1058;
t1405 = pkin(4) * t998;
t1344 = t1304 * t1185 - t1306 * t1314;
t1106 = t1266 * t1407 + t1344;
t1030 = -t1306 * t1106 + t1304 * t1107;
t1404 = pkin(2) * t1030;
t1195 = -t1224 + t1395;
t1197 = t1257 + t1225;
t1130 = t1304 * t1195 - t1306 * t1197;
t1403 = pkin(2) * t1130;
t1402 = pkin(3) * t1304;
t1401 = t1303 * t888;
t1400 = t1305 * t888;
t1072 = qJDD(4) - t1410 * qJ(4) - qJDD(2) * pkin(3) + (t1407 + t1214) * t1266 + t1344;
t1399 = -pkin(3) * t1072 + qJ(4) * t948;
t1393 = t1259 * t1179;
t1392 = t1259 * t1307;
t1391 = t1259 * t1310;
t1390 = t1264 * t1303;
t1389 = t1264 * t1304;
t1388 = t1264 * t1305;
t1387 = t1264 * t1306;
t1386 = t1266 * t1304;
t1385 = t1266 * t1306;
t1300 = t1308 ^ 2;
t1384 = t1300 * t1313;
t1068 = t1303 * t1072;
t1151 = t1190 + t1224;
t1383 = t1303 * t1151;
t1382 = t1304 * t1191;
t1217 = qJDD(2) + t1223;
t1381 = t1304 * t1217;
t1380 = t1304 * t1224;
t1069 = t1305 * t1072;
t1379 = t1305 * t1151;
t1378 = t1306 * t1191;
t1377 = t1306 * t1217;
t1018 = -pkin(4) * t1331 - t1409 * pkin(7) + t1235 * t1200 + t1072;
t1376 = t1307 * t1018;
t1097 = t1128 + t1221;
t1375 = t1307 * t1097;
t1374 = t1308 * t1030;
t1268 = t1313 * pkin(6) + t1328;
t1373 = t1308 * t1268;
t1292 = t1311 * t1370;
t1281 = qJDD(2) + t1292;
t1372 = t1308 * t1281;
t1282 = qJDD(2) - t1292;
t1371 = t1308 * t1282;
t1369 = t1310 * t1018;
t1368 = t1310 * t1097;
t1367 = t1311 * t1030;
t1366 = t1311 * t1268;
t1275 = -0.2e1 * t1350 + t1359;
t1236 = t1311 * t1275;
t1365 = t1311 * t1282;
t1362 = t1300 + t1301;
t1360 = t1309 * qJDD(1);
t1358 = t1312 * qJDD(1);
t1357 = t1312 * qJDD(2);
t1356 = -pkin(3) * t1306 - pkin(2);
t1355 = t1234 * t1390;
t1354 = t1304 * t1128;
t1353 = t1304 * t1190;
t1352 = t1306 * t1128;
t1351 = t1306 * t1190;
t1349 = t1309 * t1223;
t1348 = t1312 * t1223;
t1170 = -t1233 - t1408;
t1095 = -t1303 * t1170 - t1379;
t1347 = pkin(3) * t1149 + qJ(4) * t1095 + t1068;
t1163 = -t1408 - t1409;
t1087 = t1305 * t1163 - t1421;
t1145 = -t1212 + t1331;
t1346 = pkin(3) * t1145 + qJ(4) * t1087 - t1069;
t889 = t1307 * t926 + t1310 * t927;
t1031 = t1304 * t1106 + t1306 * t1107;
t1246 = t1311 * g(3) + t1316;
t1188 = t1308 * t1246 + t1311 * t1247;
t1342 = -t1309 * t1283 - t1312 * t1284;
t1076 = -t1173 - t1174;
t1000 = -t1055 * t1310 + t1307 * t1058;
t880 = -pkin(4) * t1076 + pkin(7) * t1000 + t889;
t882 = -pkin(7) * t998 - t888;
t938 = t1305 * t1000 - t1303 * t998;
t1341 = -pkin(3) * t1076 + qJ(4) * t938 + t1303 * t882 + t1305 * t880;
t1054 = (qJD(5) + t1259) * t1179 + t1343;
t1108 = -t1258 - t1173;
t1029 = t1310 * t1108 - t1417;
t940 = -pkin(4) * t1054 + pkin(7) * t1029 - t1369;
t1028 = t1307 * t1108 + t1416;
t958 = -pkin(7) * t1028 + t1376;
t967 = -t1303 * t1028 + t1305 * t1029;
t1340 = -pkin(3) * t1054 + qJ(4) * t967 + t1303 * t958 + t1305 * t940;
t1133 = -t1174 - t1258;
t1037 = -t1307 * t1133 - t1368;
t945 = -pkin(4) * t1414 + pkin(7) * t1037 + t1376;
t1036 = t1310 * t1133 - t1375;
t965 = -pkin(7) * t1036 + t1369;
t980 = -t1303 * t1036 + t1305 * t1037;
t1339 = -pkin(3) * t1414 + qJ(4) * t980 + t1303 * t965 + t1305 * t945;
t931 = -t1306 * t1072 + t1304 * t948;
t1338 = pkin(2) * t931 + t1399;
t1337 = t1309 * t1292;
t1336 = t1312 * t1292;
t1081 = t1305 * t1146 - t1303 * t1412;
t1156 = t1233 + t1409;
t1335 = pkin(3) * t1156 + qJ(4) * t1081 + t948;
t1252 = -t1263 - t1410;
t1166 = t1306 * t1252 - t1381;
t1334 = pkin(2) * t1166 - t1363;
t1278 = -t1309 * t1313 + t1358;
t1333 = -pkin(5) * t1278 - t1309 * g(3);
t947 = -t1305 * t1010 + t1303 * t1011;
t1187 = t1311 * t1246 - t1308 * t1247;
t1332 = t1312 * t1283 - t1309 * t1284;
t1038 = t1304 * t1087 + t1306 * t1145;
t1330 = pkin(2) * t1038 + t1346;
t1040 = t1304 * t1095 + t1306 * t1149;
t1329 = pkin(2) * t1040 + t1347;
t1327 = pkin(4) * t1028 - t926;
t942 = -t1306 * t1054 + t1304 * t967;
t1326 = pkin(2) * t942 + t1340;
t949 = t1304 * t980 - t1306 * t1414;
t1325 = pkin(2) * t949 + t1339;
t923 = -t1306 * t1076 + t1304 * t938;
t1324 = pkin(2) * t923 + t1341;
t870 = t1305 * t889 - t1401;
t885 = -pkin(4) * t1018 + pkin(7) * t889;
t1323 = -pkin(3) * t1018 - pkin(7) * t1401 + qJ(4) * t870 + t1305 * t885;
t1034 = t1304 * t1081 + t1306 * t1156;
t1321 = pkin(2) * t1034 + t1335;
t1320 = pkin(4) * t1036 - t927;
t864 = -t1306 * t1018 + t1304 * t870;
t1319 = pkin(2) * t864 + t1323;
t1215 = -t1408 - t1410;
t1154 = t1304 * t1215 + t1418;
t1318 = pkin(2) * t1154 - t1106;
t1296 = t1309 * qJDD(2);
t1290 = t1297 - t1410;
t1289 = -t1384 - t1410;
t1288 = -t1384 + t1410;
t1280 = -t1297 + t1384;
t1279 = t1297 + t1384;
t1277 = t1312 * t1313 + t1360;
t1276 = t1362 * qJDD(1);
t1273 = 0.2e1 * t1294 + t1361;
t1271 = t1311 * t1281;
t1270 = t1362 * t1398;
t1260 = -pkin(5) * t1277 + t1312 * g(3);
t1251 = -t1263 + t1410;
t1250 = t1408 - t1410;
t1249 = t1311 * t1274 - t1300 * t1398;
t1248 = -t1301 * t1398 + t1308 * t1322;
t1245 = -t1308 * t1289 - t1365;
t1244 = -t1308 * t1288 + t1271;
t1243 = t1311 * t1291 - t1372;
t1242 = t1311 * t1290 - t1371;
t1241 = t1311 * t1289 - t1371;
t1240 = t1311 * t1288 + t1372;
t1239 = t1308 * t1291 + t1271;
t1238 = t1308 * t1290 + t1365;
t1237 = (t1274 + t1294) * t1308;
t1227 = -t1308 * t1273 + t1236;
t1226 = t1311 * t1273 + t1308 * t1275;
t1220 = t1263 - t1408;
t1219 = t1306 * t1224;
t1209 = (t1386 - t1387) * qJD(2);
t1208 = (-t1385 - t1389) * qJD(2);
t1207 = -t1233 + t1408;
t1206 = -t1408 + t1409;
t1203 = t1235 * t1388;
t1202 = -pkin(6) * t1241 - t1366;
t1201 = -pkin(6) * t1239 - t1373;
t1199 = -pkin(1) * t1241 + t1247;
t1198 = -pkin(1) * t1239 + t1246;
t1192 = -t1408 - t1263;
t1189 = t1233 - t1409;
t1184 = -qJD(2) * t1386 + t1306 * t1225;
t1183 = qJD(2) * t1385 + t1304 * t1225;
t1182 = qJD(2) * t1387 + t1380;
t1181 = qJD(2) * t1389 - t1219;
t1176 = pkin(1) * t1275 + pkin(6) * t1243 + t1366;
t1175 = -pkin(1) * t1273 + pkin(6) * t1245 - t1373;
t1169 = -t1304 * t1252 - t1377;
t1168 = -t1304 * t1251 + t1418;
t1167 = t1306 * t1250 - t1381;
t1165 = t1306 * t1251 + t1420;
t1164 = t1304 * t1250 + t1377;
t1162 = pkin(1) * t1268 + pkin(6) * t1188;
t1159 = -t1174 + t1258;
t1158 = t1173 - t1258;
t1157 = pkin(1) * t1279 + pkin(6) * t1276 + t1188;
t1155 = t1306 * t1215 - t1420;
t1141 = t1305 * t1213 - t1235 * t1390;
t1140 = -t1303 * t1213 - t1203;
t1139 = t1234 * t1388 - t1303 * t1331;
t1138 = t1305 * t1331 + t1355;
t1137 = (-t1234 * t1305 + t1235 * t1303) * t1264;
t1136 = -t1203 - t1355;
t1135 = -t1308 * t1208 + t1311 * t1209;
t1134 = t1311 * t1208 + t1308 * t1209;
t1132 = t1306 * t1195 + t1304 * t1197;
t1131 = -t1306 * t1193 - t1304 * t1196;
t1129 = -t1304 * t1193 + t1306 * t1196;
t1127 = t1174 - t1173;
t1126 = -qJ(3) * t1166 - t1378;
t1125 = -t1308 * t1183 + t1311 * t1184;
t1124 = -t1308 * t1181 + t1311 * t1182;
t1123 = t1311 * t1183 + t1308 * t1184;
t1122 = t1311 * t1181 + t1308 * t1182;
t1121 = t1306 * t1137 + t1380;
t1120 = t1304 * t1137 - t1219;
t1119 = -t1308 * t1166 + t1311 * t1169;
t1118 = -t1308 * t1165 + t1311 * t1168;
t1117 = -t1308 * t1164 + t1311 * t1167;
t1116 = t1311 * t1166 + t1308 * t1169;
t1115 = t1311 * t1165 + t1308 * t1168;
t1114 = t1311 * t1164 + t1308 * t1167;
t1113 = t1305 * t1206 - t1383;
t1112 = -t1303 * t1207 + t1419;
t1111 = t1303 * t1206 + t1379;
t1110 = t1305 * t1207 + t1421;
t1109 = -qJ(3) * t1154 - t1382;
t1104 = t1306 * t1141 + t1353;
t1103 = t1306 * t1139 - t1353;
t1102 = t1304 * t1141 - t1351;
t1101 = t1304 * t1139 + t1351;
t1099 = -t1179 * qJD(5) - t1343;
t1094 = t1305 * t1170 - t1383;
t1089 = (-t1177 * t1310 + t1179 * t1307) * t1259;
t1088 = (-t1177 * t1307 - t1179 * t1310) * t1259;
t1086 = t1303 * t1163 + t1419;
t1084 = -t1308 * t1154 + t1311 * t1155;
t1083 = t1311 * t1154 + t1308 * t1155;
t1082 = -pkin(2) * t1196 + qJ(3) * t1169 - t1382;
t1080 = t1305 * t1145 + t1149 * t1303;
t1079 = t1303 * t1146 + t1305 * t1412;
t1078 = t1303 * t1145 - t1149 * t1305;
t1074 = -pkin(2) * t1193 + qJ(3) * t1155 + t1378;
t1067 = -t1308 * t1130 + t1311 * t1132;
t1066 = -t1308 * t1129 + t1311 * t1131;
t1065 = t1311 * t1130 + t1308 * t1132;
t1064 = t1311 * t1129 + t1308 * t1131;
t1063 = t1306 * t1113 + t1146 * t1304;
t1062 = t1306 * t1112 - t1304 * t1412;
t1061 = t1304 * t1113 - t1146 * t1306;
t1060 = t1304 * t1112 + t1306 * t1412;
t1051 = t1310 * t1100 - t1179 * t1392;
t1050 = t1307 * t1100 + t1179 * t1391;
t1049 = -t1307 * t1099 + t1177 * t1391;
t1048 = t1310 * t1099 + t1177 * t1392;
t1047 = t1310 * t1158 - t1375;
t1046 = -t1307 * t1159 + t1416;
t1045 = t1307 * t1158 + t1368;
t1044 = t1310 * t1159 + t1417;
t1043 = t1306 * t1080 + t1304 * t1189;
t1042 = t1304 * t1080 - t1306 * t1189;
t1041 = t1306 * t1095 - t1304 * t1149;
t1039 = t1306 * t1087 - t1304 * t1145;
t1035 = t1306 * t1081 - t1304 * t1156;
t1033 = -t1308 * t1120 + t1311 * t1121;
t1032 = t1311 * t1120 + t1308 * t1121;
t1027 = -t1308 * t1102 + t1311 * t1104;
t1026 = -t1308 * t1101 + t1311 * t1103;
t1025 = t1311 * t1102 + t1308 * t1104;
t1024 = t1311 * t1101 + t1308 * t1103;
t1023 = -t1303 * t1088 + t1305 * t1089;
t1022 = t1305 * t1088 + t1303 * t1089;
t1021 = -pkin(1) * t1065 - t1403;
t1020 = -pkin(1) * t1116 + t1253 - t1334;
t1019 = pkin(2) * t1191 + qJ(3) * t1031;
t1016 = -qJ(4) * t1094 + t1069;
t1015 = -qJ(4) * t1086 + t1068;
t1014 = t1306 * t1023 + t1304 * t1221;
t1013 = t1304 * t1023 - t1306 * t1221;
t1012 = -pkin(1) * t1083 - t1318;
t1009 = -qJ(3) * t1130 - t1030;
t1006 = -pkin(6) * t1116 - t1308 * t1082 + t1311 * t1126;
t1005 = -pkin(2) * t1192 + qJ(3) * t1132 + t1031;
t1004 = -t1308 * t1061 + t1311 * t1063;
t1003 = -t1308 * t1060 + t1311 * t1062;
t1002 = t1311 * t1061 + t1308 * t1063;
t1001 = t1311 * t1060 + t1308 * t1062;
t999 = -t1310 * t1054 - t1307 * t1414;
t997 = -t1307 * t1054 + t1310 * t1414;
t996 = -t1303 * t1050 + t1305 * t1051;
t995 = -t1303 * t1048 + t1305 * t1049;
t994 = t1305 * t1050 + t1303 * t1051;
t993 = t1305 * t1048 + t1303 * t1049;
t992 = -t1308 * t1042 + t1311 * t1043;
t991 = t1311 * t1042 + t1308 * t1043;
t990 = -t1303 * t1045 + t1305 * t1047;
t989 = -t1303 * t1044 + t1305 * t1046;
t988 = t1305 * t1045 + t1303 * t1047;
t987 = t1305 * t1044 + t1303 * t1046;
t986 = -pkin(1) * t1196 + pkin(6) * t1119 + t1311 * t1082 + t1308 * t1126;
t985 = -t1308 * t1040 + t1311 * t1041;
t984 = t1311 * t1040 + t1308 * t1041;
t983 = -pkin(6) * t1083 - t1308 * t1074 + t1311 * t1109;
t982 = -t1308 * t1038 + t1311 * t1039;
t981 = t1311 * t1038 + t1308 * t1039;
t979 = t1305 * t1036 + t1303 * t1037;
t976 = -t1308 * t1034 + t1311 * t1035;
t975 = t1311 * t1034 + t1308 * t1035;
t974 = -pkin(3) * t1094 + t1011;
t973 = -pkin(1) * t1193 + pkin(6) * t1084 + t1311 * t1074 + t1308 * t1109;
t972 = -pkin(3) * t1086 + t1010;
t969 = t1311 * t1031 - t1374;
t968 = t1308 * t1031 + t1367;
t966 = t1305 * t1028 + t1303 * t1029;
t963 = t1306 * t996 + t1354;
t962 = t1306 * t995 - t1354;
t961 = t1304 * t996 - t1352;
t960 = t1304 * t995 + t1352;
t956 = -t1304 * t1055 + t1306 * t990;
t955 = t1304 * t1058 + t1306 * t989;
t954 = t1306 * t1055 + t1304 * t990;
t953 = -t1306 * t1058 + t1304 * t989;
t952 = -t1308 * t1013 + t1311 * t1014;
t951 = t1311 * t1013 + t1308 * t1014;
t950 = t1304 * t1414 + t1306 * t980;
t943 = t1304 * t1054 + t1306 * t967;
t941 = -pkin(1) * t968 - t1404;
t937 = -t1303 * t997 + t1305 * t999;
t936 = t1303 * t1000 + t1305 * t998;
t935 = t1303 * t999 + t1305 * t997;
t933 = -qJ(4) * t1079 - t947;
t932 = t1304 * t1072 + t1306 * t948;
t930 = t1304 * t1127 + t1306 * t937;
t929 = -t1306 * t1127 + t1304 * t937;
t928 = -pkin(6) * t1065 - t1308 * t1005 + t1311 * t1009;
t924 = t1304 * t1076 + t1306 * t938;
t922 = -pkin(1) * t1192 + pkin(6) * t1067 + t1311 * t1005 + t1308 * t1009;
t921 = -qJ(3) * t1040 + t1306 * t1016 - t1304 * t974;
t920 = -qJ(3) * t1038 + t1306 * t1015 - t1304 * t972;
t919 = -pkin(1) * t984 - t1329;
t918 = -pkin(6) * t968 - qJ(3) * t1367 - t1308 * t1019;
t917 = -t1308 * t961 + t1311 * t963;
t916 = -t1308 * t960 + t1311 * t962;
t915 = t1308 * t963 + t1311 * t961;
t914 = t1308 * t962 + t1311 * t960;
t913 = -pkin(1) * t981 - t1330;
t912 = pkin(1) * t1191 + pkin(6) * t969 - qJ(3) * t1374 + t1311 * t1019;
t911 = -pkin(2) * t1094 + qJ(3) * t1041 + t1304 * t1016 + t1306 * t974;
t910 = -pkin(2) * t1086 + qJ(3) * t1039 + t1304 * t1015 + t1306 * t972;
t909 = -pkin(3) * t936 - t1405;
t908 = -qJ(3) * t1034 + t1079 * t1402 + t1306 * t933;
t907 = -t1308 * t954 + t1311 * t956;
t906 = -t1308 * t953 + t1311 * t955;
t905 = t1308 * t956 + t1311 * t954;
t904 = t1308 * t955 + t1311 * t953;
t903 = -t1308 * t949 + t1311 * t950;
t902 = t1308 * t950 + t1311 * t949;
t901 = qJ(3) * t1035 + t1079 * t1356 + t1304 * t933;
t900 = -t1308 * t942 + t1311 * t943;
t899 = t1308 * t943 + t1311 * t942;
t898 = -pkin(1) * t975 - t1321;
t897 = -pkin(3) * t979 - t1320;
t896 = -qJ(4) * t979 - t1303 * t945 + t1305 * t965;
t895 = -pkin(3) * t966 - t1327;
t894 = -t1308 * t931 + t1311 * t932;
t893 = t1308 * t932 + t1311 * t931;
t892 = -qJ(4) * t966 - t1303 * t940 + t1305 * t958;
t891 = -t1308 * t929 + t1311 * t930;
t890 = t1308 * t930 + t1311 * t929;
t887 = -t1308 * t923 + t1311 * t924;
t886 = t1308 * t924 + t1311 * t923;
t883 = -qJ(3) * t931 + (-qJ(4) * t1306 + t1402) * t947;
t878 = -pkin(6) * t984 - t1308 * t911 + t1311 * t921;
t877 = -pkin(6) * t981 - t1308 * t910 + t1311 * t920;
t876 = -pkin(1) * t1094 + pkin(6) * t985 + t1308 * t921 + t1311 * t911;
t875 = -pkin(1) * t1086 + pkin(6) * t982 + t1308 * t920 + t1311 * t910;
t874 = qJ(3) * t932 + (-qJ(4) * t1304 + t1356) * t947;
t873 = -pkin(6) * t975 - t1308 * t901 + t1311 * t908;
t872 = -pkin(1) * t1079 + pkin(6) * t976 + t1308 * t908 + t1311 * t901;
t871 = -pkin(1) * t893 - t1338;
t869 = t1303 * t889 + t1400;
t867 = -pkin(1) * t902 - t1325;
t866 = -qJ(3) * t949 - t1304 * t897 + t1306 * t896;
t865 = t1304 * t1018 + t1306 * t870;
t863 = -pkin(1) * t899 - t1326;
t862 = -qJ(3) * t942 - t1304 * t895 + t1306 * t892;
t861 = -pkin(2) * t979 + qJ(3) * t950 + t1304 * t896 + t1306 * t897;
t860 = -pkin(2) * t966 + qJ(3) * t943 + t1304 * t892 + t1306 * t895;
t859 = -qJ(4) * t936 - t1303 * t880 + t1305 * t882;
t858 = -pkin(3) * t869 - t1406;
t857 = -pkin(6) * t893 - t1308 * t874 + t1311 * t883;
t856 = -pkin(1) * t947 + pkin(6) * t894 + t1308 * t883 + t1311 * t874;
t855 = -qJ(3) * t923 - t1304 * t909 + t1306 * t859;
t854 = -pkin(7) * t1400 - qJ(4) * t869 - t1303 * t885;
t853 = -t1308 * t864 + t1311 * t865;
t852 = t1308 * t865 + t1311 * t864;
t851 = -pkin(2) * t936 + qJ(3) * t924 + t1304 * t859 + t1306 * t909;
t850 = -pkin(1) * t886 - t1324;
t849 = -pkin(6) * t902 - t1308 * t861 + t1311 * t866;
t848 = -pkin(1) * t979 + pkin(6) * t903 + t1308 * t866 + t1311 * t861;
t847 = -pkin(6) * t899 - t1308 * t860 + t1311 * t862;
t846 = -pkin(1) * t966 + pkin(6) * t900 + t1308 * t862 + t1311 * t860;
t845 = -qJ(3) * t864 - t1304 * t858 + t1306 * t854;
t844 = -pkin(6) * t886 - t1308 * t851 + t1311 * t855;
t843 = -pkin(1) * t936 + pkin(6) * t887 + t1308 * t855 + t1311 * t851;
t842 = -pkin(1) * t852 - t1319;
t841 = -pkin(2) * t869 + qJ(3) * t865 + t1304 * t854 + t1306 * t858;
t840 = -pkin(6) * t852 - t1308 * t841 + t1311 * t845;
t839 = -pkin(1) * t869 + pkin(6) * t853 + t1308 * t845 + t1311 * t841;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1278, 0, -t1277, 0, t1333, -t1260, -t1332, -pkin(5) * t1332, t1312 * t1249 - t1337, t1312 * t1227 + t1309 * t1280, t1312 * t1244 + t1308 * t1360, t1312 * t1248 + t1337, t1312 * t1242 + t1309 * t1359, t1312 * t1270 + t1296, t1312 * t1201 - t1309 * t1198 - pkin(5) * (t1309 * t1243 + t1312 * t1275), t1312 * t1202 - t1309 * t1199 - pkin(5) * (t1309 * t1245 - t1312 * t1273), t1312 * t1187 - pkin(5) * (t1309 * t1276 + t1312 * t1279), -pkin(5) * (t1309 * t1188 + t1312 * t1268) - (t1309 * pkin(1) - t1312 * pkin(6)) * t1187, t1312 * t1125 + t1349, t1312 * t1066 + t1309 * t1220, t1312 * t1118 + t1309 * t1197, t1312 * t1124 - t1349, t1312 * t1117 + t1195 * t1309, t1312 * t1135 + t1296, t1312 * t983 - t1309 * t1012 - pkin(5) * (t1309 * t1084 - t1312 * t1193), t1312 * t1006 - t1309 * t1020 - pkin(5) * (t1309 * t1119 - t1312 * t1196), t1312 * t928 - t1309 * t1021 - pkin(5) * (t1309 * t1067 - t1312 * t1192), t1312 * t918 - t1309 * t941 - pkin(5) * (t1312 * t1191 + t1309 * t969), t1312 * t1027 - t1309 * t1140, t1309 * t1078 + t1312 * t992, t1312 * t1003 + t1309 * t1110, t1312 * t1026 + t1309 * t1138, t1312 * t1004 + t1309 * t1111, t1312 * t1033 + t1309 * t1136, t1312 * t877 - t1309 * t913 - pkin(5) * (-t1312 * t1086 + t1309 * t982), t1312 * t878 - t1309 * t919 - pkin(5) * (-t1312 * t1094 + t1309 * t985), t1312 * t873 - t1309 * t898 - pkin(5) * (-t1312 * t1079 + t1309 * t976), t1312 * t857 - t1309 * t871 - pkin(5) * (t1309 * t894 - t1312 * t947), t1309 * t994 + t1312 * t917, t1309 * t935 + t1312 * t891, t1309 * t987 + t1312 * t906, t1309 * t993 + t1312 * t916, t1309 * t988 + t1312 * t907, t1309 * t1022 + t1312 * t952, t1312 * t847 - t1309 * t863 - pkin(5) * (t1309 * t900 - t1312 * t966), t1312 * t849 - t1309 * t867 - pkin(5) * (t1309 * t903 - t1312 * t979), t1312 * t844 - t1309 * t850 - pkin(5) * (t1309 * t887 - t1312 * t936), t1312 * t840 - t1309 * t842 - pkin(5) * (t1309 * t853 - t1312 * t869); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1277, 0, t1278, 0, t1260, t1333, t1342, pkin(5) * t1342, t1309 * t1249 + t1336, t1309 * t1227 - t1312 * t1280, t1309 * t1244 - t1308 * t1358, t1309 * t1248 - t1336, t1309 * t1242 - t1311 * t1358, t1309 * t1270 - t1357, t1309 * t1201 + t1312 * t1198 + pkin(5) * (t1312 * t1243 - t1309 * t1275), t1309 * t1202 + t1312 * t1199 + pkin(5) * (t1312 * t1245 + t1309 * t1273), t1309 * t1187 + pkin(5) * (t1312 * t1276 - t1309 * t1279), pkin(5) * (t1312 * t1188 - t1309 * t1268) - (-t1312 * pkin(1) - t1309 * pkin(6)) * t1187, t1309 * t1125 - t1348, t1309 * t1066 - t1312 * t1220, t1309 * t1118 - t1312 * t1197, t1309 * t1124 + t1348, t1309 * t1117 - t1195 * t1312, t1309 * t1135 - t1357, t1309 * t983 + t1312 * t1012 + pkin(5) * (t1312 * t1084 + t1309 * t1193), t1309 * t1006 + t1312 * t1020 + pkin(5) * (t1312 * t1119 + t1309 * t1196), t1309 * t928 + t1312 * t1021 + pkin(5) * (t1312 * t1067 + t1309 * t1192), t1309 * t918 + t1312 * t941 + pkin(5) * (-t1309 * t1191 + t1312 * t969), t1309 * t1027 + t1312 * t1140, -t1312 * t1078 + t1309 * t992, t1309 * t1003 - t1312 * t1110, t1309 * t1026 - t1312 * t1138, t1309 * t1004 - t1312 * t1111, t1309 * t1033 - t1312 * t1136, t1309 * t877 + t1312 * t913 + pkin(5) * (t1309 * t1086 + t1312 * t982), t1309 * t878 + t1312 * t919 + pkin(5) * (t1309 * t1094 + t1312 * t985), t1309 * t873 + t1312 * t898 + pkin(5) * (t1309 * t1079 + t1312 * t976), t1309 * t857 + t1312 * t871 + pkin(5) * (t1309 * t947 + t1312 * t894), t1309 * t917 - t1312 * t994, t1309 * t891 - t1312 * t935, t1309 * t906 - t1312 * t987, t1309 * t916 - t1312 * t993, t1309 * t907 - t1312 * t988, -t1312 * t1022 + t1309 * t952, t1309 * t847 + t1312 * t863 + pkin(5) * (t1309 * t966 + t1312 * t900), t1309 * t849 + t1312 * t867 + pkin(5) * (t1309 * t979 + t1312 * t903), t1309 * t844 + t1312 * t850 + pkin(5) * (t1309 * t936 + t1312 * t887), t1309 * t840 + t1312 * t842 + pkin(5) * (t1309 * t869 + t1312 * t853); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1283, t1284, 0, 0, t1237, t1226, t1240, t1236, t1238, 0, t1176, t1175, t1157, t1162, t1123, t1064, t1115, t1122, t1114, t1134, t973, t986, t922, t912, t1025, t991, t1001, t1024, t1002, t1032, t875, t876, t872, t856, t915, t890, t904, t914, t905, t951, t846, t848, t843, t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1313, 0, 0, -g(3), -t1283, 0, t1249, t1227, t1244, t1248, t1242, t1270, t1201, t1202, t1187, pkin(6) * t1187, t1125, t1066, t1118, t1124, t1117, t1135, t983, t1006, t928, t918, t1027, t992, t1003, t1026, t1004, t1033, t877, t878, t873, t857, t917, t891, t906, t916, t907, t952, t847, t849, t844, t840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1313, 0, qJDD(1), 0, g(3), 0, -t1284, 0, t1292, -t1280, -t1361, -t1292, -t1359, -qJDD(2), t1198, t1199, 0, pkin(1) * t1187, -t1223, -t1220, -t1197, t1223, -t1195, -qJDD(2), t1012, t1020, t1021, t941, t1140, -t1078, -t1110, -t1138, -t1111, -t1136, t913, t919, t898, t871, -t994, -t935, -t987, -t993, -t988, -t1022, t863, t867, t850, t842; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1283, t1284, 0, 0, t1237, t1226, t1240, t1236, t1238, 0, t1176, t1175, t1157, t1162, t1123, t1064, t1115, t1122, t1114, t1134, t973, t986, t922, t912, t1025, t991, t1001, t1024, t1002, t1032, t875, t876, t872, t856, t915, t890, t904, t914, t905, t951, t846, t848, t843, t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1274, t1275, t1281, -t1294, t1290, t1294, 0, -t1268, t1246, 0, t1184, t1131, t1168, t1182, t1167, t1209, t1109, t1126, t1009, -qJ(3) * t1030, t1104, t1043, t1062, t1103, t1063, t1121, t920, t921, t908, t883, t963, t930, t955, t962, t956, t1014, t862, t866, t855, t845; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1350, t1273, t1288, -t1322, t1282, -t1350, t1268, 0, t1247, 0, t1183, t1129, t1165, t1181, t1164, t1208, t1074, t1082, t1005, t1019, t1102, t1042, t1060, t1101, t1061, t1120, t910, t911, t901, t874, t961, t929, t953, t960, t954, t1013, t860, t861, t851, t841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1292, t1280, t1361, t1292, t1359, qJDD(2), -t1246, -t1247, 0, 0, t1223, t1220, t1197, -t1223, t1195, qJDD(2), t1318, t1334 + 0.2e1 * t1394, t1403, t1404, -t1140, t1078, t1110, t1138, t1111, t1136, t1330, t1329, t1321, t1338, t994, t935, t987, t993, t988, t1022, t1326, t1325, t1324, t1319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1225, -t1193, t1411, t1257, t1250, -t1257, 0, -t1191, t1106, 0, t1141, t1080, t1112, t1139, t1113, t1137, t1015, t1016, t933, -qJ(4) * t947, t996, t937, t989, t995, t990, t1023, t892, t896, t859, t854; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1395, t1196, t1251, -t1224, t1217, -t1395, t1191, 0, t1107, 0, -t1190, -t1189, t1412, t1190, -t1146, -t1224, t972, t974, -pkin(3) * t1079, -pkin(3) * t947, -t1128, -t1127, -t1058, t1128, t1055, -t1221, t895, t897, t909, t858; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1223, t1220, t1197, -t1223, t1195, qJDD(2), -t1106, -t1107, 0, 0, -t1140, t1078, t1110, t1138, t1111, t1136, t1346, t1347, t1335, t1399, t994, t935, t987, t993, t988, t1022, t1340, t1339, t1341, t1323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1213, t1145, t1413, t1211, t1206, -t1211, 0, t1072, t1010, 0, t1051, t999, t1046, t1049, t1047, t1089, t958, t965, t882, -pkin(7) * t888; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1212, -t1149, t1207, t1331, t1151, -t1212, -t1072, 0, t1011, 0, t1050, t997, t1044, t1048, t1045, t1088, t940, t945, t880, t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1190, t1189, -t1412, -t1190, t1146, t1224, -t1010, -t1011, 0, 0, t1128, t1127, t1058, -t1128, -t1055, t1221, t1327, t1320, t1405, t1406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1100, -t1054, t1415, t1160, t1158, -t1160, 0, t1018, t926, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1393, t1414, t1159, t1099, t1097, -t1393, -t1018, 0, t927, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1128, t1127, t1058, -t1128, -t1055, t1221, -t926, -t927, 0, 0;];
m_new_reg = t1;
