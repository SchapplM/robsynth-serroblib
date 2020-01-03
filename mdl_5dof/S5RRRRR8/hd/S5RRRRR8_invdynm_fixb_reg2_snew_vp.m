% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRRR8
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRRR8_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:26:25
% EndTime: 2019-12-31 22:26:44
% DurationCPUTime: 20.33s
% Computational Cost: add. (179918->741), mult. (363621->1042), div. (0->0), fcn. (262823->10), ass. (0->504)
t1344 = sin(qJ(5));
t1346 = sin(qJ(3));
t1347 = sin(qJ(2));
t1351 = cos(qJ(3));
t1352 = cos(qJ(2));
t1307 = (t1346 * t1352 + t1347 * t1351) * qJD(1);
t1341 = qJD(2) + qJD(3);
t1345 = sin(qJ(4));
t1350 = cos(qJ(4));
t1276 = t1345 * t1307 - t1350 * t1341;
t1277 = t1350 * t1307 + t1345 * t1341;
t1349 = cos(qJ(5));
t1232 = t1349 * t1276 + t1344 * t1277;
t1234 = -t1344 * t1276 + t1349 * t1277;
t1186 = t1234 * t1232;
t1436 = qJD(1) * t1352;
t1335 = qJD(2) * t1436;
t1404 = t1347 * qJDD(1);
t1315 = t1335 + t1404;
t1338 = t1352 * qJDD(1);
t1437 = qJD(1) * t1347;
t1396 = qJD(2) * t1437;
t1316 = t1338 - t1396;
t1387 = t1346 * t1315 - t1351 * t1316;
t1250 = -t1307 * qJD(3) - t1387;
t1247 = qJDD(4) - t1250;
t1246 = qJDD(5) + t1247;
t1451 = -t1186 + t1246;
t1457 = t1344 * t1451;
t1245 = t1277 * t1276;
t1449 = -t1245 + t1247;
t1456 = t1345 * t1449;
t1305 = t1346 * t1437 - t1351 * t1436;
t1269 = t1307 * t1305;
t1401 = qJDD(2) + qJDD(3);
t1448 = -t1269 + t1401;
t1455 = t1346 * t1448;
t1454 = t1349 * t1451;
t1453 = t1350 * t1449;
t1452 = t1351 * t1448;
t1251 = -t1305 * qJD(3) + t1351 * t1315 + t1346 * t1316;
t1296 = t1341 * t1305;
t1221 = t1251 - t1296;
t1343 = t1352 ^ 2;
t1355 = qJD(1) ^ 2;
t1348 = sin(qJ(1));
t1353 = cos(qJ(1));
t1325 = t1348 * g(1) - t1353 * g(2);
t1372 = qJDD(1) * pkin(1) + t1325;
t1373 = qJD(2) * pkin(2) - pkin(7) * t1437;
t1253 = t1316 * pkin(2) + (pkin(7) * t1343 + pkin(6)) * t1355 - t1373 * t1437 + t1372;
t1430 = t1341 * t1307;
t1128 = -t1221 * pkin(8) + (-t1250 + t1430) * pkin(3) - t1253;
t1326 = t1353 * g(1) + t1348 * g(2);
t1359 = -t1355 * pkin(1) + qJDD(1) * pkin(6) - t1326;
t1290 = -t1347 * g(3) + t1352 * t1359;
t1339 = t1343 * t1355;
t1244 = -pkin(2) * t1339 + t1316 * pkin(7) - qJD(2) * t1373 + t1290;
t1357 = t1347 * t1359;
t1417 = t1347 * t1355;
t1438 = qJD(1) * qJD(2);
t1356 = -t1357 - t1315 * pkin(7) + qJDD(2) * pkin(2) + (pkin(2) * t1417 + pkin(7) * t1438 - g(3)) * t1352;
t1189 = t1351 * t1244 + t1346 * t1356;
t1267 = t1305 * pkin(3) - t1307 * pkin(8);
t1446 = t1341 ^ 2;
t1139 = -pkin(3) * t1446 + pkin(8) * t1401 - t1305 * t1267 + t1189;
t1062 = -t1350 * t1128 + t1345 * t1139;
t1063 = t1345 * t1128 + t1350 * t1139;
t1004 = t1345 * t1062 + t1350 * t1063;
t1363 = -t1350 * t1251 - t1345 * t1401;
t1201 = -t1276 * qJD(4) - t1363;
t1364 = -t1345 * t1251 + t1350 * t1401;
t1358 = t1277 * qJD(4) - t1364;
t1117 = -t1232 * qJD(5) + t1349 * t1201 - t1344 * t1358;
t1301 = qJD(4) + t1305;
t1299 = qJD(5) + t1301;
t1208 = t1299 * t1232;
t1450 = -t1208 + t1117;
t1261 = t1301 * t1276;
t1182 = t1261 + t1201;
t1389 = t1344 * t1201 + t1349 * t1358;
t1083 = (qJD(5) - t1299) * t1234 + t1389;
t1218 = (qJD(3) - t1341) * t1307 + t1387;
t1230 = t1232 ^ 2;
t1231 = t1234 ^ 2;
t1447 = t1276 ^ 2;
t1274 = t1277 ^ 2;
t1297 = t1299 ^ 2;
t1300 = t1301 ^ 2;
t1303 = t1305 ^ 2;
t1304 = t1307 ^ 2;
t1030 = pkin(4) * t1449 - t1182 * pkin(9) - t1062;
t1254 = t1301 * pkin(4) - t1277 * pkin(9);
t1040 = -pkin(4) * t1447 - pkin(9) * t1358 - t1301 * t1254 + t1063;
t975 = -t1349 * t1030 + t1344 * t1040;
t976 = t1344 * t1030 + t1349 * t1040;
t936 = t1344 * t976 - t1349 * t975;
t1445 = pkin(4) * t936;
t1188 = t1346 * t1244 - t1351 * t1356;
t1118 = -t1351 * t1188 + t1346 * t1189;
t1444 = pkin(2) * t1118;
t1222 = t1251 + t1296;
t1161 = -t1218 * t1346 - t1351 * t1222;
t1443 = pkin(2) * t1161;
t1442 = pkin(3) * t1346;
t1086 = t1208 + t1117;
t1018 = -t1083 * t1344 - t1349 * t1086;
t1441 = pkin(4) * t1018;
t1440 = t1345 * t936;
t1439 = t1350 * t936;
t1435 = t1299 * t1234;
t1434 = t1299 * t1344;
t1433 = t1299 * t1349;
t1432 = t1301 * t1345;
t1431 = t1301 * t1350;
t1429 = t1341 * t1346;
t1428 = t1341 * t1351;
t1342 = t1347 ^ 2;
t1427 = t1342 * t1355;
t1138 = -t1401 * pkin(3) - t1446 * pkin(8) + t1307 * t1267 + t1188;
t1065 = pkin(4) * t1358 - pkin(9) * t1447 + t1277 * t1254 + t1138;
t1426 = t1344 * t1065;
t1144 = t1186 + t1246;
t1425 = t1344 * t1144;
t1134 = t1345 * t1138;
t1191 = t1245 + t1247;
t1424 = t1345 * t1191;
t1423 = t1346 * t1253;
t1265 = t1269 + t1401;
t1422 = t1346 * t1265;
t1421 = t1347 * t1118;
t1308 = t1355 * pkin(6) + t1372;
t1420 = t1347 * t1308;
t1332 = t1352 * t1417;
t1323 = qJDD(2) + t1332;
t1419 = t1347 * t1323;
t1324 = qJDD(2) - t1332;
t1418 = t1347 * t1324;
t1416 = t1349 * t1065;
t1415 = t1349 * t1144;
t1135 = t1350 * t1138;
t1414 = t1350 * t1191;
t1413 = t1351 * t1253;
t1412 = t1351 * t1265;
t1411 = t1352 * t1118;
t1410 = t1352 * t1308;
t1409 = t1352 * t1324;
t1406 = -pkin(3) * t1138 + pkin(8) * t1004;
t1405 = t1342 + t1343;
t1403 = t1348 * qJDD(1);
t1402 = t1353 * qJDD(1);
t1400 = -pkin(3) * t1351 - pkin(2);
t1398 = t1346 * t1186;
t1397 = t1346 * t1245;
t1395 = t1348 * t1269;
t1394 = t1351 * t1186;
t1393 = t1351 * t1245;
t1392 = t1353 * t1269;
t1229 = -t1274 - t1300;
t1142 = -t1345 * t1229 - t1414;
t1183 = (qJD(4) + t1301) * t1276 + t1363;
t1391 = pkin(3) * t1183 + pkin(8) * t1142 + t1134;
t1216 = -t1300 - t1447;
t1132 = t1350 * t1216 - t1456;
t1262 = t1301 * t1277;
t1179 = -t1262 - t1358;
t1390 = pkin(3) * t1179 + pkin(8) * t1132 - t1135;
t937 = t1344 * t975 + t1349 * t976;
t1119 = t1346 * t1188 + t1351 * t1189;
t1289 = t1352 * g(3) + t1357;
t1243 = t1347 * t1289 + t1352 * t1290;
t1386 = -t1348 * t1325 - t1353 * t1326;
t1133 = -t1230 - t1231;
t1020 = -t1083 * t1349 + t1344 * t1086;
t925 = -pkin(4) * t1133 + pkin(9) * t1020 + t937;
t929 = -pkin(9) * t1018 - t936;
t971 = -t1345 * t1018 + t1350 * t1020;
t1385 = -pkin(3) * t1133 + pkin(8) * t971 + t1345 * t929 + t1350 * t925;
t1166 = -t1297 - t1230;
t1090 = t1344 * t1166 + t1454;
t1009 = -pkin(9) * t1090 + t1426;
t1091 = t1349 * t1166 - t1457;
t1026 = -t1345 * t1090 + t1350 * t1091;
t1082 = (qJD(5) + t1299) * t1234 + t1389;
t983 = -pkin(4) * t1082 + pkin(9) * t1091 - t1416;
t1384 = -pkin(3) * t1082 + pkin(8) * t1026 + t1345 * t1009 + t1350 * t983;
t1193 = -t1231 - t1297;
t1105 = t1349 * t1193 - t1425;
t1021 = -pkin(9) * t1105 + t1416;
t1106 = -t1344 * t1193 - t1415;
t1038 = -t1345 * t1105 + t1350 * t1106;
t985 = -pkin(4) * t1450 + pkin(9) * t1106 + t1426;
t1383 = -pkin(3) * t1450 + pkin(8) * t1038 + t1345 * t1021 + t1350 * t985;
t1382 = t1348 * t1332;
t1381 = t1353 * t1332;
t1180 = (-qJD(4) + t1301) * t1277 + t1364;
t1111 = t1350 * t1180 + t1345 * t1182;
t1203 = t1274 + t1447;
t1380 = pkin(3) * t1203 + pkin(8) * t1111 + t1004;
t988 = t1346 * t1004 - t1351 * t1138;
t1379 = pkin(2) * t988 + t1406;
t1288 = -t1304 - t1446;
t1223 = t1351 * t1288 - t1422;
t1378 = pkin(2) * t1223 - t1189;
t1320 = -t1348 * t1355 + t1402;
t1377 = -pkin(5) * t1320 - t1348 * g(3);
t1003 = -t1350 * t1062 + t1345 * t1063;
t1242 = t1352 * t1289 - t1347 * t1290;
t1376 = t1353 * t1325 - t1348 * t1326;
t1088 = t1346 * t1132 + t1351 * t1179;
t1375 = pkin(2) * t1088 + t1390;
t1092 = t1346 * t1142 + t1351 * t1183;
t1374 = pkin(2) * t1092 + t1391;
t1263 = -t1446 - t1303;
t1206 = t1346 * t1263 + t1452;
t1371 = pkin(2) * t1206 - t1188;
t1370 = pkin(4) * t1090 - t975;
t959 = -t1351 * t1133 + t1346 * t971;
t1369 = pkin(2) * t959 + t1385;
t915 = t1350 * t937 - t1440;
t933 = -pkin(4) * t1065 + pkin(9) * t937;
t1368 = -pkin(3) * t1065 + pkin(8) * t915 - pkin(9) * t1440 + t1350 * t933;
t986 = t1346 * t1026 - t1351 * t1082;
t1367 = pkin(2) * t986 + t1384;
t990 = t1346 * t1038 - t1351 * t1450;
t1366 = pkin(2) * t990 + t1383;
t1068 = t1346 * t1111 + t1351 * t1203;
t1365 = pkin(2) * t1068 + t1380;
t1362 = pkin(4) * t1105 - t976;
t911 = -t1351 * t1065 + t1346 * t915;
t1361 = pkin(2) * t911 + t1368;
t1354 = qJD(2) ^ 2;
t1330 = -t1339 - t1354;
t1329 = t1339 - t1354;
t1328 = -t1354 - t1427;
t1327 = t1354 - t1427;
t1322 = -t1339 + t1427;
t1321 = t1339 + t1427;
t1319 = t1353 * t1355 + t1403;
t1318 = t1405 * qJDD(1);
t1317 = t1338 - 0.2e1 * t1396;
t1314 = 0.2e1 * t1335 + t1404;
t1312 = t1352 * t1323;
t1311 = t1405 * t1438;
t1302 = -pkin(5) * t1319 + t1353 * g(3);
t1294 = -t1304 + t1446;
t1293 = t1303 - t1446;
t1292 = t1352 * t1315 - t1342 * t1438;
t1291 = -t1347 * t1316 - t1343 * t1438;
t1287 = -t1347 * t1328 - t1409;
t1286 = -t1347 * t1327 + t1312;
t1285 = t1352 * t1330 - t1419;
t1284 = t1352 * t1329 - t1418;
t1283 = t1352 * t1328 - t1418;
t1282 = t1352 * t1327 + t1419;
t1281 = t1347 * t1330 + t1312;
t1280 = t1347 * t1329 + t1409;
t1279 = (t1315 + t1335) * t1347;
t1278 = (t1316 - t1396) * t1352;
t1271 = -t1347 * t1314 + t1352 * t1317;
t1270 = t1352 * t1314 + t1347 * t1317;
t1268 = t1304 - t1303;
t1260 = -t1274 + t1300;
t1259 = -t1300 + t1447;
t1258 = (-t1305 * t1351 + t1307 * t1346) * t1341;
t1257 = (-t1305 * t1346 - t1307 * t1351) * t1341;
t1256 = -pkin(6) * t1283 - t1410;
t1255 = -pkin(6) * t1281 - t1420;
t1252 = -t1303 - t1304;
t1249 = -pkin(1) * t1283 + t1290;
t1248 = -pkin(1) * t1281 + t1289;
t1240 = t1274 - t1447;
t1236 = pkin(1) * t1317 + pkin(6) * t1285 + t1410;
t1235 = -pkin(1) * t1314 + pkin(6) * t1287 - t1420;
t1228 = t1351 * t1293 - t1422;
t1227 = -t1346 * t1294 + t1452;
t1226 = t1346 * t1293 + t1412;
t1225 = t1351 * t1294 + t1455;
t1224 = -t1346 * t1288 - t1412;
t1217 = (qJD(3) + t1341) * t1307 + t1387;
t1215 = pkin(1) * t1308 + pkin(6) * t1243;
t1214 = t1351 * t1251 - t1307 * t1429;
t1213 = t1346 * t1251 + t1307 * t1428;
t1212 = -t1346 * t1250 + t1305 * t1428;
t1211 = t1351 * t1250 + t1305 * t1429;
t1210 = pkin(1) * t1321 + pkin(6) * t1318 + t1243;
t1207 = t1351 * t1263 - t1455;
t1205 = -t1231 + t1297;
t1204 = t1230 - t1297;
t1197 = (-t1276 * t1350 + t1277 * t1345) * t1301;
t1196 = (-t1276 * t1345 - t1277 * t1350) * t1301;
t1195 = -t1347 * t1257 + t1352 * t1258;
t1194 = t1352 * t1257 + t1347 * t1258;
t1185 = -pkin(7) * t1223 - t1413;
t1184 = t1231 - t1230;
t1181 = -t1261 + t1201;
t1178 = -t1262 + t1358;
t1177 = -pkin(7) * t1206 - t1423;
t1174 = -t1347 * t1226 + t1352 * t1228;
t1173 = -t1347 * t1225 + t1352 * t1227;
t1172 = t1352 * t1226 + t1347 * t1228;
t1171 = t1352 * t1225 + t1347 * t1227;
t1170 = t1350 * t1201 - t1277 * t1432;
t1169 = t1345 * t1201 + t1277 * t1431;
t1168 = t1276 * t1431 + t1345 * t1358;
t1167 = -t1276 * t1432 + t1350 * t1358;
t1165 = -t1347 * t1223 + t1352 * t1224;
t1164 = t1352 * t1223 + t1347 * t1224;
t1163 = -t1218 * t1351 + t1346 * t1222;
t1162 = -t1351 * t1217 - t1346 * t1221;
t1160 = -t1346 * t1217 + t1351 * t1221;
t1159 = t1351 * t1197 + t1346 * t1247;
t1158 = t1346 * t1197 - t1351 * t1247;
t1157 = t1350 * t1259 - t1424;
t1156 = -t1345 * t1260 + t1453;
t1155 = t1345 * t1259 + t1414;
t1154 = t1350 * t1260 + t1456;
t1153 = -t1347 * t1213 + t1352 * t1214;
t1152 = -t1347 * t1211 + t1352 * t1212;
t1151 = t1352 * t1213 + t1347 * t1214;
t1150 = t1352 * t1211 + t1347 * t1212;
t1149 = (-t1232 * t1349 + t1234 * t1344) * t1299;
t1148 = (-t1232 * t1344 - t1234 * t1349) * t1299;
t1147 = -t1347 * t1206 + t1352 * t1207;
t1146 = t1352 * t1206 + t1347 * t1207;
t1141 = t1350 * t1229 - t1424;
t1131 = t1345 * t1216 + t1453;
t1127 = -pkin(2) * t1221 + pkin(7) * t1224 - t1423;
t1124 = t1351 * t1170 + t1397;
t1123 = t1351 * t1168 - t1397;
t1122 = t1346 * t1170 - t1393;
t1121 = t1346 * t1168 + t1393;
t1120 = -pkin(2) * t1217 + pkin(7) * t1207 + t1413;
t1116 = -t1234 * qJD(5) - t1389;
t1115 = t1349 * t1204 - t1425;
t1114 = -t1344 * t1205 + t1454;
t1113 = t1344 * t1204 + t1415;
t1112 = t1349 * t1205 + t1457;
t1110 = t1350 * t1179 - t1345 * t1181;
t1109 = t1345 * t1180 - t1350 * t1182;
t1108 = t1345 * t1179 + t1350 * t1181;
t1104 = pkin(2) * t1253 + pkin(7) * t1119;
t1103 = -t1347 * t1161 + t1352 * t1163;
t1102 = -t1347 * t1160 + t1352 * t1162;
t1101 = t1352 * t1161 + t1347 * t1163;
t1100 = t1352 * t1160 + t1347 * t1162;
t1099 = t1351 * t1157 - t1346 * t1178;
t1098 = t1351 * t1156 + t1346 * t1182;
t1097 = t1346 * t1157 + t1351 * t1178;
t1096 = t1346 * t1156 - t1351 * t1182;
t1095 = -t1347 * t1158 + t1352 * t1159;
t1094 = t1352 * t1158 + t1347 * t1159;
t1093 = t1351 * t1142 - t1346 * t1183;
t1089 = t1351 * t1132 - t1346 * t1179;
t1079 = -pkin(1) * t1164 - t1378;
t1078 = t1349 * t1117 - t1234 * t1434;
t1077 = t1344 * t1117 + t1234 * t1433;
t1076 = -t1344 * t1116 + t1232 * t1433;
t1075 = t1349 * t1116 + t1232 * t1434;
t1074 = t1351 * t1110 + t1346 * t1240;
t1073 = t1346 * t1110 - t1351 * t1240;
t1072 = -t1345 * t1148 + t1350 * t1149;
t1071 = t1350 * t1148 + t1345 * t1149;
t1070 = -pkin(8) * t1141 + t1135;
t1069 = t1351 * t1111 - t1346 * t1203;
t1067 = -pkin(1) * t1146 - t1371;
t1066 = -pkin(8) * t1131 + t1134;
t1059 = -pkin(7) * t1161 - t1118;
t1058 = t1351 * t1072 + t1346 * t1246;
t1057 = t1346 * t1072 - t1351 * t1246;
t1056 = -t1347 * t1122 + t1352 * t1124;
t1055 = -t1347 * t1121 + t1352 * t1123;
t1054 = t1352 * t1122 + t1347 * t1124;
t1053 = t1352 * t1121 + t1347 * t1123;
t1052 = -pkin(2) * t1252 + pkin(7) * t1163 + t1119;
t1051 = -pkin(1) * t1101 - t1443;
t1050 = -pkin(6) * t1164 - t1347 * t1127 + t1352 * t1185;
t1049 = t1352 * t1119 - t1421;
t1048 = t1347 * t1119 + t1411;
t1047 = -pkin(6) * t1146 - t1347 * t1120 + t1352 * t1177;
t1046 = -t1345 * t1113 + t1350 * t1115;
t1045 = -t1345 * t1112 + t1350 * t1114;
t1044 = t1350 * t1113 + t1345 * t1115;
t1043 = t1350 * t1112 + t1345 * t1114;
t1042 = -pkin(1) * t1221 + pkin(6) * t1165 + t1352 * t1127 + t1347 * t1185;
t1041 = -pkin(3) * t1141 + t1063;
t1039 = -pkin(3) * t1131 + t1062;
t1037 = t1350 * t1105 + t1345 * t1106;
t1035 = -pkin(1) * t1217 + pkin(6) * t1147 + t1352 * t1120 + t1347 * t1177;
t1034 = -t1347 * t1097 + t1352 * t1099;
t1033 = -t1347 * t1096 + t1352 * t1098;
t1032 = t1352 * t1097 + t1347 * t1099;
t1031 = t1352 * t1096 + t1347 * t1098;
t1028 = -t1347 * t1092 + t1352 * t1093;
t1027 = t1352 * t1092 + t1347 * t1093;
t1025 = t1350 * t1090 + t1345 * t1091;
t1023 = -t1347 * t1088 + t1352 * t1089;
t1022 = t1352 * t1088 + t1347 * t1089;
t1019 = -t1349 * t1082 - t1344 * t1450;
t1017 = -t1344 * t1082 + t1349 * t1450;
t1015 = -t1345 * t1077 + t1350 * t1078;
t1014 = -t1345 * t1075 + t1350 * t1076;
t1013 = t1350 * t1077 + t1345 * t1078;
t1012 = t1350 * t1075 + t1345 * t1076;
t1011 = -t1347 * t1073 + t1352 * t1074;
t1010 = t1352 * t1073 + t1347 * t1074;
t1007 = -t1347 * t1068 + t1352 * t1069;
t1006 = t1352 * t1068 + t1347 * t1069;
t1005 = -pkin(1) * t1048 - t1444;
t1001 = -t1347 * t1057 + t1352 * t1058;
t1000 = t1352 * t1057 + t1347 * t1058;
t999 = t1351 * t1015 + t1398;
t998 = t1351 * t1014 - t1398;
t997 = t1346 * t1015 - t1394;
t996 = t1346 * t1014 + t1394;
t995 = t1351 * t1046 - t1346 * t1083;
t994 = t1351 * t1045 + t1346 * t1086;
t993 = t1346 * t1046 + t1351 * t1083;
t992 = t1346 * t1045 - t1351 * t1086;
t991 = t1351 * t1038 + t1346 * t1450;
t989 = t1351 * t1004 + t1346 * t1138;
t987 = t1351 * t1026 + t1346 * t1082;
t981 = -pkin(8) * t1109 - t1003;
t980 = -pkin(6) * t1048 - pkin(7) * t1411 - t1347 * t1104;
t979 = pkin(1) * t1253 + pkin(6) * t1049 - pkin(7) * t1421 + t1352 * t1104;
t978 = -pkin(6) * t1101 - t1347 * t1052 + t1352 * t1059;
t977 = -pkin(1) * t1252 + pkin(6) * t1103 + t1352 * t1052 + t1347 * t1059;
t973 = -pkin(7) * t1092 - t1346 * t1041 + t1351 * t1070;
t972 = -pkin(7) * t1088 - t1346 * t1039 + t1351 * t1066;
t970 = -t1345 * t1017 + t1350 * t1019;
t969 = t1350 * t1018 + t1345 * t1020;
t968 = t1350 * t1017 + t1345 * t1019;
t966 = -pkin(2) * t1141 + pkin(7) * t1093 + t1351 * t1041 + t1346 * t1070;
t965 = -pkin(1) * t1027 - t1374;
t964 = -pkin(2) * t1131 + pkin(7) * t1089 + t1351 * t1039 + t1346 * t1066;
t963 = t1346 * t1184 + t1351 * t970;
t962 = -t1351 * t1184 + t1346 * t970;
t961 = -pkin(1) * t1022 - t1375;
t960 = t1346 * t1133 + t1351 * t971;
t958 = -t1347 * t997 + t1352 * t999;
t957 = -t1347 * t996 + t1352 * t998;
t956 = t1347 * t999 + t1352 * t997;
t955 = t1347 * t998 + t1352 * t996;
t954 = -pkin(7) * t1068 + t1109 * t1442 + t1351 * t981;
t953 = -t1347 * t993 + t1352 * t995;
t952 = -t1347 * t992 + t1352 * t994;
t951 = t1347 * t995 + t1352 * t993;
t950 = t1347 * t994 + t1352 * t992;
t949 = -t1347 * t990 + t1352 * t991;
t948 = t1347 * t991 + t1352 * t990;
t947 = -t1347 * t988 + t1352 * t989;
t946 = t1347 * t989 + t1352 * t988;
t945 = pkin(7) * t1069 + t1109 * t1400 + t1346 * t981;
t944 = -t1347 * t986 + t1352 * t987;
t943 = t1347 * t987 + t1352 * t986;
t942 = -pkin(3) * t969 - t1441;
t941 = -pkin(3) * t1037 - t1362;
t940 = -pkin(3) * t1025 - t1370;
t939 = -pkin(1) * t1006 - t1365;
t938 = -pkin(8) * t1037 + t1350 * t1021 - t1345 * t985;
t935 = -pkin(8) * t1025 + t1350 * t1009 - t1345 * t983;
t934 = -pkin(7) * t988 + (-pkin(8) * t1351 + t1442) * t1003;
t931 = -t1347 * t962 + t1352 * t963;
t930 = t1347 * t963 + t1352 * t962;
t927 = -t1347 * t959 + t1352 * t960;
t926 = t1347 * t960 + t1352 * t959;
t923 = pkin(7) * t989 + (-pkin(8) * t1346 + t1400) * t1003;
t922 = -pkin(6) * t1027 - t1347 * t966 + t1352 * t973;
t921 = -pkin(6) * t1022 - t1347 * t964 + t1352 * t972;
t920 = -pkin(1) * t1141 + pkin(6) * t1028 + t1347 * t973 + t1352 * t966;
t919 = -pkin(1) * t1131 + pkin(6) * t1023 + t1347 * t972 + t1352 * t964;
t918 = -pkin(1) * t946 - t1379;
t917 = -pkin(6) * t1006 - t1347 * t945 + t1352 * t954;
t916 = -pkin(1) * t1109 + pkin(6) * t1007 + t1347 * t954 + t1352 * t945;
t914 = t1345 * t937 + t1439;
t912 = t1346 * t1065 + t1351 * t915;
t910 = -pkin(7) * t990 - t1346 * t941 + t1351 * t938;
t909 = -pkin(1) * t948 - t1366;
t908 = -pkin(7) * t986 - t1346 * t940 + t1351 * t935;
t907 = -pkin(1) * t943 - t1367;
t906 = -pkin(2) * t1037 + pkin(7) * t991 + t1346 * t938 + t1351 * t941;
t905 = -pkin(2) * t1025 + pkin(7) * t987 + t1346 * t935 + t1351 * t940;
t904 = -pkin(3) * t914 - t1445;
t903 = -pkin(8) * t969 - t1345 * t925 + t1350 * t929;
t902 = -pkin(6) * t946 - t1347 * t923 + t1352 * t934;
t901 = -pkin(1) * t1003 + pkin(6) * t947 + t1347 * t934 + t1352 * t923;
t900 = -pkin(8) * t914 - pkin(9) * t1439 - t1345 * t933;
t899 = -t1347 * t911 + t1352 * t912;
t898 = t1347 * t912 + t1352 * t911;
t897 = -pkin(7) * t959 - t1346 * t942 + t1351 * t903;
t896 = -pkin(1) * t926 - t1369;
t895 = -pkin(2) * t969 + pkin(7) * t960 + t1346 * t903 + t1351 * t942;
t894 = -pkin(6) * t948 - t1347 * t906 + t1352 * t910;
t893 = -pkin(1) * t1037 + pkin(6) * t949 + t1347 * t910 + t1352 * t906;
t892 = -pkin(6) * t943 - t1347 * t905 + t1352 * t908;
t891 = -pkin(1) * t1025 + pkin(6) * t944 + t1347 * t908 + t1352 * t905;
t890 = -pkin(7) * t911 - t1346 * t904 + t1351 * t900;
t889 = -pkin(1) * t898 - t1361;
t888 = -pkin(6) * t926 - t1347 * t895 + t1352 * t897;
t887 = -pkin(2) * t914 + pkin(7) * t912 + t1346 * t900 + t1351 * t904;
t886 = -pkin(1) * t969 + pkin(6) * t927 + t1347 * t897 + t1352 * t895;
t885 = -pkin(6) * t898 - t1347 * t887 + t1352 * t890;
t884 = -pkin(1) * t914 + pkin(6) * t899 + t1347 * t890 + t1352 * t887;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1320, 0, -t1319, 0, t1377, -t1302, -t1376, -pkin(5) * t1376, t1353 * t1292 - t1382, t1353 * t1271 + t1348 * t1322, t1353 * t1286 + t1347 * t1403, t1353 * t1291 + t1382, t1353 * t1284 + t1338 * t1348, t1348 * qJDD(2) + t1353 * t1311, t1353 * t1255 - t1348 * t1248 - pkin(5) * (t1348 * t1285 + t1353 * t1317), t1353 * t1256 - t1348 * t1249 - pkin(5) * (t1348 * t1287 - t1353 * t1314), t1353 * t1242 - pkin(5) * (t1348 * t1318 + t1353 * t1321), -pkin(5) * (t1348 * t1243 + t1353 * t1308) - (t1348 * pkin(1) - t1353 * pkin(6)) * t1242, t1353 * t1153 + t1395, t1353 * t1102 + t1348 * t1268, t1353 * t1173 + t1348 * t1222, t1353 * t1152 - t1395, t1353 * t1174 - t1348 * t1218, t1353 * t1195 + t1348 * t1401, t1353 * t1047 - t1348 * t1067 - pkin(5) * (t1348 * t1147 - t1353 * t1217), t1353 * t1050 - t1348 * t1079 - pkin(5) * (t1348 * t1165 - t1221 * t1353), t1353 * t978 - t1348 * t1051 - pkin(5) * (t1348 * t1103 - t1353 * t1252), t1353 * t980 - t1348 * t1005 - pkin(5) * (t1348 * t1049 + t1353 * t1253), t1353 * t1056 + t1348 * t1169, t1353 * t1011 + t1348 * t1108, t1353 * t1033 + t1348 * t1154, t1353 * t1055 - t1348 * t1167, t1353 * t1034 + t1348 * t1155, t1353 * t1095 + t1348 * t1196, t1353 * t921 - t1348 * t961 - pkin(5) * (t1348 * t1023 - t1353 * t1131), t1353 * t922 - t1348 * t965 - pkin(5) * (t1348 * t1028 - t1353 * t1141), t1353 * t917 - t1348 * t939 - pkin(5) * (t1348 * t1007 - t1353 * t1109), t1353 * t902 - t1348 * t918 - pkin(5) * (-t1353 * t1003 + t1348 * t947), t1348 * t1013 + t1353 * t958, t1348 * t968 + t1353 * t931, t1348 * t1043 + t1353 * t952, t1348 * t1012 + t1353 * t957, t1348 * t1044 + t1353 * t953, t1353 * t1001 + t1348 * t1071, t1353 * t892 - t1348 * t907 - pkin(5) * (-t1353 * t1025 + t1348 * t944), t1353 * t894 - t1348 * t909 - pkin(5) * (-t1353 * t1037 + t1348 * t949), t1353 * t888 - t1348 * t896 - pkin(5) * (t1348 * t927 - t1353 * t969), t1353 * t885 - t1348 * t889 - pkin(5) * (t1348 * t899 - t1353 * t914); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1319, 0, t1320, 0, t1302, t1377, t1386, pkin(5) * t1386, t1348 * t1292 + t1381, t1348 * t1271 - t1353 * t1322, t1348 * t1286 - t1347 * t1402, t1348 * t1291 - t1381, t1348 * t1284 - t1352 * t1402, -t1353 * qJDD(2) + t1348 * t1311, t1348 * t1255 + t1353 * t1248 + pkin(5) * (t1353 * t1285 - t1348 * t1317), t1348 * t1256 + t1353 * t1249 + pkin(5) * (t1353 * t1287 + t1348 * t1314), t1348 * t1242 + pkin(5) * (t1353 * t1318 - t1348 * t1321), pkin(5) * (t1353 * t1243 - t1348 * t1308) - (-t1353 * pkin(1) - t1348 * pkin(6)) * t1242, t1348 * t1153 - t1392, t1348 * t1102 - t1353 * t1268, t1348 * t1173 - t1353 * t1222, t1348 * t1152 + t1392, t1348 * t1174 + t1353 * t1218, t1348 * t1195 - t1353 * t1401, t1348 * t1047 + t1353 * t1067 + pkin(5) * (t1353 * t1147 + t1348 * t1217), t1348 * t1050 + t1353 * t1079 + pkin(5) * (t1353 * t1165 + t1221 * t1348), t1348 * t978 + t1353 * t1051 + pkin(5) * (t1353 * t1103 + t1348 * t1252), t1348 * t980 + t1353 * t1005 + pkin(5) * (t1353 * t1049 - t1348 * t1253), t1348 * t1056 - t1353 * t1169, t1348 * t1011 - t1353 * t1108, t1348 * t1033 - t1353 * t1154, t1348 * t1055 + t1353 * t1167, t1348 * t1034 - t1353 * t1155, t1348 * t1095 - t1353 * t1196, t1348 * t921 + t1353 * t961 + pkin(5) * (t1353 * t1023 + t1348 * t1131), t1348 * t922 + t1353 * t965 + pkin(5) * (t1353 * t1028 + t1348 * t1141), t1348 * t917 + t1353 * t939 + pkin(5) * (t1353 * t1007 + t1348 * t1109), t1348 * t902 + t1353 * t918 + pkin(5) * (t1348 * t1003 + t1353 * t947), -t1353 * t1013 + t1348 * t958, t1348 * t931 - t1353 * t968, -t1353 * t1043 + t1348 * t952, -t1353 * t1012 + t1348 * t957, -t1353 * t1044 + t1348 * t953, t1348 * t1001 - t1353 * t1071, t1348 * t892 + t1353 * t907 + pkin(5) * (t1348 * t1025 + t1353 * t944), t1348 * t894 + t1353 * t909 + pkin(5) * (t1348 * t1037 + t1353 * t949), t1348 * t888 + t1353 * t896 + pkin(5) * (t1348 * t969 + t1353 * t927), t1348 * t885 + t1353 * t889 + pkin(5) * (t1348 * t914 + t1353 * t899); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1325, t1326, 0, 0, t1279, t1270, t1282, t1278, t1280, 0, t1236, t1235, t1210, t1215, t1151, t1100, t1171, t1150, t1172, t1194, t1035, t1042, t977, t979, t1054, t1010, t1031, t1053, t1032, t1094, t919, t920, t916, t901, t956, t930, t950, t955, t951, t1000, t891, t893, t886, t884; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1355, 0, 0, -g(3), -t1325, 0, t1292, t1271, t1286, t1291, t1284, t1311, t1255, t1256, t1242, pkin(6) * t1242, t1153, t1102, t1173, t1152, t1174, t1195, t1047, t1050, t978, t980, t1056, t1011, t1033, t1055, t1034, t1095, t921, t922, t917, t902, t958, t931, t952, t957, t953, t1001, t892, t894, t888, t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1355, 0, qJDD(1), 0, g(3), 0, -t1326, 0, t1332, -t1322, -t1404, -t1332, -t1338, -qJDD(2), t1248, t1249, 0, pkin(1) * t1242, -t1269, -t1268, -t1222, t1269, t1218, -t1401, t1067, t1079, t1051, t1005, -t1169, -t1108, -t1154, t1167, -t1155, -t1196, t961, t965, t939, t918, -t1013, -t968, -t1043, -t1012, -t1044, -t1071, t907, t909, t896, t889; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1325, t1326, 0, 0, t1279, t1270, t1282, t1278, t1280, 0, t1236, t1235, t1210, t1215, t1151, t1100, t1171, t1150, t1172, t1194, t1035, t1042, t977, t979, t1054, t1010, t1031, t1053, t1032, t1094, t919, t920, t916, t901, t956, t930, t950, t955, t951, t1000, t891, t893, t886, t884; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1315, t1317, t1323, -t1335, t1329, t1335, 0, -t1308, t1289, 0, t1214, t1162, t1227, t1212, t1228, t1258, t1177, t1185, t1059, -pkin(7) * t1118, t1124, t1074, t1098, t1123, t1099, t1159, t972, t973, t954, t934, t999, t963, t994, t998, t995, t1058, t908, t910, t897, t890; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1396, t1314, t1327, t1316, t1324, -t1396, t1308, 0, t1290, 0, t1213, t1160, t1225, t1211, t1226, t1257, t1120, t1127, t1052, t1104, t1122, t1073, t1096, t1121, t1097, t1158, t964, t966, t945, t923, t997, t962, t992, t996, t993, t1057, t905, t906, t895, t887; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1332, t1322, t1404, t1332, t1338, qJDD(2), -t1289, -t1290, 0, 0, t1269, t1268, t1222, -t1269, -t1218, t1401, t1371, t1378, t1443, t1444, t1169, t1108, t1154, -t1167, t1155, t1196, t1375, t1374, t1365, t1379, t1013, t968, t1043, t1012, t1044, t1071, t1367, t1366, t1369, t1361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1251, -t1217, t1448, t1296, t1293, -t1296, 0, -t1253, t1188, 0, t1170, t1110, t1156, t1168, t1157, t1197, t1066, t1070, t981, -pkin(8) * t1003, t1015, t970, t1045, t1014, t1046, t1072, t935, t938, t903, t900; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1430, t1221, t1294, t1250, t1265, -t1430, t1253, 0, t1189, 0, -t1245, -t1240, -t1182, t1245, t1178, -t1247, t1039, t1041, -pkin(3) * t1109, -pkin(3) * t1003, -t1186, -t1184, -t1086, t1186, t1083, -t1246, t940, t941, t942, t904; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1269, t1268, t1222, -t1269, -t1218, t1401, -t1188, -t1189, 0, 0, t1169, t1108, t1154, -t1167, t1155, t1196, t1390, t1391, t1380, t1406, t1013, t968, t1043, t1012, t1044, t1071, t1384, t1383, t1385, t1368; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1201, t1179, t1449, t1261, t1259, -t1261, 0, t1138, t1062, 0, t1078, t1019, t1114, t1076, t1115, t1149, t1009, t1021, t929, -pkin(9) * t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, t1181, t1260, -t1358, t1191, -t1262, -t1138, 0, t1063, 0, t1077, t1017, t1112, t1075, t1113, t1148, t983, t985, t925, t933; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1245, t1240, t1182, -t1245, -t1178, t1247, -t1062, -t1063, 0, 0, t1186, t1184, t1086, -t1186, -t1083, t1246, t1370, t1362, t1441, t1445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1117, -t1082, t1451, t1208, t1204, -t1208, 0, t1065, t975, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1435, t1450, t1205, t1116, t1144, -t1435, -t1065, 0, t976, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1186, t1184, t1086, -t1186, -t1083, t1246, -t975, -t976, 0, 0;];
m_new_reg = t1;
