% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPRP12_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:34:43
% EndTime: 2019-05-07 09:36:27
% DurationCPUTime: 65.44s
% Computational Cost: add. (145396->872), mult. (314263->1284), div. (0->0), fcn. (244691->10), ass. (0->668)
t1170 = sin(qJ(1));
t1174 = cos(qJ(1));
t1165 = sin(pkin(6));
t1166 = cos(pkin(6));
t1169 = sin(qJ(2));
t1173 = cos(qJ(2));
t1168 = sin(qJ(3));
t1172 = cos(qJ(3));
t1341 = qJD(1) * t1166;
t1289 = qJD(2) + t1341;
t1315 = t1165 * t1169;
t1302 = qJD(1) * t1315;
t1118 = t1168 * t1302 - t1172 * t1289;
t1305 = t1169 * qJDD(1);
t1340 = qJD(1) * t1173;
t1129 = (qJD(2) * t1340 + t1305) * t1165;
t1157 = qJDD(1) * t1166 + qJDD(2);
t1061 = -t1118 * qJD(3) + t1172 * t1129 + t1168 * t1157;
t1314 = t1165 * t1173;
t1150 = qJD(1) * t1314 - qJD(3);
t1325 = t1118 * t1150;
t1393 = t1061 + t1325;
t1146 = t1150 ^ 2;
t1120 = t1168 * t1289 + t1172 * t1302;
t1369 = t1120 ^ 2;
t1071 = t1369 + t1146;
t1306 = qJDD(1) * t1165;
t1130 = -qJD(2) * t1302 + t1173 * t1306;
t1124 = -qJDD(3) + t1130;
t1326 = t1118 * t1120;
t1191 = t1124 - t1326;
t1396 = t1191 * t1172;
t972 = t1071 * t1168 + t1396;
t1248 = t1169 * t972 - t1173 * t1393;
t1397 = t1191 * t1168;
t971 = t1071 * t1172 - t1397;
t836 = t1165 * t971 + t1166 * t1248;
t886 = t1169 * t1393 + t1173 * t972;
t764 = t1170 * t886 + t1174 * t836;
t1579 = pkin(7) * t764;
t767 = t1170 * t836 - t1174 * t886;
t1578 = pkin(7) * t767;
t1167 = sin(qJ(5));
t1171 = cos(qJ(5));
t1082 = t1118 * t1167 - t1150 * t1171;
t1286 = -t1168 * t1129 + t1172 * t1157;
t1060 = qJD(3) * t1120 - t1286;
t1288 = t1171 * t1060 + t1167 * t1124;
t1113 = qJD(5) + t1120;
t1309 = qJD(5) - t1113;
t1196 = -t1082 * t1309 + t1288;
t1110 = t1113 ^ 2;
t1080 = -t1171 * t1118 - t1150 * t1167;
t1371 = t1080 ^ 2;
t1039 = t1371 - t1110;
t1026 = t1082 * t1080;
t1057 = qJDD(5) + t1061;
t1394 = t1026 + t1057;
t1346 = t1171 * t1394;
t891 = t1039 * t1167 + t1346;
t812 = t1168 * t891 + t1172 * t1196;
t1353 = t1167 * t1394;
t895 = t1039 * t1171 - t1353;
t1260 = t1169 * t812 - t1173 * t895;
t809 = -t1168 * t1196 + t1172 * t891;
t719 = t1165 * t809 + t1166 * t1260;
t754 = t1169 * t895 + t1173 * t812;
t1577 = t1170 * t719 - t1174 * t754;
t1576 = t1170 * t754 + t1174 * t719;
t834 = t1165 * t1248 - t1166 * t971;
t1575 = pkin(8) * (t1165 * t834 + t1166 * t836);
t1574 = pkin(1) * t834;
t1573 = pkin(1) * t836;
t1392 = t1325 - t1061;
t1087 = -t1369 + t1146;
t1190 = t1124 + t1326;
t1401 = t1172 * t1190;
t979 = t1087 * t1168 + t1401;
t1249 = t1169 * t979 - t1173 * t1392;
t1405 = t1168 * t1190;
t976 = t1087 * t1172 - t1405;
t844 = t1165 * t976 + t1166 * t1249;
t900 = t1169 * t1392 + t1173 * t979;
t1570 = t1170 * t844 - t1174 * t900;
t1569 = t1170 * t900 + t1174 * t844;
t1568 = t1165 * t1260 - t1166 * t809;
t1566 = pkin(8) * t886;
t1390 = -t1371 - t1110;
t1395 = -t1026 + t1057;
t946 = t1167 * t1395;
t1424 = t1171 * t1390 - t946;
t1308 = qJD(5) + t1113;
t1197 = t1082 * t1308 - t1288;
t1193 = t1171 * t1395;
t1419 = t1167 * t1390 + t1193;
t1457 = t1168 * t1419 + t1172 * t1197;
t1488 = t1169 * t1424 + t1173 * t1457;
t1454 = t1168 * t1197 - t1172 * t1419;
t1492 = t1169 * t1457 - t1173 * t1424;
t1508 = -t1165 * t1454 + t1166 * t1492;
t1547 = t1170 * t1488 + t1174 * t1508;
t1565 = pkin(7) * t1547;
t1240 = t1167 * t1060 - t1171 * t1124;
t1188 = t1080 * t1309 - t1240;
t1387 = -t1167 * t1188 + t1171 * t1196;
t1384 = t1167 * t1196 + t1171 * t1188;
t1079 = t1082 ^ 2;
t967 = -t1371 - t1079;
t1462 = t1168 * t1384 + t1172 * t967;
t1489 = t1169 * t1387 + t1173 * t1462;
t1463 = t1168 * t967 - t1172 * t1384;
t1493 = t1169 * t1462 - t1173 * t1387;
t1507 = -t1165 * t1463 + t1166 * t1493;
t1548 = t1170 * t1489 + t1174 * t1507;
t1564 = pkin(7) * t1548;
t1549 = -t1170 * t1508 + t1174 * t1488;
t1563 = pkin(7) * t1549;
t1550 = -t1170 * t1507 + t1174 * t1489;
t1562 = pkin(7) * t1550;
t1555 = t1165 * t1249 - t1166 * t976;
t1554 = pkin(1) * t1507;
t1553 = pkin(1) * t1508;
t1509 = t1165 * t1493 + t1166 * t1463;
t1552 = pkin(1) * t1509;
t1510 = t1165 * t1492 + t1166 * t1454;
t1551 = pkin(1) * t1510;
t1546 = (-t1165 * t1509 - t1166 * t1507) * pkin(8);
t1545 = (-t1165 * t1510 - t1166 * t1508) * pkin(8);
t1303 = t1110 + t1079;
t878 = t1167 * t1303 - t1346;
t1544 = pkin(2) * t878;
t1370 = t1118 ^ 2;
t1064 = -t1146 - t1370;
t959 = t1064 * t1168 - t1401;
t1543 = pkin(2) * t959;
t1542 = pkin(2) * t971;
t1541 = pkin(9) * t959;
t961 = t1064 * t1172 + t1405;
t1540 = pkin(9) * t961;
t1539 = pkin(9) * t971;
t1538 = pkin(9) * t972;
t1537 = pkin(8) * t1488;
t1536 = pkin(8) * t1489;
t1048 = t1080 * t1113;
t965 = -qJD(5) * t1080 + t1240;
t1482 = t1048 - t965;
t1535 = qJ(6) * t1482;
t1534 = t1165 * t959;
t1086 = t1370 - t1146;
t977 = t1086 * t1168 - t1396;
t1531 = t1165 * t977;
t1530 = t1166 * t959;
t1527 = t1166 * t977;
t869 = t1171 * t1303 + t1353;
t1526 = t1168 * t869;
t1525 = t1169 * t878;
t1523 = t1169 * t961;
t981 = t1086 * t1172 + t1397;
t1522 = t1169 * t981;
t1521 = t1172 * t869;
t1520 = t1173 * t878;
t1518 = t1173 * t961;
t1517 = t1173 * t981;
t1366 = pkin(3) + pkin(10);
t1516 = t1366 * t869;
t1515 = t1366 * t878;
t1506 = -pkin(4) * t869 - qJ(4) * t878;
t1368 = -2 * qJD(4);
t1505 = pkin(9) * t1454;
t1504 = pkin(9) * t1463;
t1040 = -t1079 + t1110;
t1466 = t1167 * t1040 - t1193;
t1499 = t1169 * t1466;
t1498 = t1173 * t1466;
t1497 = -pkin(2) * t1454 + t1366 * t1419;
t1496 = -pkin(2) * t1463 - qJ(4) * t967 + t1366 * t1384;
t1495 = -pkin(2) * t1387 + pkin(9) * t1462;
t1494 = -pkin(2) * t1424 + pkin(9) * t1457;
t1327 = t1113 * t1171;
t964 = -qJD(5) * t1082 + t1288;
t1281 = -t1080 * t1327 + t1167 * t964;
t1328 = t1113 * t1167;
t1299 = t1080 * t1328;
t1227 = -t1171 * t964 - t1299;
t1300 = t1172 * t1026;
t1383 = -t1168 * t1227 - t1300;
t1416 = -t1169 * t1281 + t1173 * t1383;
t1301 = t1168 * t1026;
t1380 = t1172 * t1227 - t1301;
t1418 = t1169 * t1383 + t1173 * t1281;
t1458 = -t1165 * t1380 + t1166 * t1418;
t1491 = -t1170 * t1458 + t1174 * t1416;
t1218 = (t1080 * t1171 - t1082 * t1167) * t1113;
t1037 = t1082 * t1327;
t1280 = t1037 + t1299;
t1385 = t1172 * t1057 - t1168 * t1280;
t1415 = -t1169 * t1218 + t1173 * t1385;
t1386 = t1168 * t1057 + t1172 * t1280;
t1417 = t1169 * t1385 + t1173 * t1218;
t1459 = -t1165 * t1386 + t1166 * t1417;
t1490 = -t1170 * t1459 + t1174 * t1415;
t1487 = t1170 * t1416 + t1174 * t1458;
t1486 = t1170 * t1415 + t1174 * t1459;
t1036 = t1370 + t1369;
t1485 = pkin(2) * t1036;
t1483 = qJ(4) * t1393;
t1023 = t1079 - t1371;
t1481 = t1023 * t1168;
t1480 = t1023 * t1172;
t1479 = t1036 * t1169;
t1478 = t1036 * t1173;
t1449 = -t1040 * t1171 - t946;
t1477 = t1168 * t1449;
t1473 = t1172 * t1449;
t1468 = t1366 * t1424;
t1364 = pkin(2) * t1173;
t1284 = -pkin(9) * t1169 - t1364;
t1342 = qJD(1) * t1165;
t1128 = t1284 * t1342;
t1285 = t1289 ^ 2;
t1153 = g(1) * t1174 + g(2) * t1170;
t1175 = qJD(1) ^ 2;
t1125 = -pkin(1) * t1175 + pkin(8) * t1306 - t1153;
t1152 = t1170 * g(1) - t1174 * g(2);
t1359 = pkin(8) * t1165;
t1189 = qJDD(1) * pkin(1) + t1175 * t1359 + t1152;
t1183 = t1166 * t1189;
t1287 = t1169 * t1125 - t1173 * t1183;
t991 = (qJD(1) * t1128 * t1169 + g(3) * t1173) * t1165 - t1157 * pkin(2) - t1285 * pkin(9) + t1287;
t1178 = t1060 * pkin(3) - t1483 + t991;
t1467 = t1120 * t1368 + t1178;
t1465 = pkin(4) * t967 - t1366 * t1387;
t1461 = t1165 * t1417 + t1166 * t1386;
t1460 = t1165 * t1418 + t1166 * t1380;
t1318 = t1150 * t1172;
t1219 = t1060 * t1168 - t1118 * t1318;
t1298 = t1169 * t1326;
t1378 = t1173 * t1219 - t1298;
t1319 = t1150 * t1168;
t1278 = -t1172 * t1060 - t1118 * t1319;
t1297 = t1173 * t1326;
t1382 = t1169 * t1219 + t1297;
t1420 = -t1165 * t1278 + t1166 * t1382;
t1456 = -t1170 * t1420 + t1174 * t1378;
t1277 = t1172 * t1061 + t1120 * t1319;
t1376 = t1173 * t1277 + t1298;
t1279 = t1168 * t1061 - t1120 * t1318;
t1381 = t1169 * t1277 - t1297;
t1421 = -t1165 * t1279 + t1166 * t1381;
t1455 = -t1170 * t1421 + t1174 * t1376;
t1453 = t1170 * t1378 + t1174 * t1420;
t1452 = t1170 * t1376 + t1174 * t1421;
t1451 = pkin(4) * t1419 - qJ(4) * t1424;
t1443 = t1168 * t1392;
t1442 = t1168 * t1393;
t1391 = -t1369 + t1370;
t1440 = t1169 * t1391;
t1434 = t1172 * t1392;
t1433 = t1172 * t1393;
t1431 = t1173 * t1391;
t1423 = t1165 * t1381 + t1166 * t1279;
t1422 = t1165 * t1382 + t1166 * t1278;
t738 = pkin(4) * t1384 - qJ(4) * t1387;
t1216 = (t1118 * t1168 + t1120 * t1172) * t1150;
t1217 = (t1118 * t1172 - t1120 * t1168) * t1150;
t1312 = t1166 * t1173;
t1313 = t1166 * t1169;
t1372 = t1124 * t1312 - t1165 * t1216 + t1217 * t1313;
t1379 = -t1124 * t1169 + t1173 * t1217;
t1414 = t1170 * t1379 + t1174 * t1372;
t1413 = -t1170 * t1372 + t1174 * t1379;
t915 = -t1167 * t965 - t1037;
t1228 = -t1168 * t915 + t1300;
t1282 = t1172 * t915 + t1301;
t916 = t1082 * t1328 - t1171 * t965;
t1375 = -t1165 * t1282 + t1228 * t1313 + t916 * t1312;
t1377 = -t1169 * t916 + t1173 * t1228;
t1412 = t1170 * t1377 + t1174 * t1375;
t1411 = -t1170 * t1375 + t1174 * t1377;
t1162 = t1165 ^ 2;
t1316 = t1162 * t1175;
t1241 = t1289 * qJD(1);
t1410 = t1162 * (-t1166 * t1175 + t1241);
t1389 = -t964 * pkin(5) + t1535;
t1068 = pkin(3) * t1118 - qJ(4) * t1120;
t1179 = -g(3) * t1315 + t1169 * t1183;
t992 = t1157 * pkin(9) - t1285 * pkin(2) + (t1128 * t1342 + t1125) * t1173 + t1179;
t1222 = t1169 * t1241;
t1223 = t1173 * t1241;
t1358 = t1166 * g(3);
t993 = -t1130 * pkin(2) - t1129 * pkin(9) - t1358 + (pkin(2) * t1222 - pkin(9) * t1223 - t1189) * t1165;
t1357 = -t1168 * t992 + t1172 * t993;
t1229 = t1124 * pkin(3) - t1146 * qJ(4) + t1120 * t1068 + qJDD(4) - t1357;
t786 = -pkin(4) * t1392 + t1190 * pkin(10) + t1229;
t1085 = pkin(4) * t1120 + pkin(10) * t1150;
t1292 = -pkin(3) * t1150 + t1368;
t794 = -t1370 * pkin(4) + t1060 * pkin(10) + (-t1085 + t1292) * t1120 + t1178;
t1290 = -t1167 * t794 + t1171 * t786;
t728 = t1167 * t786 + t1171 * t794;
t672 = -t1167 * t1290 + t1171 * t728;
t671 = t1167 * t728 + t1171 * t1290;
t1374 = t1166 * t1282 + t1228 * t1315 + t916 * t1314;
t1373 = t1124 * t1314 + t1166 * t1216 + t1217 * t1315;
t1367 = 2 * qJD(6);
t1365 = pkin(2) * t1169;
t1363 = pkin(3) * t1168;
t1362 = pkin(3) * t1172;
t1361 = pkin(5) * t1167;
t1360 = pkin(5) * t1171;
t899 = t1168 * t993 + t1172 * t992;
t1232 = t1146 * pkin(3) + t1118 * t1068 - t899;
t1195 = t1124 * qJ(4) + t1232;
t1184 = -t1060 * pkin(4) - pkin(10) * t1370 - t1195;
t793 = (t1368 - t1085) * t1150 + t1184;
t1355 = t1167 * t793;
t1354 = t1167 * t1197;
t1351 = t1168 * t991;
t1348 = t1171 * t793;
t1347 = t1171 * t1197;
t1344 = t1172 * t991;
t1343 = -t1110 - t967;
t1338 = qJD(4) * t1150;
t1098 = t1150 * t1120;
t1007 = t1060 - t1098;
t1337 = t1007 * t1172;
t1329 = t1082 * t1113;
t1149 = t1173 * t1169 * t1316;
t1126 = -t1149 + t1157;
t1323 = t1126 * t1169;
t1322 = t1126 * t1173;
t1127 = t1149 + t1157;
t1321 = t1127 * t1169;
t1320 = t1127 * t1173;
t1317 = t1157 * t1165;
t1100 = t1165 * t1189 + t1358;
t1311 = t1169 * t1100;
t1310 = t1173 * t1100;
t1163 = t1169 ^ 2;
t1164 = t1173 ^ 2;
t1307 = t1163 + t1164;
t1296 = t1163 * t1316;
t1295 = t1164 * t1316;
t1294 = qJ(4) * t1168 + pkin(2);
t1293 = qJ(6) * t1167 + pkin(4);
t1291 = qJ(6) * t1171 - qJ(4);
t806 = -t1168 * t1357 + t1172 * t899;
t1103 = -t1152 * t1170 - t1174 * t1153;
t1148 = qJDD(1) * t1174 - t1170 * t1175;
t1283 = -pkin(7) * t1148 - g(3) * t1170;
t1022 = pkin(5) * t1080 - qJ(6) * t1082;
t1276 = t1057 * qJ(6) - t1080 * t1022 + t1113 * t1367 + t728;
t1111 = -t1296 - t1285;
t1067 = -t1111 * t1169 - t1322;
t1275 = pkin(8) * t1067 - t1311;
t1134 = -t1285 - t1295;
t1075 = t1134 * t1173 - t1321;
t1274 = pkin(8) * t1075 + t1310;
t805 = t1168 * t899 + t1172 * t1357;
t1200 = t1082 * t1022 + qJDD(6) - t1290;
t1185 = -t1057 * pkin(5) + t1200;
t1180 = t1110 * qJ(6) - t1185;
t703 = -t1110 * pkin(5) + t1276;
t659 = t1167 * t703 + t1171 * t1180;
t731 = (pkin(5) * t1113 - (2 * qJD(6))) * t1082 + t793 + t1389;
t640 = t1168 * t659 + t1172 * t731;
t660 = -t1167 * t1180 + t1171 * t703;
t1273 = t1169 * t640 - t1173 * t660;
t658 = t1168 * t671 + t1172 * t793;
t1272 = t1169 * t658 - t1173 * t672;
t842 = -t1195 - 0.2e1 * t1338;
t751 = t1168 * t1229 + t1172 * t842;
t850 = t1120 * t1292 + t1178;
t1271 = t1169 * t751 - t1173 * t850;
t1187 = -t1080 * t1308 + t1240;
t819 = -t1171 * t1187 + t1354;
t783 = -t1168 * t819 + t1480;
t823 = t1167 * t1187 + t1347;
t1268 = t1169 * t783 + t1173 * t823;
t818 = -t1171 * t1482 - t1354;
t784 = -t1168 * t818 - t1480;
t822 = t1167 * t1482 - t1347;
t1267 = t1169 * t784 + t1173 * t822;
t798 = t1172 * t1482 + t1526;
t1265 = t1169 * t798 + t1520;
t804 = t1172 * t1187 - t1526;
t1263 = t1169 * t804 - t1520;
t1262 = t1169 * t806 - t1173 * t991;
t811 = -t1172 * t1188 - t1477;
t1261 = t1169 * t811 + t1498;
t925 = t1048 + t965;
t814 = t1172 * t925 - t1477;
t1258 = t1169 * t814 + t1498;
t1253 = t1007 * t1173 - t1523;
t1008 = (-qJD(3) - t1150) * t1120 + t1286;
t1252 = -t1008 * t1173 + t1522;
t1009 = -t1060 - t1098;
t1251 = t1009 * t1173 - t1522;
t1016 = (-qJD(3) + t1150) * t1120 + t1286;
t1246 = t1016 * t1173 + t1523;
t930 = t1008 * t1172 - t1443;
t1245 = t1169 * t930 + t1478;
t931 = t1009 * t1172 - t1443;
t1244 = t1169 * t931 + t1478;
t932 = t1016 * t1172 - t1442;
t1243 = t1169 * t932 + t1431;
t933 = -t1337 - t1442;
t1242 = t1169 * t933 + t1431;
t1065 = g(3) * t1314 + t1287;
t1066 = t1173 * t1125 + t1179;
t1239 = -t1173 * t1065 + t1169 * t1066;
t969 = t1065 * t1169 + t1066 * t1173;
t1138 = t1165 * t1223;
t1091 = t1138 + t1129;
t1137 = t1165 * t1222;
t1095 = t1130 - t1137;
t1238 = t1091 * t1173 + t1095 * t1169;
t1093 = -t1138 + t1129;
t1094 = t1130 + t1137;
t1237 = -t1093 * t1173 + t1094 * t1169;
t1236 = t1111 * t1173 - t1323;
t1133 = -t1285 + t1295;
t1235 = t1133 * t1169 + t1322;
t1132 = t1285 - t1296;
t1234 = t1132 * t1173 + t1321;
t1233 = t1134 * t1169 + t1320;
t1102 = t1152 * t1174 - t1153 * t1170;
t1224 = t1165 * t1241;
t613 = -t1366 * t660 + (t1293 + t1360) * t731;
t615 = pkin(4) * t659 + pkin(5) * t1180 - qJ(4) * t660 + qJ(6) * t703;
t639 = t1168 * t731 - t1172 * t659;
t593 = -pkin(9) * t639 - t1168 * t613 + t1172 * t615;
t604 = -pkin(2) * t639 + t1366 * t659 + (t1291 - t1361) * t731;
t623 = t1169 * t660 + t1173 * t640;
t1214 = pkin(8) * t623 + t1169 * t593 + t1173 * t604;
t631 = pkin(4) * t793 - t1366 * t672;
t634 = pkin(4) * t671 - qJ(4) * t672;
t657 = t1168 * t793 - t1172 * t671;
t605 = -pkin(9) * t657 - t1168 * t631 + t1172 * t634;
t620 = -pkin(2) * t657 - qJ(4) * t793 + t1366 * t671;
t630 = t1169 * t672 + t1173 * t658;
t1213 = pkin(8) * t630 + t1169 * t605 + t1173 * t620;
t698 = pkin(5) * t1343 + t1276;
t701 = qJ(6) * t1343 + t1185;
t635 = -t1167 * t701 - t1171 * t698 + t1465;
t708 = pkin(5) * t1188 + qJ(6) * t1196 + t738;
t624 = -t1168 * t635 + t1172 * t708 - t1504;
t633 = t1167 * t698 - t1171 * t701 + t1496;
t1212 = t1169 * t624 + t1173 * t633 + t1536;
t1139 = 0.2e1 * t1338;
t1177 = t1082 * t1367 + t1150 * t1085 + t1139 - t1184 - t1389;
t716 = -pkin(5) * t1329 + t1177 - t1535;
t663 = -t1167 * t716 - (-pkin(4) - t1360) * t1482 + t1515;
t666 = pkin(5) * t1079 + qJ(6) * t1394 + t1276 - t1506;
t796 = t1168 * t1482 - t1521;
t628 = -pkin(9) * t796 - t1168 * t663 + t1172 * t666;
t651 = -pkin(2) * t796 - t1171 * t716 - (qJ(4) + t1361) * t1482 + t1516;
t746 = t1173 * t798 - t1525;
t1211 = pkin(8) * t746 + t1169 * t628 + t1173 * t651;
t717 = (-t1197 - t1329) * pkin(5) + t1177;
t665 = -t1171 * t717 + t1197 * t1293 - t1468;
t674 = (t1390 + t1110) * qJ(6) + (t1057 + t1395) * pkin(5) - t1200 + t1451;
t629 = -t1168 * t665 + t1172 * t674 - t1505;
t652 = t1167 * t717 + t1197 * t1291 + t1497;
t1210 = t1169 * t629 + t1173 * t652 + t1537;
t646 = t1465 - t672;
t632 = -t1168 * t646 + t1172 * t738 - t1504;
t637 = t1496 + t671;
t1209 = t1169 * t632 + t1173 * t637 + t1536;
t692 = t1290 + t1451;
t709 = pkin(4) * t1197 + t1348 - t1468;
t642 = -t1168 * t709 + t1172 * t692 - t1505;
t682 = -qJ(4) * t1197 - t1355 + t1497;
t1208 = t1169 * t642 + t1173 * t682 + t1537;
t694 = t1506 - t728;
t711 = pkin(4) * t1187 - t1355 - t1515;
t802 = t1168 * t1187 + t1521;
t643 = -pkin(9) * t802 - t1168 * t711 + t1172 * t694;
t684 = -pkin(2) * t802 - qJ(4) * t1187 - t1348 - t1516;
t749 = t1173 * t804 + t1525;
t1207 = pkin(8) * t749 + t1169 * t643 + t1173 * t684;
t750 = t1168 * t842 - t1172 * t1229;
t699 = -pkin(2) * t750 + pkin(3) * t1229 - qJ(4) * t842;
t700 = -pkin(9) * t750 + (-qJ(4) * t1172 + t1363) * t850;
t724 = t1169 * t850 + t1173 * t751;
t1206 = pkin(8) * t724 + t1169 * t700 + t1173 * t699;
t817 = pkin(3) * t1036 + t842;
t826 = qJ(4) * t1036 + t1229;
t926 = t1008 * t1168 + t1434;
t723 = -pkin(9) * t926 - t1168 * t817 + t1172 * t826;
t827 = -pkin(2) * t926 - pkin(3) * t1392 - qJ(4) * t1008;
t863 = t1173 * t930 - t1479;
t1205 = pkin(8) * t863 + t1169 * t723 + t1173 * t827;
t816 = (t1007 - t1098) * pkin(3) + t1467;
t752 = qJ(4) * t1337 - t1168 * t816 + t1541;
t763 = -pkin(3) * t1190 + qJ(4) * t1064 - t1229 + t1543;
t884 = -t1007 * t1169 - t1518;
t1204 = pkin(8) * t884 + t1169 * t752 + t1173 * t763;
t815 = pkin(3) * t1098 - t1467 + t1483;
t757 = t1172 * t815 - t1363 * t1393 - t1539;
t768 = -t1542 - pkin(3) * t1071 + t1139 + (t1191 + t1124) * qJ(4) + t1232;
t1203 = t1169 * t757 + t1173 * t768 - t1566;
t838 = -t1357 - t1543;
t880 = t1351 - t1541;
t883 = -t1016 * t1169 + t1518;
t1202 = pkin(8) * t883 + t1169 * t880 + t1173 * t838;
t841 = t899 + t1542;
t885 = t1344 + t1539;
t1201 = t1169 * t885 + t1173 * t841 + t1566;
t1027 = t1093 * t1169 + t1094 * t1173;
t1199 = pkin(8) * t1027 + t969;
t927 = t1009 * t1168 + t1434;
t758 = -pkin(9) * t927 - t805;
t864 = t1173 * t931 - t1479;
t1198 = pkin(8) * t864 + t1169 * t758 - t1364 * t927;
t771 = t1169 * t991 + t1173 * t806;
t1192 = pkin(8) * t771 + t1284 * t805;
t1182 = t1165 * t1316 + t1166 * t1224;
t1147 = qJDD(1) * t1170 + t1174 * t1175;
t1136 = t1307 * t1316;
t1135 = (t1163 - t1164) * t1316;
t1131 = -pkin(7) * t1147 + g(3) * t1174;
t1099 = t1289 * t1307 * t1342;
t1092 = (t1305 + (0.2e1 * qJD(2) + t1341) * t1340) * t1165;
t1084 = t1173 * t1129 - t1163 * t1224;
t1083 = -t1169 * t1130 - t1164 * t1224;
t1074 = t1133 * t1173 - t1323;
t1073 = -t1132 * t1169 + t1320;
t1063 = (t1166 * t1129 + t1173 * t1182) * t1169;
t1062 = (t1166 * t1130 - t1169 * t1182) * t1173;
t1028 = -t1091 * t1169 + t1095 * t1173;
t1021 = t1165 * t1095 + t1166 * t1233;
t1020 = -t1165 * t1094 + t1166 * t1235;
t1019 = -t1165 * t1093 + t1166 * t1234;
t1018 = -t1166 * t1095 + t1165 * t1233;
t1005 = -t1165 * t1092 + t1166 * t1236;
t1004 = t1166 * t1092 + t1165 * t1236;
t990 = -t1165 * t1135 + t1166 * t1238;
t989 = t1165 * t1136 + t1166 * t1237;
t988 = -t1166 * t1136 + t1165 * t1237;
t952 = -t1021 * t1170 + t1075 * t1174;
t951 = t1021 * t1174 + t1075 * t1170;
t941 = -t1005 * t1170 + t1067 * t1174;
t940 = t1005 * t1174 + t1067 * t1170;
t939 = t1165 * t1100 + t1166 * t1239;
t938 = -t1166 * t1100 + t1165 * t1239;
t935 = t1027 * t1174 - t1170 * t989;
t934 = t1027 * t1170 + t1174 * t989;
t929 = -t1007 * t1168 + t1433;
t928 = t1016 * t1168 + t1433;
t903 = -t1009 * t1169 - t1517;
t902 = t1008 * t1169 + t1517;
t888 = -t1311 + (-t1018 * t1165 - t1021 * t1166) * pkin(8);
t882 = -t1310 + (-t1004 * t1165 - t1005 * t1166) * pkin(8);
t881 = -pkin(1) * t1018 + t1165 * t1065 + t1166 * t1274;
t875 = t1173 * t933 - t1440;
t874 = t1173 * t932 - t1440;
t873 = -pkin(1) * t1004 + t1165 * t1066 + t1166 * t1275;
t862 = pkin(8) * t1166 * t969 - pkin(1) * t938;
t861 = -t1170 * t939 + t1174 * t969;
t860 = t1170 * t969 + t1174 * t939;
t851 = -pkin(1) * t988 + t1166 * t1199;
t848 = -pkin(2) * t1393 + t1351 + t1538;
t847 = t1166 * t1251 + t1531;
t846 = t1166 * t1252 - t1531;
t843 = (-t1165 * t938 - t1166 * t939) * pkin(8);
t840 = pkin(2) * t1016 - t1344 + t1540;
t839 = (-t1165 * t988 - t1166 * t989) * pkin(8) - t1239;
t831 = t1166 * t1253 + t1534;
t830 = t1166 * t1246 - t1534;
t829 = t1165 * t1253 - t1530;
t828 = t1165 * t1246 + t1530;
t810 = t1168 * t925 + t1473;
t807 = -t1168 * t1188 + t1473;
t800 = -t1165 * t929 + t1166 * t1242;
t799 = -t1165 * t928 + t1166 * t1243;
t790 = -t1165 * t927 + t1166 * t1244;
t789 = -t1165 * t926 + t1166 * t1245;
t788 = t1165 * t1244 + t1166 * t927;
t787 = t1165 * t1245 + t1166 * t926;
t782 = t1172 * t818 - t1481;
t781 = t1172 * t819 + t1481;
t772 = -pkin(2) * t991 + pkin(9) * t806;
t762 = -t1170 * t831 + t1174 * t884;
t761 = -t1170 * t830 + t1174 * t883;
t760 = t1170 * t884 + t1174 * t831;
t759 = t1170 * t883 + t1174 * t830;
t756 = t1173 * t814 - t1499;
t753 = t1173 * t811 - t1499;
t747 = pkin(9) * t931 + t1485 + t806;
t744 = -t1538 + t1168 * t815 + (pkin(2) + t1362) * t1393;
t743 = t1007 * t1294 + t1172 * t816 - t1540;
t742 = -t1170 * t790 + t1174 * t864;
t741 = -t1170 * t789 + t1174 * t863;
t740 = t1170 * t864 + t1174 * t790;
t739 = t1170 * t863 + t1174 * t789;
t733 = -t1169 * t822 + t1173 * t784;
t732 = -t1169 * t823 + t1173 * t783;
t726 = -t1165 * t805 + t1166 * t1262;
t725 = t1165 * t1262 + t1166 * t805;
t722 = pkin(9) * t930 + t1168 * t826 + t1172 * t817 + t1485;
t721 = -t1165 * t810 + t1166 * t1258;
t718 = -t1165 * t807 + t1166 * t1261;
t715 = -t1165 * t802 + t1166 * t1263;
t713 = t1165 * t1263 + t1166 * t802;
t707 = -t1165 * t796 + t1166 * t1265;
t705 = t1165 * t1265 + t1166 * t796;
t702 = -t1169 * t841 + t1173 * t885 - t1575;
t697 = -t1169 * t838 + t1173 * t880 + (-t1165 * t828 - t1166 * t830) * pkin(8);
t696 = -t1165 * t782 + t1166 * t1267;
t695 = -t1165 * t781 + t1166 * t1268;
t693 = -t1165 * t848 + t1166 * t1201 - t1574;
t687 = -t1170 * t726 + t1174 * t771;
t686 = t1170 * t771 + t1174 * t726;
t685 = -pkin(1) * t828 - t1165 * t840 + t1166 * t1202;
t683 = pkin(9) * t751 + (-t1294 - t1362) * t850;
t681 = -t1165 * t750 + t1166 * t1271;
t680 = t1165 * t1271 + t1166 * t750;
t679 = t927 * t1365 + t1173 * t758 + (-t1165 * t788 - t1166 * t790) * pkin(8);
t678 = -t1170 * t715 + t1174 * t749;
t676 = t1170 * t749 + t1174 * t715;
t673 = -t1169 * t768 + t1173 * t757 + t1575;
t670 = -t1170 * t707 + t1174 * t746;
t668 = t1170 * t746 + t1174 * t707;
t664 = -t1169 * t763 + t1173 * t752 + (-t1165 * t829 - t1166 * t831) * pkin(8);
t662 = -pkin(1) * t788 - t1165 * t747 + t1166 * t1198;
t661 = -t1169 * t827 + t1173 * t723 + (-t1165 * t787 - t1166 * t789) * pkin(8);
t650 = -t1165 * t744 + t1166 * t1203 + t1574;
t649 = -pkin(1) * t829 - t1165 * t743 + t1166 * t1204;
t648 = -t1170 * t681 + t1174 * t724;
t647 = t1170 * t724 + t1174 * t681;
t645 = (-pkin(9) * t1173 + t1365) * t805 + (-t1165 * t725 - t1166 * t726) * pkin(8);
t644 = -pkin(1) * t725 - t1165 * t772 + t1166 * t1192;
t641 = -pkin(1) * t787 - t1165 * t722 + t1166 * t1205;
t638 = pkin(9) * t804 + t1168 * t694 + t1172 * t711 - t1544;
t636 = t1168 * t692 + t1172 * t709 + t1494;
t627 = t1168 * t674 + t1172 * t665 + t1494;
t626 = t1168 * t738 + t1172 * t646 + t1495;
t625 = pkin(9) * t798 + t1168 * t666 + t1172 * t663 + t1544;
t622 = t1168 * t708 + t1172 * t635 + t1495;
t621 = -t1169 * t699 + t1173 * t700 + (-t1165 * t680 - t1166 * t681) * pkin(8);
t619 = -t1165 * t657 + t1166 * t1272;
t618 = t1165 * t1272 + t1166 * t657;
t617 = -pkin(1) * t680 - t1165 * t683 + t1166 * t1206;
t616 = -t1169 * t684 + t1173 * t643 + (-t1165 * t713 - t1166 * t715) * pkin(8);
t614 = -t1169 * t682 + t1173 * t642 + t1545;
t612 = -t1165 * t639 + t1166 * t1273;
t611 = t1165 * t1273 + t1166 * t639;
t610 = -t1169 * t652 + t1173 * t629 + t1545;
t609 = -t1169 * t651 + t1173 * t628 + (-t1165 * t705 - t1166 * t707) * pkin(8);
t608 = -pkin(1) * t713 - t1165 * t638 + t1166 * t1207;
t607 = -t1165 * t636 + t1166 * t1208 - t1551;
t606 = -t1169 * t637 + t1173 * t632 + t1546;
t603 = -t1170 * t619 + t1174 * t630;
t602 = t1170 * t630 + t1174 * t619;
t601 = -t1169 * t633 + t1173 * t624 + t1546;
t600 = -pkin(2) * t672 + pkin(9) * t658 + t1168 * t634 + t1172 * t631;
t599 = -t1165 * t627 + t1166 * t1210 - t1551;
t598 = -pkin(1) * t705 - t1165 * t625 + t1166 * t1211;
t597 = -t1165 * t626 + t1166 * t1209 - t1552;
t596 = -t1170 * t612 + t1174 * t623;
t595 = t1170 * t623 + t1174 * t612;
t594 = -t1165 * t622 + t1166 * t1212 - t1552;
t592 = -pkin(2) * t660 + pkin(9) * t640 + t1168 * t615 + t1172 * t613;
t591 = -t1169 * t620 + t1173 * t605 + (-t1165 * t618 - t1166 * t619) * pkin(8);
t590 = -pkin(1) * t618 - t1165 * t600 + t1166 * t1213;
t589 = -t1169 * t604 + t1173 * t593 + (-t1165 * t611 - t1166 * t612) * pkin(8);
t588 = -pkin(1) * t611 - t1165 * t592 + t1166 * t1214;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1147, -t1148, 0, t1103, 0, 0, 0, 0, 0, 0, t952, t941, t935, t861, 0, 0, 0, 0, 0, 0, t761, -t767, t742, t687, 0, 0, 0, 0, 0, 0, t741, t762, t767, t648, 0, 0, 0, 0, 0, 0, t1549, t678, t1550, t603, 0, 0, 0, 0, 0, 0, t1549, t1550, t670, t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1148, -t1147, 0, t1102, 0, 0, 0, 0, 0, 0, t951, t940, t934, t860, 0, 0, 0, 0, 0, 0, t759, t764, t740, t686, 0, 0, 0, 0, 0, 0, t739, t760, -t764, t647, 0, 0, 0, 0, 0, 0, t1547, t676, t1548, t602, 0, 0, 0, 0, 0, 0, t1547, t1548, t668, t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1018, t1004, t988, t938, 0, 0, 0, 0, 0, 0, t828, t834, t788, t725, 0, 0, 0, 0, 0, 0, t787, t829, -t834, t680, 0, 0, 0, 0, 0, 0, t1510, t713, t1509, t618, 0, 0, 0, 0, 0, 0, t1510, t1509, t705, t611; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1148, 0, -t1147, 0, t1283, -t1131, -t1102, -pkin(7) * t1102, -t1063 * t1170 + t1084 * t1174, t1028 * t1174 - t1170 * t990, -t1019 * t1170 + t1073 * t1174, -t1062 * t1170 + t1083 * t1174, -t1020 * t1170 + t1074 * t1174, t1099 * t1174 + t1170 * t1317, -pkin(7) * t951 - t1170 * t881 + t1174 * t888, -pkin(7) * t940 - t1170 * t873 + t1174 * t882, -pkin(7) * t934 - t1170 * t851 + t1174 * t839, -pkin(7) * t860 - t1170 * t862 + t1174 * t843, t1455, -t1170 * t799 + t1174 * t874, t1570, t1456, -t1170 * t846 + t1174 * t902, t1413, -pkin(7) * t759 - t1170 * t685 + t1174 * t697, -t1170 * t693 + t1174 * t702 - t1579, -pkin(7) * t740 - t1170 * t662 + t1174 * t679, -pkin(7) * t686 - t1170 * t644 + t1174 * t645, t1413, -t1570, -t1170 * t847 + t1174 * t903, t1455, -t1170 * t800 + t1174 * t875, t1456, -pkin(7) * t739 - t1170 * t641 + t1174 * t661, -pkin(7) * t760 - t1170 * t649 + t1174 * t664, -t1170 * t650 + t1174 * t673 + t1579, -pkin(7) * t647 - t1170 * t617 + t1174 * t621, t1411, -t1170 * t695 + t1174 * t732, -t1170 * t718 + t1174 * t753, t1491, -t1577, t1490, -t1170 * t607 + t1174 * t614 - t1565, -pkin(7) * t676 - t1170 * t608 + t1174 * t616, -t1170 * t597 + t1174 * t606 - t1564, -pkin(7) * t602 - t1170 * t590 + t1174 * t591, t1411, -t1170 * t721 + t1174 * t756, -t1170 * t696 + t1174 * t733, t1490, t1577, t1491, -t1170 * t599 + t1174 * t610 - t1565, -t1170 * t594 + t1174 * t601 - t1564, -pkin(7) * t668 - t1170 * t598 + t1174 * t609, -pkin(7) * t595 - t1170 * t588 + t1174 * t589; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1147, 0, t1148, 0, t1131, t1283, t1103, pkin(7) * t1103, t1063 * t1174 + t1084 * t1170, t1028 * t1170 + t1174 * t990, t1019 * t1174 + t1073 * t1170, t1062 * t1174 + t1083 * t1170, t1020 * t1174 + t1074 * t1170, t1099 * t1170 - t1174 * t1317, pkin(7) * t952 + t1170 * t888 + t1174 * t881, pkin(7) * t941 + t1170 * t882 + t1174 * t873, pkin(7) * t935 + t1170 * t839 + t1174 * t851, pkin(7) * t861 + t1170 * t843 + t1174 * t862, t1452, t1170 * t874 + t1174 * t799, -t1569, t1453, t1170 * t902 + t1174 * t846, t1414, pkin(7) * t761 + t1170 * t697 + t1174 * t685, t1170 * t702 + t1174 * t693 - t1578, pkin(7) * t742 + t1170 * t679 + t1174 * t662, pkin(7) * t687 + t1170 * t645 + t1174 * t644, t1414, t1569, t1170 * t903 + t1174 * t847, t1452, t1170 * t875 + t1174 * t800, t1453, pkin(7) * t741 + t1170 * t661 + t1174 * t641, pkin(7) * t762 + t1170 * t664 + t1174 * t649, t1170 * t673 + t1174 * t650 + t1578, pkin(7) * t648 + t1170 * t621 + t1174 * t617, t1412, t1170 * t732 + t1174 * t695, t1170 * t753 + t1174 * t718, t1487, t1576, t1486, t1170 * t614 + t1174 * t607 + t1563, pkin(7) * t678 + t1170 * t616 + t1174 * t608, t1170 * t606 + t1174 * t597 + t1562, pkin(7) * t603 + t1170 * t591 + t1174 * t590, t1412, t1170 * t756 + t1174 * t721, t1170 * t733 + t1174 * t696, t1486, -t1576, t1487, t1170 * t610 + t1174 * t599 + t1563, t1170 * t601 + t1174 * t594 + t1562, pkin(7) * t670 + t1170 * t609 + t1174 * t598, pkin(7) * t596 + t1170 * t589 + t1174 * t588; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1152, t1153, 0, 0, (t1165 * t1129 + t1173 * t1410) * t1169, t1166 * t1135 + t1165 * t1238, t1166 * t1093 + t1165 * t1234, (t1165 * t1130 - t1169 * t1410) * t1173, t1166 * t1094 + t1165 * t1235, t1166 * t1157, pkin(1) * t1021 - t1166 * t1065 + t1165 * t1274, pkin(1) * t1005 - t1166 * t1066 + t1165 * t1275, pkin(1) * t989 + t1165 * t1199, pkin(1) * t939 + t1359 * t969, t1423, t1165 * t1243 + t1166 * t928, -t1555, t1422, t1165 * t1252 + t1527, t1373, pkin(1) * t830 + t1165 * t1202 + t1166 * t840, t1165 * t1201 + t1166 * t848 + t1573, pkin(1) * t790 + t1165 * t1198 + t1166 * t747, pkin(1) * t726 + t1165 * t1192 + t1166 * t772, t1373, t1555, t1165 * t1251 - t1527, t1423, t1165 * t1242 + t1166 * t929, t1422, pkin(1) * t789 + t1165 * t1205 + t1166 * t722, pkin(1) * t831 + t1165 * t1204 + t1166 * t743, t1165 * t1203 + t1166 * t744 - t1573, pkin(1) * t681 + t1165 * t1206 + t1166 * t683, t1374, t1165 * t1268 + t1166 * t781, t1165 * t1261 + t1166 * t807, t1460, t1568, t1461, t1165 * t1208 + t1166 * t636 + t1553, pkin(1) * t715 + t1165 * t1207 + t1166 * t638, t1165 * t1209 + t1166 * t626 + t1554, pkin(1) * t619 + t1165 * t1213 + t1166 * t600, t1374, t1165 * t1258 + t1166 * t810, t1165 * t1267 + t1166 * t782, t1461, -t1568, t1460, t1165 * t1210 + t1166 * t627 + t1553, t1165 * t1212 + t1166 * t622 + t1554, pkin(1) * t707 + t1165 * t1211 + t1166 * t625, pkin(1) * t612 + t1165 * t1214 + t1166 * t592;];
tauB_reg  = t1;