% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRPPR6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:09:29
% EndTime: 2019-05-05 17:10:00
% DurationCPUTime: 32.46s
% Computational Cost: add. (211830->866), mult. (485588->1123), div. (0->0), fcn. (333958->10), ass. (0->569)
t1389 = sin(pkin(9));
t1391 = cos(pkin(9));
t1393 = sin(qJ(3));
t1396 = cos(qJ(3));
t1349 = (t1389 * t1396 + t1391 * t1393) * qJD(1);
t1537 = qJD(1) * t1396;
t1538 = qJD(1) * t1393;
t1351 = -t1389 * t1538 + t1391 * t1537;
t1510 = t1351 * t1349;
t1561 = qJDD(3) - t1510;
t1563 = t1389 * t1561;
t1562 = t1391 * t1561;
t1379 = t1396 * qJDD(1);
t1481 = qJD(3) * t1538;
t1358 = t1379 - t1481;
t1480 = qJD(3) * t1537;
t1492 = t1393 * qJDD(1);
t1428 = -t1480 - t1492;
t1301 = t1391 * t1358 + t1389 * t1428;
t1388 = sin(pkin(10));
t1390 = cos(pkin(10));
t1436 = qJDD(3) * t1390 - t1301 * t1388;
t1312 = qJD(3) * t1388 + t1351 * t1390;
t1515 = t1349 * t1312;
t1223 = -t1515 - t1436;
t1287 = qJDD(3) * t1388 + t1301 * t1390;
t1311 = -t1390 * qJD(3) + t1351 * t1388;
t1516 = t1349 * t1311;
t1226 = -t1516 + t1287;
t1227 = t1516 + t1287;
t1535 = qJD(3) * t1349;
t1275 = -t1535 + t1301;
t1270 = t1312 * t1311;
t1300 = t1358 * t1389 - t1391 * t1428;
t1554 = -t1270 + t1300;
t1560 = t1388 * t1554;
t1559 = t1390 * t1554;
t1392 = sin(qJ(6));
t1395 = cos(qJ(6));
t1261 = t1395 * t1311 + t1312 * t1392;
t1263 = -t1311 * t1392 + t1312 * t1395;
t1201 = t1263 * t1261;
t1296 = qJDD(6) + t1300;
t1556 = -t1201 + t1296;
t1558 = t1392 * t1556;
t1557 = t1395 * t1556;
t1534 = qJD(4) * t1349;
t1333 = -0.2e1 * t1534;
t1398 = qJD(1) ^ 2;
t1394 = sin(qJ(1));
t1397 = cos(qJ(1));
t1367 = t1394 * g(1) - t1397 * g(2);
t1432 = qJDD(2) - t1367;
t1407 = -t1398 * qJ(2) + t1432;
t1549 = pkin(7) + pkin(1);
t1403 = -qJDD(1) * t1549 + t1407;
t1306 = t1396 * g(3) - t1393 * t1403;
t1385 = t1393 ^ 2;
t1501 = t1385 * t1398;
t1552 = qJD(3) ^ 2;
t1371 = -t1501 - t1552;
t1258 = pkin(3) * t1371 - qJ(4) * t1492 - t1306;
t1402 = t1396 * t1403;
t1497 = t1396 * t1398;
t1539 = qJD(1) * qJD(3);
t1399 = t1402 - t1358 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t1497 - qJ(4) * t1539 + g(3)) * t1393;
t1494 = t1391 * t1258 + t1389 * t1399;
t1177 = t1333 + t1494;
t1289 = pkin(4) * t1349 - qJ(5) * t1351;
t1141 = -pkin(4) * t1552 + qJDD(3) * qJ(5) - t1289 * t1349 + t1177;
t1337 = qJD(3) * t1351;
t1272 = t1300 + t1337;
t1536 = qJD(2) * qJD(1);
t1382 = 0.2e1 * t1536;
t1368 = t1397 * g(1) + t1394 * g(2);
t1383 = qJDD(1) * qJ(2);
t1433 = t1368 - t1383;
t1553 = -(qJ(4) * t1385 + t1549) * t1398 - pkin(3) * t1428 + qJDD(4) + (qJD(3) * pkin(3) - qJ(4) * t1537) * t1537 - t1433;
t1156 = t1272 * pkin(4) - qJ(5) * t1275 + t1382 + t1553;
t1077 = 0.2e1 * qJD(5) * t1312 + t1388 * t1141 - t1390 * t1156;
t1078 = -0.2e1 * qJD(5) * t1311 + t1390 * t1141 + t1388 * t1156;
t1010 = t1388 * t1077 + t1390 * t1078;
t1171 = -t1261 * qJD(6) + t1395 * t1287 + t1392 * t1436;
t1339 = qJD(6) + t1349;
t1241 = t1339 * t1261;
t1555 = -t1241 + t1171;
t1305 = t1393 * g(3) + t1402;
t1247 = t1396 * t1305 - t1393 * t1306;
t1273 = t1300 - t1337;
t1462 = t1392 * t1287 - t1395 * t1436;
t1123 = (qJD(6) - t1339) * t1263 + t1462;
t1259 = t1261 ^ 2;
t1260 = t1263 ^ 2;
t1551 = t1311 ^ 2;
t1310 = t1312 ^ 2;
t1338 = t1339 ^ 2;
t1550 = t1349 ^ 2;
t1347 = t1351 ^ 2;
t1038 = t1554 * pkin(5) - pkin(8) * t1227 - t1077;
t1277 = pkin(5) * t1349 - pkin(8) * t1312;
t1043 = -pkin(5) * t1551 + pkin(8) * t1436 - t1349 * t1277 + t1078;
t985 = -t1395 * t1038 + t1043 * t1392;
t986 = t1038 * t1392 + t1043 * t1395;
t943 = t1392 * t986 - t1395 * t985;
t1548 = pkin(5) * t943;
t1547 = pkin(2) * t1247;
t1491 = -0.2e1 * t1536;
t1406 = t1433 + t1491;
t1326 = t1398 * t1549 + t1406;
t1546 = pkin(2) * t1326;
t1386 = t1396 ^ 2;
t1493 = t1385 + t1386;
t1360 = t1493 * qJDD(1);
t1545 = pkin(2) * t1360;
t1544 = pkin(4) * t1389;
t1126 = t1241 + t1171;
t1066 = -t1123 * t1392 - t1126 * t1395;
t1543 = pkin(5) * t1066;
t1542 = qJDD(1) * pkin(1);
t1541 = t1388 * t943;
t1540 = t1390 * t943;
t1533 = qJD(4) * t1351;
t1463 = t1389 * t1258 - t1391 * t1399;
t1140 = qJDD(5) - t1552 * qJ(5) - qJDD(3) * pkin(4) + (0.2e1 * qJD(4) + t1289) * t1351 + t1463;
t1085 = -pkin(5) * t1436 - pkin(8) * t1551 + t1312 * t1277 + t1140;
t1532 = t1085 * t1392;
t1531 = t1085 * t1395;
t1490 = 0.2e1 * t1533;
t1176 = t1463 + t1490;
t1093 = -t1176 * t1391 + t1177 * t1389;
t1530 = t1093 * t1393;
t1529 = t1093 * t1396;
t1167 = t1201 + t1296;
t1528 = t1167 * t1392;
t1527 = t1167 * t1395;
t1230 = t1270 + t1300;
t1526 = t1230 * t1388;
t1525 = t1230 * t1390;
t1264 = t1491 - t1553;
t1524 = t1264 * t1389;
t1523 = t1264 * t1391;
t1292 = qJDD(3) + t1510;
t1522 = t1292 * t1389;
t1521 = t1292 * t1391;
t1520 = t1300 * t1389;
t1519 = t1339 * t1263;
t1518 = t1339 * t1392;
t1517 = t1339 * t1395;
t1514 = t1349 * t1388;
t1513 = t1349 * t1389;
t1512 = t1349 * t1390;
t1511 = t1349 * t1391;
t1509 = t1351 * t1389;
t1508 = t1351 * t1391;
t1357 = 0.2e1 * t1480 + t1492;
t1314 = t1357 * t1393;
t1507 = t1360 * t1394;
t1506 = t1360 * t1397;
t1374 = t1393 * t1497;
t1365 = qJDD(3) + t1374;
t1505 = t1365 * t1393;
t1504 = t1365 * t1396;
t1366 = qJDD(3) - t1374;
t1503 = t1366 * t1393;
t1502 = t1366 * t1396;
t1500 = t1386 * t1398;
t1136 = t1388 * t1140;
t1137 = t1390 * t1140;
t1498 = t1393 * t1326;
t1313 = t1396 * t1326;
t1495 = -pkin(4) * t1140 + qJ(5) * t1010;
t1489 = -pkin(4) * t1391 - pkin(3);
t993 = t1010 * t1389 - t1140 * t1391;
t1488 = pkin(3) * t993 + t1495;
t1487 = t1389 * t1201;
t1486 = t1391 * t1201;
t1485 = t1389 * t1270;
t1484 = t1391 * t1270;
t1483 = t1394 * t1510;
t1482 = t1397 * t1510;
t1257 = -t1310 - t1550;
t1165 = -t1257 * t1388 - t1525;
t1479 = -pkin(4) * t1226 + qJ(5) * t1165 + t1136;
t1244 = -t1550 - t1551;
t1158 = t1244 * t1390 - t1560;
t1224 = -t1515 + t1436;
t1478 = pkin(4) * t1224 + qJ(5) * t1158 - t1137;
t1331 = -t1347 - t1552;
t1253 = t1331 * t1391 - t1522;
t1477 = pkin(3) * t1253 - t1494;
t1094 = t1176 * t1389 + t1391 * t1177;
t1027 = t1094 * t1393 + t1529;
t1092 = pkin(3) * t1093;
t1476 = -pkin(2) * t1027 - t1092;
t1276 = t1535 + t1301;
t1204 = -t1273 * t1389 - t1276 * t1391;
t1206 = -t1273 * t1391 + t1276 * t1389;
t1133 = t1204 * t1396 + t1206 * t1393;
t1202 = pkin(3) * t1204;
t1475 = -pkin(2) * t1133 - t1202;
t944 = t1392 * t985 + t1395 * t986;
t925 = t1388 * t944 + t1540;
t941 = -pkin(5) * t1085 + pkin(8) * t944;
t907 = -pkin(8) * t1540 - qJ(5) * t925 - t1388 * t941;
t910 = -pkin(4) * t925 - t1548;
t926 = t1390 * t944 - t1541;
t921 = t1085 * t1389 + t1391 * t926;
t892 = -pkin(3) * t925 + qJ(4) * t921 + t1389 * t907 + t1391 * t910;
t920 = -t1085 * t1391 + t1389 * t926;
t896 = -qJ(4) * t920 - t1389 * t910 + t1391 * t907;
t1474 = -t1393 * t892 + t1396 * t896;
t1068 = -t1123 * t1395 + t1126 * t1392;
t1145 = -t1259 - t1260;
t935 = -pkin(5) * t1145 + pkin(8) * t1068 + t944;
t937 = -pkin(8) * t1066 - t943;
t998 = t1066 * t1390 + t1068 * t1388;
t912 = -qJ(5) * t998 - t1388 * t935 + t1390 * t937;
t971 = -pkin(4) * t998 - t1543;
t1000 = -t1066 * t1388 + t1068 * t1390;
t989 = t1000 * t1391 + t1145 * t1389;
t904 = -pkin(3) * t998 + qJ(4) * t989 + t1389 * t912 + t1391 * t971;
t988 = t1000 * t1389 - t1145 * t1391;
t909 = -qJ(4) * t988 - t1389 * t971 + t1391 * t912;
t1473 = -t1393 * t904 + t1396 * t909;
t1180 = -t1338 - t1259;
t1095 = t1180 * t1392 + t1557;
t1096 = t1180 * t1395 - t1558;
t1036 = -t1095 * t1388 + t1096 * t1390;
t1122 = (qJD(6) + t1339) * t1263 + t1462;
t1007 = t1036 * t1391 + t1122 * t1389;
t1035 = t1095 * t1390 + t1096 * t1388;
t1002 = -pkin(5) * t1122 + pkin(8) * t1096 - t1531;
t1023 = -pkin(8) * t1095 + t1532;
t947 = -qJ(5) * t1035 - t1002 * t1388 + t1023 * t1390;
t1431 = pkin(5) * t1095 - t985;
t952 = -pkin(4) * t1035 - t1431;
t915 = -pkin(3) * t1035 + qJ(4) * t1007 + t1389 * t947 + t1391 * t952;
t1006 = t1036 * t1389 - t1122 * t1391;
t918 = -qJ(4) * t1006 - t1389 * t952 + t1391 * t947;
t1472 = -t1393 * t915 + t1396 * t918;
t1208 = -t1260 - t1338;
t1103 = t1208 * t1395 - t1528;
t1104 = -t1208 * t1392 - t1527;
t1048 = -t1103 * t1388 + t1104 * t1390;
t1013 = t1048 * t1391 + t1389 * t1555;
t1047 = t1103 * t1390 + t1104 * t1388;
t1004 = -pkin(5) * t1555 + pkin(8) * t1104 + t1532;
t1029 = -pkin(8) * t1103 + t1531;
t953 = -qJ(5) * t1047 - t1004 * t1388 + t1029 * t1390;
t1409 = pkin(5) * t1103 - t986;
t956 = -pkin(4) * t1047 - t1409;
t916 = -pkin(3) * t1047 + qJ(4) * t1013 + t1389 * t953 + t1391 * t956;
t1012 = t1048 * t1389 - t1391 * t1555;
t923 = -qJ(4) * t1012 - t1389 * t956 + t1391 * t953;
t1471 = -t1393 * t916 + t1396 * t923;
t1009 = -t1077 * t1390 + t1078 * t1388;
t994 = t1010 * t1391 + t1140 * t1389;
t933 = qJ(4) * t994 + (-qJ(5) * t1389 + t1489) * t1009;
t939 = -qJ(4) * t993 + (-qJ(5) * t1391 + t1544) * t1009;
t1470 = -t1393 * t933 + t1396 * t939;
t1150 = -t1223 * t1390 + t1227 * t1388;
t1238 = t1310 + t1551;
t1101 = t1150 * t1391 - t1238 * t1389;
t1148 = -t1223 * t1388 - t1227 * t1390;
t995 = -qJ(5) * t1148 - t1009;
t961 = qJ(4) * t1101 + t1148 * t1489 + t1389 * t995;
t1100 = t1150 * t1389 + t1238 * t1391;
t966 = -qJ(4) * t1100 + t1148 * t1544 + t1391 * t995;
t1469 = -t1393 * t961 + t1396 * t966;
t1157 = t1244 * t1388 + t1559;
t1039 = -pkin(4) * t1157 + t1077;
t1079 = -qJ(5) * t1157 + t1136;
t1106 = t1158 * t1391 - t1224 * t1389;
t974 = -pkin(3) * t1157 + qJ(4) * t1106 + t1039 * t1391 + t1079 * t1389;
t1105 = t1158 * t1389 + t1224 * t1391;
t981 = -qJ(4) * t1105 - t1039 * t1389 + t1079 * t1391;
t1468 = -t1393 * t974 + t1396 * t981;
t1164 = t1257 * t1390 - t1526;
t1042 = -pkin(4) * t1164 + t1078;
t1082 = -qJ(5) * t1164 + t1137;
t1109 = t1165 * t1391 + t1226 * t1389;
t975 = -pkin(3) * t1164 + qJ(4) * t1109 + t1042 * t1391 + t1082 * t1389;
t1108 = t1165 * t1389 - t1226 * t1391;
t983 = -qJ(4) * t1108 - t1042 * t1389 + t1082 * t1391;
t1467 = -t1393 * t975 + t1396 * t983;
t1271 = -t1550 - t1347;
t1064 = -pkin(3) * t1271 + qJ(4) * t1206 + t1094;
t1074 = -qJ(4) * t1204 - t1093;
t1466 = -t1393 * t1064 + t1396 * t1074;
t1290 = -t1550 - t1552;
t1237 = t1290 * t1391 - t1563;
t1142 = -pkin(3) * t1272 + qJ(4) * t1237 + t1523;
t1236 = t1290 * t1389 + t1562;
t1179 = -qJ(4) * t1236 - t1524;
t1465 = -t1393 * t1142 + t1396 * t1179;
t1256 = -t1331 * t1389 - t1521;
t1144 = -pkin(3) * t1275 + qJ(4) * t1256 - t1524;
t1195 = -qJ(4) * t1253 - t1523;
t1464 = -t1393 * t1144 + t1396 * t1195;
t1332 = t1398 * pkin(1) + t1406;
t1340 = -t1407 + t1542;
t1460 = -t1397 * t1332 - t1340 * t1394;
t1459 = -t1367 * t1394 - t1397 * t1368;
t1458 = -pkin(4) * t1145 + qJ(5) * t1000 + t1388 * t937 + t1390 * t935;
t1456 = t1394 * t1374;
t1455 = t1397 * t1374;
t1454 = -pkin(4) * t1122 + qJ(5) * t1036 + t1390 * t1002 + t1388 * t1023;
t1453 = -pkin(4) * t1555 + qJ(5) * t1048 + t1390 * t1004 + t1388 * t1029;
t1452 = pkin(4) * t1238 + qJ(5) * t1150 + t1010;
t1451 = pkin(3) * t1108 + t1479;
t1450 = pkin(3) * t1105 + t1478;
t1361 = qJDD(1) * t1397 - t1394 * t1398;
t1449 = pkin(6) * t1361 + g(3) * t1394;
t1362 = qJDD(1) * t1394 + t1397 * t1398;
t1448 = -pkin(6) * t1362 + g(3) * t1397;
t1447 = pkin(3) * t988 + t1458;
t1446 = pkin(3) * t1236 - t1463;
t1445 = pkin(2) * t1357 - t1313;
t1359 = t1379 - 0.2e1 * t1481;
t1444 = pkin(2) * t1359 + t1498;
t950 = t1393 * t994 + t1396 * t993;
t1443 = -pkin(2) * t950 - t1488;
t1442 = pkin(3) * t1006 + t1454;
t1441 = pkin(3) * t1012 + t1453;
t1440 = pkin(3) * t1100 + t1452;
t1248 = -t1305 * t1393 - t1306 * t1396;
t1439 = t1332 * t1394 - t1340 * t1397;
t1438 = t1367 * t1397 - t1368 * t1394;
t1187 = t1253 * t1396 + t1256 * t1393;
t1435 = -pkin(2) * t1187 - t1477;
t1373 = -t1500 - t1552;
t1318 = t1373 * t1396 - t1505;
t1434 = -pkin(2) * t1318 - t1306;
t1430 = -pkin(4) * t1085 - pkin(8) * t1541 + qJ(5) * t926 + t1390 * t941;
t1083 = pkin(3) * t1264 + qJ(4) * t1094;
t1429 = -qJ(4) * t1529 - t1393 * t1083;
t1049 = t1105 * t1396 + t1106 * t1393;
t1427 = -pkin(2) * t1049 - t1450;
t1051 = t1108 * t1396 + t1109 * t1393;
t1426 = -pkin(2) * t1051 - t1451;
t945 = t1393 * t989 + t1396 * t988;
t1425 = -pkin(2) * t945 - t1447;
t1153 = t1236 * t1396 + t1237 * t1393;
t1424 = -pkin(2) * t1153 - t1446;
t1423 = pkin(3) * t920 + t1430;
t1422 = pkin(2) * t925 - t1393 * t896 - t1396 * t892;
t1421 = pkin(2) * t998 - t1393 * t909 - t1396 * t904;
t1420 = pkin(2) * t1009 - t1393 * t939 - t1396 * t933;
t1419 = pkin(2) * t1035 - t1393 * t918 - t1396 * t915;
t1418 = pkin(2) * t1047 - t1393 * t923 - t1396 * t916;
t1417 = pkin(2) * t1148 - t1393 * t966 - t1396 * t961;
t1416 = pkin(2) * t1157 - t1393 * t981 - t1396 * t974;
t1415 = pkin(2) * t1164 - t1393 * t983 - t1396 * t975;
t959 = t1006 * t1396 + t1007 * t1393;
t1414 = -pkin(2) * t959 - t1442;
t963 = t1012 * t1396 + t1013 * t1393;
t1413 = -pkin(2) * t963 - t1441;
t1412 = pkin(2) * t1271 - t1396 * t1064 - t1393 * t1074;
t1411 = pkin(2) * t1272 - t1396 * t1142 - t1393 * t1179;
t1410 = pkin(2) * t1275 - t1396 * t1144 - t1393 * t1195;
t1044 = t1100 * t1396 + t1101 * t1393;
t1408 = -pkin(2) * t1044 - t1440;
t1405 = -pkin(2) * t1264 + qJ(4) * t1530 - t1396 * t1083;
t905 = t1393 * t921 + t1396 * t920;
t1404 = -pkin(2) * t905 - t1423;
t1316 = t1371 * t1393 + t1502;
t1400 = -pkin(2) * t1316 - t1305;
t1380 = t1397 * qJDD(3);
t1378 = t1394 * qJDD(3);
t1372 = -t1500 + t1552;
t1370 = t1501 - t1552;
t1364 = (-t1385 + t1386) * t1398;
t1363 = t1493 * t1398;
t1355 = t1493 * t1539;
t1354 = t1432 - 0.2e1 * t1542;
t1348 = -t1368 + t1382 + 0.2e1 * t1383;
t1335 = -0.2e1 * t1533;
t1334 = 0.2e1 * t1534;
t1330 = -t1347 + t1552;
t1329 = t1550 - t1552;
t1328 = t1358 * t1393 + t1386 * t1539;
t1327 = t1385 * t1539 + t1396 * t1428;
t1323 = -t1373 * t1393 - t1504;
t1322 = -t1372 * t1393 + t1502;
t1321 = (t1358 - t1481) * t1396;
t1320 = t1371 * t1396 - t1503;
t1319 = t1370 * t1396 - t1505;
t1317 = t1372 * t1396 + t1503;
t1315 = t1370 * t1393 + t1504;
t1303 = -t1357 * t1396 - t1359 * t1393;
t1302 = t1359 * t1396 - t1314;
t1295 = t1347 - t1550;
t1294 = t1391 * t1300;
t1288 = pkin(1) * t1340 - qJ(2) * t1332;
t1283 = (t1509 - t1511) * qJD(3);
t1282 = (-t1508 - t1513) * qJD(3);
t1281 = -t1310 + t1550;
t1280 = -t1550 + t1551;
t1269 = t1310 - t1551;
t1268 = -qJD(3) * t1509 + t1301 * t1391;
t1267 = qJD(3) * t1508 + t1301 * t1389;
t1266 = qJD(3) * t1511 + t1520;
t1265 = qJD(3) * t1513 - t1294;
t1255 = -t1330 * t1389 + t1562;
t1254 = t1329 * t1391 - t1522;
t1252 = t1330 * t1391 + t1563;
t1251 = t1329 * t1389 + t1521;
t1243 = pkin(2) * t1363 + t1248;
t1240 = -t1260 + t1338;
t1239 = t1259 - t1338;
t1233 = -qJ(2) * t1323 - t1434;
t1232 = -qJ(2) * t1320 - t1400;
t1220 = t1287 * t1390 - t1312 * t1514;
t1219 = t1287 * t1388 + t1312 * t1512;
t1218 = t1311 * t1512 - t1388 * t1436;
t1217 = -t1311 * t1514 - t1390 * t1436;
t1216 = (-t1311 * t1390 + t1312 * t1388) * t1349;
t1215 = (-t1311 * t1388 - t1312 * t1390) * t1349;
t1214 = -t1320 * t1549 + t1445;
t1213 = -t1323 * t1549 + t1444;
t1212 = qJ(2) * t1359 - t1318 * t1549 - t1313;
t1211 = qJ(2) * t1357 - t1316 * t1549 - t1498;
t1210 = -t1282 * t1393 + t1283 * t1396;
t1209 = t1282 * t1396 + t1283 * t1393;
t1207 = -qJ(2) * t1363 + t1360 * t1549 - t1247;
t1205 = -t1272 * t1391 - t1275 * t1389;
t1203 = -t1272 * t1389 + t1275 * t1391;
t1200 = t1260 - t1259;
t1199 = -t1267 * t1393 + t1268 * t1396;
t1198 = -t1265 * t1393 + t1266 * t1396;
t1197 = t1267 * t1396 + t1268 * t1393;
t1196 = t1265 * t1396 + t1266 * t1393;
t1193 = t1216 * t1391 + t1520;
t1192 = t1216 * t1389 - t1294;
t1191 = -qJ(2) * t1248 + t1547;
t1190 = -t1253 * t1393 + t1256 * t1396;
t1189 = -t1252 * t1393 + t1255 * t1396;
t1188 = -t1251 * t1393 + t1254 * t1396;
t1186 = t1252 * t1396 + t1255 * t1393;
t1185 = t1251 * t1396 + t1254 * t1393;
t1184 = t1280 * t1390 - t1526;
t1183 = -t1281 * t1388 + t1559;
t1182 = t1280 * t1388 + t1525;
t1181 = t1281 * t1390 + t1560;
t1175 = t1220 * t1391 + t1485;
t1174 = t1218 * t1391 - t1485;
t1173 = t1220 * t1389 - t1484;
t1172 = t1218 * t1389 + t1484;
t1170 = -qJD(6) * t1263 - t1462;
t1162 = -t1248 * t1549 - t1546;
t1161 = -qJ(2) * t1326 - t1247 * t1549;
t1160 = (-t1261 * t1395 + t1263 * t1392) * t1339;
t1159 = (-t1261 * t1392 - t1263 * t1395) * t1339;
t1154 = -t1236 * t1393 + t1237 * t1396;
t1149 = t1224 * t1390 - t1226 * t1388;
t1147 = t1224 * t1388 + t1226 * t1390;
t1135 = -t1204 * t1393 + t1206 * t1396;
t1134 = -t1203 * t1393 + t1205 * t1396;
t1132 = t1203 * t1396 + t1205 * t1393;
t1131 = t1184 * t1391 - t1223 * t1389;
t1130 = t1183 * t1391 + t1227 * t1389;
t1129 = t1184 * t1389 + t1223 * t1391;
t1128 = t1183 * t1389 - t1227 * t1391;
t1119 = t1171 * t1395 - t1263 * t1518;
t1118 = t1171 * t1392 + t1263 * t1517;
t1117 = -t1170 * t1392 + t1261 * t1517;
t1116 = t1170 * t1395 + t1261 * t1518;
t1115 = t1239 * t1395 - t1528;
t1114 = -t1240 * t1392 + t1557;
t1113 = t1239 * t1392 + t1527;
t1112 = t1240 * t1395 + t1558;
t1111 = t1149 * t1391 + t1269 * t1389;
t1110 = t1149 * t1389 - t1269 * t1391;
t1098 = -t1192 * t1393 + t1193 * t1396;
t1097 = t1192 * t1396 + t1193 * t1393;
t1091 = -t1173 * t1393 + t1175 * t1396;
t1090 = -t1172 * t1393 + t1174 * t1396;
t1089 = t1173 * t1396 + t1175 * t1393;
t1088 = t1172 * t1396 + t1174 * t1393;
t1087 = -t1159 * t1388 + t1160 * t1390;
t1086 = t1159 * t1390 + t1160 * t1388;
t1081 = t1087 * t1391 + t1296 * t1389;
t1080 = t1087 * t1389 - t1296 * t1391;
t1072 = -t1129 * t1393 + t1131 * t1396;
t1071 = -t1128 * t1393 + t1130 * t1396;
t1070 = t1129 * t1396 + t1131 * t1393;
t1069 = t1128 * t1396 + t1130 * t1393;
t1067 = -t1122 * t1395 - t1392 * t1555;
t1065 = -t1122 * t1392 + t1395 * t1555;
t1063 = -t1118 * t1388 + t1119 * t1390;
t1062 = -t1116 * t1388 + t1117 * t1390;
t1061 = t1118 * t1390 + t1119 * t1388;
t1060 = t1116 * t1390 + t1117 * t1388;
t1059 = -t1110 * t1393 + t1111 * t1396;
t1058 = t1110 * t1396 + t1111 * t1393;
t1057 = -t1113 * t1388 + t1115 * t1390;
t1056 = -t1112 * t1388 + t1114 * t1390;
t1055 = t1113 * t1390 + t1115 * t1388;
t1054 = t1112 * t1390 + t1114 * t1388;
t1053 = -qJ(2) * t1190 + t1334 - t1435;
t1052 = -t1108 * t1393 + t1109 * t1396;
t1050 = -t1105 * t1393 + t1106 * t1396;
t1045 = -t1100 * t1393 + t1101 * t1396;
t1041 = -qJ(2) * t1135 - t1475;
t1040 = -qJ(2) * t1154 + t1335 - t1424;
t1033 = t1063 * t1391 + t1487;
t1032 = t1062 * t1391 - t1487;
t1031 = t1063 * t1389 - t1486;
t1030 = t1062 * t1389 + t1486;
t1028 = t1094 * t1396 - t1530;
t1025 = -t1190 * t1549 + t1410;
t1024 = qJ(2) * t1275 - t1187 * t1549 + t1464;
t1021 = t1057 * t1391 - t1123 * t1389;
t1020 = t1056 * t1391 + t1126 * t1389;
t1019 = t1057 * t1389 + t1123 * t1391;
t1018 = t1056 * t1389 - t1126 * t1391;
t1017 = -t1154 * t1549 + t1411;
t1016 = qJ(2) * t1272 - t1153 * t1549 + t1465;
t1015 = -t1080 * t1393 + t1081 * t1396;
t1014 = t1080 * t1396 + t1081 * t1393;
t999 = -t1065 * t1388 + t1067 * t1390;
t997 = t1065 * t1390 + t1067 * t1388;
t991 = t1200 * t1389 + t1391 * t999;
t990 = -t1200 * t1391 + t1389 * t999;
t979 = -t1031 * t1393 + t1033 * t1396;
t978 = -t1030 * t1393 + t1032 * t1396;
t977 = t1031 * t1396 + t1033 * t1393;
t976 = t1030 * t1396 + t1032 * t1393;
t973 = -t1135 * t1549 + t1412;
t972 = qJ(2) * t1271 - t1133 * t1549 + t1466;
t970 = -t1019 * t1393 + t1021 * t1396;
t969 = -t1018 * t1393 + t1020 * t1396;
t968 = t1019 * t1396 + t1021 * t1393;
t967 = t1018 * t1396 + t1020 * t1393;
t964 = -t1012 * t1393 + t1013 * t1396;
t962 = -qJ(2) * t1028 - t1476;
t960 = -t1006 * t1393 + t1007 * t1396;
t958 = -qJ(2) * t1052 - t1426;
t957 = -qJ(2) * t1050 - t1427;
t955 = -t1028 * t1549 + t1405;
t954 = -qJ(2) * t1264 - t1027 * t1549 + t1429;
t951 = -t1393 * t993 + t1396 * t994;
t949 = -t1393 * t990 + t1396 * t991;
t948 = t1393 * t991 + t1396 * t990;
t946 = -t1393 * t988 + t1396 * t989;
t942 = -qJ(2) * t1045 - t1408;
t932 = -t1052 * t1549 + t1415;
t931 = qJ(2) * t1164 - t1051 * t1549 + t1467;
t930 = -t1050 * t1549 + t1416;
t929 = qJ(2) * t1157 - t1049 * t1549 + t1468;
t928 = -t1045 * t1549 + t1417;
t927 = qJ(2) * t1148 - t1044 * t1549 + t1469;
t914 = -qJ(2) * t964 - t1413;
t913 = -qJ(2) * t951 - t1443;
t911 = -qJ(2) * t960 - t1414;
t906 = -t1393 * t920 + t1396 * t921;
t903 = -t1549 * t951 + t1420;
t902 = qJ(2) * t1009 - t1549 * t950 + t1470;
t901 = -qJ(2) * t946 - t1425;
t900 = -t1549 * t964 + t1418;
t899 = qJ(2) * t1047 - t1549 * t963 + t1471;
t898 = -t1549 * t960 + t1419;
t897 = qJ(2) * t1035 - t1549 * t959 + t1472;
t894 = -t1549 * t946 + t1421;
t893 = qJ(2) * t998 - t1549 * t945 + t1473;
t891 = -qJ(2) * t906 - t1404;
t890 = -t1549 * t906 + t1422;
t889 = qJ(2) * t925 - t1549 * t905 + t1474;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1361, 0, -t1362, 0, -t1449, -t1448, -t1438, -pkin(6) * t1438, 0, -t1361, t1362, 0, 0, 0, t1439, t1449, t1448, pkin(6) * t1439 + (-pkin(1) * t1394 + qJ(2) * t1397) * g(3), t1328 * t1394 + t1455, t1302 * t1394 + t1364 * t1397, t1317 * t1394 + t1379 * t1397, t1327 * t1394 - t1455, t1315 * t1394 - t1397 * t1492, -t1355 * t1394 + t1380, t1397 * t1232 - t1394 * t1214 - pkin(6) * (-t1316 * t1397 + t1357 * t1394), t1397 * t1233 - t1394 * t1213 - pkin(6) * (-t1318 * t1397 + t1359 * t1394), -pkin(2) * t1506 + t1394 * t1243 - pkin(6) * (-t1363 * t1394 + t1506), t1397 * t1191 - t1394 * t1162 - pkin(6) * (-t1247 * t1397 - t1326 * t1394), t1197 * t1394 + t1482, t1132 * t1394 + t1295 * t1397, t1186 * t1394 + t1276 * t1397, t1196 * t1394 - t1482, t1185 * t1394 - t1273 * t1397, t1209 * t1394 + t1380, t1397 * t1040 - t1394 * t1017 - pkin(6) * (-t1153 * t1397 + t1272 * t1394), t1397 * t1053 - t1394 * t1025 - pkin(6) * (-t1187 * t1397 + t1275 * t1394), t1397 * t1041 - t1394 * t973 - pkin(6) * (-t1133 * t1397 + t1271 * t1394), t1397 * t962 - t1394 * t955 - pkin(6) * (-t1027 * t1397 - t1264 * t1394), t1089 * t1394 + t1219 * t1397, t1058 * t1394 + t1147 * t1397, t1069 * t1394 + t1181 * t1397, t1088 * t1394 - t1217 * t1397, t1070 * t1394 + t1182 * t1397, t1097 * t1394 + t1215 * t1397, t1397 * t957 - t1394 * t930 - pkin(6) * (-t1049 * t1397 + t1157 * t1394), t1397 * t958 - t1394 * t932 - pkin(6) * (-t1051 * t1397 + t1164 * t1394), t1397 * t942 - t1394 * t928 - pkin(6) * (-t1044 * t1397 + t1148 * t1394), t1397 * t913 - t1394 * t903 - pkin(6) * (t1009 * t1394 - t1397 * t950), t1061 * t1397 + t1394 * t977, t1394 * t948 + t1397 * t997, t1054 * t1397 + t1394 * t967, t1060 * t1397 + t1394 * t976, t1055 * t1397 + t1394 * t968, t1014 * t1394 + t1086 * t1397, t1397 * t911 - t1394 * t898 - pkin(6) * (t1035 * t1394 - t1397 * t959), t1397 * t914 - t1394 * t900 - pkin(6) * (t1047 * t1394 - t1397 * t963), t1397 * t901 - t1394 * t894 - pkin(6) * (t1394 * t998 - t1397 * t945), t1397 * t891 - t1394 * t890 - pkin(6) * (t1394 * t925 - t1397 * t905); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1362, 0, t1361, 0, t1448, -t1449, t1459, pkin(6) * t1459, 0, -t1362, -t1361, 0, 0, 0, t1460, -t1448, t1449, pkin(6) * t1460 + (pkin(1) * t1397 + qJ(2) * t1394) * g(3), -t1328 * t1397 + t1456, -t1302 * t1397 + t1364 * t1394, -t1317 * t1397 + t1379 * t1394, -t1327 * t1397 - t1456, -t1315 * t1397 - t1394 * t1492, t1355 * t1397 + t1378, t1394 * t1232 + t1397 * t1214 + pkin(6) * (t1316 * t1394 + t1357 * t1397), t1394 * t1233 + t1397 * t1213 + pkin(6) * (t1318 * t1394 + t1359 * t1397), -pkin(2) * t1507 - t1397 * t1243 + pkin(6) * (-t1363 * t1397 - t1507), t1394 * t1191 + t1397 * t1162 + pkin(6) * (t1247 * t1394 - t1326 * t1397), -t1197 * t1397 + t1483, -t1132 * t1397 + t1295 * t1394, -t1186 * t1397 + t1276 * t1394, -t1196 * t1397 - t1483, -t1185 * t1397 - t1273 * t1394, -t1209 * t1397 + t1378, t1394 * t1040 + t1397 * t1017 + pkin(6) * (t1153 * t1394 + t1272 * t1397), t1394 * t1053 + t1397 * t1025 + pkin(6) * (t1187 * t1394 + t1275 * t1397), t1394 * t1041 + t1397 * t973 + pkin(6) * (t1133 * t1394 + t1271 * t1397), t1394 * t962 + t1397 * t955 + pkin(6) * (t1027 * t1394 - t1264 * t1397), -t1089 * t1397 + t1219 * t1394, -t1058 * t1397 + t1147 * t1394, -t1069 * t1397 + t1181 * t1394, -t1088 * t1397 - t1217 * t1394, -t1070 * t1397 + t1182 * t1394, -t1097 * t1397 + t1215 * t1394, t1394 * t957 + t1397 * t930 + pkin(6) * (t1049 * t1394 + t1157 * t1397), t1394 * t958 + t1397 * t932 + pkin(6) * (t1051 * t1394 + t1164 * t1397), t1394 * t942 + t1397 * t928 + pkin(6) * (t1044 * t1394 + t1148 * t1397), t1394 * t913 + t1397 * t903 + pkin(6) * (t1009 * t1397 + t1394 * t950), t1061 * t1394 - t1397 * t977, t1394 * t997 - t1397 * t948, t1054 * t1394 - t1397 * t967, t1060 * t1394 - t1397 * t976, t1055 * t1394 - t1397 * t968, -t1014 * t1397 + t1086 * t1394, t1394 * t911 + t1397 * t898 + pkin(6) * (t1035 * t1397 + t1394 * t959), t1394 * t914 + t1397 * t900 + pkin(6) * (t1047 * t1397 + t1394 * t963), t1394 * t901 + t1397 * t894 + pkin(6) * (t1394 * t945 + t1397 * t998), t1394 * t891 + t1397 * t890 + pkin(6) * (t1394 * t905 + t1397 * t925); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1367, t1368, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1354, t1348, t1288, t1321, t1303, t1322, t1314, t1319, 0, t1211, t1212, t1207, t1161, t1199, t1134, t1189, t1198, t1188, t1210, t1016, t1024, t972, t954, t1091, t1059, t1071, t1090, t1072, t1098, t929, t931, t927, t902, t979, t949, t969, t978, t970, t1015, t897, t899, t893, t889; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1398, 0, 0, -g(3), -t1367, 0, 0, -qJDD(1), t1398, 0, 0, 0, -t1340, 0, g(3), qJ(2) * g(3), t1374, t1364, t1379, -t1374, -t1492, qJDD(3), t1232, t1233, -t1545, t1191, t1510, t1295, t1276, -t1510, -t1273, qJDD(3), t1040, t1053, t1041, t962, t1219, t1147, t1181, -t1217, t1182, t1215, t957, t958, t942, t913, t1061, t997, t1054, t1060, t1055, t1086, t911, t914, t901, t891; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1398, 0, qJDD(1), 0, g(3), 0, -t1368, 0, 0, -t1398, -qJDD(1), 0, 0, 0, -t1332, -g(3), 0, pkin(1) * g(3), -t1328, -t1302, -t1317, -t1327, -t1315, t1355, t1214, t1213, -t1243, t1162, -t1197, -t1132, -t1186, -t1196, -t1185, -t1209, t1017, t1025, t973, t955, -t1089, -t1058, -t1069, -t1088, -t1070, -t1097, t930, t932, t928, t903, -t977, -t948, -t967, -t976, -t968, -t1014, t898, t900, t894, t890; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1367, t1368, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1354, t1348, t1288, t1321, t1303, t1322, t1314, t1319, 0, t1211, t1212, t1207, t1161, t1199, t1134, t1189, t1198, t1188, t1210, t1016, t1024, t972, t954, t1091, t1059, t1071, t1090, t1072, t1098, t929, t931, t927, t902, t979, t949, t969, t978, t970, t1015, t897, t899, t893, t889; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t1340, -t1332, 0, t1321, t1303, t1322, t1314, t1319, 0, -pkin(7) * t1316 - t1498, -pkin(7) * t1318 - t1313, pkin(7) * t1360 - t1247, -pkin(7) * t1247, t1199, t1134, t1189, t1198, t1188, t1210, -pkin(7) * t1153 + t1465, -pkin(7) * t1187 + t1464, -pkin(7) * t1133 + t1466, -pkin(7) * t1027 + t1429, t1091, t1059, t1071, t1090, t1072, t1098, -pkin(7) * t1049 + t1468, -pkin(7) * t1051 + t1467, -pkin(7) * t1044 + t1469, -pkin(7) * t950 + t1470, t979, t949, t969, t978, t970, t1015, -pkin(7) * t959 + t1472, -pkin(7) * t963 + t1471, -pkin(7) * t945 + t1473, -pkin(7) * t905 + t1474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1398, 0, 0, 0, t1340, 0, -g(3), 0, -t1374, -t1364, -t1379, t1374, t1492, -qJDD(3), t1400, t1434, t1545, -t1547, -t1510, -t1295, -t1276, t1510, t1273, -qJDD(3), t1424 + t1490, t1333 + t1435, t1475, t1476, -t1219, -t1147, -t1181, t1217, -t1182, -t1215, t1427, t1426, t1408, t1443, -t1061, -t997, -t1054, -t1060, -t1055, -t1086, t1414, t1413, t1425, t1404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1398, qJDD(1), 0, 0, 0, t1332, g(3), 0, 0, t1328, t1302, t1317, t1327, t1315, -t1355, pkin(7) * t1320 - t1445, pkin(7) * t1323 - t1444, t1243, pkin(7) * t1248 + t1546, t1197, t1132, t1186, t1196, t1185, t1209, pkin(7) * t1154 - t1411, pkin(7) * t1190 - t1410, pkin(7) * t1135 - t1412, pkin(7) * t1028 - t1405, t1089, t1058, t1069, t1088, t1070, t1097, pkin(7) * t1050 - t1416, pkin(7) * t1052 - t1415, pkin(7) * t1045 - t1417, pkin(7) * t951 - t1420, t977, t948, t967, t976, t968, t1014, pkin(7) * t960 - t1419, pkin(7) * t964 - t1418, pkin(7) * t946 - t1421, pkin(7) * t906 - t1422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1358, -t1357, t1366, t1481, t1370, -t1481, 0, -t1326, -t1305, 0, t1268, t1205, t1255, t1266, t1254, t1283, t1179, t1195, t1074, -qJ(4) * t1093, t1175, t1111, t1130, t1174, t1131, t1193, t981, t983, t966, t939, t1033, t991, t1020, t1032, t1021, t1081, t918, t923, t909, t896; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1480, t1359, t1372, t1428, t1365, -t1480, t1326, 0, -t1306, 0, t1267, t1203, t1252, t1265, t1251, t1282, t1142, t1144, t1064, t1083, t1173, t1110, t1128, t1172, t1129, t1192, t974, t975, t961, t933, t1031, t990, t1018, t1030, t1019, t1080, t915, t916, t904, t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1374, t1364, t1379, -t1374, -t1492, qJDD(3), t1305, t1306, 0, 0, t1510, t1295, t1276, -t1510, -t1273, qJDD(3), t1335 + t1446, t1334 + t1477, t1202, t1092, t1219, t1147, t1181, -t1217, t1182, t1215, t1450, t1451, t1440, t1488, t1061, t997, t1054, t1060, t1055, t1086, t1442, t1441, t1447, t1423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1301, -t1272, t1561, t1535, t1329, -t1535, 0, -t1264, t1176, 0, t1220, t1149, t1183, t1218, t1184, t1216, t1079, t1082, t995, -qJ(5) * t1009, t1063, t999, t1056, t1062, t1057, t1087, t947, t953, t912, t907; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1337, t1275, t1330, -t1300, t1292, -t1337, t1264, 0, t1177, 0, -t1270, -t1269, -t1227, t1270, t1223, -t1300, t1039, t1042, -pkin(4) * t1148, -pkin(4) * t1009, -t1201, -t1200, -t1126, t1201, t1123, -t1296, t952, t956, t971, t910; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1510, t1295, t1276, -t1510, -t1273, qJDD(3), -t1176, -t1177, 0, 0, t1219, t1147, t1181, -t1217, t1182, t1215, t1478, t1479, t1452, t1495, t1061, t997, t1054, t1060, t1055, t1086, t1454, t1453, t1458, t1430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1287, t1224, t1554, t1516, t1280, -t1516, 0, t1140, t1077, 0, t1119, t1067, t1114, t1117, t1115, t1160, t1023, t1029, t937, -pkin(8) * t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1515, t1226, t1281, t1436, t1230, -t1515, -t1140, 0, t1078, 0, t1118, t1065, t1112, t1116, t1113, t1159, t1002, t1004, t935, t941; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1270, t1269, t1227, -t1270, -t1223, t1300, -t1077, -t1078, 0, 0, t1201, t1200, t1126, -t1201, -t1123, t1296, t1431, t1409, t1543, t1548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1171, -t1122, t1556, t1241, t1239, -t1241, 0, t1085, t985, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1519, t1555, t1240, t1170, t1167, -t1519, -t1085, 0, t986, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1201, t1200, t1126, -t1201, -t1123, t1296, -t985, -t986, 0, 0;];
m_new_reg  = t1;