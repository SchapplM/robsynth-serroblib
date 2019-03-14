% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:38
% EndTime: 2019-03-09 07:03:47
% DurationCPUTime: 59.49s
% Computational Cost: add. (42645->1036), mult. (85342->1347), div. (0->0), fcn. (90961->10), ass. (0->764)
t964 = sin(qJ(4));
t966 = cos(qJ(5));
t1296 = t966 * t964;
t963 = sin(qJ(5));
t967 = cos(qJ(4));
t1303 = t963 * t967;
t899 = t1296 + t1303;
t965 = sin(qJ(3));
t858 = t899 * t965;
t962 = sin(qJ(6));
t1312 = t962 * t858;
t1401 = cos(qJ(6));
t1295 = t966 * t967;
t1301 = t964 * t965;
t861 = t1295 * t965 - t1301 * t963;
t839 = t1401 * t861;
t1473 = t839 - t1312;
t1436 = -t1473 / 0.2e1;
t1500 = t1473 / 0.2e1;
t1556 = -t1436 - t1500;
t1125 = t1401 * t858 + t962 * t861;
t1438 = t1125 / 0.2e1;
t1439 = -t1125 / 0.2e1;
t1555 = -t1438 - t1439;
t1445 = -pkin(9) - pkin(8);
t1192 = t1445 * t964;
t1474 = t1445 * t967;
t1510 = t966 * t1192 + t963 * t1474;
t1521 = -t899 * pkin(10) + t1510;
t1533 = t962 * t1521;
t1138 = -t1533 / 0.2e1;
t1117 = t963 * t1192;
t1494 = t966 * t1474;
t787 = t1494 - t1117;
t895 = t963 * t964 - t1295;
t679 = t895 * pkin(10) + t787;
t1542 = t1401 * t679;
t651 = t1542 / 0.2e1;
t1551 = t1138 + t651;
t1554 = 0.2e1 * t1551;
t1309 = t962 * t895;
t881 = t1401 * t899;
t1472 = t881 - t1309;
t1429 = t1472 / 0.2e1;
t1397 = t858 * pkin(10);
t968 = cos(qJ(3));
t1100 = -t968 * pkin(3) - t965 * pkin(8);
t945 = -cos(pkin(11)) * pkin(1) - pkin(2);
t1001 = t1100 + t945;
t944 = sin(pkin(11)) * pkin(1) + pkin(7);
t915 = t968 * t944;
t1188 = t967 * t915;
t771 = t1001 * t964 + t1188;
t714 = -pkin(9) * t1301 + t771;
t677 = t966 * t714;
t1298 = t965 * t967;
t1202 = pkin(9) * t1298;
t1325 = t944 * t964;
t863 = t967 * t1001;
t988 = -t1202 + t863 + (-pkin(4) - t1325) * t968;
t434 = t963 * t988 + t677;
t349 = t434 - t1397;
t1321 = t962 * t349;
t1396 = t861 * pkin(10);
t1305 = t963 * t714;
t637 = t966 * t988;
t433 = -t637 + t1305;
t1060 = -t433 - t1396;
t1390 = t968 * pkin(5);
t1009 = t1060 - t1390;
t328 = t1401 * t1009;
t201 = -t328 + t1321;
t1444 = t201 / 0.2e1;
t1532 = t1401 * t1521;
t1543 = t962 * t679;
t1544 = t1532 + t1543;
t1178 = t1401 * t349;
t990 = t962 * t1009;
t202 = t1178 + t990;
t428 = t1542 - t1533;
t748 = t1401 * t895 + t962 * t899;
t1552 = t202 * t1429 + t1444 * t748 + t1555 * t1544 + t1556 * t428;
t1204 = qJD(4) + qJD(5);
t1434 = -t1532 / 0.2e1;
t1522 = t1434 - t1543 / 0.2e1;
t1548 = 0.2e1 * t1522;
t1403 = -t968 / 0.2e1;
t1547 = -t1542 / 0.2e1;
t1106 = t1532 / 0.2e1;
t1540 = t1434 + t1106;
t1546 = qJD(3) * t1540;
t1545 = qJD(6) * t1540;
t1426 = -t1510 / 0.2e1;
t1541 = qJD(5) + qJD(6);
t1414 = -t899 / 0.2e1;
t955 = -pkin(4) * t967 - pkin(3);
t1501 = -t955 / 0.2e1;
t914 = t965 * t944;
t872 = pkin(4) * t1301 + t914;
t1539 = t1403 * t787 - t1414 * t872 - t1501 * t861;
t1538 = qJD(4) * t1554;
t1537 = qJD(4) * t1548;
t1300 = t964 * t968;
t1189 = t944 * t1298;
t1389 = t968 * pkin(8);
t1394 = t965 * pkin(3);
t916 = -t1389 + t1394;
t901 = t964 * t916;
t792 = t901 - t1189;
t763 = -pkin(9) * t1300 + t792;
t1304 = t963 * t763;
t1293 = t967 * t968;
t1393 = t965 * pkin(4);
t904 = t967 * t916;
t906 = t964 * t914;
t791 = t906 + t904;
t728 = -pkin(9) * t1293 + t1393 + t791;
t684 = t966 * t728;
t1097 = t684 / 0.2e1 - t1304 / 0.2e1;
t256 = t1097 - t1539;
t1415 = -t895 / 0.2e1;
t1536 = -t1415 * t872 + t1426 * t968 - t1501 * t858;
t259 = t1097 + t1539;
t1190 = t964 * t915;
t770 = -t863 + t1190;
t713 = -t770 - t1202;
t1306 = t963 * t713;
t450 = -t677 - t1306;
t1516 = t450 + t434;
t1517 = qJD(6) + t1204;
t1256 = qJD(1) * t968;
t1488 = t1473 * t1256;
t919 = t963 * t1300;
t862 = t1293 * t966 - t919;
t1310 = t962 * t862;
t860 = t899 * t968;
t838 = t1401 * t860;
t1096 = -t838 / 0.2e1 - t1310 / 0.2e1;
t1493 = t968 * t1472;
t1503 = t1493 / 0.2e1 + t1096;
t1520 = qJD(3) * t1503;
t1531 = t1488 + t1520;
t1504 = -t1493 / 0.2e1 + t1096;
t1519 = qJD(3) * t1504;
t1530 = t1519 - t1488;
t683 = t963 * t728;
t743 = t966 * t763;
t1133 = -t683 / 0.2e1 - t743 / 0.2e1;
t257 = t1133 + t1536;
t258 = t1133 - t1536;
t1527 = t1204 * t1510;
t1297 = t966 * t713;
t451 = t1297 - t1305;
t1526 = t451 / 0.2e1;
t1525 = t1510 / 0.2e1;
t1230 = t1503 * qJD(2);
t1518 = t1504 * qJD(1);
t1502 = -t748 / 0.2e1;
t726 = pkin(5) * t858 + t872;
t844 = pkin(5) * t895 + t955;
t1478 = t1439 * t844 + t1502 * t726;
t1092 = t1494 / 0.2e1;
t1515 = t451 + t433;
t1430 = -t1472 / 0.2e1;
t367 = t450 + t1397;
t1177 = t1401 * t367;
t368 = t451 - t1396;
t1319 = t962 * t368;
t213 = t1177 - t1319;
t1176 = t1401 * t368;
t1320 = t962 * t367;
t214 = t1176 + t1320;
t1513 = -t1430 * t213 - t1502 * t214 + t1552;
t1036 = t962 * t1060;
t205 = -t1036 - t1178;
t1008 = t1401 * t1060;
t206 = t1008 - t1321;
t1512 = -t1430 * t205 - t1502 * t206 + t1552;
t1413 = t899 / 0.2e1;
t530 = -t1413 * t858 + t1415 * t861;
t1511 = t1204 * t530;
t1484 = t1125 ^ 2 - t1473 ^ 2;
t1509 = qJD(1) * t1484;
t1485 = -t1472 ^ 2 + t748 ^ 2;
t1508 = qJD(3) * t1485;
t1507 = t1125 * t1502;
t1506 = t1473 * t1502;
t1212 = t748 * qJD(6);
t1505 = -t1204 * t748 - t1212;
t1466 = t1125 * qJD(6);
t1491 = t1204 * t1125;
t1487 = t1491 + t1466;
t1423 = -t858 / 0.2e1;
t1422 = -t861 / 0.2e1;
t1359 = t1125 * t726;
t1499 = t1204 * t861;
t274 = t1125 * t1430 + t1506;
t1498 = t274 * qJD(6);
t466 = t726 * t1473;
t1251 = qJD(3) * t1472;
t1157 = t748 * t1251;
t1123 = t1204 * t899;
t1497 = t895 * t1123;
t1492 = t1125 * t1256;
t1260 = qJD(1) * t1473;
t1162 = t1125 * t1260;
t1490 = t1204 * t1473;
t1486 = t1544 * t1556 + t1555 * t428;
t1249 = qJD(3) * t899;
t1483 = -qJD(1) * t530 + t895 * t1249;
t1482 = -qJD(1) * t274 + t1157;
t1257 = qJD(1) * t861;
t1481 = qJD(3) * t530 - t858 * t1257;
t1480 = -qJD(3) * t274 + t1162;
t1479 = 0.2e1 * t964;
t810 = t962 * pkin(5);
t1477 = qJD(6) * t810;
t1476 = t1204 * t1472;
t1475 = t1204 * t858;
t1105 = -t839 / 0.2e1;
t639 = t1105 + t839 / 0.2e1;
t1218 = t639 * qJD(2);
t987 = t1533 / 0.2e1;
t243 = t1547 + t987 + t1551;
t1032 = -t1319 / 0.2e1 + t1177 / 0.2e1;
t1109 = t1178 / 0.2e1;
t1179 = t966 * pkin(4) + pkin(5);
t1101 = t962 * t1179;
t1165 = t1401 * t963;
t875 = pkin(4) * t1165 + t1101;
t1330 = t875 * t968;
t66 = -t1330 / 0.2e1 + t1109 + t990 / 0.2e1 + t1032;
t1464 = qJD(1) * t66 + qJD(3) * t243 + qJD(4) * t875 - t1218;
t242 = t1138 + t987;
t1194 = t1390 / 0.2e1;
t1110 = -t1178 / 0.2e1;
t1289 = -t990 / 0.2e1 + t1110;
t1061 = t1194 * t962 + t1289;
t973 = t1036 / 0.2e1 + t1109;
t72 = t973 + t1061;
t1463 = qJD(1) * t72 - qJD(3) * t242 - qJD(4) * t810 + t1218;
t734 = t919 / 0.2e1 + (-t1295 / 0.2e1 + t1415) * t968;
t1213 = t734 * qJD(2);
t1248 = qJD(3) * t955;
t1462 = qJD(1) * t257 + t1248 * t895 + t1213;
t1049 = -t1303 / 0.2e1 - t1296 / 0.2e1;
t731 = (t1413 + t1049) * t968;
t1215 = t731 * qJD(2);
t1461 = qJD(1) * t256 - t1248 * t899 + t1215;
t1302 = t964 * t858;
t1404 = -t967 / 0.2e1;
t1405 = t966 / 0.2e1;
t210 = (-t1302 / 0.2e1 + (t1404 * t895 + t1405) * t965) * pkin(4) + t256;
t1400 = pkin(4) * t964;
t724 = t1400 * t895 + t899 * t955;
t1460 = qJD(1) * t210 - qJD(3) * t724 + t1215;
t1409 = -t964 / 0.2e1;
t1411 = -t963 / 0.2e1;
t209 = (t861 * t1409 + (t1404 * t899 + t1411) * t965) * pkin(4) + t257;
t725 = t1400 * t899 - t895 * t955;
t1459 = qJD(1) * t209 - qJD(3) * t725 + t1213;
t1402 = t968 / 0.2e1;
t1144 = t748 * t1402;
t1166 = t1401 * t862;
t1311 = t962 * t860;
t1273 = t1311 / 0.2e1 - t1166 / 0.2e1;
t461 = -t1144 + t1273;
t1223 = t461 * qJD(2);
t1398 = pkin(5) * t899;
t567 = t844 * t748;
t376 = -t1398 * t1472 + t567;
t1371 = t1544 * t968;
t1408 = -t965 / 0.2e1;
t497 = t743 + t683;
t394 = -pkin(10) * t860 + t497;
t1174 = t1401 * t394;
t1107 = -t1174 / 0.2e1;
t1392 = t965 * pkin(5);
t496 = t684 - t1304;
t380 = -pkin(10) * t862 + t1392 + t496;
t1318 = t962 * t380;
t1031 = -t1318 / 0.2e1 + t1107;
t971 = t1031 - t1478;
t52 = -t1371 / 0.2e1 + (t1408 * t962 + t1414 * t1473 + t1422 * t1472) * pkin(5) + t971;
t1458 = qJD(1) * t52 + qJD(3) * t376 + t1223;
t566 = t844 * t1472;
t375 = -t1398 * t748 - t566;
t1200 = t1401 / 0.2e1;
t1103 = t965 * t1200;
t1373 = t428 * t968;
t1175 = t1401 * t380;
t1317 = t962 * t394;
t1030 = -t1317 / 0.2e1 + t1175 / 0.2e1;
t1448 = t1430 * t726 + t1436 * t844;
t972 = t1030 + t1448;
t53 = t1373 / 0.2e1 + (t1125 * t1414 + t1422 * t748 + t1103) * pkin(5) + t972;
t1457 = qJD(1) * t53 + qJD(3) * t375 + t1230;
t849 = t1398 + t1400;
t348 = t1472 * t849 - t567;
t1142 = t875 * t1408;
t1424 = t849 / 0.2e1;
t1203 = pkin(4) * t1298;
t1399 = pkin(5) * t861;
t780 = t1203 + t1399;
t1428 = t780 / 0.2e1;
t995 = t1402 * t1544 + t1424 * t1473 + t1428 * t1472;
t49 = t1142 + t971 - t995;
t1456 = qJD(1) * t49 - qJD(3) * t348 + t1223;
t347 = t748 * t849 + t566;
t1308 = t962 * t963;
t917 = t1401 * t1179;
t874 = pkin(4) * t1308 - t917;
t1143 = t874 * t1408;
t997 = t1125 * t1424 + t1403 * t428 + t1428 * t748;
t48 = t1143 + t972 - t997;
t1455 = qJD(1) * t48 - qJD(3) * t347 + t1230;
t1051 = t1429 * t1473 + t1507;
t1052 = t1436 * t1472 - t1507;
t130 = t1051 + t1052;
t1239 = t130 * qJD(2);
t1201 = -t1401 / 0.2e1;
t657 = t838 + t1310;
t1437 = -t657 / 0.2e1;
t660 = t1166 - t1311;
t1002 = (t1201 * t660 + t1437 * t962) * pkin(5);
t8 = t1002 + t1512;
t1454 = -t8 * qJD(1) + t1239;
t1419 = -t875 / 0.2e1;
t1420 = t874 / 0.2e1;
t1050 = t1419 * t657 + t1420 * t660;
t6 = t1050 + t1513;
t1453 = -t6 * qJD(1) + t1239;
t1443 = t202 / 0.2e1;
t1057 = -t1125 * t1443 + t1500 * t201;
t20 = t1403 * t780 + t1439 * t213 + t1500 * t214 + t1057;
t1452 = t20 * qJD(1);
t19 = t1390 * t1422 + t1439 * t205 + t1500 * t206 + t1057;
t1451 = t19 * qJD(1);
t1098 = -t677 / 0.2e1 - t1306 / 0.2e1;
t1447 = t1517 * t1503;
t1446 = t1517 * t1504;
t1442 = t434 / 0.2e1;
t1441 = -t637 / 0.2e1;
t1427 = t787 / 0.2e1;
t1425 = t844 / 0.2e1;
t1140 = t860 / 0.2e1;
t1421 = t861 / 0.2e1;
t1104 = -t881 / 0.2e1;
t1307 = t962 * t966;
t883 = (t1165 + t1307) * pkin(4);
t1418 = -t883 / 0.2e1;
t1417 = t883 / 0.2e1;
t1164 = t1401 * t966;
t884 = (t1164 - t1308) * pkin(4);
t1416 = -t884 / 0.2e1;
t1412 = t962 / 0.2e1;
t1410 = t963 / 0.2e1;
t1407 = t965 / 0.2e1;
t1406 = -t966 / 0.2e1;
t1391 = t968 * pkin(4);
t1388 = pkin(4) * qJD(4);
t1387 = pkin(4) * qJD(5);
t1386 = pkin(5) * qJD(5);
t97 = -t1125 * t1399 + t205 * t968 - t466;
t1384 = qJD(1) * t97;
t98 = t1399 * t1473 + t206 * t968 - t1359;
t1383 = qJD(1) * t98;
t99 = -t1125 * t780 + t213 * t968 - t466;
t1382 = qJD(1) * t99;
t1087 = -t1125 * t201 - t1473 * t202;
t23 = -t1125 * t206 - t1473 * t205 + t1087;
t1376 = t23 * qJD(1);
t24 = -t1125 * t214 - t1473 * t213 + t1087;
t1375 = t24 * qJD(1);
t221 = t1175 - t1317;
t222 = t1174 + t1318;
t25 = -t1125 * t222 - t1473 * t221 + t201 * t660 - t202 * t657;
t1374 = t25 * qJD(1);
t943 = pkin(4) * t1300;
t873 = t915 + t943;
t727 = pkin(5) * t860 + t873;
t45 = t1473 * t727 - t202 * t965 + t222 * t968 + t726 * t660;
t1367 = t45 * qJD(1);
t1365 = t1473 * t874;
t1364 = t1473 * t962;
t1363 = t657 * t1472;
t1361 = t1125 * t875;
t1360 = t1125 * t962;
t1356 = t1472 * t962;
t1354 = t770 * t968;
t1353 = t771 * t968;
t1345 = t791 * t965;
t1344 = t792 * t965;
t1341 = t860 * t899;
t1134 = t1526 + t433 / 0.2e1;
t1135 = t1442 + t450 / 0.2e1;
t1136 = -t1298 / 0.2e1;
t87 = t1134 * t861 - t1135 * t858 + t1136 * t1391;
t1339 = t87 * qJD(1);
t1338 = t872 * t858;
t1335 = t874 * t1125;
t1334 = t874 * t748;
t1333 = t874 * t968;
t1332 = t875 * t1473;
t1331 = t875 * t1472;
t1329 = t883 * t968;
t1328 = t884 * t968;
t92 = -t1515 * t858 - t1516 * t861;
t1327 = t92 * qJD(1);
t93 = t433 * t862 - t434 * t860 - t496 * t861 - t497 * t858;
t1326 = t93 * qJD(1);
t958 = t964 ^ 2;
t1322 = t958 * t968;
t951 = t965 * t968;
t959 = t965 ^ 2;
t1294 = t967 * t959;
t1139 = t1321 / 0.2e1;
t1288 = -t328 / 0.2e1 + t1139;
t1287 = -t1036 / 0.2e1 + t1110;
t1286 = t1139 - t1008 / 0.2e1;
t391 = t858 * t895 - t899 * t861;
t1285 = t1204 * t391;
t1282 = t1204 * t639;
t1276 = t1204 * t734;
t733 = -t919 / 0.2e1 + (t1295 / 0.2e1 + t1415) * t968;
t1275 = t1204 * t733;
t960 = t967 ^ 2;
t935 = t960 - t958;
t961 = t968 ^ 2;
t936 = t961 - t959;
t100 = t1473 * t780 + t214 * t968 - t1359;
t1271 = qJD(1) * t100;
t134 = -t201 * t968 - t1359;
t1270 = qJD(1) * t134;
t135 = -t202 * t968 - t466;
t1269 = qJD(1) * t135;
t739 = t872 * t861;
t276 = -t1203 * t858 + t450 * t968 - t739;
t1268 = qJD(1) * t276;
t277 = t1203 * t861 + t451 * t968 - t1338;
t1267 = qJD(1) * t277;
t309 = -t433 * t968 - t1338;
t1266 = qJD(1) * t309;
t310 = -t434 * t968 - t739;
t1265 = qJD(1) * t310;
t402 = -t1473 * t965 + t660 * t968;
t1264 = qJD(1) * t402;
t602 = -t1325 * t959 - t1354;
t1263 = qJD(1) * t602;
t603 = -t1294 * t944 - t1353;
t1262 = qJD(1) * t603;
t1261 = qJD(1) * t1125;
t662 = -t861 * t965 + t862 * t968;
t1259 = qJD(1) * t662;
t1258 = qJD(1) * t858;
t1254 = qJD(3) * t130;
t1252 = qJD(3) * t748;
t1250 = qJD(3) * t844;
t1247 = qJD(3) * t964;
t1246 = qJD(3) * t965;
t1245 = qJD(3) * t967;
t1244 = qJD(3) * t968;
t1243 = qJD(4) * t964;
t1242 = qJD(4) * t967;
t1241 = qJD(4) * t968;
t1240 = qJD(5) * t955;
t183 = -t434 * t965 + t497 * t968 + t873 * t861 + t872 * t862;
t1238 = t183 * qJD(1);
t1195 = t1391 / 0.2e1;
t234 = t1441 + (t1195 + t713 / 0.2e1) * t966;
t1236 = t234 * qJD(1);
t255 = -t1125 * t660 - t1473 * t657;
t1234 = t255 * qJD(1);
t280 = (t1345 - t1354) * t967 + (t1344 + t1353) * t964;
t1233 = t280 * qJD(1);
t401 = t1125 * t965 - t657 * t968;
t1232 = t401 * qJD(1);
t1145 = t748 * t1403;
t458 = t1145 - t1273;
t1227 = t458 * qJD(1);
t1226 = t458 * qJD(3);
t459 = -t1144 - t1273;
t315 = t459 * qJD(1);
t1225 = t459 * qJD(3);
t460 = t1145 + t1273;
t1224 = t460 * qJD(3);
t303 = t461 * qJD(3);
t472 = t770 * t965 + (t791 - 0.2e1 * t906) * t968;
t1222 = t472 * qJD(1);
t473 = t792 * t968 + (-t771 + 0.2e1 * t1188) * t965;
t1221 = t473 * qJD(1);
t504 = -t862 * t858 - t861 * t860;
t1220 = t504 * qJD(1);
t636 = t1473 * qJD(6);
t661 = t965 * t858 - t860 * t968;
t1216 = t661 * qJD(1);
t732 = (t1414 + t1049) * t968;
t696 = t732 * qJD(1);
t697 = t733 * qJD(1);
t1214 = t733 * qJD(3);
t698 = t734 * qJD(3);
t1211 = t1472 * qJD(6);
t886 = (t958 / 0.2e1 - t960 / 0.2e1) * t965;
t1210 = t886 * qJD(4);
t898 = t936 * t964;
t1209 = t898 * qJD(1);
t900 = t961 * t967 - t1294;
t1208 = t900 * qJD(1);
t1207 = t936 * qJD(1);
t1206 = t965 * qJD(1);
t1205 = t965 * qJD(4);
t1199 = t1400 / 0.2e1;
t1198 = t1399 / 0.2e1;
t1197 = t1398 / 0.2e1;
t1196 = t1393 / 0.2e1;
t1187 = t1204 * t130;
t131 = -t1051 + t1052;
t192 = t1125 * t748 - t1472 * t1473;
t1186 = t192 * qJD(6) + t1204 * t131;
t271 = t1438 * t1472 + t1500 * t748;
t1185 = t1204 * t271 - t1498;
t272 = -t1125 * t1429 + t1506;
t1184 = t1204 * t272 + t1498;
t1183 = -t458 * qJD(6) - t1204 * t459;
t1182 = -t460 * qJD(6) - t1204 * t461;
t524 = 0.2e1 * t1105 + t1312;
t1181 = t524 * qJD(6) - t1490;
t1180 = t1204 * t524 - t636;
t1173 = t1401 * t1473;
t1172 = t1401 * t1125;
t1167 = t1401 * t748;
t1159 = t861 * t1256;
t1158 = t967 * t1206;
t1156 = t748 * t1246;
t1154 = t895 * t1246;
t1153 = t964 * t1245;
t1152 = t964 * t1241;
t1151 = t967 * t1241;
t1150 = t945 * t1206;
t1149 = t945 * t1256;
t1148 = t964 * t1242;
t1147 = t965 * t1244;
t940 = t968 * t1206;
t1146 = t965 * t1245;
t1141 = -t860 / 0.2e1;
t1132 = t1525 + t1426;
t1131 = -t787 / 0.2e1 + t1427;
t1130 = t1401 * qJD(5);
t1129 = t1401 * qJD(6);
t1128 = -qJD(4) / 0.2e1 - qJD(5) / 0.2e1;
t1127 = pkin(4) * t1204;
t1126 = (t958 + t960) * t968;
t1121 = t1204 * t968;
t1120 = t967 * t1196;
t1119 = -qJD(4) + t1256;
t1112 = t964 * t1146;
t1111 = t959 * t1148;
t1108 = -t1176 / 0.2e1;
t1102 = -t1164 / 0.2e1;
t1093 = -qJD(5) + t1119;
t18 = t1403 * t727 + t1407 * t726 + t1439 * t221 + t1443 * t660 + t1444 * t657 + t1500 * t222;
t30 = -t201 * t221 + t202 * t222 + t726 * t727;
t1090 = t30 * qJD(1) + t18 * qJD(2);
t28 = t1399 * t726 - t201 * t205 + t202 * t206;
t1089 = t28 * qJD(1) + t19 * qJD(2);
t29 = -t201 * t213 + t202 * t214 + t726 * t780;
t1088 = t29 * qJD(1) + t20 * qJD(2);
t1085 = -t791 * t964 + t792 * t967;
t101 = t1203 * t872 - t433 * t450 + t434 * t451;
t1084 = t101 * qJD(1) + t87 * qJD(2);
t102 = -t433 * t496 + t434 * t497 + t872 * t873;
t82 = t1140 * t433 + t1403 * t873 + t1407 * t872 + t1421 * t497 + t1423 * t496 + t1442 * t862;
t1083 = t102 * qJD(1) + t82 * qJD(2);
t230 = t1125 * t657 + t1473 * t660 - t951;
t1082 = t18 * qJD(1) + t230 * qJD(2);
t44 = -t1125 * t727 + t201 * t965 + t221 * t968 - t726 * t657;
t1081 = t44 * qJD(1);
t470 = t858 * t860 + t861 * t862 - t951;
t1080 = t82 * qJD(1) + t470 * qJD(2);
t1025 = (t1406 * t862 + t1411 * t860) * pkin(4);
t40 = t1131 * t861 + t1132 * t858 + t1134 * t895 + t1135 * t899 + t1025;
t1079 = t40 * qJD(1);
t1074 = t1119 * t965;
t182 = t433 * t965 + t496 * t968 - t873 * t858 - t872 * t860;
t1073 = t182 * qJD(1);
t260 = (t1344 / 0.2e1 + t1353 / 0.2e1) * t967 + (-t1345 / 0.2e1 + t1354 / 0.2e1) * t964 + (t959 / 0.2e1 - t961 / 0.2e1) * t944;
t311 = t944 ^ 2 * t951 - t770 * t791 + t771 * t792;
t1072 = t311 * qJD(1) + t260 * qJD(2);
t797 = t1126 * t965 - t951;
t1071 = -t260 * qJD(1) - t797 * qJD(2);
t1005 = (-t1364 / 0.2e1 + t1172 / 0.2e1) * pkin(5);
t975 = -t1335 / 0.2e1 - t1332 / 0.2e1 + t1473 * t1417 + t1125 * t1416;
t110 = t1005 - t975;
t1003 = (-t1356 / 0.2e1 + t1167 / 0.2e1) * pkin(5);
t974 = -t1334 / 0.2e1 - t1331 / 0.2e1 + t1472 * t1417 + t748 * t1416;
t140 = t1003 - t974;
t1070 = qJD(1) * t110 + qJD(3) * t140;
t58 = qJD(3) * t131 + t1509;
t80 = qJD(1) * t131 + t1508;
t1069 = qJD(3) * t192 + t1509;
t1068 = qJD(1) * t192 + t1508;
t970 = t1391 * t963 + t1098;
t232 = t970 - t1098;
t782 = -t1494 / 0.2e1 + t1092;
t1065 = t232 * qJD(1) + t782 * qJD(3);
t563 = t858 ^ 2 - t861 ^ 2;
t269 = qJD(1) * t563 + qJD(3) * t391;
t678 = t895 ^ 2 - t899 ^ 2;
t278 = qJD(1) * t391 + qJD(3) * t678;
t745 = t1104 + t881 / 0.2e1;
t438 = qJD(1) * t639 + qJD(3) * t745;
t1063 = -t917 / 0.2e1 + pkin(5) * t1201;
t1062 = t1389 / 0.2e1 - t1394 / 0.2e1;
t1038 = t1062 * t964;
t778 = t901 / 0.2e1 - t1038;
t1059 = pkin(3) * t1245 - qJD(1) * t778;
t1037 = t1062 * t967;
t779 = -t904 / 0.2e1 + t1037;
t1058 = pkin(3) * t1247 - qJD(1) * t779;
t1056 = t1419 * t222 + t1420 * t221;
t1055 = t1405 * t496 + t1410 * t497;
t1053 = t1419 * t660 + t1437 * t874;
t1048 = t967 * t1074;
t169 = qJD(3) * t271 + t1261 * t1473;
t197 = qJD(1) * t271 + t1252 * t1472;
t171 = qJD(3) * t272 - t1162;
t199 = qJD(1) * t272 - t1157;
t796 = -qJD(1) * t886 + t1153;
t1040 = t1194 * t1401 + t1288;
t1039 = t636 + t1490;
t1035 = t1200 * t221 + t1412 * t222;
t1034 = t1201 * t657 + t1412 * t660;
t1033 = -t1320 / 0.2e1 + t1108;
t783 = qJD(1) * t1294 * t964 + qJD(3) * t886;
t897 = t935 * t959;
t1028 = qJD(1) * t897 + 0.2e1 * t1112;
t1027 = -qJD(3) * t935 + t1158 * t1479;
t1026 = (t1406 * t860 + t1410 * t862) * pkin(4);
t969 = t726 * t1424 + t780 * t1425 + (t1443 + t213 / 0.2e1) * t1544 - (t201 + t214) * t428 / 0.2e1;
t3 = t969 + t1056;
t984 = t1403 * t849;
t33 = t984 + t1053;
t94 = t844 * t849;
t1024 = t3 * qJD(1) + t33 * qJD(2) + t94 * qJD(3);
t38 = (t1140 + t1034) * pkin(5) - t1486;
t983 = -(-t1444 - t206 / 0.2e1) * t428 - (t202 + t205) * t1544 / 0.2e1;
t4 = (t1414 * t726 + t1422 * t844 + t1035) * pkin(5) + t983;
t95 = t1398 * t844;
t1023 = -t4 * qJD(1) - t38 * qJD(2) + t95 * qJD(3);
t1022 = t1049 * t968;
t190 = t943 / 0.2e1 - t1132 * t861 + t1131 * t858 + t1026;
t320 = t1400 * t955;
t977 = t1426 * t1516 + t433 * t1427 + t787 * t1526;
t36 = (t1136 * t955 + t1409 * t872 + t1055) * pkin(4) + t977;
t1016 = -t36 * qJD(1) - t190 * qJD(2) + t320 * qJD(3);
t240 = t651 + t1547;
t76 = -t1329 / 0.2e1 + t973 + t1032;
t1015 = qJD(1) * t76 - qJD(3) * t240 + qJD(4) * t883;
t989 = t1008 / 0.2e1;
t74 = t989 - t1321 / 0.2e1 + t1040;
t811 = pkin(4) * t1102 - t1063;
t1012 = -t74 * qJD(1) + t811 * qJD(4) - t1546;
t991 = (-t367 / 0.2e1 - t349 / 0.2e1) * t962 + t1108;
t67 = t1333 / 0.2e1 + t328 / 0.2e1 + t991;
t1011 = qJD(1) * t67 - qJD(4) * t874 - t1546;
t248 = t1106 + t1543 / 0.2e1 + t1522;
t77 = -t1328 / 0.2e1 + t989 + t991;
t1010 = qJD(1) * t77 - qJD(3) * t248 + qJD(4) * t884;
t1007 = (t1200 * t213 + t1412 * t214) * pkin(5);
t1006 = (t1200 * t428 + t1412 * t1544) * pkin(5);
t1004 = (-t1360 / 0.2e1 - t1173 / 0.2e1) * pkin(5);
t994 = -t1402 * t428 + t1425 * t1473 + t1429 * t726;
t83 = -t994 + t1030;
t1000 = qJD(1) * t83 - t1250 * t1472 + t1230;
t996 = -t1403 * t1544 + t1478;
t84 = -t996 + t1031;
t999 = qJD(1) * t84 + qJD(2) * t460 + t1250 * t748;
t993 = t1030 - t1448;
t976 = -t1365 / 0.2e1 + t1125 * t1418 + t1361 / 0.2e1 + t1473 * t1416;
t111 = t1004 + t976;
t981 = t1416 * t202 + t1418 * t201 + t1419 * t206 + t1420 * t205;
t12 = t1007 + t981;
t978 = -(t1416 - t1420) * t428 + (-t1418 + t1419) * t1544;
t42 = t1006 + t978;
t533 = t874 * t883 + t875 * t884;
t986 = -t12 * qJD(1) - t111 * qJD(2) - t42 * qJD(3) + t533 * qJD(4);
t953 = t960 * t968;
t950 = t1246 / 0.2e1;
t949 = -t1206 / 0.2e1;
t948 = t1206 / 0.2e1;
t939 = t964 * t1246;
t894 = t940 - t1205 / 0.2e1;
t877 = t884 * qJD(5);
t876 = t883 * qJD(5);
t871 = t899 * t1246;
t854 = t875 * qJD(6);
t853 = t874 * qJD(6);
t852 = t1128 * t965 + t940;
t794 = t940 + (-qJD(6) / 0.2e1 + t1128) * t965;
t777 = (t1102 + t1308) * pkin(4) + t1063;
t776 = -t810 / 0.2e1 - t1101 / 0.2e1 + (-t1165 - t1307 / 0.2e1) * pkin(4);
t769 = t862 * t895;
t736 = t1140 + t1022;
t735 = t1141 + t1022;
t721 = t1472 * t1246;
t686 = 0.2e1 * t1092 - t1117;
t673 = t906 + t904 / 0.2e1 + t1037;
t672 = t1189 - t901 / 0.2e1 - t1038;
t606 = t639 * qJD(6);
t576 = 0.2e1 * t1104 + t1309;
t555 = t732 * qJD(3) - t1159;
t554 = -t1256 * t858 + t1214;
t515 = -t1123 - t696;
t514 = -t1204 * t895 - t697;
t494 = t660 * t748;
t471 = -t1214 + t1475;
t435 = qJD(3) * t731;
t370 = t736 * qJD(3) + t1159 - t1499;
t369 = t1093 * t858 - t698;
t314 = qJD(6) * t745 + t1518;
t312 = qJD(3) * t735 - t1499;
t237 = t1225 - t1492;
t236 = t606 + t1530;
t235 = t966 * t1195 + t1305 + t1441 - t1297 / 0.2e1;
t233 = t970 + t1098;
t218 = -t315 + t1505;
t217 = t576 * qJD(6) - t1476 - t1518;
t212 = t899 * t1120 + t861 * t1199 + t1393 * t1411 + t258;
t211 = pkin(4) * t1302 / 0.2e1 + t895 * t1120 + t966 * t1196 + t259;
t195 = t606 + t1520;
t191 = t1510 * t1421 - t943 / 0.2e1 + t1026 + t861 * t1426;
t186 = qJD(3) * t260;
t181 = -t1225 + t1487;
t141 = t1003 + t974;
t139 = t1093 * t1125 - t1466 - t303;
t138 = t1181 + t1531;
t113 = t1004 - t976;
t112 = t1005 + t975;
t109 = t1181 + t1519;
t86 = t994 + t1030;
t85 = t996 + t1031;
t79 = t1328 / 0.2e1 + t1033 + t1286;
t78 = t1329 / 0.2e1 + t1032 + t1287;
t75 = t1040 + t1286;
t73 = t1061 + t1287;
t69 = t1330 / 0.2e1 + t1032 + t1289;
t68 = -t1333 / 0.2e1 + t1033 + t1288;
t55 = t1371 / 0.2e1 + t1473 * t1197 + t1472 * t1198 + t1107 + (-t1392 / 0.2e1 - t380 / 0.2e1) * t962 + t1478;
t54 = -t1373 / 0.2e1 + t1125 * t1197 + t748 * t1198 + pkin(5) * t1103 + t993;
t51 = t1142 + t995 + t1031 + t1478;
t50 = t1143 + t993 + t997;
t43 = t1006 - t978;
t41 = t1516 * t1414 + t1515 * t1415 + t1423 * t1510 + t1525 * t858 + t1025;
t39 = (t1034 + t1141) * pkin(5) + t1486;
t37 = pkin(4) * t1055 + t1120 * t955 + t1199 * t872 - t977;
t32 = t984 - t1053;
t15 = qJD(3) * t82 + qJD(4) * t87;
t13 = t1007 - t981;
t9 = t1002 - t1512;
t7 = t1050 - t1513;
t5 = pkin(5) * t1035 + t1197 * t726 + t1198 * t844 - t983;
t2 = t969 - t1056;
t1 = qJD(3) * t18 + qJD(4) * t20 + qJD(5) * t19;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1147, t936 * qJD(3), 0, -t1147, 0, 0, t945 * t1246, t945 * t1244, 0, 0, t1147 * t960 - t1111, -t897 * qJD(4) - 0.2e1 * t1112 * t968, -t900 * qJD(3) + t1152 * t965, t1147 * t958 + t1111, t898 * qJD(3) + t1151 * t965, -t1147, -qJD(3) * t472 - qJD(4) * t603, qJD(3) * t473 + qJD(4) * t602, -qJD(3) * t280, qJD(3) * t311 (qJD(3) * t862 - t1475) * t861, t504 * qJD(3) + t1204 * t563, -t662 * qJD(3) + t1121 * t858 (qJD(3) * t860 + t1499) * t858, -t661 * qJD(3) + t1121 * t861, -t1147, -qJD(3) * t182 - qJD(4) * t276 - qJD(5) * t310, qJD(3) * t183 + qJD(4) * t277 + qJD(5) * t309, qJD(3) * t93 + qJD(4) * t92, qJD(3) * t102 + qJD(4) * t101 (qJD(3) * t660 - t1487) * t1473, qJD(3) * t255 + t1484 * t1517, -t402 * qJD(3) + t1487 * t968 (qJD(3) * t657 + t1039) * t1125, -t401 * qJD(3) + t1039 * t968, -t1147, -qJD(3) * t44 - qJD(4) * t99 - qJD(5) * t97 - qJD(6) * t135, qJD(3) * t45 + qJD(4) * t100 + qJD(5) * t98 + qJD(6) * t134, qJD(3) * t25 + qJD(4) * t24 + qJD(5) * t23, qJD(3) * t30 + qJD(4) * t29 + qJD(5) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t940, t1207, t1244, -t940, -t1246, 0, -t1244 * t944 + t1150, t1246 * t944 + t1149, 0, 0, -t1210 + (t1206 * t960 + t1153) * t968 (t953 - t1322) * qJD(3) + (-qJD(4) - t1256) * t1298 * t1479, t939 - t1208, t1210 + (t1206 * t958 - t1153) * t968, t1146 + t1209, -t894, -t1222 + (t1100 * t964 - t1188) * qJD(3) + t673 * qJD(4), t1221 + (t1100 * t967 + t1190) * qJD(3) + t672 * qJD(4), qJD(3) * t1085 - t1233 (-pkin(3) * t915 + pkin(8) * t1085) * qJD(3) + t1072 (t1249 + t1257) * t862 + t1511, t1220 + (-t769 - t1341) * qJD(3) + t1285, t871 - t1259 - t1276 (qJD(3) * t895 + t1258) * t860 - t1511, t1204 * t736 - t1154 - t1216, -t852 (t1510 * t965 + t860 * t955 + t873 * t895) * qJD(3) + t211 * qJD(4) + t259 * qJD(5) - t1073, t1238 + (t787 * t965 + t862 * t955 + t873 * t899) * qJD(3) + t212 * qJD(4) + t258 * qJD(5), t1326 + (-t1510 * t862 - t496 * t899 - t497 * t895 + t787 * t860) * qJD(3) + t41 * qJD(4) (t1510 * t496 - t497 * t787 + t873 * t955) * qJD(3) + t37 * qJD(4) + t1083 (t1251 + t1260) * t660 + t1184, t1234 + (-t494 - t1363) * qJD(3) + t1186, t1182 + t721 - t1264 (t1252 + t1261) * t657 + t1185, -t1156 - t1232 + t1447, -t794 (t1544 * t965 + t657 * t844 + t727 * t748) * qJD(3) + t50 * qJD(4) + t54 * qJD(5) + t86 * qJD(6) - t1081, t1367 + (t1472 * t727 + t428 * t965 + t660 * t844) * qJD(3) + t51 * qJD(4) + t55 * qJD(5) + t85 * qJD(6), t1374 + (-t1472 * t221 - t1544 * t660 - t222 * t748 + t428 * t657) * qJD(3) + t7 * qJD(4) + t9 * qJD(5) (t1544 * t221 - t222 * t428 + t727 * t844) * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t1090; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, -t1028, t964 * t1074, t783, t1048, t950, qJD(3) * t673 - qJD(4) * t771 - t1262, qJD(3) * t672 + qJD(4) * t770 + t1263, 0, 0, t1481, t269, t369, -t1481, t370, t950, qJD(3) * t211 + qJD(4) * t450 + qJD(5) * t233 - t1268, qJD(3) * t212 - qJD(4) * t451 + qJD(5) * t235 + t1267, t1327 + t41 * qJD(3) + (t858 * t966 - t861 * t963) * t1388, t37 * qJD(3) + (t450 * t966 + t451 * t963) * t1388 + t1084, t171, t58, t139, t169, t138, t950, qJD(3) * t50 + qJD(4) * t213 + qJD(5) * t78 + qJD(6) * t69 - t1382, qJD(3) * t51 - qJD(4) * t214 + qJD(5) * t79 + qJD(6) * t68 + t1271, t1375 + t7 * qJD(3) + (-t1332 - t1335) * qJD(4) + t112 * qJD(5), t2 * qJD(3) + (-t213 * t874 + t214 * t875) * qJD(4) + t13 * qJD(5) + t1088; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1481, t269, t369, -t1481, t370, t950, qJD(3) * t259 + qJD(4) * t233 - qJD(5) * t434 - t1265, qJD(3) * t258 + qJD(4) * t235 + qJD(5) * t433 + t1266, 0, 0, t171, t58, t139, t169, t138, t950, qJD(3) * t54 + qJD(4) * t78 + qJD(5) * t205 + qJD(6) * t73 - t1384, qJD(3) * t55 + qJD(4) * t79 - qJD(5) * t206 + qJD(6) * t75 + t1383, t1376 + t9 * qJD(3) + t112 * qJD(4) + (t1172 - t1364) * t1386, t5 * qJD(3) + t13 * qJD(4) + (t1401 * t205 + t206 * t962) * t1386 + t1089; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1480, t1069, -t1224 + (-qJD(6) + t1256) * t1125 - t1491, t1480, t1180 + t1531, t950, qJD(3) * t86 + qJD(4) * t69 + qJD(5) * t73 - qJD(6) * t202 - t1269, qJD(3) * t85 + qJD(4) * t68 + qJD(5) * t75 + qJD(6) * t201 + t1270, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t797 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t470, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1246, -t1244, 0, 0, 0, 0, 0, 0, 0, 0, -t1146 - t1152, t939 - t1151 (t953 + t1322) * qJD(3) (pkin(8) * t1126 - t1394) * qJD(3) - t1071, 0, 0, 0, 0, 0, 0, t1204 * t735 + t1154, t871 - t1275 (-t769 + t1341) * qJD(3) (-t1510 * t860 - t787 * t862 + t955 * t965) * qJD(3) + t191 * qJD(4) + t1080, 0, 0, 0, 0, 0, 0, t1156 + t1446, t721 + t1183 (-t494 + t1363) * qJD(3) + t1187 (-t1544 * t657 - t428 * t660 + t844 * t965) * qJD(3) + t32 * qJD(4) + t39 * qJD(5) + t1082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1205 * t967 - t1244 * t964, t1205 * t964 - t1244 * t967, 0, 0, 0, 0, 0, 0, 0, 0, t312, t471, 0, t1339 + t191 * qJD(3) + (-t858 * t963 - t861 * t966) * t1388, 0, 0, 0, 0, 0, 0, t109, t181, t1254, t32 * qJD(3) + (-t1361 + t1365) * qJD(4) + t113 * qJD(5) + t1452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, t471, 0, 0, 0, 0, 0, 0, 0, 0, t109, t181, t1254, t39 * qJD(3) + t113 * qJD(4) + (-t1173 - t1360) * t1386 + t1451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1180 + t1519, -t1226 + t1487, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t940, -t1207, 0, t940, 0, 0, -t1150, -t1149, 0, 0, -t940 * t960 - t1210, t1048 * t1479, -t1151 + t1208, -t940 * t958 + t1210, t1152 - t1209, t894, qJD(4) * t779 + t1222, qJD(4) * t778 - t1221, t1233, -t1072, -t1257 * t862 + t1511, -t1220 + t1285, t1259 - t1275, -t1258 * t860 - t1511, -t1204 * t732 + t1216, t852, -qJD(4) * t210 - qJD(5) * t256 + t1073, -qJD(4) * t209 - qJD(5) * t257 - t1238, -qJD(4) * t40 - t1326, -qJD(4) * t36 - t1083, -t1260 * t660 + t1184, t1186 - t1234, t1183 + t1264, -t1261 * t657 + t1185, t1232 - t1446, t794, -qJD(4) * t48 - qJD(5) * t53 - qJD(6) * t83 + t1081, -qJD(4) * t49 - qJD(5) * t52 - qJD(6) * t84 - t1367, -qJD(4) * t6 - qJD(5) * t8 - t1374, qJD(4) * t3 - qJD(5) * t4 - t1090; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1071, 0, 0, 0, 0, 0, 0, -t1204 * t731, -t1276, 0, -qJD(4) * t190 - t1080, 0, 0, 0, 0, 0, 0, -t1447, t1182, t1187, qJD(4) * t33 - qJD(5) * t38 - t1082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1148, t935 * qJD(4), 0, -t1148, 0, 0, -pkin(3) * t1243, -pkin(3) * t1242, 0, 0, -t1497, t1204 * t678, 0, t1497, 0, 0, qJD(4) * t724 + t1240 * t899, qJD(4) * t725 - t1240 * t895, 0, qJD(4) * t320, t1505 * t1472, t1517 * t1485, 0 (t1211 + t1476) * t748, 0, 0, qJD(4) * t347 - qJD(5) * t375 + t1211 * t844, qJD(4) * t348 - qJD(5) * t376 - t1212 * t844, 0, qJD(4) * t94 + qJD(5) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t796, -t1027, -t1119 * t967, -t796, t1119 * t964, t949, -pkin(8) * t1242 - t1058, pkin(8) * t1243 - t1059, 0, 0, -t1483, t278, t514, t1483, t515, t949, qJD(4) * t787 + qJD(5) * t686 - t1460, -t1459 - t1527 (t895 * t966 - t899 * t963) * t1388 - t1079 (t1510 * t963 + t787 * t966) * t1388 + t1016, t199, t80, t218, t197, t217, t949, qJD(4) * t428 + t1541 * t1554 - t1455, -qJD(4) * t1544 + t1541 * t1548 - t1456 (-t1331 - t1334) * qJD(4) + t141 * qJD(5) + t1453 (t1544 * t875 - t428 * t874) * qJD(4) + t43 * qJD(5) + t1024; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1483, t278, t514, t1483, t515, t949, qJD(4) * t686 + qJD(5) * t787 - t1461, -t1462 - t1527, 0, 0, t199, t80, t218, t197, t217, t949, qJD(5) * t428 + qJD(6) * t1554 - t1457 + t1538, -qJD(5) * t1544 + qJD(6) * t1548 - t1458 + t1537, t141 * qJD(4) + (t1167 - t1356) * t1386 + t1454, t43 * qJD(4) + (t1401 * t428 + t1544 * t962) * t1386 + t1023; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1482, t1068, -t1227 + t1505, t1482, t1204 * t576 - t1211 - t1518, t949, qJD(5) * t1554 + qJD(6) * t428 - t1000 + t1538, qJD(5) * t1548 - qJD(6) * t1544 + t1537 - t999, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, t1028 (-t1206 * t964 + t1245) * t968, -t783 (-t1158 - t1247) * t968, t950, -qJD(3) * t779 + t1262, -qJD(3) * t778 - t1263, 0, 0, -t1481, -t269, t554, t1481, t555, t950, qJD(3) * t210 + qJD(5) * t232 + t1268, qJD(3) * t209 + qJD(5) * t234 - t1267, qJD(3) * t40 - t1327, qJD(3) * t36 - t1084, -t171, -t58, t237, -t169, t236, t950, qJD(3) * t48 - qJD(5) * t76 - qJD(6) * t66 + t1382, qJD(3) * t49 - qJD(5) * t77 - qJD(6) * t67 - t1271, qJD(3) * t6 - qJD(5) * t110 - t1375, -qJD(3) * t3 - qJD(5) * t12 - t1088; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t698, 0, qJD(3) * t190 - t1339, 0, 0, 0, 0, 0, 0, t195, t303, -t1254, -qJD(3) * t33 - qJD(5) * t111 - t1452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t796, t1027, t967 * t1256, t796, -t964 * t1256, t948, t1058, t1059, 0, 0, t1483, -t278, t697, -t1483, t696, t948, qJD(5) * t782 + t1460, t1459, t1079, -t1016, -t199, -t80, t315, -t197, t314, t948, qJD(5) * t240 - qJD(6) * t243 + t1455, qJD(5) * t248 + t1456 + t1545, -qJD(5) * t140 - t1453, -qJD(5) * t42 - t1024; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t963 * t1387, -t966 * t1387, 0, 0, 0, 0, 0, 0, 0, 0, -t876 - t854, -t877 + t853, 0, qJD(5) * t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1127 * t963 + t1065, -t1127 * t966 + t1236, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6) * t776 - t1015 - t876, qJD(6) * t777 - t1010 - t877, -t1070 (-t1401 * t883 + t884 * t962) * t1386 + t986; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, 0, qJD(5) * t776 - t1464 - t854, qJD(5) * t777 - t1011 + t853, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1481, -t269, t554, t1481, t555, t950, qJD(3) * t256 - qJD(4) * t232 + t1265, qJD(3) * t257 - qJD(4) * t234 - t1266, 0, 0, -t171, -t58, t237, -t169, t236, t950, qJD(3) * t53 + qJD(4) * t76 + qJD(6) * t72 + t1384, qJD(3) * t52 + qJD(4) * t77 + qJD(6) * t74 - t1383, qJD(3) * t8 + qJD(4) * t110 - t1376, qJD(3) * t4 + qJD(4) * t12 - t1089; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t698, 0, 0, 0, 0, 0, 0, 0, 0, t195, t303, -t1254, qJD(3) * t38 + qJD(4) * t111 - t1451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1483, -t278, t697, -t1483, t696, t948, -qJD(4) * t782 + t1461, t1462, 0, 0, -t199, -t80, t315, -t197, t314, t948, -qJD(4) * t240 - qJD(6) * t242 + t1457, -qJD(4) * t248 + t1458 + t1545, qJD(4) * t140 - t1454, qJD(4) * t42 - t1023; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1388 * t963 - t1065, t1388 * t966 - t1236, 0, 0, 0, 0, 0, 0, 0, 0, t1015 - t1477, -qJD(6) * t811 + t1010, t1070, -t986; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1477, -pkin(5) * t1129, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t438, 0, -t1541 * t810 + t1463 (-t1130 - t1129) * pkin(5) - t1012, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1480, -t1069, t1226 - t1492, -t1480, -t1282 + t1530, t950, qJD(3) * t83 + qJD(4) * t66 - qJD(5) * t72 + t1269, qJD(3) * t84 + qJD(4) * t67 - qJD(5) * t74 - t1270, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1282 + t1520, t1224, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1482, -t1068, t1227, -t1482, -t1204 * t745 + t1518, t948, qJD(4) * t243 + qJD(5) * t242 + t1000, -t1204 * t1540 + t999, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438, 0, qJD(5) * t810 + t1464, qJD(5) * t811 + t1011, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438, 0, t1386 * t962 - t1463, pkin(5) * t1130 + t1012, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t10;