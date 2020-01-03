% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:13:13
% EndTime: 2019-12-05 18:13:17
% DurationCPUTime: 4.45s
% Computational Cost: add. (31700->241), mult. (79098->342), div. (0->0), fcn. (61146->10), ass. (0->193)
t1526 = qJD(1) ^ 2;
t1515 = sin(pkin(9));
t1512 = t1515 ^ 2;
t1516 = cos(pkin(9));
t1513 = t1516 ^ 2;
t1544 = t1512 + t1513;
t1499 = t1544 * t1526;
t1514 = qJD(3) + qJD(4);
t1510 = qJD(5) + t1514;
t1564 = qJD(5) + t1510;
t1519 = sin(qJ(3));
t1523 = cos(qJ(3));
t1530 = t1515 * t1523 + t1516 * t1519;
t1563 = t1530 * qJDD(1);
t1555 = qJD(1) * t1516;
t1556 = qJD(1) * t1515;
t1492 = t1519 * t1556 - t1523 * t1555;
t1494 = t1530 * qJD(1);
t1518 = sin(qJ(4));
t1522 = cos(qJ(4));
t1469 = t1522 * t1492 + t1518 * t1494;
t1471 = -t1518 * t1492 + t1522 * t1494;
t1517 = sin(qJ(5));
t1521 = cos(qJ(5));
t1444 = t1521 * t1469 + t1517 * t1471;
t1562 = t1444 ^ 2;
t1446 = -t1517 * t1469 + t1521 * t1471;
t1561 = t1446 ^ 2;
t1468 = t1469 ^ 2;
t1560 = t1471 ^ 2;
t1490 = t1492 ^ 2;
t1559 = t1494 ^ 2;
t1558 = t1510 ^ 2;
t1557 = t1514 ^ 2;
t1554 = t1446 * t1444;
t1553 = t1471 * t1469;
t1552 = t1492 * qJD(3);
t1551 = t1494 * qJD(3);
t1550 = t1494 * t1492;
t1549 = t1513 * t1526;
t1548 = t1514 * t1469;
t1547 = t1516 * t1526;
t1546 = qJD(4) - t1514;
t1545 = qJD(5) - t1510;
t1520 = sin(qJ(1));
t1524 = cos(qJ(1));
t1503 = -t1524 * g(1) - t1520 * g(2);
t1495 = -t1526 * pkin(1) + qJDD(1) * qJ(2) + t1503;
t1539 = -t1516 * g(3) - 0.2e1 * qJD(2) * t1556;
t1464 = (pkin(2) * t1547 - pkin(6) * qJDD(1) - t1495) * t1515 + t1539;
t1483 = -t1515 * g(3) + 0.2e1 * qJD(2) * t1555 + t1516 * t1495;
t1509 = t1516 * qJDD(1);
t1466 = -pkin(2) * t1549 + pkin(6) * t1509 + t1483;
t1440 = t1523 * t1464 - t1519 * t1466;
t1476 = qJDD(3) - t1550;
t1481 = t1563 - t1552;
t1414 = (-t1481 - t1552) * pkin(7) + t1476 * pkin(3) + t1440;
t1441 = t1519 * t1464 + t1523 * t1466;
t1543 = t1515 * qJDD(1);
t1458 = t1523 * t1509 - t1519 * t1543;
t1479 = t1458 - t1551;
t1535 = qJD(3) * pkin(3) - t1494 * pkin(7);
t1421 = -t1490 * pkin(3) + t1479 * pkin(7) - qJD(3) * t1535 + t1441;
t1398 = t1518 * t1414 + t1522 * t1421;
t1542 = t1520 * qJDD(1);
t1541 = t1524 * qJDD(1);
t1540 = qJDD(3) + qJDD(4);
t1502 = t1520 * g(1) - t1524 * g(2);
t1538 = -qJDD(5) - t1540;
t1397 = t1522 * t1414 - t1518 * t1421;
t1536 = -t1522 * t1479 + t1518 * t1481;
t1437 = -t1471 * qJD(4) - t1536;
t1531 = -t1518 * t1479 - t1522 * t1481;
t1438 = -t1469 * qJD(4) - t1531;
t1537 = t1521 * t1437 - t1517 * t1438;
t1534 = t1514 * pkin(4) - t1471 * pkin(8);
t1533 = -qJDD(2) + t1502;
t1448 = t1540 - t1553;
t1532 = -t1517 * t1437 - t1521 * t1438;
t1475 = (pkin(2) * t1516 + pkin(1)) * qJDD(1) + (t1544 * pkin(6) + qJ(2)) * t1526 + t1533;
t1434 = t1479 * pkin(3) + t1490 * pkin(7) - t1494 * t1535 + t1475;
t1525 = qJD(3) ^ 2;
t1504 = t1515 * t1547;
t1501 = -t1524 * t1526 - t1542;
t1500 = -t1520 * t1526 + t1541;
t1498 = t1544 * qJDD(1);
t1497 = t1516 * t1499;
t1496 = t1515 * t1499;
t1489 = qJDD(1) * pkin(1) + t1526 * qJ(2) + t1533;
t1484 = -t1525 - t1559;
t1482 = -t1515 * t1495 + t1539;
t1480 = t1563 - 0.2e1 * t1552;
t1478 = -t1458 + 0.2e1 * t1551;
t1477 = -qJDD(3) - t1550;
t1474 = -t1525 - t1490;
t1462 = -t1557 - t1560;
t1459 = -t1490 - t1559;
t1457 = t1523 * t1477 - t1519 * t1484;
t1456 = t1519 * t1477 + t1523 * t1484;
t1455 = -t1515 * t1482 + t1516 * t1483;
t1454 = t1516 * t1482 + t1515 * t1483;
t1453 = t1523 * t1458 + t1519 * t1563;
t1452 = t1519 * t1458 - t1523 * t1563;
t1451 = t1523 * t1474 - t1519 * t1476;
t1450 = t1519 * t1474 + t1523 * t1476;
t1449 = -t1540 - t1553;
t1447 = -t1557 - t1468;
t1439 = -t1558 - t1561;
t1436 = -t1468 - t1560;
t1433 = -t1515 * t1456 + t1516 * t1457;
t1432 = t1516 * t1456 + t1515 * t1457;
t1431 = t1522 * t1449 - t1518 * t1462;
t1430 = t1518 * t1449 + t1522 * t1462;
t1429 = -t1515 * t1452 + t1516 * t1453;
t1428 = t1516 * t1452 + t1515 * t1453;
t1427 = t1546 * t1469 + t1531;
t1426 = t1438 - t1548;
t1425 = -t1546 * t1471 - t1536;
t1424 = (qJD(4) + t1514) * t1471 + t1536;
t1423 = -t1515 * t1450 + t1516 * t1451;
t1422 = t1516 * t1450 + t1515 * t1451;
t1420 = t1522 * t1447 - t1518 * t1448;
t1419 = t1518 * t1447 + t1522 * t1448;
t1417 = t1538 - t1554;
t1416 = -t1538 - t1554;
t1415 = -t1558 - t1562;
t1411 = -t1519 * t1440 + t1523 * t1441;
t1410 = t1523 * t1440 + t1519 * t1441;
t1409 = -t1561 - t1562;
t1408 = t1521 * t1417 - t1517 * t1439;
t1407 = t1517 * t1417 + t1521 * t1439;
t1406 = -t1519 * t1430 + t1523 * t1431;
t1405 = t1523 * t1430 + t1519 * t1431;
t1404 = t1522 * t1425 - t1518 * t1427;
t1403 = t1518 * t1425 + t1522 * t1427;
t1402 = -t1519 * t1419 + t1523 * t1420;
t1401 = t1523 * t1419 + t1519 * t1420;
t1400 = t1521 * t1415 - t1517 * t1416;
t1399 = t1517 * t1415 + t1521 * t1416;
t1396 = t1437 * pkin(4) + t1468 * pkin(8) - t1471 * t1534 + t1434;
t1395 = -t1515 * t1410 + t1516 * t1411;
t1394 = t1516 * t1410 + t1515 * t1411;
t1393 = t1545 * t1444 + t1532;
t1392 = -t1564 * t1444 - t1532;
t1391 = -t1545 * t1446 + t1537;
t1390 = t1564 * t1446 - t1537;
t1389 = -t1518 * t1407 + t1522 * t1408;
t1388 = t1522 * t1407 + t1518 * t1408;
t1387 = -t1515 * t1405 + t1516 * t1406;
t1386 = t1516 * t1405 + t1515 * t1406;
t1385 = -t1468 * pkin(4) + t1437 * pkin(8) - t1514 * t1534 + t1398;
t1384 = (-t1438 - t1548) * pkin(8) + t1448 * pkin(4) + t1397;
t1383 = -t1519 * t1403 + t1523 * t1404;
t1382 = t1523 * t1403 + t1519 * t1404;
t1381 = -t1515 * t1401 + t1516 * t1402;
t1380 = t1516 * t1401 + t1515 * t1402;
t1379 = -t1518 * t1399 + t1522 * t1400;
t1378 = t1522 * t1399 + t1518 * t1400;
t1377 = -t1518 * t1397 + t1522 * t1398;
t1376 = t1522 * t1397 + t1518 * t1398;
t1375 = t1521 * t1391 - t1517 * t1393;
t1374 = t1517 * t1391 + t1521 * t1393;
t1373 = -t1519 * t1388 + t1523 * t1389;
t1372 = t1523 * t1388 + t1519 * t1389;
t1371 = t1517 * t1384 + t1521 * t1385;
t1370 = t1521 * t1384 - t1517 * t1385;
t1369 = -t1515 * t1382 + t1516 * t1383;
t1368 = t1516 * t1382 + t1515 * t1383;
t1367 = -t1519 * t1378 + t1523 * t1379;
t1366 = t1523 * t1378 + t1519 * t1379;
t1365 = -t1519 * t1376 + t1523 * t1377;
t1364 = t1523 * t1376 + t1519 * t1377;
t1363 = -t1518 * t1374 + t1522 * t1375;
t1362 = t1522 * t1374 + t1518 * t1375;
t1361 = -t1515 * t1372 + t1516 * t1373;
t1360 = t1516 * t1372 + t1515 * t1373;
t1359 = -t1517 * t1370 + t1521 * t1371;
t1358 = t1521 * t1370 + t1517 * t1371;
t1357 = -t1515 * t1366 + t1516 * t1367;
t1356 = t1516 * t1366 + t1515 * t1367;
t1355 = -t1515 * t1364 + t1516 * t1365;
t1354 = t1516 * t1364 + t1515 * t1365;
t1353 = -t1519 * t1362 + t1523 * t1363;
t1352 = t1523 * t1362 + t1519 * t1363;
t1351 = -t1518 * t1358 + t1522 * t1359;
t1350 = t1522 * t1358 + t1518 * t1359;
t1349 = -t1515 * t1352 + t1516 * t1353;
t1348 = t1516 * t1352 + t1515 * t1353;
t1347 = -t1519 * t1350 + t1523 * t1351;
t1346 = t1523 * t1350 + t1519 * t1351;
t1345 = -t1515 * t1346 + t1516 * t1347;
t1344 = t1516 * t1346 + t1515 * t1347;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1501, -t1500, 0, -t1520 * t1502 + t1524 * t1503, 0, 0, 0, 0, 0, 0, -t1524 * t1497 - t1516 * t1542, t1524 * t1496 + t1515 * t1542, t1524 * t1498 - t1520 * t1499, t1524 * t1455 - t1520 * t1489, 0, 0, 0, 0, 0, 0, t1524 * t1423 + t1520 * t1478, t1524 * t1433 + t1520 * t1480, t1524 * t1429 + t1520 * t1459, t1524 * t1395 - t1520 * t1475, 0, 0, 0, 0, 0, 0, t1524 * t1381 + t1520 * t1424, t1524 * t1387 + t1520 * t1426, t1524 * t1369 + t1520 * t1436, t1524 * t1355 - t1520 * t1434, 0, 0, 0, 0, 0, 0, t1524 * t1357 + t1520 * t1390, t1524 * t1361 + t1520 * t1392, t1524 * t1349 + t1520 * t1409, t1524 * t1345 - t1520 * t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1500, t1501, 0, t1524 * t1502 + t1520 * t1503, 0, 0, 0, 0, 0, 0, -t1520 * t1497 + t1516 * t1541, t1520 * t1496 - t1515 * t1541, t1520 * t1498 + t1524 * t1499, t1520 * t1455 + t1524 * t1489, 0, 0, 0, 0, 0, 0, t1520 * t1423 - t1524 * t1478, t1520 * t1433 - t1524 * t1480, t1520 * t1429 - t1524 * t1459, t1520 * t1395 + t1524 * t1475, 0, 0, 0, 0, 0, 0, t1520 * t1381 - t1524 * t1424, t1520 * t1387 - t1524 * t1426, t1520 * t1369 - t1524 * t1436, t1520 * t1355 + t1524 * t1434, 0, 0, 0, 0, 0, 0, t1520 * t1357 - t1524 * t1390, t1520 * t1361 - t1524 * t1392, t1520 * t1349 - t1524 * t1409, t1520 * t1345 + t1524 * t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1454, 0, 0, 0, 0, 0, 0, t1422, t1432, t1428, t1394, 0, 0, 0, 0, 0, 0, t1380, t1386, t1368, t1354, 0, 0, 0, 0, 0, 0, t1356, t1360, t1348, t1344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1526, -qJDD(1), 0, t1503, 0, 0, 0, 0, 0, 0, -t1497, t1496, t1498, t1455, 0, 0, 0, 0, 0, 0, t1423, t1433, t1429, t1395, 0, 0, 0, 0, 0, 0, t1381, t1387, t1369, t1355, 0, 0, 0, 0, 0, 0, t1357, t1361, t1349, t1345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1526, 0, t1502, 0, 0, 0, 0, 0, 0, t1509, -t1543, t1499, t1489, 0, 0, 0, 0, 0, 0, -t1478, -t1480, -t1459, t1475, 0, 0, 0, 0, 0, 0, -t1424, -t1426, -t1436, t1434, 0, 0, 0, 0, 0, 0, -t1390, -t1392, -t1409, t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1454, 0, 0, 0, 0, 0, 0, t1422, t1432, t1428, t1394, 0, 0, 0, 0, 0, 0, t1380, t1386, t1368, t1354, 0, 0, 0, 0, 0, 0, t1356, t1360, t1348, t1344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1549, t1504, t1509, t1483, 0, 0, 0, 0, 0, 0, t1451, t1457, t1453, t1411, 0, 0, 0, 0, 0, 0, t1402, t1406, t1383, t1365, 0, 0, 0, 0, 0, 0, t1367, t1373, t1353, t1347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1504, -t1512 * t1526, -t1543, t1482, 0, 0, 0, 0, 0, 0, t1450, t1456, t1452, t1410, 0, 0, 0, 0, 0, 0, t1401, t1405, t1382, t1364, 0, 0, 0, 0, 0, 0, t1366, t1372, t1352, t1346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1509, t1543, -t1499, -t1489, 0, 0, 0, 0, 0, 0, t1478, t1480, t1459, -t1475, 0, 0, 0, 0, 0, 0, t1424, t1426, t1436, -t1434, 0, 0, 0, 0, 0, 0, t1390, t1392, t1409, -t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1474, t1477, t1458, t1441, 0, 0, 0, 0, 0, 0, t1420, t1431, t1404, t1377, 0, 0, 0, 0, 0, 0, t1379, t1389, t1363, t1351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1476, t1484, -t1563, t1440, 0, 0, 0, 0, 0, 0, t1419, t1430, t1403, t1376, 0, 0, 0, 0, 0, 0, t1378, t1388, t1362, t1350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1478, t1480, t1459, -t1475, 0, 0, 0, 0, 0, 0, t1424, t1426, t1436, -t1434, 0, 0, 0, 0, 0, 0, t1390, t1392, t1409, -t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1447, t1449, t1425, t1398, 0, 0, 0, 0, 0, 0, t1400, t1408, t1375, t1359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1448, t1462, t1427, t1397, 0, 0, 0, 0, 0, 0, t1399, t1407, t1374, t1358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1424, t1426, t1436, -t1434, 0, 0, 0, 0, 0, 0, t1390, t1392, t1409, -t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1415, t1417, t1391, t1371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1416, t1439, t1393, t1370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1390, t1392, t1409, -t1396;];
f_new_reg = t1;