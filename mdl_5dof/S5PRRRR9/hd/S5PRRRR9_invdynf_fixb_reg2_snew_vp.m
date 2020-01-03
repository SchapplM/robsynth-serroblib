% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5PRRRR9_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:22:05
% EndTime: 2019-12-05 17:22:08
% DurationCPUTime: 3.54s
% Computational Cost: add. (19265->243), mult. (37759->369), div. (0->0), fcn. (27714->12), ass. (0->208)
t1558 = sin(pkin(10));
t1560 = cos(pkin(10));
t1535 = t1558 * g(1) - t1560 * g(2);
t1561 = cos(pkin(5));
t1602 = t1561 * t1535;
t1555 = -g(3) + qJDD(1);
t1559 = sin(pkin(5));
t1603 = t1559 * t1555;
t1616 = t1602 + t1603;
t1570 = cos(qJ(3));
t1601 = t1570 * qJD(2);
t1546 = -qJD(4) + t1601;
t1541 = -qJD(5) + t1546;
t1615 = qJD(5) - t1541;
t1614 = qJD(3) ^ 2;
t1565 = sin(qJ(4));
t1569 = cos(qJ(4));
t1566 = sin(qJ(3));
t1607 = qJD(2) * t1566;
t1524 = -t1569 * qJD(3) + t1565 * t1607;
t1526 = t1565 * qJD(3) + t1569 * t1607;
t1564 = sin(qJ(5));
t1568 = cos(qJ(5));
t1502 = t1568 * t1524 + t1564 * t1526;
t1613 = t1502 ^ 2;
t1504 = -t1564 * t1524 + t1568 * t1526;
t1612 = t1504 ^ 2;
t1611 = t1524 ^ 2;
t1610 = t1526 ^ 2;
t1609 = t1541 ^ 2;
t1608 = t1546 ^ 2;
t1606 = t1504 * t1502;
t1605 = t1524 * t1546;
t1604 = t1526 * t1524;
t1600 = qJD(4) + t1546;
t1599 = qJD(5) + t1541;
t1536 = -t1560 * g(1) - t1558 * g(2);
t1567 = sin(qJ(2));
t1571 = cos(qJ(2));
t1495 = t1571 * t1536 + t1616 * t1567;
t1572 = qJD(2) ^ 2;
t1487 = -t1572 * pkin(2) + qJDD(2) * pkin(7) + t1495;
t1577 = -t1559 * t1535 + t1561 * t1555;
t1479 = t1570 * t1487 + t1566 * t1577;
t1527 = (-pkin(3) * t1570 - pkin(8) * t1566) * qJD(2);
t1462 = -t1614 * pkin(3) + qJDD(3) * pkin(8) + t1527 * t1601 + t1479;
t1591 = t1567 * t1536 - t1616 * t1571;
t1486 = -qJDD(2) * pkin(2) - t1572 * pkin(7) + t1591;
t1596 = qJD(3) * t1601;
t1597 = t1566 * qJDD(2);
t1529 = t1596 + t1597;
t1549 = qJD(3) * t1607;
t1551 = t1570 * qJDD(2);
t1530 = t1551 - 0.2e1 * t1549;
t1465 = (-t1529 - t1596) * pkin(8) - t1530 * pkin(3) + t1486;
t1437 = t1569 * t1462 + t1565 * t1465;
t1553 = t1566 ^ 2;
t1554 = t1570 ^ 2;
t1598 = t1553 + t1554;
t1595 = t1551 - qJDD(4) - t1549;
t1436 = -t1565 * t1462 + t1569 * t1465;
t1575 = -t1565 * qJDD(3) - t1569 * t1529;
t1498 = -t1524 * qJD(4) - t1575;
t1593 = -t1569 * qJDD(3) + t1565 * t1529;
t1574 = t1526 * qJD(4) + t1593;
t1594 = -t1564 * t1498 - t1568 * t1574;
t1592 = -qJDD(5) + t1595;
t1493 = -t1595 - t1604;
t1429 = (-t1498 + t1605) * pkin(9) + t1493 * pkin(4) + t1436;
t1513 = -t1546 * pkin(4) - t1526 * pkin(9);
t1432 = -t1611 * pkin(4) - t1574 * pkin(9) + t1546 * t1513 + t1437;
t1406 = t1568 * t1429 - t1564 * t1432;
t1407 = t1564 * t1429 + t1568 * t1432;
t1397 = t1568 * t1406 + t1564 * t1407;
t1398 = -t1564 * t1406 + t1568 * t1407;
t1387 = -t1565 * t1397 + t1569 * t1398;
t1512 = t1570 * t1577;
t1461 = -t1512 - qJDD(3) * pkin(3) - t1614 * pkin(8) + (qJD(2) * t1527 + t1487) * t1566;
t1438 = pkin(4) * t1574 - t1611 * pkin(9) + t1526 * t1513 + t1461;
t1385 = t1570 * t1387 + t1566 * t1438;
t1386 = t1569 * t1397 + t1565 * t1398;
t1590 = t1385 * t1567 - t1386 * t1571;
t1444 = -t1599 * t1504 + t1594;
t1573 = -t1568 * t1498 + t1564 * t1574;
t1446 = t1599 * t1502 + t1573;
t1426 = t1564 * t1444 + t1568 * t1446;
t1427 = t1568 * t1444 - t1564 * t1446;
t1405 = -t1565 * t1426 + t1569 * t1427;
t1459 = -t1612 - t1613;
t1402 = t1570 * t1405 + t1566 * t1459;
t1404 = t1569 * t1426 + t1565 * t1427;
t1589 = t1402 * t1567 - t1404 * t1571;
t1468 = -t1592 - t1606;
t1475 = -t1609 - t1613;
t1441 = t1568 * t1468 + t1564 * t1475;
t1442 = -t1564 * t1468 + t1568 * t1475;
t1425 = -t1565 * t1441 + t1569 * t1442;
t1443 = t1615 * t1504 - t1594;
t1409 = t1570 * t1425 + t1566 * t1443;
t1424 = t1569 * t1441 + t1565 * t1442;
t1588 = t1409 * t1567 - t1424 * t1571;
t1419 = -t1565 * t1436 + t1569 * t1437;
t1411 = t1570 * t1419 + t1566 * t1461;
t1418 = t1569 * t1436 + t1565 * t1437;
t1587 = t1411 * t1567 - t1418 * t1571;
t1467 = t1592 - t1606;
t1484 = -t1609 - t1612;
t1451 = t1564 * t1467 + t1568 * t1484;
t1452 = t1568 * t1467 - t1564 * t1484;
t1431 = -t1565 * t1451 + t1569 * t1452;
t1445 = -t1615 * t1502 - t1573;
t1415 = t1570 * t1431 + t1566 * t1445;
t1430 = t1569 * t1451 + t1565 * t1452;
t1586 = t1415 * t1567 - t1430 * t1571;
t1481 = -t1600 * t1526 - t1593;
t1483 = t1600 * t1524 + t1575;
t1456 = t1569 * t1481 - t1565 * t1483;
t1491 = -t1610 - t1611;
t1440 = t1570 * t1456 + t1566 * t1491;
t1455 = t1565 * t1481 + t1569 * t1483;
t1585 = t1440 * t1567 - t1455 * t1571;
t1478 = -t1566 * t1487 + t1512;
t1448 = -t1566 * t1478 + t1570 * t1479;
t1584 = t1448 * t1567 - t1486 * t1571;
t1499 = -t1608 - t1611;
t1470 = -t1565 * t1493 + t1569 * t1499;
t1480 = (qJD(4) - t1546) * t1526 + t1593;
t1450 = t1570 * t1470 + t1566 * t1480;
t1469 = t1569 * t1493 + t1565 * t1499;
t1583 = t1450 * t1567 - t1469 * t1571;
t1492 = t1595 - t1604;
t1505 = -t1608 - t1610;
t1477 = t1569 * t1492 - t1565 * t1505;
t1482 = t1498 + t1605;
t1454 = t1570 * t1477 + t1566 * t1482;
t1476 = t1565 * t1492 + t1569 * t1505;
t1582 = t1454 * t1567 - t1476 * t1571;
t1581 = t1567 * t1495 - t1571 * t1591;
t1545 = t1566 * t1572 * t1570;
t1537 = qJDD(3) + t1545;
t1544 = -t1554 * t1572 - t1614;
t1509 = -t1566 * t1537 + t1570 * t1544;
t1580 = t1509 * t1567 + t1530 * t1571;
t1538 = -qJDD(3) + t1545;
t1543 = -t1553 * t1572 - t1614;
t1510 = t1570 * t1538 - t1566 * t1543;
t1528 = 0.2e1 * t1596 + t1597;
t1579 = t1510 * t1567 - t1528 * t1571;
t1531 = t1598 * qJDD(2);
t1534 = t1598 * t1572;
t1578 = t1531 * t1567 + t1534 * t1571;
t1576 = t1571 * qJDD(2) - t1567 * t1572;
t1533 = -t1567 * qJDD(2) - t1571 * t1572;
t1518 = t1576 * t1561;
t1517 = t1533 * t1561;
t1516 = t1576 * t1559;
t1515 = t1533 * t1559;
t1508 = t1566 * t1538 + t1570 * t1543;
t1507 = t1570 * t1537 + t1566 * t1544;
t1506 = t1571 * t1531 - t1567 * t1534;
t1501 = t1578 * t1561;
t1500 = t1578 * t1559;
t1489 = t1571 * t1510 + t1567 * t1528;
t1488 = t1571 * t1509 - t1567 * t1530;
t1474 = -t1559 * t1508 + t1561 * t1579;
t1473 = -t1559 * t1507 + t1561 * t1580;
t1472 = t1561 * t1508 + t1559 * t1579;
t1471 = t1561 * t1507 + t1559 * t1580;
t1466 = t1571 * t1495 + t1567 * t1591;
t1458 = t1559 ^ 2 * t1535 + (t1581 - t1603) * t1561;
t1457 = t1561 ^ 2 * t1555 + (t1581 - t1602) * t1559;
t1453 = t1566 * t1477 - t1570 * t1482;
t1449 = t1566 * t1470 - t1570 * t1480;
t1447 = t1570 * t1478 + t1566 * t1479;
t1439 = t1566 * t1456 - t1570 * t1491;
t1435 = t1571 * t1448 + t1567 * t1486;
t1434 = t1571 * t1454 + t1567 * t1476;
t1433 = t1571 * t1450 + t1567 * t1469;
t1428 = t1571 * t1440 + t1567 * t1455;
t1423 = -t1559 * t1453 + t1561 * t1582;
t1422 = t1561 * t1453 + t1559 * t1582;
t1421 = -t1559 * t1447 + t1561 * t1584;
t1420 = t1561 * t1447 + t1559 * t1584;
t1417 = -t1559 * t1449 + t1561 * t1583;
t1416 = t1561 * t1449 + t1559 * t1583;
t1414 = t1566 * t1431 - t1570 * t1445;
t1413 = -t1559 * t1439 + t1561 * t1585;
t1412 = t1561 * t1439 + t1559 * t1585;
t1410 = t1566 * t1419 - t1570 * t1461;
t1408 = t1566 * t1425 - t1570 * t1443;
t1403 = t1571 * t1415 + t1567 * t1430;
t1401 = t1566 * t1405 - t1570 * t1459;
t1400 = t1571 * t1409 + t1567 * t1424;
t1399 = t1571 * t1411 + t1567 * t1418;
t1396 = -t1559 * t1414 + t1561 * t1586;
t1395 = t1561 * t1414 + t1559 * t1586;
t1394 = -t1559 * t1408 + t1561 * t1588;
t1393 = t1561 * t1408 + t1559 * t1588;
t1392 = -t1559 * t1410 + t1561 * t1587;
t1391 = t1561 * t1410 + t1559 * t1587;
t1390 = t1571 * t1402 + t1567 * t1404;
t1389 = -t1559 * t1401 + t1561 * t1589;
t1388 = t1561 * t1401 + t1559 * t1589;
t1384 = t1566 * t1387 - t1570 * t1438;
t1383 = t1571 * t1385 + t1567 * t1386;
t1382 = -t1559 * t1384 + t1561 * t1590;
t1381 = t1561 * t1384 + t1559 * t1590;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1558 * t1535 + t1560 * t1536, 0, 0, 0, 0, 0, 0, -t1558 * t1518 + t1560 * t1533, -t1558 * t1517 - t1560 * t1576, 0, -t1558 * t1458 + t1560 * t1466, 0, 0, 0, 0, 0, 0, -t1558 * t1473 + t1560 * t1488, -t1558 * t1474 + t1560 * t1489, -t1558 * t1501 + t1560 * t1506, -t1558 * t1421 + t1560 * t1435, 0, 0, 0, 0, 0, 0, -t1558 * t1417 + t1560 * t1433, -t1558 * t1423 + t1560 * t1434, -t1558 * t1413 + t1560 * t1428, -t1558 * t1392 + t1560 * t1399, 0, 0, 0, 0, 0, 0, -t1558 * t1394 + t1560 * t1400, -t1558 * t1396 + t1560 * t1403, -t1558 * t1389 + t1560 * t1390, -t1558 * t1382 + t1560 * t1383; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1560 * t1535 + t1558 * t1536, 0, 0, 0, 0, 0, 0, t1560 * t1518 + t1558 * t1533, t1560 * t1517 - t1558 * t1576, 0, t1560 * t1458 + t1558 * t1466, 0, 0, 0, 0, 0, 0, t1560 * t1473 + t1558 * t1488, t1560 * t1474 + t1558 * t1489, t1560 * t1501 + t1558 * t1506, t1560 * t1421 + t1558 * t1435, 0, 0, 0, 0, 0, 0, t1560 * t1417 + t1558 * t1433, t1560 * t1423 + t1558 * t1434, t1560 * t1413 + t1558 * t1428, t1560 * t1392 + t1558 * t1399, 0, 0, 0, 0, 0, 0, t1560 * t1394 + t1558 * t1400, t1560 * t1396 + t1558 * t1403, t1560 * t1389 + t1558 * t1390, t1560 * t1382 + t1558 * t1383; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1555, 0, 0, 0, 0, 0, 0, t1516, t1515, 0, t1457, 0, 0, 0, 0, 0, 0, t1471, t1472, t1500, t1420, 0, 0, 0, 0, 0, 0, t1416, t1422, t1412, t1391, 0, 0, 0, 0, 0, 0, t1393, t1395, t1388, t1381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1536, 0, 0, 0, 0, 0, 0, t1533, -t1576, 0, t1466, 0, 0, 0, 0, 0, 0, t1488, t1489, t1506, t1435, 0, 0, 0, 0, 0, 0, t1433, t1434, t1428, t1399, 0, 0, 0, 0, 0, 0, t1400, t1403, t1390, t1383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1535, 0, 0, 0, 0, 0, 0, t1518, t1517, 0, t1458, 0, 0, 0, 0, 0, 0, t1473, t1474, t1501, t1421, 0, 0, 0, 0, 0, 0, t1417, t1423, t1413, t1392, 0, 0, 0, 0, 0, 0, t1394, t1396, t1389, t1382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1555, 0, 0, 0, 0, 0, 0, t1516, t1515, 0, t1457, 0, 0, 0, 0, 0, 0, t1471, t1472, t1500, t1420, 0, 0, 0, 0, 0, 0, t1416, t1422, t1412, t1391, 0, 0, 0, 0, 0, 0, t1393, t1395, t1388, t1381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1572, -qJDD(2), 0, t1495, 0, 0, 0, 0, 0, 0, t1509, t1510, t1531, t1448, 0, 0, 0, 0, 0, 0, t1450, t1454, t1440, t1411, 0, 0, 0, 0, 0, 0, t1409, t1415, t1402, t1385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1572, 0, -t1591, 0, 0, 0, 0, 0, 0, t1530, -t1528, t1534, -t1486, 0, 0, 0, 0, 0, 0, -t1469, -t1476, -t1455, -t1418, 0, 0, 0, 0, 0, 0, -t1424, -t1430, -t1404, -t1386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1577, 0, 0, 0, 0, 0, 0, t1507, t1508, 0, t1447, 0, 0, 0, 0, 0, 0, t1449, t1453, t1439, t1410, 0, 0, 0, 0, 0, 0, t1408, t1414, t1401, t1384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1544, t1538, t1551, t1479, 0, 0, 0, 0, 0, 0, t1470, t1477, t1456, t1419, 0, 0, 0, 0, 0, 0, t1425, t1431, t1405, t1387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1543, -t1597, t1478, 0, 0, 0, 0, 0, 0, -t1480, -t1482, -t1491, -t1461, 0, 0, 0, 0, 0, 0, -t1443, -t1445, -t1459, -t1438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1530, t1528, -t1534, t1486, 0, 0, 0, 0, 0, 0, t1469, t1476, t1455, t1418, 0, 0, 0, 0, 0, 0, t1424, t1430, t1404, t1386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1499, t1492, t1481, t1437, 0, 0, 0, 0, 0, 0, t1442, t1452, t1427, t1398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1493, t1505, t1483, t1436, 0, 0, 0, 0, 0, 0, t1441, t1451, t1426, t1397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1480, t1482, t1491, t1461, 0, 0, 0, 0, 0, 0, t1443, t1445, t1459, t1438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1475, t1467, t1444, t1407; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1468, t1484, t1446, t1406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1443, t1445, t1459, t1438;];
f_new_reg = t1;