% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRRR7
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
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRRR7_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:23:24
% EndTime: 2019-12-31 22:23:29
% DurationCPUTime: 4.45s
% Computational Cost: add. (28354->259), mult. (61938->358), div. (0->0), fcn. (45995->10), ass. (0->206)
t1619 = cos(qJ(2));
t1606 = t1619 * qJDD(1);
t1614 = sin(qJ(2));
t1655 = qJD(1) * t1614;
t1640 = qJD(2) * t1655;
t1587 = t1606 - t1640;
t1610 = t1619 ^ 2;
t1622 = qJD(1) ^ 2;
t1615 = sin(qJ(1));
t1620 = cos(qJ(1));
t1595 = t1615 * g(1) - t1620 * g(2);
t1628 = qJDD(1) * pkin(1) + t1595;
t1629 = qJD(2) * pkin(2) - pkin(7) * t1655;
t1557 = t1587 * pkin(2) - t1629 * t1655 + t1628 + (t1610 * pkin(7) + pkin(6)) * t1622;
t1613 = sin(qJ(3));
t1618 = cos(qJ(3));
t1580 = (t1613 * t1619 + t1614 * t1618) * qJD(1);
t1656 = qJD(1) * qJD(2);
t1639 = t1619 * t1656;
t1643 = t1614 * qJDD(1);
t1586 = t1639 + t1643;
t1634 = t1613 * t1586 - t1618 * t1587;
t1554 = -t1580 * qJD(3) - t1634;
t1578 = (-t1613 * t1614 + t1618 * t1619) * qJD(1);
t1577 = t1578 ^ 2;
t1608 = qJD(2) + qJD(3);
t1632 = t1608 * pkin(3) - t1580 * pkin(8);
t1516 = t1554 * pkin(3) + t1577 * pkin(8) - t1580 * t1632 + t1557;
t1612 = sin(qJ(4));
t1617 = cos(qJ(4));
t1561 = -t1617 * t1578 + t1612 * t1580;
t1559 = qJD(5) + t1561;
t1665 = qJD(5) + t1559;
t1563 = t1612 * t1578 + t1617 * t1580;
t1605 = qJD(4) + t1608;
t1611 = sin(qJ(5));
t1616 = cos(qJ(5));
t1550 = t1611 * t1563 - t1616 * t1605;
t1664 = t1550 ^ 2;
t1552 = t1616 * t1563 + t1611 * t1605;
t1663 = t1552 ^ 2;
t1662 = t1559 ^ 2;
t1661 = t1561 ^ 2;
t1660 = t1563 ^ 2;
t1659 = t1580 ^ 2;
t1658 = t1605 ^ 2;
t1657 = t1608 ^ 2;
t1654 = t1552 * t1550;
t1653 = t1563 * t1561;
t1652 = t1580 * t1578;
t1651 = t1608 * t1578;
t1650 = t1610 * t1622;
t1596 = -t1620 * g(1) - t1615 * g(2);
t1582 = -t1622 * pkin(1) + qJDD(1) * pkin(6) + t1596;
t1649 = t1614 * t1582;
t1648 = t1614 * t1622;
t1647 = qJD(3) - t1608;
t1646 = qJD(4) - t1605;
t1645 = qJD(5) - t1559;
t1545 = qJDD(2) * pkin(2) - t1586 * pkin(7) - t1649 + (pkin(2) * t1648 + pkin(7) * t1656 - g(3)) * t1619;
t1573 = -t1614 * g(3) + t1619 * t1582;
t1548 = -pkin(2) * t1650 + t1587 * pkin(7) - qJD(2) * t1629 + t1573;
t1529 = t1613 * t1545 + t1618 * t1548;
t1507 = -t1577 * pkin(3) + t1554 * pkin(8) - t1608 * t1632 + t1529;
t1528 = t1618 * t1545 - t1613 * t1548;
t1630 = -t1618 * t1586 - t1613 * t1587;
t1555 = t1578 * qJD(3) - t1630;
t1642 = qJDD(2) + qJDD(3);
t1565 = t1642 + t1652;
t1623 = (-t1555 + t1651) * pkin(8) + t1565 * pkin(3) + t1528;
t1483 = t1617 * t1507 + t1612 * t1623;
t1609 = t1614 ^ 2;
t1644 = t1609 + t1610;
t1638 = qJDD(4) + t1642;
t1482 = -t1612 * t1507 + t1617 * t1623;
t1631 = -t1612 * t1554 - t1617 * t1555;
t1518 = -t1561 * qJD(4) - t1631;
t1637 = t1605 * t1561 - t1518;
t1636 = -t1611 * t1518 + t1616 * t1638;
t1635 = -t1617 * t1554 + t1612 * t1555;
t1627 = -t1563 * qJD(4) - qJDD(5) - t1635;
t1508 = (qJD(4) + t1605) * t1563 + t1635;
t1626 = -t1616 * t1518 - t1611 * t1638;
t1621 = qJD(2) ^ 2;
t1601 = t1619 * t1648;
t1600 = -t1621 - t1650;
t1599 = -t1609 * t1622 - t1621;
t1594 = -qJDD(2) + t1601;
t1593 = qJDD(2) + t1601;
t1592 = t1644 * t1622;
t1591 = -t1615 * qJDD(1) - t1620 * t1622;
t1590 = t1620 * qJDD(1) - t1615 * t1622;
t1589 = t1644 * qJDD(1);
t1588 = t1606 - 0.2e1 * t1640;
t1585 = 0.2e1 * t1639 + t1643;
t1581 = t1622 * pkin(6) + t1628;
t1572 = -t1619 * g(3) - t1649;
t1571 = -t1657 - t1659;
t1570 = t1619 * t1594 - t1614 * t1599;
t1569 = -t1614 * t1593 + t1619 * t1600;
t1568 = t1614 * t1594 + t1619 * t1599;
t1567 = t1619 * t1593 + t1614 * t1600;
t1566 = -t1642 + t1652;
t1564 = -t1657 - t1577;
t1556 = -t1577 - t1659;
t1553 = -t1658 - t1660;
t1547 = -t1614 * t1572 + t1619 * t1573;
t1546 = t1619 * t1572 + t1614 * t1573;
t1542 = t1618 * t1566 - t1613 * t1571;
t1541 = t1613 * t1566 + t1618 * t1571;
t1540 = -t1647 * t1578 + t1630;
t1539 = t1555 + t1651;
t1538 = -t1647 * t1580 - t1634;
t1537 = (qJD(3) + t1608) * t1580 + t1634;
t1536 = t1618 * t1564 - t1613 * t1565;
t1535 = t1613 * t1564 + t1618 * t1565;
t1534 = t1561 * pkin(4) - t1563 * pkin(9);
t1533 = -t1638 - t1653;
t1532 = t1638 - t1653;
t1531 = -t1658 - t1661;
t1527 = -t1660 - t1661;
t1526 = -t1662 - t1663;
t1525 = t1617 * t1533 - t1612 * t1553;
t1524 = t1612 * t1533 + t1617 * t1553;
t1523 = -t1614 * t1541 + t1619 * t1542;
t1522 = t1619 * t1541 + t1614 * t1542;
t1521 = -t1662 - t1664;
t1520 = t1618 * t1538 - t1613 * t1540;
t1519 = t1613 * t1538 + t1618 * t1540;
t1517 = -t1663 - t1664;
t1515 = -t1614 * t1535 + t1619 * t1536;
t1514 = t1619 * t1535 + t1614 * t1536;
t1513 = t1617 * t1531 - t1612 * t1532;
t1512 = t1612 * t1531 + t1617 * t1532;
t1511 = t1646 * t1561 + t1631;
t1509 = -t1646 * t1563 - t1635;
t1503 = -t1613 * t1528 + t1618 * t1529;
t1502 = t1618 * t1528 + t1613 * t1529;
t1501 = t1627 - t1654;
t1500 = -t1627 - t1654;
t1499 = -t1613 * t1524 + t1618 * t1525;
t1498 = t1618 * t1524 + t1613 * t1525;
t1497 = -t1614 * t1519 + t1619 * t1520;
t1496 = t1619 * t1519 + t1614 * t1520;
t1495 = t1645 * t1550 + t1626;
t1494 = -t1665 * t1550 - t1626;
t1493 = -t1645 * t1552 + t1636;
t1492 = t1665 * t1552 - t1636;
t1491 = -t1613 * t1512 + t1618 * t1513;
t1490 = t1618 * t1512 + t1613 * t1513;
t1489 = t1617 * t1509 - t1612 * t1511;
t1488 = t1612 * t1509 + t1617 * t1511;
t1487 = t1616 * t1501 - t1611 * t1526;
t1486 = t1611 * t1501 + t1616 * t1526;
t1485 = -t1611 * t1500 + t1616 * t1521;
t1484 = t1616 * t1500 + t1611 * t1521;
t1481 = -t1614 * t1502 + t1619 * t1503;
t1480 = t1619 * t1502 + t1614 * t1503;
t1479 = pkin(4) * t1508 + pkin(9) * t1637 - t1516;
t1478 = -t1658 * pkin(4) + pkin(9) * t1638 - t1561 * t1534 + t1483;
t1477 = -pkin(4) * t1638 - t1658 * pkin(9) + t1563 * t1534 - t1482;
t1476 = -t1614 * t1498 + t1619 * t1499;
t1475 = t1619 * t1498 + t1614 * t1499;
t1474 = t1616 * t1493 - t1611 * t1495;
t1473 = t1611 * t1493 + t1616 * t1495;
t1472 = t1617 * t1487 + t1612 * t1494;
t1471 = t1612 * t1487 - t1617 * t1494;
t1470 = t1617 * t1485 + t1612 * t1492;
t1469 = t1612 * t1485 - t1617 * t1492;
t1468 = -t1614 * t1490 + t1619 * t1491;
t1467 = t1619 * t1490 + t1614 * t1491;
t1466 = t1617 * t1474 + t1612 * t1517;
t1465 = t1612 * t1474 - t1617 * t1517;
t1464 = -t1613 * t1488 + t1618 * t1489;
t1463 = t1618 * t1488 + t1613 * t1489;
t1462 = -t1612 * t1482 + t1617 * t1483;
t1461 = t1617 * t1482 + t1612 * t1483;
t1460 = t1616 * t1478 + t1611 * t1479;
t1459 = -t1611 * t1478 + t1616 * t1479;
t1458 = -t1613 * t1471 + t1618 * t1472;
t1457 = t1618 * t1471 + t1613 * t1472;
t1456 = -t1613 * t1469 + t1618 * t1470;
t1455 = t1618 * t1469 + t1613 * t1470;
t1454 = -t1613 * t1465 + t1618 * t1466;
t1453 = t1618 * t1465 + t1613 * t1466;
t1452 = -t1614 * t1463 + t1619 * t1464;
t1451 = t1619 * t1463 + t1614 * t1464;
t1450 = -t1613 * t1461 + t1618 * t1462;
t1449 = t1618 * t1461 + t1613 * t1462;
t1448 = -t1611 * t1459 + t1616 * t1460;
t1447 = t1616 * t1459 + t1611 * t1460;
t1446 = -t1614 * t1457 + t1619 * t1458;
t1445 = t1619 * t1457 + t1614 * t1458;
t1444 = t1617 * t1448 + t1612 * t1477;
t1443 = t1612 * t1448 - t1617 * t1477;
t1442 = -t1614 * t1455 + t1619 * t1456;
t1441 = t1619 * t1455 + t1614 * t1456;
t1440 = -t1614 * t1453 + t1619 * t1454;
t1439 = t1619 * t1453 + t1614 * t1454;
t1438 = -t1614 * t1449 + t1619 * t1450;
t1437 = t1619 * t1449 + t1614 * t1450;
t1436 = -t1613 * t1443 + t1618 * t1444;
t1435 = t1618 * t1443 + t1613 * t1444;
t1434 = -t1614 * t1435 + t1619 * t1436;
t1433 = t1619 * t1435 + t1614 * t1436;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1591, -t1590, 0, -t1615 * t1595 + t1620 * t1596, 0, 0, 0, 0, 0, 0, t1620 * t1569 - t1615 * t1588, t1620 * t1570 + t1615 * t1585, t1620 * t1589 - t1615 * t1592, t1620 * t1547 - t1615 * t1581, 0, 0, 0, 0, 0, 0, t1620 * t1515 + t1615 * t1537, t1620 * t1523 + t1615 * t1539, t1620 * t1497 + t1615 * t1556, t1620 * t1481 - t1615 * t1557, 0, 0, 0, 0, 0, 0, t1620 * t1468 + t1615 * t1508, t1620 * t1476 - t1615 * t1637, t1620 * t1452 + t1615 * t1527, t1620 * t1438 - t1615 * t1516, 0, 0, 0, 0, 0, 0, t1620 * t1442 + t1615 * t1484, t1620 * t1446 + t1615 * t1486, t1620 * t1440 + t1615 * t1473, t1620 * t1434 + t1615 * t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1590, t1591, 0, t1620 * t1595 + t1615 * t1596, 0, 0, 0, 0, 0, 0, t1615 * t1569 + t1620 * t1588, t1615 * t1570 - t1620 * t1585, t1615 * t1589 + t1620 * t1592, t1615 * t1547 + t1620 * t1581, 0, 0, 0, 0, 0, 0, t1615 * t1515 - t1620 * t1537, t1615 * t1523 - t1620 * t1539, t1615 * t1497 - t1620 * t1556, t1615 * t1481 + t1620 * t1557, 0, 0, 0, 0, 0, 0, t1615 * t1468 - t1620 * t1508, t1615 * t1476 + t1620 * t1637, t1615 * t1452 - t1620 * t1527, t1615 * t1438 + t1620 * t1516, 0, 0, 0, 0, 0, 0, t1615 * t1442 - t1620 * t1484, t1615 * t1446 - t1620 * t1486, t1615 * t1440 - t1620 * t1473, t1615 * t1434 - t1620 * t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1567, t1568, 0, t1546, 0, 0, 0, 0, 0, 0, t1514, t1522, t1496, t1480, 0, 0, 0, 0, 0, 0, t1467, t1475, t1451, t1437, 0, 0, 0, 0, 0, 0, t1441, t1445, t1439, t1433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1622, -qJDD(1), 0, t1596, 0, 0, 0, 0, 0, 0, t1569, t1570, t1589, t1547, 0, 0, 0, 0, 0, 0, t1515, t1523, t1497, t1481, 0, 0, 0, 0, 0, 0, t1468, t1476, t1452, t1438, 0, 0, 0, 0, 0, 0, t1442, t1446, t1440, t1434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1622, 0, t1595, 0, 0, 0, 0, 0, 0, t1588, -t1585, t1592, t1581, 0, 0, 0, 0, 0, 0, -t1537, -t1539, -t1556, t1557, 0, 0, 0, 0, 0, 0, -t1508, t1637, -t1527, t1516, 0, 0, 0, 0, 0, 0, -t1484, -t1486, -t1473, -t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1567, t1568, 0, t1546, 0, 0, 0, 0, 0, 0, t1514, t1522, t1496, t1480, 0, 0, 0, 0, 0, 0, t1467, t1475, t1451, t1437, 0, 0, 0, 0, 0, 0, t1441, t1445, t1439, t1433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1600, t1594, t1606, t1573, 0, 0, 0, 0, 0, 0, t1536, t1542, t1520, t1503, 0, 0, 0, 0, 0, 0, t1491, t1499, t1464, t1450, 0, 0, 0, 0, 0, 0, t1456, t1458, t1454, t1436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1593, t1599, -t1643, t1572, 0, 0, 0, 0, 0, 0, t1535, t1541, t1519, t1502, 0, 0, 0, 0, 0, 0, t1490, t1498, t1463, t1449, 0, 0, 0, 0, 0, 0, t1455, t1457, t1453, t1435; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1588, t1585, -t1592, -t1581, 0, 0, 0, 0, 0, 0, t1537, t1539, t1556, -t1557, 0, 0, 0, 0, 0, 0, t1508, -t1637, t1527, -t1516, 0, 0, 0, 0, 0, 0, t1484, t1486, t1473, t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1564, t1566, t1538, t1529, 0, 0, 0, 0, 0, 0, t1513, t1525, t1489, t1462, 0, 0, 0, 0, 0, 0, t1470, t1472, t1466, t1444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1565, t1571, t1540, t1528, 0, 0, 0, 0, 0, 0, t1512, t1524, t1488, t1461, 0, 0, 0, 0, 0, 0, t1469, t1471, t1465, t1443; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1539, t1556, -t1557, 0, 0, 0, 0, 0, 0, t1508, -t1637, t1527, -t1516, 0, 0, 0, 0, 0, 0, t1484, t1486, t1473, t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1531, t1533, t1509, t1483, 0, 0, 0, 0, 0, 0, t1485, t1487, t1474, t1448; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1532, t1553, t1511, t1482, 0, 0, 0, 0, 0, 0, -t1492, -t1494, -t1517, -t1477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1508, -t1637, t1527, -t1516, 0, 0, 0, 0, 0, 0, t1484, t1486, t1473, t1447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1521, t1501, t1493, t1460; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1500, t1526, t1495, t1459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1492, t1494, t1517, t1477;];
f_new_reg = t1;
