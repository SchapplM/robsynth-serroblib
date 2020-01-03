% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRPR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:21:16
% EndTime: 2019-12-31 21:21:20
% DurationCPUTime: 3.54s
% Computational Cost: add. (9555->248), mult. (20155->287), div. (0->0), fcn. (13840->8), ass. (0->173)
t1604 = qJDD(2) + qJDD(3);
t1609 = sin(qJ(3));
t1613 = cos(qJ(3));
t1614 = cos(qJ(2));
t1654 = qJD(1) * t1614;
t1610 = sin(qJ(2));
t1655 = qJD(1) * t1610;
t1572 = t1609 * t1655 - t1613 * t1654;
t1574 = (t1609 * t1614 + t1610 * t1613) * qJD(1);
t1650 = t1574 * t1572;
t1549 = -t1604 - t1650;
t1605 = qJD(2) + qJD(3);
t1603 = t1605 ^ 2;
t1657 = t1574 ^ 2;
t1633 = -t1603 - t1657;
t1530 = t1609 * t1549 + t1613 * t1633;
t1532 = t1613 * t1549 - t1609 * t1633;
t1505 = t1610 * t1530 - t1614 * t1532;
t1611 = sin(qJ(1));
t1673 = t1611 * t1505;
t1615 = cos(qJ(1));
t1672 = t1615 * t1505;
t1571 = t1572 ^ 2;
t1545 = -t1603 - t1571;
t1629 = t1604 - t1650;
t1516 = t1609 * t1545 + t1613 * t1629;
t1519 = -t1613 * t1545 + t1609 * t1629;
t1498 = t1610 * t1516 + t1614 * t1519;
t1671 = t1611 * t1498;
t1670 = t1615 * t1498;
t1504 = t1614 * t1530 + t1610 * t1532;
t1635 = qJD(2) * t1654;
t1638 = t1610 * qJDD(1);
t1581 = t1635 + t1638;
t1601 = t1614 * qJDD(1);
t1636 = qJD(2) * t1655;
t1582 = t1601 - t1636;
t1626 = t1613 * t1581 + t1609 * t1582;
t1643 = qJD(3) - t1605;
t1528 = t1643 * t1572 - t1626;
t1669 = t1609 * t1528;
t1668 = t1613 * t1528;
t1495 = t1614 * t1516 - t1610 * t1519;
t1662 = -t1571 - t1657;
t1667 = t1611 * t1662;
t1666 = t1615 * t1662;
t1607 = t1614 ^ 2;
t1617 = qJD(1) ^ 2;
t1590 = t1611 * g(1) - t1615 * g(2);
t1623 = qJDD(1) * pkin(1) + t1590;
t1624 = qJD(2) * pkin(2) - pkin(7) * t1655;
t1544 = t1582 * pkin(2) + (pkin(7) * t1607 + pkin(6)) * t1617 - t1624 * t1655 + t1623;
t1541 = -t1572 * qJD(3) + t1626;
t1651 = t1572 * t1605;
t1631 = -t1541 + t1651;
t1649 = t1605 * t1574;
t1665 = pkin(3) * t1649 + t1631 * qJ(4) - 0.2e1 * qJD(4) * t1574 - t1544;
t1608 = sin(qJ(5));
t1612 = cos(qJ(5));
t1552 = -t1612 * t1572 + t1608 * t1605;
t1660 = t1552 ^ 2;
t1554 = t1608 * t1572 + t1612 * t1605;
t1659 = t1554 ^ 2;
t1570 = qJD(5) + t1574;
t1658 = t1570 ^ 2;
t1656 = 0.2e1 * qJD(4);
t1652 = t1554 * t1552;
t1648 = t1607 * t1617;
t1591 = -t1615 * g(1) - t1611 * g(2);
t1577 = -t1617 * pkin(1) + qJDD(1) * pkin(6) + t1591;
t1646 = t1610 * t1577;
t1645 = t1610 * t1617;
t1642 = qJD(3) + t1605;
t1641 = qJD(5) - t1570;
t1640 = qJD(5) + t1570;
t1606 = t1610 ^ 2;
t1639 = t1606 + t1607;
t1634 = -t1658 - t1659;
t1562 = -t1610 * g(3) + t1614 * t1577;
t1535 = qJDD(2) * pkin(2) - t1581 * pkin(7) - t1646 + (qJD(2) * pkin(7) * qJD(1) + pkin(2) * t1645 - g(3)) * t1614;
t1538 = -pkin(2) * t1648 + t1582 * pkin(7) - qJD(2) * t1624 + t1562;
t1511 = t1613 * t1535 - t1609 * t1538;
t1630 = t1609 * t1581 - t1613 * t1582;
t1540 = t1574 * qJD(3) + t1630;
t1632 = t1612 * t1540 - t1608 * t1604;
t1563 = t1574 * pkin(4) - t1605 * pkin(8);
t1475 = -t1571 * pkin(4) - t1574 * t1563 + (pkin(3) + pkin(8)) * t1540 + t1665;
t1550 = t1572 * pkin(3) - t1574 * qJ(4);
t1492 = -t1604 * pkin(3) - t1603 * qJ(4) + t1574 * t1550 + qJDD(4) - t1511;
t1482 = -t1629 * pkin(8) + (t1541 + t1651) * pkin(4) + t1492;
t1628 = -t1608 * t1475 + t1612 * t1482;
t1512 = t1609 * t1535 + t1613 * t1538;
t1627 = -t1608 * t1540 - t1612 * t1604;
t1621 = qJDD(5) + t1541;
t1620 = t1641 * t1552 + t1627;
t1619 = -t1603 * pkin(3) + t1604 * qJ(4) - t1572 * t1550 + t1512;
t1618 = t1621 - t1652;
t1616 = qJD(2) ^ 2;
t1595 = t1614 * t1645;
t1593 = -t1616 - t1648;
t1592 = -t1606 * t1617 - t1616;
t1589 = -qJDD(2) + t1595;
t1588 = qJDD(2) + t1595;
t1587 = t1639 * t1617;
t1586 = -t1611 * qJDD(1) - t1615 * t1617;
t1585 = t1615 * qJDD(1) - t1611 * t1617;
t1584 = t1639 * qJDD(1);
t1583 = t1601 - 0.2e1 * t1636;
t1580 = 0.2e1 * t1635 + t1638;
t1576 = t1617 * pkin(6) + t1623;
t1561 = -t1614 * g(3) - t1646;
t1558 = t1614 * t1589 - t1610 * t1592;
t1557 = -t1610 * t1588 + t1614 * t1593;
t1556 = t1610 * t1589 + t1614 * t1592;
t1555 = t1614 * t1588 + t1610 * t1593;
t1537 = -t1610 * t1561 + t1614 * t1562;
t1536 = t1614 * t1561 + t1610 * t1562;
t1526 = -t1642 * t1572 + t1626;
t1524 = -t1540 + t1649;
t1523 = t1540 + t1649;
t1522 = t1643 * t1574 + t1630;
t1521 = t1642 * t1574 + t1630;
t1520 = -t1658 - t1660;
t1515 = -t1659 - t1660;
t1514 = -t1621 - t1652;
t1510 = -t1640 * t1552 - t1627;
t1509 = -t1641 * t1554 + t1632;
t1508 = t1640 * t1554 - t1632;
t1502 = t1613 * t1524 - t1669;
t1501 = -t1613 * t1522 - t1669;
t1500 = t1609 * t1524 + t1668;
t1499 = -t1609 * t1522 + t1668;
t1494 = t1612 * t1514 - t1608 * t1634;
t1493 = t1608 * t1514 + t1612 * t1634;
t1491 = t1612 * t1520 - t1608 * t1618;
t1490 = t1608 * t1520 + t1612 * t1618;
t1489 = t1605 * t1656 + t1619;
t1488 = -t1540 * pkin(3) - t1665;
t1487 = -t1609 * t1511 + t1613 * t1512;
t1486 = t1613 * t1511 + t1609 * t1512;
t1485 = t1612 * t1509 - t1608 * t1620;
t1484 = t1608 * t1509 + t1612 * t1620;
t1483 = -t1540 * pkin(4) - t1571 * pkin(8) + (t1656 + t1563) * t1605 + t1619;
t1481 = -t1610 * t1500 + t1614 * t1502;
t1480 = -t1610 * t1499 + t1614 * t1501;
t1479 = t1614 * t1500 + t1610 * t1502;
t1478 = t1614 * t1499 + t1610 * t1501;
t1477 = t1609 * t1493 + t1613 * t1510;
t1476 = -t1613 * t1493 + t1609 * t1510;
t1474 = t1609 * t1490 + t1613 * t1508;
t1473 = -t1613 * t1490 + t1609 * t1508;
t1472 = t1609 * t1484 + t1613 * t1515;
t1471 = -t1613 * t1484 + t1609 * t1515;
t1470 = t1613 * t1489 + t1609 * t1492;
t1469 = t1609 * t1489 - t1613 * t1492;
t1468 = -t1610 * t1486 + t1614 * t1487;
t1467 = t1614 * t1486 + t1610 * t1487;
t1466 = t1612 * t1475 + t1608 * t1482;
t1464 = -t1610 * t1476 + t1614 * t1477;
t1463 = t1614 * t1476 + t1610 * t1477;
t1462 = -t1610 * t1473 + t1614 * t1474;
t1461 = t1614 * t1473 + t1610 * t1474;
t1460 = -t1610 * t1471 + t1614 * t1472;
t1459 = t1614 * t1471 + t1610 * t1472;
t1458 = -t1610 * t1469 + t1614 * t1470;
t1457 = t1614 * t1469 + t1610 * t1470;
t1456 = t1612 * t1466 - t1608 * t1628;
t1455 = t1608 * t1466 + t1612 * t1628;
t1454 = t1609 * t1455 + t1613 * t1483;
t1453 = -t1613 * t1455 + t1609 * t1483;
t1452 = -t1610 * t1453 + t1614 * t1454;
t1451 = t1614 * t1453 + t1610 * t1454;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1586, -t1585, 0, -t1611 * t1590 + t1615 * t1591, 0, 0, 0, 0, 0, 0, t1615 * t1557 - t1611 * t1583, t1615 * t1558 + t1611 * t1580, t1615 * t1584 - t1611 * t1587, t1615 * t1537 - t1611 * t1576, 0, 0, 0, 0, 0, 0, t1611 * t1521 - t1670, -t1611 * t1631 - t1672, t1615 * t1481 + t1667, t1615 * t1468 - t1611 * t1544, 0, 0, 0, 0, 0, 0, t1615 * t1480 + t1667, -t1611 * t1523 + t1670, -t1611 * t1526 + t1672, t1615 * t1458 - t1611 * t1488, 0, 0, 0, 0, 0, 0, t1615 * t1462 + t1611 * t1491, t1615 * t1464 + t1611 * t1494, t1615 * t1460 + t1611 * t1485, t1615 * t1452 + t1611 * t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1585, t1586, 0, t1615 * t1590 + t1611 * t1591, 0, 0, 0, 0, 0, 0, t1611 * t1557 + t1615 * t1583, t1611 * t1558 - t1615 * t1580, t1611 * t1584 + t1615 * t1587, t1611 * t1537 + t1615 * t1576, 0, 0, 0, 0, 0, 0, -t1615 * t1521 - t1671, t1615 * t1631 - t1673, t1611 * t1481 - t1666, t1611 * t1468 + t1615 * t1544, 0, 0, 0, 0, 0, 0, t1611 * t1480 - t1666, t1615 * t1523 + t1671, t1615 * t1526 + t1673, t1611 * t1458 + t1615 * t1488, 0, 0, 0, 0, 0, 0, t1611 * t1462 - t1615 * t1491, t1611 * t1464 - t1615 * t1494, t1611 * t1460 - t1615 * t1485, t1611 * t1452 - t1615 * t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1555, t1556, 0, t1536, 0, 0, 0, 0, 0, 0, t1495, t1504, t1479, t1467, 0, 0, 0, 0, 0, 0, t1478, -t1495, -t1504, t1457, 0, 0, 0, 0, 0, 0, t1461, t1463, t1459, t1451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1617, -qJDD(1), 0, t1591, 0, 0, 0, 0, 0, 0, t1557, t1558, t1584, t1537, 0, 0, 0, 0, 0, 0, -t1498, -t1505, t1481, t1468, 0, 0, 0, 0, 0, 0, t1480, t1498, t1505, t1458, 0, 0, 0, 0, 0, 0, t1462, t1464, t1460, t1452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1617, 0, t1590, 0, 0, 0, 0, 0, 0, t1583, -t1580, t1587, t1576, 0, 0, 0, 0, 0, 0, -t1521, t1631, -t1662, t1544, 0, 0, 0, 0, 0, 0, -t1662, t1523, t1526, t1488, 0, 0, 0, 0, 0, 0, -t1491, -t1494, -t1485, -t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1555, t1556, 0, t1536, 0, 0, 0, 0, 0, 0, t1495, t1504, t1479, t1467, 0, 0, 0, 0, 0, 0, t1478, -t1495, -t1504, t1457, 0, 0, 0, 0, 0, 0, t1461, t1463, t1459, t1451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1593, t1589, t1601, t1562, 0, 0, 0, 0, 0, 0, -t1519, t1532, t1502, t1487, 0, 0, 0, 0, 0, 0, t1501, t1519, -t1532, t1470, 0, 0, 0, 0, 0, 0, t1474, t1477, t1472, t1454; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1588, t1592, -t1638, t1561, 0, 0, 0, 0, 0, 0, t1516, t1530, t1500, t1486, 0, 0, 0, 0, 0, 0, t1499, -t1516, -t1530, t1469, 0, 0, 0, 0, 0, 0, t1473, t1476, t1471, t1453; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1583, t1580, -t1587, -t1576, 0, 0, 0, 0, 0, 0, t1521, -t1631, t1662, -t1544, 0, 0, 0, 0, 0, 0, t1662, -t1523, -t1526, -t1488, 0, 0, 0, 0, 0, 0, t1491, t1494, t1485, t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1545, t1549, t1524, t1512, 0, 0, 0, 0, 0, 0, -t1522, -t1545, -t1549, t1489, 0, 0, 0, 0, 0, 0, t1508, t1510, t1515, t1483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1629, t1633, t1528, t1511, 0, 0, 0, 0, 0, 0, t1528, -t1629, -t1633, -t1492, 0, 0, 0, 0, 0, 0, -t1490, -t1493, -t1484, -t1455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1521, -t1631, t1662, -t1544, 0, 0, 0, 0, 0, 0, t1662, -t1523, -t1526, -t1488, 0, 0, 0, 0, 0, 0, t1491, t1494, t1485, t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1662, -t1523, -t1526, -t1488, 0, 0, 0, 0, 0, 0, t1491, t1494, t1485, t1456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1522, t1545, t1549, -t1489, 0, 0, 0, 0, 0, 0, -t1508, -t1510, -t1515, -t1483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1528, t1629, t1633, t1492, 0, 0, 0, 0, 0, 0, t1490, t1493, t1484, t1455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1520, t1514, t1509, t1466; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1618, t1634, t1620, t1628; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1508, t1510, t1515, t1483;];
f_new_reg = t1;