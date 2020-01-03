% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPRR16_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:48:23
% EndTime: 2019-12-31 20:48:27
% DurationCPUTime: 3.78s
% Computational Cost: add. (15487->266), mult. (35628->342), div. (0->0), fcn. (25927->10), ass. (0->197)
t1671 = sin(pkin(5));
t1672 = cos(pkin(5));
t1675 = sin(qJ(2));
t1679 = cos(qJ(2));
t1640 = (qJD(1) * qJD(2) * t1679 + qJDD(1) * t1675) * t1671;
t1666 = t1672 * qJD(1) + qJD(2);
t1716 = t1671 * t1679;
t1709 = qJD(1) * t1716;
t1704 = t1666 * t1709;
t1689 = -t1640 - t1704;
t1681 = qJD(1) ^ 2;
t1720 = t1671 ^ 2 * t1681;
t1663 = t1675 ^ 2 * t1720;
t1729 = t1666 ^ 2;
t1625 = -t1663 - t1729;
t1653 = t1679 * t1675 * t1720;
t1665 = t1672 * qJDD(1) + qJDD(2);
t1639 = -t1653 + t1665;
t1692 = t1625 * t1679 - t1639 * t1675;
t1590 = t1671 * t1689 + t1692 * t1672;
t1609 = t1675 * t1625 + t1679 * t1639;
t1676 = sin(qJ(1));
t1680 = cos(qJ(1));
t1746 = t1676 * t1590 + t1680 * t1609;
t1745 = t1680 * t1590 - t1676 * t1609;
t1717 = t1671 * t1675;
t1708 = qJD(1) * t1717;
t1710 = qJDD(1) * t1671;
t1641 = -qJD(2) * t1708 + t1679 * t1710;
t1646 = t1666 * t1708;
t1621 = t1641 - t1646;
t1664 = t1679 ^ 2 * t1720;
t1642 = -t1664 - t1729;
t1736 = -t1665 - t1653;
t1691 = t1642 * t1675 - t1679 * t1736;
t1598 = t1671 * t1621 + t1691 * t1672;
t1614 = -t1679 * t1642 - t1675 * t1736;
t1744 = t1676 * t1598 + t1680 * t1614;
t1743 = t1680 * t1598 - t1676 * t1614;
t1685 = qJDD(4) + t1640;
t1674 = sin(qJ(4));
t1678 = cos(qJ(4));
t1628 = t1674 * t1666 + t1678 * t1709;
t1626 = qJD(5) + t1628;
t1738 = qJD(5) + t1626;
t1737 = pkin(2) * t1646 - 0.2e1 * qJD(3) * t1708;
t1596 = -t1672 * t1621 + t1691 * t1671;
t1588 = t1692 * t1671 - t1672 * t1689;
t1630 = t1678 * t1666 - t1674 * t1709;
t1655 = qJD(4) + t1708;
t1673 = sin(qJ(5));
t1677 = cos(qJ(5));
t1615 = t1673 * t1630 - t1677 * t1655;
t1735 = t1615 ^ 2;
t1617 = t1677 * t1630 + t1673 * t1655;
t1734 = t1617 ^ 2;
t1733 = t1626 ^ 2;
t1732 = t1628 ^ 2;
t1731 = t1630 ^ 2;
t1730 = t1655 ^ 2;
t1728 = 0.2e1 * qJD(3);
t1727 = g(3) * t1679;
t1726 = t1672 * g(3);
t1725 = t1617 * t1615;
t1724 = t1630 * t1628;
t1657 = t1676 * g(1) - t1680 * g(2);
t1634 = t1681 * t1671 * pkin(7) + qJDD(1) * pkin(1) + t1657;
t1723 = t1634 * t1672;
t1636 = (-pkin(2) * t1679 - qJ(3) * t1675) * t1671 * qJD(1);
t1722 = t1636 * t1675;
t1721 = t1666 * t1679;
t1713 = qJD(4) - t1655;
t1712 = qJD(4) + t1655;
t1711 = qJD(5) - t1626;
t1637 = pkin(3) * t1708 - t1666 * pkin(8);
t1568 = -pkin(3) * t1664 - t1726 - t1640 * qJ(3) + (-pkin(2) - pkin(8)) * t1641 + (-t1634 + (-qJ(3) * t1721 - t1637 * t1675) * qJD(1)) * t1671 + t1737;
t1658 = -t1680 * g(1) - t1676 * g(2);
t1635 = -t1681 * pkin(1) + pkin(7) * t1710 + t1658;
t1706 = t1675 * t1635 - t1679 * t1723;
t1687 = -t1665 * pkin(2) - t1729 * qJ(3) + qJDD(3) + t1706;
t1682 = t1640 * pkin(3) + t1736 * pkin(8) + (t1727 + (-pkin(3) * t1721 + t1722) * qJD(1)) * t1671 + t1687;
t1549 = t1678 * t1568 + t1674 * t1682;
t1608 = -g(3) * t1717 + t1679 * t1635 + t1675 * t1723;
t1548 = -t1674 * t1568 + t1678 * t1682;
t1690 = t1674 * t1641 - t1678 * t1665;
t1605 = -t1628 * qJD(4) - t1690;
t1707 = -t1673 * t1605 + t1677 * t1685;
t1705 = t1678 * t1641 + t1674 * t1665;
t1622 = -t1671 * t1634 - t1726;
t1611 = t1628 * pkin(4) - t1630 * pkin(9);
t1536 = -t1730 * pkin(4) + t1685 * pkin(9) - t1628 * t1611 + t1549;
t1684 = -t1729 * pkin(2) + t1665 * qJ(3) + t1636 * t1709 + t1608;
t1567 = -pkin(8) * t1664 + t1641 * pkin(3) + (t1728 + t1637) * t1666 + t1684;
t1591 = t1712 * t1630 + t1705;
t1543 = (t1655 * t1628 - t1605) * pkin(9) + t1591 * pkin(4) + t1567;
t1523 = -t1673 * t1536 + t1677 * t1543;
t1524 = t1677 * t1536 + t1673 * t1543;
t1512 = -t1673 * t1523 + t1677 * t1524;
t1535 = -t1685 * pkin(4) - t1730 * pkin(9) + t1630 * t1611 - t1548;
t1509 = t1674 * t1512 - t1678 * t1535;
t1511 = t1677 * t1523 + t1673 * t1524;
t1703 = -t1509 * t1679 + t1511 * t1675;
t1527 = t1678 * t1548 + t1674 * t1549;
t1702 = -t1527 * t1679 + t1567 * t1675;
t1559 = -t1711 * t1617 + t1707;
t1683 = -t1677 * t1605 - t1673 * t1685;
t1561 = t1711 * t1615 + t1683;
t1540 = t1677 * t1559 - t1673 * t1561;
t1577 = -t1734 - t1735;
t1529 = t1674 * t1540 - t1678 * t1577;
t1539 = t1673 * t1559 + t1677 * t1561;
t1701 = -t1529 * t1679 + t1539 * t1675;
t1686 = -t1630 * qJD(4) - qJDD(5) - t1705;
t1573 = -t1686 - t1725;
t1583 = -t1733 - t1735;
t1552 = -t1673 * t1573 + t1677 * t1583;
t1558 = t1738 * t1617 - t1707;
t1533 = t1674 * t1552 - t1678 * t1558;
t1551 = t1677 * t1573 + t1673 * t1583;
t1700 = -t1533 * t1679 + t1551 * t1675;
t1574 = t1686 - t1725;
t1595 = -t1733 - t1734;
t1554 = t1677 * t1574 - t1673 * t1595;
t1560 = -t1738 * t1615 - t1683;
t1537 = t1674 * t1554 - t1678 * t1560;
t1553 = t1673 * t1574 + t1677 * t1595;
t1699 = -t1537 * t1679 + t1553 * t1675;
t1592 = -t1713 * t1630 - t1705;
t1594 = t1713 * t1628 + t1690;
t1562 = t1674 * t1592 + t1678 * t1594;
t1601 = -t1731 - t1732;
t1698 = -t1562 * t1679 + t1601 * t1675;
t1602 = t1685 - t1724;
t1606 = -t1730 - t1732;
t1575 = t1678 * t1602 + t1674 * t1606;
t1697 = -t1575 * t1679 + t1591 * t1675;
t1603 = -t1724 - t1685;
t1612 = -t1730 - t1731;
t1579 = t1674 * t1603 + t1678 * t1612;
t1593 = -t1712 * t1628 - t1690;
t1696 = -t1579 * t1679 + t1593 * t1675;
t1581 = t1666 * t1728 + t1684;
t1584 = (qJD(1) * t1722 + t1727) * t1671 + t1687;
t1695 = t1581 * t1675 - t1584 * t1679;
t1607 = -g(3) * t1716 - t1706;
t1694 = t1607 * t1679 + t1608 * t1675;
t1619 = t1640 - t1704;
t1620 = t1641 + t1646;
t1693 = -t1619 * t1679 + t1620 * t1675;
t1651 = -t1676 * qJDD(1) - t1680 * t1681;
t1650 = t1680 * qJDD(1) - t1676 * t1681;
t1643 = -t1663 - t1664;
t1600 = t1675 * t1619 + t1679 * t1620;
t1586 = -t1671 * t1643 + t1693 * t1672;
t1585 = t1672 * t1643 + t1693 * t1671;
t1582 = -t1641 * pkin(2) + t1689 * qJ(3) + t1622 + t1737;
t1580 = t1678 * t1603 - t1674 * t1612;
t1578 = -t1675 * t1607 + t1679 * t1608;
t1576 = -t1674 * t1602 + t1678 * t1606;
t1572 = -t1671 * t1622 + t1694 * t1672;
t1571 = t1672 * t1622 + t1694 * t1671;
t1565 = -t1676 * t1586 + t1680 * t1600;
t1564 = t1680 * t1586 + t1676 * t1600;
t1563 = t1678 * t1592 - t1674 * t1594;
t1557 = t1675 * t1579 + t1679 * t1593;
t1556 = t1679 * t1581 + t1675 * t1584;
t1555 = t1675 * t1575 + t1679 * t1591;
t1550 = t1675 * t1562 + t1679 * t1601;
t1547 = -t1671 * t1582 + t1695 * t1672;
t1546 = t1672 * t1582 + t1695 * t1671;
t1545 = -t1671 * t1580 + t1696 * t1672;
t1544 = t1672 * t1580 + t1696 * t1671;
t1542 = -t1671 * t1576 + t1697 * t1672;
t1541 = t1672 * t1576 + t1697 * t1671;
t1538 = t1678 * t1554 + t1674 * t1560;
t1534 = t1678 * t1552 + t1674 * t1558;
t1532 = -t1671 * t1563 + t1698 * t1672;
t1531 = t1672 * t1563 + t1698 * t1671;
t1530 = t1678 * t1540 + t1674 * t1577;
t1528 = -t1674 * t1548 + t1678 * t1549;
t1526 = t1675 * t1537 + t1679 * t1553;
t1525 = t1675 * t1533 + t1679 * t1551;
t1522 = t1675 * t1527 + t1679 * t1567;
t1521 = t1675 * t1529 + t1679 * t1539;
t1520 = -t1671 * t1538 + t1699 * t1672;
t1519 = t1672 * t1538 + t1699 * t1671;
t1518 = -t1671 * t1534 + t1700 * t1672;
t1517 = t1672 * t1534 + t1700 * t1671;
t1516 = -t1671 * t1530 + t1701 * t1672;
t1515 = t1672 * t1530 + t1701 * t1671;
t1514 = -t1671 * t1528 + t1702 * t1672;
t1513 = t1672 * t1528 + t1702 * t1671;
t1510 = t1678 * t1512 + t1674 * t1535;
t1508 = t1675 * t1509 + t1679 * t1511;
t1507 = -t1671 * t1510 + t1703 * t1672;
t1506 = t1672 * t1510 + t1703 * t1671;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1651, -t1650, 0, -t1676 * t1657 + t1680 * t1658, 0, 0, 0, 0, 0, 0, -t1744, -t1746, t1565, -t1676 * t1572 + t1680 * t1578, 0, 0, 0, 0, 0, 0, t1565, t1744, t1746, -t1676 * t1547 + t1680 * t1556, 0, 0, 0, 0, 0, 0, -t1676 * t1542 + t1680 * t1555, -t1676 * t1545 + t1680 * t1557, -t1676 * t1532 + t1680 * t1550, -t1676 * t1514 + t1680 * t1522, 0, 0, 0, 0, 0, 0, -t1676 * t1518 + t1680 * t1525, -t1676 * t1520 + t1680 * t1526, -t1676 * t1516 + t1680 * t1521, -t1676 * t1507 + t1680 * t1508; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1650, t1651, 0, t1680 * t1657 + t1676 * t1658, 0, 0, 0, 0, 0, 0, t1743, t1745, t1564, t1680 * t1572 + t1676 * t1578, 0, 0, 0, 0, 0, 0, t1564, -t1743, -t1745, t1680 * t1547 + t1676 * t1556, 0, 0, 0, 0, 0, 0, t1680 * t1542 + t1676 * t1555, t1680 * t1545 + t1676 * t1557, t1680 * t1532 + t1676 * t1550, t1680 * t1514 + t1676 * t1522, 0, 0, 0, 0, 0, 0, t1680 * t1518 + t1676 * t1525, t1680 * t1520 + t1676 * t1526, t1680 * t1516 + t1676 * t1521, t1680 * t1507 + t1676 * t1508; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1596, t1588, t1585, t1571, 0, 0, 0, 0, 0, 0, t1585, -t1596, -t1588, t1546, 0, 0, 0, 0, 0, 0, t1541, t1544, t1531, t1513, 0, 0, 0, 0, 0, 0, t1517, t1519, t1515, t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1681, -qJDD(1), 0, t1658, 0, 0, 0, 0, 0, 0, -t1614, -t1609, t1600, t1578, 0, 0, 0, 0, 0, 0, t1600, t1614, t1609, t1556, 0, 0, 0, 0, 0, 0, t1555, t1557, t1550, t1522, 0, 0, 0, 0, 0, 0, t1525, t1526, t1521, t1508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1681, 0, t1657, 0, 0, 0, 0, 0, 0, t1598, t1590, t1586, t1572, 0, 0, 0, 0, 0, 0, t1586, -t1598, -t1590, t1547, 0, 0, 0, 0, 0, 0, t1542, t1545, t1532, t1514, 0, 0, 0, 0, 0, 0, t1518, t1520, t1516, t1507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1596, t1588, t1585, t1571, 0, 0, 0, 0, 0, 0, t1585, -t1596, -t1588, t1546, 0, 0, 0, 0, 0, 0, t1541, t1544, t1531, t1513, 0, 0, 0, 0, 0, 0, t1517, t1519, t1515, t1506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1642, -t1639, t1620, t1608, 0, 0, 0, 0, 0, 0, t1620, -t1642, t1639, t1581, 0, 0, 0, 0, 0, 0, t1591, t1593, t1601, t1567, 0, 0, 0, 0, 0, 0, t1551, t1553, t1539, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1736, t1625, -t1619, t1607, 0, 0, 0, 0, 0, 0, -t1619, t1736, -t1625, -t1584, 0, 0, 0, 0, 0, 0, -t1575, -t1579, -t1562, -t1527, 0, 0, 0, 0, 0, 0, -t1533, -t1537, -t1529, -t1509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1621, -t1689, t1643, t1622, 0, 0, 0, 0, 0, 0, t1643, t1621, t1689, t1582, 0, 0, 0, 0, 0, 0, t1576, t1580, t1563, t1528, 0, 0, 0, 0, 0, 0, t1534, t1538, t1530, t1510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1643, t1621, t1689, t1582, 0, 0, 0, 0, 0, 0, t1576, t1580, t1563, t1528, 0, 0, 0, 0, 0, 0, t1534, t1538, t1530, t1510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1620, t1642, -t1639, -t1581, 0, 0, 0, 0, 0, 0, -t1591, -t1593, -t1601, -t1567, 0, 0, 0, 0, 0, 0, -t1551, -t1553, -t1539, -t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1619, -t1736, t1625, t1584, 0, 0, 0, 0, 0, 0, t1575, t1579, t1562, t1527, 0, 0, 0, 0, 0, 0, t1533, t1537, t1529, t1509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1606, t1603, t1592, t1549, 0, 0, 0, 0, 0, 0, t1552, t1554, t1540, t1512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1602, t1612, t1594, t1548, 0, 0, 0, 0, 0, 0, -t1558, -t1560, -t1577, -t1535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1591, t1593, t1601, t1567, 0, 0, 0, 0, 0, 0, t1551, t1553, t1539, t1511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1583, t1574, t1559, t1524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1573, t1595, t1561, t1523; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1558, t1560, t1577, t1535;];
f_new_reg = t1;