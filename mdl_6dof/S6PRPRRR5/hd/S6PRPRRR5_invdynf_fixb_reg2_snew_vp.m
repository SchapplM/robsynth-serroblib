% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6PRPRRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:19:38
% EndTime: 2019-05-05 01:19:43
% DurationCPUTime: 4.75s
% Computational Cost: add. (21409->281), mult. (41341->387), div. (0->0), fcn. (29972->12), ass. (0->223)
t1762 = sin(pkin(11));
t1764 = cos(pkin(11));
t1738 = g(1) * t1762 - g(2) * t1764;
t1760 = -g(3) + qJDD(1);
t1763 = sin(pkin(6));
t1765 = cos(pkin(6));
t1834 = t1738 * t1765 + t1760 * t1763;
t1769 = sin(qJ(5));
t1770 = sin(qJ(4));
t1773 = cos(qJ(5));
t1774 = cos(qJ(4));
t1725 = (t1769 * t1774 + t1770 * t1773) * qJD(2);
t1718 = qJD(6) + t1725;
t1833 = qJD(6) + t1718;
t1771 = sin(qJ(2));
t1775 = cos(qJ(2));
t1777 = qJD(2) ^ 2;
t1735 = qJDD(2) * t1775 - t1771 * t1777;
t1824 = qJD(2) * t1774;
t1825 = qJD(2) * t1770;
t1727 = -t1769 * t1825 + t1773 * t1824;
t1756 = qJD(4) + qJD(5);
t1768 = sin(qJ(6));
t1772 = cos(qJ(6));
t1707 = t1727 * t1768 - t1772 * t1756;
t1832 = t1707 ^ 2;
t1709 = t1727 * t1772 + t1756 * t1768;
t1831 = t1709 ^ 2;
t1830 = t1718 ^ 2;
t1829 = t1725 ^ 2;
t1828 = t1727 ^ 2;
t1827 = t1756 ^ 2;
t1826 = 2 * qJD(3);
t1823 = t1707 * t1709;
t1822 = t1725 * t1727;
t1758 = t1770 ^ 2;
t1820 = t1758 * t1777;
t1717 = -t1738 * t1763 + t1760 * t1765;
t1818 = t1763 * t1717;
t1739 = -g(1) * t1764 - g(2) * t1762;
t1698 = -t1739 * t1771 + t1775 * t1834;
t1694 = -qJDD(2) * pkin(2) - t1777 * qJ(3) + qJDD(3) - t1698;
t1690 = -qJDD(2) * pkin(8) + t1694;
t1816 = t1774 * t1690;
t1815 = t1774 * t1777;
t1814 = qJD(5) - t1756;
t1813 = qJD(5) + t1756;
t1812 = qJD(6) - t1718;
t1675 = t1770 * t1690 + t1774 * t1717;
t1805 = qJD(4) * t1824;
t1809 = t1770 * qJDD(2);
t1731 = -t1805 - t1809;
t1744 = qJD(4) * pkin(4) - pkin(9) * t1824;
t1665 = -pkin(4) * t1820 + pkin(9) * t1731 - qJD(4) * t1744 + t1675;
t1806 = qJD(4) * t1825;
t1808 = t1774 * qJDD(2);
t1732 = -t1806 + t1808;
t1778 = qJDD(4) * pkin(4) - t1732 * pkin(9) + t1816 + (-pkin(9) * qJD(2) * qJD(4) - pkin(4) * t1815 - t1717) * t1770;
t1638 = t1773 * t1665 + t1769 * t1778;
t1759 = t1774 ^ 2;
t1811 = t1758 + t1759;
t1807 = qJDD(4) + qJDD(5);
t1804 = t1770 * t1815;
t1699 = t1775 * t1739 + t1771 * t1834;
t1637 = -t1769 * t1665 + t1773 * t1778;
t1784 = t1769 * t1731 + t1773 * t1732;
t1692 = -qJD(5) * t1725 + t1784;
t1803 = -t1768 * t1692 + t1772 * t1807;
t1802 = -t1773 * t1731 + t1769 * t1732;
t1703 = pkin(5) * t1725 - pkin(10) * t1727;
t1633 = -pkin(5) * t1827 + pkin(10) * t1807 - t1725 * t1703 + t1638;
t1766 = t1777 * pkin(8);
t1781 = -t1777 * pkin(2) + qJDD(2) * qJ(3) + t1699;
t1670 = -pkin(9) * t1820 - t1731 * pkin(4) - t1766 + (t1744 * t1774 + t1826) * qJD(2) + t1781;
t1681 = t1727 * t1813 + t1802;
t1639 = (t1725 * t1756 - t1692) * pkin(10) + t1681 * pkin(5) + t1670;
t1613 = -t1633 * t1768 + t1639 * t1772;
t1614 = t1633 * t1772 + t1639 * t1768;
t1598 = -t1613 * t1768 + t1614 * t1772;
t1632 = -pkin(5) * t1807 - pkin(10) * t1827 + t1727 * t1703 - t1637;
t1595 = t1598 * t1769 - t1632 * t1773;
t1596 = t1598 * t1773 + t1632 * t1769;
t1585 = t1595 * t1774 + t1596 * t1770;
t1597 = t1613 * t1772 + t1614 * t1768;
t1801 = -t1585 * t1775 + t1597 * t1771;
t1615 = t1637 * t1773 + t1638 * t1769;
t1616 = -t1637 * t1769 + t1638 * t1773;
t1601 = t1615 * t1774 + t1616 * t1770;
t1800 = -t1601 * t1775 + t1670 * t1771;
t1657 = -t1709 * t1812 + t1803;
t1780 = -t1772 * t1692 - t1768 * t1807;
t1659 = t1707 * t1812 + t1780;
t1635 = t1657 * t1772 - t1659 * t1768;
t1671 = -t1831 - t1832;
t1623 = t1635 * t1769 - t1671 * t1773;
t1624 = t1635 * t1773 + t1671 * t1769;
t1605 = t1623 * t1774 + t1624 * t1770;
t1634 = t1657 * t1768 + t1659 * t1772;
t1799 = -t1605 * t1775 + t1634 * t1771;
t1779 = -qJD(5) * t1727 - qJDD(6) - t1802;
t1667 = -t1779 - t1823;
t1680 = -t1830 - t1832;
t1643 = -t1667 * t1768 + t1680 * t1772;
t1656 = t1709 * t1833 - t1803;
t1626 = t1643 * t1769 - t1656 * t1773;
t1627 = t1643 * t1773 + t1656 * t1769;
t1609 = t1626 * t1774 + t1627 * t1770;
t1642 = t1667 * t1772 + t1680 * t1768;
t1798 = -t1609 * t1775 + t1642 * t1771;
t1668 = t1779 - t1823;
t1687 = -t1830 - t1831;
t1645 = t1668 * t1772 - t1687 * t1768;
t1658 = -t1707 * t1833 - t1780;
t1628 = t1645 * t1769 - t1658 * t1773;
t1629 = t1645 * t1773 + t1658 * t1769;
t1611 = t1628 * t1774 + t1629 * t1770;
t1644 = t1668 * t1768 + t1687 * t1772;
t1797 = -t1611 * t1775 + t1644 * t1771;
t1682 = -t1727 * t1814 - t1802;
t1684 = t1725 * t1814 - t1784;
t1650 = t1682 * t1769 + t1684 * t1773;
t1651 = t1682 * t1773 - t1684 * t1769;
t1630 = t1650 * t1774 + t1651 * t1770;
t1693 = -t1828 - t1829;
t1796 = -t1630 * t1775 + t1693 * t1771;
t1700 = -t1827 - t1829;
t1701 = t1807 - t1822;
t1672 = t1700 * t1769 + t1701 * t1773;
t1673 = t1700 * t1773 - t1701 * t1769;
t1646 = t1672 * t1774 + t1673 * t1770;
t1795 = -t1646 * t1775 + t1681 * t1771;
t1674 = -t1770 * t1717 + t1816;
t1648 = t1674 * t1774 + t1675 * t1770;
t1691 = qJD(2) * t1826 + t1781;
t1689 = t1691 - t1766;
t1794 = -t1648 * t1775 + t1689 * t1771;
t1702 = -t1807 - t1822;
t1715 = -t1827 - t1828;
t1685 = t1702 * t1769 + t1715 * t1773;
t1686 = t1702 * t1773 - t1715 * t1769;
t1652 = t1685 * t1774 + t1686 * t1770;
t1683 = -t1725 * t1813 + t1784;
t1793 = -t1652 * t1775 + t1683 * t1771;
t1792 = t1691 * t1771 - t1694 * t1775;
t1791 = t1698 * t1775 + t1699 * t1771;
t1740 = qJDD(4) - t1804;
t1776 = qJD(4) ^ 2;
t1747 = -t1776 - t1820;
t1711 = t1740 * t1774 + t1747 * t1770;
t1730 = 0.2e1 * t1805 + t1809;
t1790 = -t1711 * t1775 + t1730 * t1771;
t1741 = -qJDD(4) - t1804;
t1748 = -t1759 * t1777 - t1776;
t1712 = t1741 * t1770 + t1748 * t1774;
t1733 = -0.2e1 * t1806 + t1808;
t1789 = -t1712 * t1775 + t1733 * t1771;
t1736 = qJDD(2) * t1771 + t1775 * t1777;
t1721 = t1736 * t1765;
t1788 = t1721 * t1764 + t1735 * t1762;
t1787 = t1721 * t1762 - t1735 * t1764;
t1722 = t1735 * t1765;
t1786 = -t1722 * t1764 + t1736 * t1762;
t1785 = -t1722 * t1762 - t1736 * t1764;
t1734 = t1811 * qJDD(2);
t1737 = t1811 * t1777;
t1783 = t1734 * t1775 - t1737 * t1771;
t1720 = t1735 * t1763;
t1719 = t1736 * t1763;
t1714 = t1741 * t1774 - t1748 * t1770;
t1713 = -t1740 * t1770 + t1747 * t1774;
t1710 = t1765 * t1717;
t1706 = -t1734 * t1771 - t1737 * t1775;
t1705 = t1783 * t1765;
t1704 = t1783 * t1763;
t1697 = t1712 * t1771 + t1733 * t1775;
t1696 = t1711 * t1771 + t1730 * t1775;
t1679 = -t1763 * t1714 + t1765 * t1789;
t1678 = -t1763 * t1713 + t1765 * t1790;
t1677 = t1765 * t1714 + t1763 * t1789;
t1676 = t1765 * t1713 + t1763 * t1790;
t1669 = -t1698 * t1771 + t1699 * t1775;
t1666 = t1691 * t1775 + t1694 * t1771;
t1661 = t1765 * t1791 - t1818;
t1660 = t1763 * t1791 + t1710;
t1655 = t1765 * t1792 - t1818;
t1654 = t1763 * t1792 + t1710;
t1653 = -t1685 * t1770 + t1686 * t1774;
t1649 = -t1674 * t1770 + t1675 * t1774;
t1647 = -t1672 * t1770 + t1673 * t1774;
t1641 = t1652 * t1771 + t1683 * t1775;
t1640 = t1648 * t1771 + t1689 * t1775;
t1636 = t1646 * t1771 + t1681 * t1775;
t1631 = -t1650 * t1770 + t1651 * t1774;
t1625 = t1630 * t1771 + t1693 * t1775;
t1622 = -t1763 * t1653 + t1765 * t1793;
t1621 = t1765 * t1653 + t1763 * t1793;
t1620 = -t1763 * t1649 + t1765 * t1794;
t1619 = t1765 * t1649 + t1763 * t1794;
t1618 = -t1763 * t1647 + t1765 * t1795;
t1617 = t1765 * t1647 + t1763 * t1795;
t1612 = -t1628 * t1770 + t1629 * t1774;
t1610 = -t1626 * t1770 + t1627 * t1774;
t1608 = -t1763 * t1631 + t1765 * t1796;
t1607 = t1765 * t1631 + t1763 * t1796;
t1606 = -t1623 * t1770 + t1624 * t1774;
t1604 = t1611 * t1771 + t1644 * t1775;
t1603 = t1609 * t1771 + t1642 * t1775;
t1602 = -t1615 * t1770 + t1616 * t1774;
t1600 = t1605 * t1771 + t1634 * t1775;
t1599 = t1601 * t1771 + t1670 * t1775;
t1594 = -t1763 * t1612 + t1765 * t1797;
t1593 = t1765 * t1612 + t1763 * t1797;
t1592 = -t1763 * t1610 + t1765 * t1798;
t1591 = t1765 * t1610 + t1763 * t1798;
t1590 = -t1763 * t1606 + t1765 * t1799;
t1589 = t1765 * t1606 + t1763 * t1799;
t1588 = -t1763 * t1602 + t1765 * t1800;
t1587 = t1765 * t1602 + t1763 * t1800;
t1586 = -t1595 * t1770 + t1596 * t1774;
t1584 = t1585 * t1771 + t1597 * t1775;
t1583 = -t1763 * t1586 + t1765 * t1801;
t1582 = t1765 * t1586 + t1763 * t1801;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1738 * t1762 + t1739 * t1764, 0, 0, 0, 0, 0, 0, t1785, t1787, 0, -t1661 * t1762 + t1669 * t1764, 0, 0, 0, 0, 0, 0, 0, -t1785, -t1787, -t1655 * t1762 + t1666 * t1764, 0, 0, 0, 0, 0, 0, -t1678 * t1762 + t1696 * t1764, -t1679 * t1762 + t1697 * t1764, -t1705 * t1762 + t1706 * t1764, -t1620 * t1762 + t1640 * t1764, 0, 0, 0, 0, 0, 0, -t1618 * t1762 + t1636 * t1764, -t1622 * t1762 + t1641 * t1764, -t1608 * t1762 + t1625 * t1764, -t1588 * t1762 + t1599 * t1764, 0, 0, 0, 0, 0, 0, -t1592 * t1762 + t1603 * t1764, -t1594 * t1762 + t1604 * t1764, -t1590 * t1762 + t1600 * t1764, -t1583 * t1762 + t1584 * t1764; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1738 * t1764 + t1739 * t1762, 0, 0, 0, 0, 0, 0, -t1786, -t1788, 0, t1661 * t1764 + t1669 * t1762, 0, 0, 0, 0, 0, 0, 0, t1786, t1788, t1655 * t1764 + t1666 * t1762, 0, 0, 0, 0, 0, 0, t1678 * t1764 + t1696 * t1762, t1679 * t1764 + t1697 * t1762, t1705 * t1764 + t1706 * t1762, t1620 * t1764 + t1640 * t1762, 0, 0, 0, 0, 0, 0, t1618 * t1764 + t1636 * t1762, t1622 * t1764 + t1641 * t1762, t1608 * t1764 + t1625 * t1762, t1588 * t1764 + t1599 * t1762, 0, 0, 0, 0, 0, 0, t1592 * t1764 + t1603 * t1762, t1594 * t1764 + t1604 * t1762, t1590 * t1764 + t1600 * t1762, t1583 * t1764 + t1584 * t1762; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1760, 0, 0, 0, 0, 0, 0, t1720, -t1719, 0, t1660, 0, 0, 0, 0, 0, 0, 0, -t1720, t1719, t1654, 0, 0, 0, 0, 0, 0, t1676, t1677, t1704, t1619, 0, 0, 0, 0, 0, 0, t1617, t1621, t1607, t1587, 0, 0, 0, 0, 0, 0, t1591, t1593, t1589, t1582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1739, 0, 0, 0, 0, 0, 0, -t1736, -t1735, 0, t1669, 0, 0, 0, 0, 0, 0, 0, t1736, t1735, t1666, 0, 0, 0, 0, 0, 0, t1696, t1697, t1706, t1640, 0, 0, 0, 0, 0, 0, t1636, t1641, t1625, t1599, 0, 0, 0, 0, 0, 0, t1603, t1604, t1600, t1584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1738, 0, 0, 0, 0, 0, 0, t1722, -t1721, 0, t1661, 0, 0, 0, 0, 0, 0, 0, -t1722, t1721, t1655, 0, 0, 0, 0, 0, 0, t1678, t1679, t1705, t1620, 0, 0, 0, 0, 0, 0, t1618, t1622, t1608, t1588, 0, 0, 0, 0, 0, 0, t1592, t1594, t1590, t1583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1760, 0, 0, 0, 0, 0, 0, t1720, -t1719, 0, t1660, 0, 0, 0, 0, 0, 0, 0, -t1720, t1719, t1654, 0, 0, 0, 0, 0, 0, t1676, t1677, t1704, t1619, 0, 0, 0, 0, 0, 0, t1617, t1621, t1607, t1587, 0, 0, 0, 0, 0, 0, t1591, t1593, t1589, t1582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1777, -qJDD(2), 0, t1699, 0, 0, 0, 0, 0, 0, 0, t1777, qJDD(2), t1691, 0, 0, 0, 0, 0, 0, t1730, t1733, -t1737, t1689, 0, 0, 0, 0, 0, 0, t1681, t1683, t1693, t1670, 0, 0, 0, 0, 0, 0, t1642, t1644, t1634, t1597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1777, 0, t1698, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t1777, -t1694, 0, 0, 0, 0, 0, 0, -t1711, -t1712, t1734, -t1648, 0, 0, 0, 0, 0, 0, -t1646, -t1652, -t1630, -t1601, 0, 0, 0, 0, 0, 0, -t1609, -t1611, -t1605, -t1585; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1717, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1717, 0, 0, 0, 0, 0, 0, t1713, t1714, 0, t1649, 0, 0, 0, 0, 0, 0, t1647, t1653, t1631, t1602, 0, 0, 0, 0, 0, 0, t1610, t1612, t1606, t1586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1717, 0, 0, 0, 0, 0, 0, t1713, t1714, 0, t1649, 0, 0, 0, 0, 0, 0, t1647, t1653, t1631, t1602, 0, 0, 0, 0, 0, 0, t1610, t1612, t1606, t1586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1777, -qJDD(2), -t1691, 0, 0, 0, 0, 0, 0, -t1730, -t1733, t1737, -t1689, 0, 0, 0, 0, 0, 0, -t1681, -t1683, -t1693, -t1670, 0, 0, 0, 0, 0, 0, -t1642, -t1644, -t1634, -t1597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1777, t1694, 0, 0, 0, 0, 0, 0, t1711, t1712, -t1734, t1648, 0, 0, 0, 0, 0, 0, t1646, t1652, t1630, t1601, 0, 0, 0, 0, 0, 0, t1609, t1611, t1605, t1585; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1747, t1741, -t1809, t1675, 0, 0, 0, 0, 0, 0, t1673, t1686, t1651, t1616, 0, 0, 0, 0, 0, 0, t1627, t1629, t1624, t1596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1740, t1748, -t1808, t1674, 0, 0, 0, 0, 0, 0, t1672, t1685, t1650, t1615, 0, 0, 0, 0, 0, 0, t1626, t1628, t1623, t1595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1730, t1733, -t1737, t1689, 0, 0, 0, 0, 0, 0, t1681, t1683, t1693, t1670, 0, 0, 0, 0, 0, 0, t1642, t1644, t1634, t1597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1700, t1702, t1682, t1638, 0, 0, 0, 0, 0, 0, t1643, t1645, t1635, t1598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1701, t1715, t1684, t1637, 0, 0, 0, 0, 0, 0, -t1656, -t1658, -t1671, -t1632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1681, t1683, t1693, t1670, 0, 0, 0, 0, 0, 0, t1642, t1644, t1634, t1597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1680, t1668, t1657, t1614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1667, t1687, t1659, t1613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1656, t1658, t1671, t1632;];
f_new_reg  = t1;
