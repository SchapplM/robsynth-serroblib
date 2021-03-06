% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6PRPRRR6
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
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6PRPRRR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:35:32
% EndTime: 2019-05-05 01:35:36
% DurationCPUTime: 4.77s
% Computational Cost: add. (22963->284), mult. (43887->384), div. (0->0), fcn. (31000->12), ass. (0->221)
t1768 = sin(pkin(11));
t1770 = cos(pkin(11));
t1743 = g(1) * t1768 - g(2) * t1770;
t1765 = -g(3) + qJDD(1);
t1769 = sin(pkin(6));
t1771 = cos(pkin(6));
t1836 = t1743 * t1771 + t1765 * t1769;
t1776 = sin(qJ(4));
t1758 = t1776 * qJD(2);
t1754 = t1758 + qJD(5);
t1749 = qJD(6) + t1754;
t1835 = qJD(6) + t1749;
t1777 = sin(qJ(2));
t1781 = cos(qJ(2));
t1782 = qJD(2) ^ 2;
t1740 = qJDD(2) * t1781 - t1777 * t1782;
t1834 = qJD(4) ^ 2;
t1775 = sin(qJ(5));
t1779 = cos(qJ(5));
t1780 = cos(qJ(4));
t1825 = qJD(2) * t1780;
t1732 = -t1779 * qJD(4) + t1775 * t1825;
t1734 = qJD(4) * t1775 + t1779 * t1825;
t1774 = sin(qJ(6));
t1778 = cos(qJ(6));
t1709 = t1778 * t1732 + t1734 * t1774;
t1833 = t1709 ^ 2;
t1711 = -t1732 * t1774 + t1734 * t1778;
t1832 = t1711 ^ 2;
t1831 = t1732 ^ 2;
t1830 = t1734 ^ 2;
t1829 = t1749 ^ 2;
t1828 = t1754 ^ 2;
t1827 = 2 * qJD(3);
t1826 = t1782 * pkin(8);
t1824 = t1709 * t1711;
t1823 = t1732 * t1734;
t1821 = t1754 * t1732;
t1721 = -t1743 * t1769 + t1765 * t1771;
t1819 = t1769 * t1721;
t1817 = qJD(5) - t1754;
t1816 = qJD(6) - t1749;
t1744 = -g(1) * t1770 - g(2) * t1768;
t1701 = -t1744 * t1777 + t1781 * t1836;
t1694 = -qJDD(2) * pkin(2) - t1782 * qJ(3) + qJDD(3) - t1701;
t1784 = -qJDD(2) * pkin(8) + t1694;
t1677 = t1780 * t1721 + t1776 * t1784;
t1735 = (pkin(4) * t1776 - pkin(9) * t1780) * qJD(2);
t1666 = -pkin(4) * t1834 + qJDD(4) * pkin(9) - t1735 * t1758 + t1677;
t1811 = qJD(4) * t1758;
t1812 = t1780 * qJDD(2);
t1737 = -t1811 + t1812;
t1702 = t1781 * t1744 + t1777 * t1836;
t1785 = -t1782 * pkin(2) + qJDD(2) * qJ(3) + t1702;
t1755 = qJD(4) * t1825;
t1756 = t1776 * qJDD(2);
t1815 = -t1756 - t1755;
t1670 = -t1826 - t1737 * pkin(9) - t1815 * pkin(4) + (t1827 + (pkin(4) * t1780 + pkin(9) * t1776) * qJD(4)) * qJD(2) + t1785;
t1639 = t1779 * t1666 + t1775 * t1670;
t1763 = t1776 ^ 2;
t1764 = t1780 ^ 2;
t1814 = t1763 + t1764;
t1810 = t1776 * t1780 * t1782;
t1809 = qJDD(5) - t1815;
t1638 = -t1775 * t1666 + t1779 * t1670;
t1676 = -t1776 * t1721 + t1780 * t1784;
t1787 = -t1775 * qJDD(4) - t1779 * t1737;
t1705 = -qJD(5) * t1732 - t1787;
t1807 = -t1779 * qJDD(4) + t1775 * t1737;
t1786 = qJD(5) * t1734 + t1807;
t1808 = -t1774 * t1705 - t1778 * t1786;
t1806 = -qJDD(6) - t1809;
t1699 = t1809 - t1823;
t1630 = (-t1705 - t1821) * pkin(10) + t1699 * pkin(5) + t1638;
t1720 = pkin(5) * t1754 - pkin(10) * t1734;
t1634 = -pkin(5) * t1831 - pkin(10) * t1786 - t1754 * t1720 + t1639;
t1608 = t1630 * t1778 - t1634 * t1774;
t1609 = t1630 * t1774 + t1634 * t1778;
t1599 = t1608 * t1778 + t1609 * t1774;
t1600 = -t1608 * t1774 + t1609 * t1778;
t1589 = -t1599 * t1775 + t1600 * t1779;
t1665 = -qJDD(4) * pkin(4) - pkin(9) * t1834 + t1735 * t1825 - t1676;
t1640 = pkin(5) * t1786 - pkin(10) * t1831 + t1734 * t1720 + t1665;
t1586 = t1589 * t1776 - t1640 * t1780;
t1588 = t1599 * t1779 + t1600 * t1775;
t1805 = -t1586 * t1781 + t1588 * t1777;
t1648 = -t1711 * t1816 + t1808;
t1783 = -t1778 * t1705 + t1774 * t1786;
t1650 = t1709 * t1816 + t1783;
t1628 = t1648 * t1774 + t1650 * t1778;
t1629 = t1648 * t1778 - t1650 * t1774;
t1607 = -t1628 * t1775 + t1629 * t1779;
t1667 = -t1832 - t1833;
t1603 = t1607 * t1776 - t1667 * t1780;
t1606 = t1628 * t1779 + t1629 * t1775;
t1804 = -t1603 * t1781 + t1606 * t1777;
t1619 = -t1638 * t1775 + t1639 * t1779;
t1610 = t1619 * t1776 - t1665 * t1780;
t1618 = t1638 * t1779 + t1639 * t1775;
t1803 = -t1610 * t1781 + t1618 * t1777;
t1672 = -t1806 - t1824;
t1682 = -t1829 - t1833;
t1645 = t1672 * t1778 + t1682 * t1774;
t1646 = -t1672 * t1774 + t1682 * t1778;
t1627 = -t1645 * t1775 + t1646 * t1779;
t1647 = t1711 * t1835 - t1808;
t1612 = t1627 * t1776 - t1647 * t1780;
t1626 = t1645 * t1779 + t1646 * t1775;
t1802 = -t1612 * t1781 + t1626 * t1777;
t1673 = t1806 - t1824;
t1692 = -t1829 - t1832;
t1653 = t1673 * t1774 + t1692 * t1778;
t1654 = t1673 * t1778 - t1692 * t1774;
t1633 = -t1653 * t1775 + t1654 * t1779;
t1649 = -t1709 * t1835 - t1783;
t1616 = t1633 * t1776 - t1649 * t1780;
t1632 = t1653 * t1779 + t1654 * t1775;
t1801 = -t1616 * t1781 + t1632 * t1777;
t1686 = -t1734 * t1817 - t1807;
t1688 = t1732 * t1817 + t1787;
t1660 = t1686 * t1779 - t1688 * t1775;
t1698 = -t1830 - t1831;
t1641 = t1660 * t1776 - t1698 * t1780;
t1659 = t1686 * t1775 + t1688 * t1779;
t1800 = -t1641 * t1781 + t1659 * t1777;
t1643 = t1676 * t1780 + t1677 * t1776;
t1693 = qJD(2) * t1827 + t1785;
t1691 = t1693 - t1826;
t1799 = -t1643 * t1781 + t1691 * t1777;
t1706 = -t1828 - t1831;
t1675 = -t1699 * t1775 + t1706 * t1779;
t1685 = (qJD(5) + t1754) * t1734 + t1807;
t1651 = t1675 * t1776 - t1685 * t1780;
t1674 = t1699 * t1779 + t1706 * t1775;
t1798 = -t1651 * t1781 + t1674 * t1777;
t1700 = -t1809 - t1823;
t1712 = -t1828 - t1830;
t1684 = t1700 * t1779 - t1712 * t1775;
t1687 = t1705 - t1821;
t1655 = t1684 * t1776 - t1687 * t1780;
t1683 = t1700 * t1775 + t1712 * t1779;
t1797 = -t1655 * t1781 + t1683 * t1777;
t1796 = t1693 * t1777 - t1694 * t1781;
t1795 = t1701 * t1781 + t1702 * t1777;
t1745 = qJDD(4) - t1810;
t1752 = -t1763 * t1782 - t1834;
t1715 = t1745 * t1780 + t1752 * t1776;
t1736 = t1756 + 0.2e1 * t1755;
t1794 = -t1715 * t1781 + t1736 * t1777;
t1746 = -qJDD(4) - t1810;
t1753 = -t1764 * t1782 - t1834;
t1716 = t1746 * t1776 + t1753 * t1780;
t1738 = -0.2e1 * t1811 + t1812;
t1793 = -t1716 * t1781 + t1738 * t1777;
t1741 = qJDD(2) * t1777 + t1781 * t1782;
t1725 = t1741 * t1771;
t1792 = t1725 * t1770 + t1740 * t1768;
t1791 = t1725 * t1768 - t1740 * t1770;
t1726 = t1740 * t1771;
t1790 = -t1726 * t1770 + t1741 * t1768;
t1789 = -t1726 * t1768 - t1741 * t1770;
t1739 = t1814 * qJDD(2);
t1742 = t1814 * t1782;
t1788 = t1739 * t1781 - t1742 * t1777;
t1724 = t1740 * t1769;
t1723 = t1741 * t1769;
t1718 = t1746 * t1780 - t1753 * t1776;
t1717 = -t1745 * t1776 + t1752 * t1780;
t1714 = t1771 * t1721;
t1713 = -t1739 * t1777 - t1742 * t1781;
t1708 = t1788 * t1771;
t1707 = t1788 * t1769;
t1696 = t1716 * t1777 + t1738 * t1781;
t1695 = t1715 * t1777 + t1736 * t1781;
t1681 = -t1769 * t1718 + t1771 * t1793;
t1680 = -t1769 * t1717 + t1771 * t1794;
t1679 = t1771 * t1718 + t1769 * t1793;
t1678 = t1771 * t1717 + t1769 * t1794;
t1671 = -t1701 * t1777 + t1702 * t1781;
t1664 = t1693 * t1781 + t1694 * t1777;
t1662 = t1771 * t1795 - t1819;
t1661 = t1769 * t1795 + t1714;
t1658 = t1771 * t1796 - t1819;
t1657 = t1769 * t1796 + t1714;
t1656 = t1684 * t1780 + t1687 * t1776;
t1652 = t1675 * t1780 + t1685 * t1776;
t1644 = -t1676 * t1776 + t1677 * t1780;
t1642 = t1660 * t1780 + t1698 * t1776;
t1637 = t1655 * t1777 + t1683 * t1781;
t1636 = t1643 * t1777 + t1691 * t1781;
t1635 = t1651 * t1777 + t1674 * t1781;
t1631 = t1641 * t1777 + t1659 * t1781;
t1625 = -t1769 * t1656 + t1771 * t1797;
t1624 = t1771 * t1656 + t1769 * t1797;
t1623 = -t1769 * t1652 + t1771 * t1798;
t1622 = t1771 * t1652 + t1769 * t1798;
t1621 = -t1769 * t1644 + t1771 * t1799;
t1620 = t1771 * t1644 + t1769 * t1799;
t1617 = t1633 * t1780 + t1649 * t1776;
t1615 = -t1769 * t1642 + t1771 * t1800;
t1614 = t1771 * t1642 + t1769 * t1800;
t1613 = t1627 * t1780 + t1647 * t1776;
t1611 = t1619 * t1780 + t1665 * t1776;
t1605 = t1616 * t1777 + t1632 * t1781;
t1604 = t1607 * t1780 + t1667 * t1776;
t1602 = t1612 * t1777 + t1626 * t1781;
t1601 = t1610 * t1777 + t1618 * t1781;
t1598 = -t1769 * t1617 + t1771 * t1801;
t1597 = t1771 * t1617 + t1769 * t1801;
t1596 = -t1769 * t1613 + t1771 * t1802;
t1595 = t1771 * t1613 + t1769 * t1802;
t1594 = t1603 * t1777 + t1606 * t1781;
t1593 = -t1769 * t1611 + t1771 * t1803;
t1592 = t1771 * t1611 + t1769 * t1803;
t1591 = -t1769 * t1604 + t1771 * t1804;
t1590 = t1771 * t1604 + t1769 * t1804;
t1587 = t1589 * t1780 + t1640 * t1776;
t1585 = t1586 * t1777 + t1588 * t1781;
t1584 = -t1769 * t1587 + t1771 * t1805;
t1583 = t1771 * t1587 + t1769 * t1805;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1743 * t1768 + t1744 * t1770, 0, 0, 0, 0, 0, 0, t1789, t1791, 0, -t1662 * t1768 + t1671 * t1770, 0, 0, 0, 0, 0, 0, 0, -t1789, -t1791, -t1658 * t1768 + t1664 * t1770, 0, 0, 0, 0, 0, 0, -t1680 * t1768 + t1695 * t1770, -t1681 * t1768 + t1696 * t1770, -t1708 * t1768 + t1713 * t1770, -t1621 * t1768 + t1636 * t1770, 0, 0, 0, 0, 0, 0, -t1623 * t1768 + t1635 * t1770, -t1625 * t1768 + t1637 * t1770, -t1615 * t1768 + t1631 * t1770, -t1593 * t1768 + t1601 * t1770, 0, 0, 0, 0, 0, 0, -t1596 * t1768 + t1602 * t1770, -t1598 * t1768 + t1605 * t1770, -t1591 * t1768 + t1594 * t1770, -t1584 * t1768 + t1585 * t1770; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1743 * t1770 + t1744 * t1768, 0, 0, 0, 0, 0, 0, -t1790, -t1792, 0, t1662 * t1770 + t1671 * t1768, 0, 0, 0, 0, 0, 0, 0, t1790, t1792, t1658 * t1770 + t1664 * t1768, 0, 0, 0, 0, 0, 0, t1680 * t1770 + t1695 * t1768, t1681 * t1770 + t1696 * t1768, t1708 * t1770 + t1713 * t1768, t1621 * t1770 + t1636 * t1768, 0, 0, 0, 0, 0, 0, t1623 * t1770 + t1635 * t1768, t1625 * t1770 + t1637 * t1768, t1615 * t1770 + t1631 * t1768, t1593 * t1770 + t1601 * t1768, 0, 0, 0, 0, 0, 0, t1596 * t1770 + t1602 * t1768, t1598 * t1770 + t1605 * t1768, t1591 * t1770 + t1594 * t1768, t1584 * t1770 + t1585 * t1768; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1765, 0, 0, 0, 0, 0, 0, t1724, -t1723, 0, t1661, 0, 0, 0, 0, 0, 0, 0, -t1724, t1723, t1657, 0, 0, 0, 0, 0, 0, t1678, t1679, t1707, t1620, 0, 0, 0, 0, 0, 0, t1622, t1624, t1614, t1592, 0, 0, 0, 0, 0, 0, t1595, t1597, t1590, t1583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1744, 0, 0, 0, 0, 0, 0, -t1741, -t1740, 0, t1671, 0, 0, 0, 0, 0, 0, 0, t1741, t1740, t1664, 0, 0, 0, 0, 0, 0, t1695, t1696, t1713, t1636, 0, 0, 0, 0, 0, 0, t1635, t1637, t1631, t1601, 0, 0, 0, 0, 0, 0, t1602, t1605, t1594, t1585; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1743, 0, 0, 0, 0, 0, 0, t1726, -t1725, 0, t1662, 0, 0, 0, 0, 0, 0, 0, -t1726, t1725, t1658, 0, 0, 0, 0, 0, 0, t1680, t1681, t1708, t1621, 0, 0, 0, 0, 0, 0, t1623, t1625, t1615, t1593, 0, 0, 0, 0, 0, 0, t1596, t1598, t1591, t1584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1765, 0, 0, 0, 0, 0, 0, t1724, -t1723, 0, t1661, 0, 0, 0, 0, 0, 0, 0, -t1724, t1723, t1657, 0, 0, 0, 0, 0, 0, t1678, t1679, t1707, t1620, 0, 0, 0, 0, 0, 0, t1622, t1624, t1614, t1592, 0, 0, 0, 0, 0, 0, t1595, t1597, t1590, t1583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1782, -qJDD(2), 0, t1702, 0, 0, 0, 0, 0, 0, 0, t1782, qJDD(2), t1693, 0, 0, 0, 0, 0, 0, t1736, t1738, -t1742, t1691, 0, 0, 0, 0, 0, 0, t1674, t1683, t1659, t1618, 0, 0, 0, 0, 0, 0, t1626, t1632, t1606, t1588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1782, 0, t1701, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t1782, -t1694, 0, 0, 0, 0, 0, 0, -t1715, -t1716, t1739, -t1643, 0, 0, 0, 0, 0, 0, -t1651, -t1655, -t1641, -t1610, 0, 0, 0, 0, 0, 0, -t1612, -t1616, -t1603, -t1586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1721, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1721, 0, 0, 0, 0, 0, 0, t1717, t1718, 0, t1644, 0, 0, 0, 0, 0, 0, t1652, t1656, t1642, t1611, 0, 0, 0, 0, 0, 0, t1613, t1617, t1604, t1587; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1721, 0, 0, 0, 0, 0, 0, t1717, t1718, 0, t1644, 0, 0, 0, 0, 0, 0, t1652, t1656, t1642, t1611, 0, 0, 0, 0, 0, 0, t1613, t1617, t1604, t1587; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1782, -qJDD(2), -t1693, 0, 0, 0, 0, 0, 0, -t1736, -t1738, t1742, -t1691, 0, 0, 0, 0, 0, 0, -t1674, -t1683, -t1659, -t1618, 0, 0, 0, 0, 0, 0, -t1626, -t1632, -t1606, -t1588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1782, t1694, 0, 0, 0, 0, 0, 0, t1715, t1716, -t1739, t1643, 0, 0, 0, 0, 0, 0, t1651, t1655, t1641, t1610, 0, 0, 0, 0, 0, 0, t1612, t1616, t1603, t1586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1752, t1746, -t1756, t1677, 0, 0, 0, 0, 0, 0, t1675, t1684, t1660, t1619, 0, 0, 0, 0, 0, 0, t1627, t1633, t1607, t1589; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1745, t1753, -t1812, t1676, 0, 0, 0, 0, 0, 0, -t1685, -t1687, -t1698, -t1665, 0, 0, 0, 0, 0, 0, -t1647, -t1649, -t1667, -t1640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1736, t1738, -t1742, t1691, 0, 0, 0, 0, 0, 0, t1674, t1683, t1659, t1618, 0, 0, 0, 0, 0, 0, t1626, t1632, t1606, t1588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1706, t1700, t1686, t1639, 0, 0, 0, 0, 0, 0, t1646, t1654, t1629, t1600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1699, t1712, t1688, t1638, 0, 0, 0, 0, 0, 0, t1645, t1653, t1628, t1599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1685, t1687, t1698, t1665, 0, 0, 0, 0, 0, 0, t1647, t1649, t1667, t1640; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1682, t1673, t1648, t1609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1672, t1692, t1650, t1608; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1647, t1649, t1667, t1640;];
f_new_reg  = t1;
