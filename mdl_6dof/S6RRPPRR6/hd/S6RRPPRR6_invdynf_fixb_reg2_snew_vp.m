% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPPRR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:56:09
% EndTime: 2019-05-06 10:56:15
% DurationCPUTime: 6.01s
% Computational Cost: add. (32923->315), mult. (75979->380), div. (0->0), fcn. (51877->10), ass. (0->217)
t1924 = sin(qJ(2));
t1928 = cos(qJ(2));
t1931 = qJD(1) ^ 2;
t1949 = t1928 * t1931;
t1900 = t1924 * t1949;
t1892 = qJDD(2) - t1900;
t1917 = t1924 ^ 2;
t1930 = qJD(2) ^ 2;
t1898 = t1917 * t1931 + t1930;
t1861 = t1892 * t1928 - t1898 * t1924;
t1956 = qJD(2) * t1928;
t1942 = qJD(1) * t1956;
t1945 = t1924 * qJDD(1);
t1882 = 0.2e1 * t1942 + t1945;
t1925 = sin(qJ(1));
t1929 = cos(qJ(1));
t1971 = t1861 * t1925 + t1882 * t1929;
t1970 = t1861 * t1929 - t1882 * t1925;
t1920 = sin(pkin(10));
t1921 = cos(pkin(10));
t1958 = qJD(1) * t1928;
t1959 = qJD(1) * t1924;
t1875 = -t1920 * t1959 - t1921 * t1958;
t1876 = -t1920 * t1958 + t1921 * t1959;
t1923 = sin(qJ(5));
t1927 = cos(qJ(5));
t1845 = -t1927 * t1875 + t1876 * t1923;
t1844 = qJD(6) + t1845;
t1969 = qJD(6) + t1844;
t1883 = t1942 + t1945;
t1907 = qJD(2) * t1959;
t1944 = t1928 * qJDD(1);
t1884 = -t1907 + t1944;
t1854 = t1883 * t1921 - t1884 * t1920;
t1870 = qJD(2) * t1875;
t1838 = -t1854 - t1870;
t1847 = t1875 * t1923 + t1876 * t1927;
t1913 = qJD(2) - qJD(5);
t1922 = sin(qJ(6));
t1926 = cos(qJ(6));
t1831 = t1847 * t1922 + t1926 * t1913;
t1968 = t1831 ^ 2;
t1833 = t1847 * t1926 - t1913 * t1922;
t1967 = t1833 ^ 2;
t1966 = t1844 ^ 2;
t1965 = t1845 ^ 2;
t1964 = t1847 ^ 2;
t1874 = t1875 ^ 2;
t1963 = t1876 ^ 2;
t1962 = t1913 ^ 2;
t1961 = 2 * qJD(3);
t1960 = t1928 * g(3);
t1957 = qJD(2) * t1876;
t1955 = t1831 * t1833;
t1954 = t1845 * t1847;
t1953 = t1875 * t1876;
t1918 = t1928 ^ 2;
t1950 = t1918 * t1931;
t1948 = qJD(5) + t1913;
t1947 = qJD(6) - t1844;
t1894 = -g(1) * t1929 - g(2) * t1925;
t1879 = -pkin(1) * t1931 + qJDD(1) * pkin(7) + t1894;
t1864 = -g(3) * t1924 + t1928 * t1879;
t1881 = (-pkin(2) * t1928 - qJ(3) * t1924) * qJD(1);
t1830 = -pkin(2) * t1930 + qJDD(2) * qJ(3) + qJD(2) * t1961 + t1881 * t1958 + t1864;
t1890 = -qJD(2) * pkin(3) - qJ(4) * t1959;
t1821 = -pkin(3) * t1950 - qJ(4) * t1884 + qJD(2) * t1890 + t1830;
t1936 = -qJDD(2) * pkin(2) - t1930 * qJ(3) + qJDD(3) + t1960;
t1941 = qJD(1) * t1881 + t1879;
t1822 = -qJDD(2) * pkin(3) + (-t1883 + t1942) * qJ(4) + (-pkin(3) * t1949 + t1941) * t1924 + t1936;
t1790 = 0.2e1 * qJD(4) * t1875 + t1921 * t1821 + t1920 * t1822;
t1853 = -t1883 * t1920 - t1921 * t1884;
t1865 = -qJD(2) * pkin(4) - pkin(8) * t1876;
t1784 = -pkin(4) * t1874 + pkin(8) * t1853 + qJD(2) * t1865 + t1790;
t1789 = -0.2e1 * qJD(4) * t1876 - t1920 * t1821 + t1921 * t1822;
t1850 = -qJDD(2) + t1953;
t1932 = t1850 * pkin(4) + pkin(8) * t1838 + t1789;
t1761 = t1927 * t1784 + t1923 * t1932;
t1946 = t1917 + t1918;
t1943 = qJDD(2) - qJDD(5);
t1893 = t1925 * g(1) - t1929 * g(2);
t1760 = -t1923 * t1784 + t1927 * t1932;
t1937 = -t1923 * t1853 - t1927 * t1854;
t1808 = -qJD(5) * t1845 - t1937;
t1940 = -t1845 * t1913 - t1808;
t1939 = -t1922 * t1808 - t1926 * t1943;
t1938 = -t1927 * t1853 + t1923 * t1854;
t1858 = t1892 * t1924 + t1898 * t1928;
t1935 = -t1926 * t1808 + t1922 * t1943;
t1934 = -qJD(5) * t1847 - qJDD(6) - t1938;
t1796 = (qJD(5) - t1913) * t1847 + t1938;
t1878 = qJDD(1) * pkin(1) + t1931 * pkin(7) + t1893;
t1933 = t1878 + (t1884 - t1907) * pkin(2);
t1813 = t1883 * qJ(3) + qJDD(4) + t1884 * pkin(3) - qJ(4) * t1950 + (qJ(3) * t1956 + (t1961 + t1890) * t1924) * qJD(1) + t1933;
t1791 = -t1853 * pkin(4) - t1874 * pkin(8) + t1876 * t1865 + t1813;
t1899 = -t1930 - t1950;
t1891 = qJDD(2) + t1900;
t1889 = t1946 * t1931;
t1888 = -qJDD(1) * t1925 - t1929 * t1931;
t1887 = qJDD(1) * t1929 - t1925 * t1931;
t1886 = t1946 * qJDD(1);
t1885 = -0.2e1 * t1907 + t1944;
t1866 = -t1930 - t1963;
t1863 = -t1924 * t1879 - t1960;
t1860 = -t1891 * t1924 + t1899 * t1928;
t1857 = t1891 * t1928 + t1899 * t1924;
t1856 = t1886 * t1929 - t1889 * t1925;
t1855 = t1886 * t1925 + t1889 * t1929;
t1849 = qJDD(2) + t1953;
t1848 = -t1930 - t1874;
t1843 = t1860 * t1929 - t1885 * t1925;
t1842 = t1860 * t1925 + t1885 * t1929;
t1840 = t1941 * t1924 + t1936;
t1839 = -t1962 - t1964;
t1837 = t1854 - t1870;
t1836 = t1853 - t1957;
t1835 = -t1853 - t1957;
t1834 = -t1874 - t1963;
t1829 = -t1863 * t1924 + t1864 * t1928;
t1828 = t1863 * t1928 + t1864 * t1924;
t1827 = t1959 * t1961 + (t1883 + t1942) * qJ(3) + t1933;
t1826 = t1849 * t1921 - t1866 * t1920;
t1825 = t1849 * t1920 + t1866 * t1921;
t1824 = t1848 * t1921 - t1850 * t1920;
t1823 = t1848 * t1920 + t1850 * t1921;
t1820 = pkin(5) * t1845 - pkin(9) * t1847;
t1819 = -t1943 - t1954;
t1818 = t1943 - t1954;
t1816 = -t1962 - t1965;
t1812 = t1830 * t1928 + t1840 * t1924;
t1811 = t1830 * t1924 - t1840 * t1928;
t1810 = t1836 * t1921 - t1838 * t1920;
t1809 = t1836 * t1920 + t1838 * t1921;
t1807 = -t1964 - t1965;
t1806 = t1825 * t1924 + t1826 * t1928;
t1805 = -t1825 * t1928 + t1826 * t1924;
t1804 = -t1966 - t1967;
t1803 = t1818 * t1927 - t1839 * t1923;
t1802 = t1818 * t1923 + t1839 * t1927;
t1801 = -t1966 - t1968;
t1800 = -t1967 - t1968;
t1799 = t1845 * t1948 + t1937;
t1797 = -t1847 * t1948 - t1938;
t1795 = t1823 * t1924 + t1824 * t1928;
t1794 = -t1823 * t1928 + t1824 * t1924;
t1793 = t1816 * t1927 - t1819 * t1923;
t1792 = t1816 * t1923 + t1819 * t1927;
t1788 = t1934 - t1955;
t1787 = -t1934 - t1955;
t1786 = t1809 * t1924 + t1810 * t1928;
t1785 = -t1809 * t1928 + t1810 * t1924;
t1782 = t1831 * t1947 + t1935;
t1781 = -t1831 * t1969 - t1935;
t1780 = -t1833 * t1947 + t1939;
t1779 = t1833 * t1969 - t1939;
t1776 = -t1802 * t1920 + t1803 * t1921;
t1775 = t1802 * t1921 + t1803 * t1920;
t1774 = t1797 * t1927 - t1799 * t1923;
t1773 = t1797 * t1923 + t1799 * t1927;
t1772 = t1788 * t1926 - t1804 * t1922;
t1771 = t1788 * t1922 + t1804 * t1926;
t1770 = -t1792 * t1920 + t1793 * t1921;
t1769 = t1792 * t1921 + t1793 * t1920;
t1768 = -t1787 * t1922 + t1801 * t1926;
t1767 = t1787 * t1926 + t1801 * t1922;
t1766 = -t1789 * t1920 + t1790 * t1921;
t1765 = t1789 * t1921 + t1790 * t1920;
t1764 = pkin(5) * t1796 + pkin(9) * t1940 + t1791;
t1763 = t1780 * t1926 - t1782 * t1922;
t1762 = t1780 * t1922 + t1782 * t1926;
t1759 = t1775 * t1924 + t1776 * t1928;
t1758 = -t1775 * t1928 + t1776 * t1924;
t1757 = t1772 * t1927 + t1781 * t1923;
t1756 = t1772 * t1923 - t1781 * t1927;
t1755 = t1768 * t1927 + t1779 * t1923;
t1754 = t1768 * t1923 - t1779 * t1927;
t1753 = -t1773 * t1920 + t1774 * t1921;
t1752 = t1773 * t1921 + t1774 * t1920;
t1751 = -pkin(5) * t1962 - pkin(9) * t1943 - t1845 * t1820 + t1761;
t1750 = pkin(5) * t1943 - pkin(9) * t1962 + t1847 * t1820 - t1760;
t1749 = t1763 * t1927 + t1800 * t1923;
t1748 = t1763 * t1923 - t1800 * t1927;
t1747 = t1769 * t1924 + t1770 * t1928;
t1746 = -t1769 * t1928 + t1770 * t1924;
t1745 = t1765 * t1924 + t1766 * t1928;
t1744 = -t1765 * t1928 + t1766 * t1924;
t1743 = -t1760 * t1923 + t1761 * t1927;
t1742 = t1760 * t1927 + t1761 * t1923;
t1741 = t1751 * t1926 + t1764 * t1922;
t1740 = -t1751 * t1922 + t1764 * t1926;
t1739 = -t1756 * t1920 + t1757 * t1921;
t1738 = t1756 * t1921 + t1757 * t1920;
t1737 = -t1754 * t1920 + t1755 * t1921;
t1736 = t1754 * t1921 + t1755 * t1920;
t1735 = t1752 * t1924 + t1753 * t1928;
t1734 = -t1752 * t1928 + t1753 * t1924;
t1733 = -t1748 * t1920 + t1749 * t1921;
t1732 = t1748 * t1921 + t1749 * t1920;
t1731 = -t1742 * t1920 + t1743 * t1921;
t1730 = t1742 * t1921 + t1743 * t1920;
t1729 = -t1740 * t1922 + t1741 * t1926;
t1728 = t1740 * t1926 + t1741 * t1922;
t1727 = t1738 * t1924 + t1739 * t1928;
t1726 = -t1738 * t1928 + t1739 * t1924;
t1725 = t1736 * t1924 + t1737 * t1928;
t1724 = -t1736 * t1928 + t1737 * t1924;
t1723 = t1732 * t1924 + t1733 * t1928;
t1722 = -t1732 * t1928 + t1733 * t1924;
t1721 = t1729 * t1927 + t1750 * t1923;
t1720 = t1729 * t1923 - t1750 * t1927;
t1719 = t1730 * t1924 + t1731 * t1928;
t1718 = -t1730 * t1928 + t1731 * t1924;
t1717 = -t1720 * t1920 + t1721 * t1921;
t1716 = t1720 * t1921 + t1721 * t1920;
t1715 = t1716 * t1924 + t1717 * t1928;
t1714 = -t1716 * t1928 + t1717 * t1924;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1888, -t1887, 0, -t1893 * t1925 + t1894 * t1929, 0, 0, 0, 0, 0, 0, t1843, -t1970, t1856, t1829 * t1929 - t1878 * t1925, 0, 0, 0, 0, 0, 0, t1843, t1856, t1970, t1812 * t1929 - t1827 * t1925, 0, 0, 0, 0, 0, 0, t1795 * t1929 - t1835 * t1925, t1806 * t1929 - t1837 * t1925, t1786 * t1929 - t1834 * t1925, t1745 * t1929 - t1813 * t1925, 0, 0, 0, 0, 0, 0, t1747 * t1929 - t1796 * t1925, t1759 * t1929 + t1925 * t1940, t1735 * t1929 - t1807 * t1925, t1719 * t1929 - t1791 * t1925, 0, 0, 0, 0, 0, 0, t1725 * t1929 - t1767 * t1925, t1727 * t1929 - t1771 * t1925, t1723 * t1929 - t1762 * t1925, t1715 * t1929 - t1728 * t1925; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1887, t1888, 0, t1893 * t1929 + t1894 * t1925, 0, 0, 0, 0, 0, 0, t1842, -t1971, t1855, t1829 * t1925 + t1878 * t1929, 0, 0, 0, 0, 0, 0, t1842, t1855, t1971, t1812 * t1925 + t1827 * t1929, 0, 0, 0, 0, 0, 0, t1795 * t1925 + t1835 * t1929, t1806 * t1925 + t1837 * t1929, t1786 * t1925 + t1834 * t1929, t1745 * t1925 + t1813 * t1929, 0, 0, 0, 0, 0, 0, t1747 * t1925 + t1796 * t1929, t1759 * t1925 - t1929 * t1940, t1735 * t1925 + t1807 * t1929, t1719 * t1925 + t1791 * t1929, 0, 0, 0, 0, 0, 0, t1725 * t1925 + t1767 * t1929, t1727 * t1925 + t1771 * t1929, t1723 * t1925 + t1762 * t1929, t1715 * t1925 + t1728 * t1929; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1857, -t1858, 0, t1828, 0, 0, 0, 0, 0, 0, t1857, 0, t1858, t1811, 0, 0, 0, 0, 0, 0, t1794, t1805, t1785, t1744, 0, 0, 0, 0, 0, 0, t1746, t1758, t1734, t1718, 0, 0, 0, 0, 0, 0, t1724, t1726, t1722, t1714; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1931, -qJDD(1), 0, t1894, 0, 0, 0, 0, 0, 0, t1860, -t1861, t1886, t1829, 0, 0, 0, 0, 0, 0, t1860, t1886, t1861, t1812, 0, 0, 0, 0, 0, 0, t1795, t1806, t1786, t1745, 0, 0, 0, 0, 0, 0, t1747, t1759, t1735, t1719, 0, 0, 0, 0, 0, 0, t1725, t1727, t1723, t1715; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1931, 0, t1893, 0, 0, 0, 0, 0, 0, t1885, -t1882, t1889, t1878, 0, 0, 0, 0, 0, 0, t1885, t1889, t1882, t1827, 0, 0, 0, 0, 0, 0, t1835, t1837, t1834, t1813, 0, 0, 0, 0, 0, 0, t1796, -t1940, t1807, t1791, 0, 0, 0, 0, 0, 0, t1767, t1771, t1762, t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1857, -t1858, 0, t1828, 0, 0, 0, 0, 0, 0, t1857, 0, t1858, t1811, 0, 0, 0, 0, 0, 0, t1794, t1805, t1785, t1744, 0, 0, 0, 0, 0, 0, t1746, t1758, t1734, t1718, 0, 0, 0, 0, 0, 0, t1724, t1726, t1722, t1714; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1899, -t1892, t1944, t1864, 0, 0, 0, 0, 0, 0, t1899, t1944, t1892, t1830, 0, 0, 0, 0, 0, 0, t1824, t1826, t1810, t1766, 0, 0, 0, 0, 0, 0, t1770, t1776, t1753, t1731, 0, 0, 0, 0, 0, 0, t1737, t1739, t1733, t1717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1891, -t1898, -t1945, t1863, 0, 0, 0, 0, 0, 0, t1891, -t1945, t1898, -t1840, 0, 0, 0, 0, 0, 0, -t1823, -t1825, -t1809, -t1765, 0, 0, 0, 0, 0, 0, -t1769, -t1775, -t1752, -t1730, 0, 0, 0, 0, 0, 0, -t1736, -t1738, -t1732, -t1716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1885, t1882, -t1889, -t1878, 0, 0, 0, 0, 0, 0, -t1885, -t1889, -t1882, -t1827, 0, 0, 0, 0, 0, 0, -t1835, -t1837, -t1834, -t1813, 0, 0, 0, 0, 0, 0, -t1796, t1940, -t1807, -t1791, 0, 0, 0, 0, 0, 0, -t1767, -t1771, -t1762, -t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1899, t1944, t1892, t1830, 0, 0, 0, 0, 0, 0, t1824, t1826, t1810, t1766, 0, 0, 0, 0, 0, 0, t1770, t1776, t1753, t1731, 0, 0, 0, 0, 0, 0, t1737, t1739, t1733, t1717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1885, -t1889, -t1882, -t1827, 0, 0, 0, 0, 0, 0, -t1835, -t1837, -t1834, -t1813, 0, 0, 0, 0, 0, 0, -t1796, t1940, -t1807, -t1791, 0, 0, 0, 0, 0, 0, -t1767, -t1771, -t1762, -t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1891, t1945, -t1898, t1840, 0, 0, 0, 0, 0, 0, t1823, t1825, t1809, t1765, 0, 0, 0, 0, 0, 0, t1769, t1775, t1752, t1730, 0, 0, 0, 0, 0, 0, t1736, t1738, t1732, t1716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1848, t1849, t1836, t1790, 0, 0, 0, 0, 0, 0, t1793, t1803, t1774, t1743, 0, 0, 0, 0, 0, 0, t1755, t1757, t1749, t1721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1850, t1866, t1838, t1789, 0, 0, 0, 0, 0, 0, t1792, t1802, t1773, t1742, 0, 0, 0, 0, 0, 0, t1754, t1756, t1748, t1720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1835, t1837, t1834, t1813, 0, 0, 0, 0, 0, 0, t1796, -t1940, t1807, t1791, 0, 0, 0, 0, 0, 0, t1767, t1771, t1762, t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1816, t1818, t1797, t1761, 0, 0, 0, 0, 0, 0, t1768, t1772, t1763, t1729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1819, t1839, t1799, t1760, 0, 0, 0, 0, 0, 0, -t1779, -t1781, -t1800, -t1750; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1796, -t1940, t1807, t1791, 0, 0, 0, 0, 0, 0, t1767, t1771, t1762, t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1801, t1788, t1780, t1741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1787, t1804, t1782, t1740; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1779, t1781, t1800, t1750;];
f_new_reg  = t1;