% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRPR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRPR10_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:31:22
% EndTime: 2019-12-31 21:31:28
% DurationCPUTime: 6.08s
% Computational Cost: add. (45942->290), mult. (102617->439), div. (0->0), fcn. (82093->12), ass. (0->247)
t1950 = sin(pkin(5));
t1961 = qJD(1) ^ 2;
t2022 = t1950 * t1961;
t1952 = cos(pkin(5));
t1945 = t1952 * qJD(1) + qJD(2);
t1954 = sin(qJ(3));
t1958 = cos(qJ(3));
t1955 = sin(qJ(2));
t2001 = t1950 * t1955;
t1994 = qJD(1) * t2001;
t1913 = -t1958 * t1945 + t1954 * t1994;
t1915 = t1954 * t1945 + t1958 * t1994;
t1949 = sin(pkin(10));
t1951 = cos(pkin(10));
t1894 = t1951 * t1913 + t1949 * t1915;
t1893 = qJD(5) + t1894;
t2021 = qJD(5) + t1893;
t1896 = -t1949 * t1913 + t1951 * t1915;
t1959 = cos(qJ(2));
t2000 = t1950 * t1959;
t1993 = qJD(1) * t2000;
t1935 = -qJD(3) + t1993;
t1953 = sin(qJ(5));
t1957 = cos(qJ(5));
t1877 = t1953 * t1896 + t1957 * t1935;
t2020 = t1877 ^ 2;
t1879 = t1957 * t1896 - t1953 * t1935;
t2019 = t1879 ^ 2;
t2018 = t1893 ^ 2;
t2017 = t1894 ^ 2;
t2016 = t1896 ^ 2;
t2015 = t1913 ^ 2;
t2014 = t1915 ^ 2;
t2013 = t1935 ^ 2;
t2012 = t1945 ^ 2;
t2011 = -2 * qJD(4);
t2010 = t1952 * g(3);
t2009 = (-pkin(2) * t1959 - pkin(8) * t1955) * t2022;
t2008 = t1879 * t1877;
t2007 = t1894 * t1935;
t2006 = t1896 * t1894;
t2005 = t1913 * t1935;
t2004 = t1915 * t1913;
t2003 = t1935 * t1896;
t2002 = t1950 ^ 2 * t1961;
t1999 = qJD(3) + t1935;
t1998 = qJD(5) - t1893;
t1985 = t1952 * qJDD(1) + qJDD(2);
t1956 = sin(qJ(1));
t1960 = cos(qJ(1));
t1938 = -t1960 * g(1) - t1956 * g(2);
t1995 = qJDD(1) * t1950;
t1919 = -t1961 * pkin(1) + pkin(7) * t1995 + t1938;
t1937 = t1956 * g(1) - t1960 * g(2);
t1964 = qJDD(1) * pkin(1) + pkin(7) * t2022 + t1937;
t1963 = t1952 * t1964;
t1997 = t1959 * t1919 + t1955 * t1963;
t1865 = t1985 * pkin(8) - t2012 * pkin(2) + (-t1955 * g(3) + t1959 * t2009) * t1950 + t1997;
t1923 = qJD(2) * t1993 + t1955 * t1995;
t1996 = -qJD(2) * t1994 + t1959 * t1995;
t1866 = -t2010 - t1923 * pkin(8) - t1996 * pkin(2) + ((pkin(2) * t1955 - pkin(8) * t1959) * t1945 * qJD(1) - t1964) * t1950;
t1837 = t1958 * t1865 + t1954 * t1866;
t1986 = t1954 * t1923 - t1958 * t1985;
t1886 = -t1915 * qJD(3) - t1986;
t1900 = -t1935 * pkin(3) - t1915 * qJ(4);
t1824 = -t2015 * pkin(3) + t1886 * qJ(4) + t1935 * t1900 + t1837;
t1836 = -t1954 * t1865 + t1958 * t1866;
t1992 = -qJDD(3) + t1996;
t1882 = -t1992 - t2004;
t1965 = -t1958 * t1923 - t1954 * t1985;
t1887 = -t1913 * qJD(3) - t1965;
t1962 = (-t1887 + t2005) * qJ(4) + t1882 * pkin(3) + t1836;
t1792 = t1951 * t1824 + t1894 * t2011 + t1949 * t1962;
t1991 = t1949 * t1824 - t1951 * t1962;
t1854 = t1949 * t1886 + t1951 * t1887;
t1990 = -t1854 - t2007;
t1989 = -t1953 * t1854 - t1957 * t1992;
t1988 = -t1951 * t1886 + t1949 * t1887;
t1987 = t1955 * t1919 - t1959 * t1963;
t1984 = t1945 * t1993;
t1983 = -qJDD(5) - t1988;
t1859 = t1894 * pkin(4) - t1896 * pkin(9);
t1779 = -pkin(4) * t2013 - t1992 * pkin(9) - t1894 * t1859 + t1792;
t1864 = -t1985 * pkin(2) - t2012 * pkin(8) + (t1959 * g(3) + t1955 * t2009) * t1950 + t1987;
t1830 = -t1886 * pkin(3) - t2015 * qJ(4) + t1915 * t1900 + qJDD(4) + t1864;
t1841 = t1988 - t2003;
t1796 = t1841 * pkin(4) + t1990 * pkin(9) + t1830;
t1770 = -t1953 * t1779 + t1957 * t1796;
t1771 = t1957 * t1779 + t1953 * t1796;
t1756 = -t1953 * t1770 + t1957 * t1771;
t1778 = t1992 * pkin(4) - t2013 * pkin(9) + ((2 * qJD(4)) + t1859) * t1896 + t1991;
t1748 = t1949 * t1756 - t1951 * t1778;
t1749 = t1951 * t1756 + t1949 * t1778;
t1743 = -t1954 * t1748 + t1958 * t1749;
t1755 = t1957 * t1770 + t1953 * t1771;
t1982 = t1743 * t1955 - t1755 * t1959;
t1791 = t1896 * t2011 - t1991;
t1772 = t1951 * t1791 + t1949 * t1792;
t1773 = -t1949 * t1791 + t1951 * t1792;
t1758 = -t1954 * t1772 + t1958 * t1773;
t1981 = t1758 * t1955 - t1830 * t1959;
t1818 = -t1998 * t1879 + t1989;
t1966 = -t1957 * t1854 + t1953 * t1992;
t1820 = t1998 * t1877 + t1966;
t1795 = t1957 * t1818 - t1953 * t1820;
t1833 = -t2019 - t2020;
t1783 = t1949 * t1795 - t1951 * t1833;
t1784 = t1951 * t1795 + t1949 * t1833;
t1763 = -t1954 * t1783 + t1958 * t1784;
t1794 = t1953 * t1818 + t1957 * t1820;
t1980 = t1763 * t1955 - t1794 * t1959;
t1825 = -t1983 - t2008;
t1838 = -t2018 - t2020;
t1802 = -t1953 * t1825 + t1957 * t1838;
t1817 = t2021 * t1879 - t1989;
t1785 = t1949 * t1802 - t1951 * t1817;
t1786 = t1951 * t1802 + t1949 * t1817;
t1767 = -t1954 * t1785 + t1958 * t1786;
t1801 = t1957 * t1825 + t1953 * t1838;
t1979 = t1767 * t1955 - t1801 * t1959;
t1826 = t1983 - t2008;
t1845 = -t2018 - t2019;
t1806 = t1957 * t1826 - t1953 * t1845;
t1819 = -t2021 * t1877 - t1966;
t1787 = t1949 * t1806 - t1951 * t1819;
t1788 = t1951 * t1806 + t1949 * t1819;
t1769 = -t1954 * t1787 + t1958 * t1788;
t1805 = t1953 * t1826 + t1957 * t1845;
t1978 = t1769 * t1955 - t1805 * t1959;
t1842 = -t1988 - t2003;
t1844 = -t1854 + t2007;
t1811 = t1949 * t1842 + t1951 * t1844;
t1812 = t1951 * t1842 - t1949 * t1844;
t1790 = -t1954 * t1811 + t1958 * t1812;
t1846 = -t2016 - t2017;
t1977 = t1790 * t1955 - t1846 * t1959;
t1852 = -t1992 - t2006;
t1856 = -t2013 - t2017;
t1828 = t1951 * t1852 + t1949 * t1856;
t1829 = -t1949 * t1852 + t1951 * t1856;
t1800 = -t1954 * t1828 + t1958 * t1829;
t1976 = t1800 * t1955 - t1841 * t1959;
t1851 = t1992 - t2006;
t1875 = -t2013 - t2016;
t1834 = t1949 * t1851 + t1951 * t1875;
t1835 = t1951 * t1851 - t1949 * t1875;
t1808 = -t1954 * t1834 + t1958 * t1835;
t1975 = t1808 * t1955 + t1959 * t1990;
t1810 = -t1954 * t1836 + t1958 * t1837;
t1974 = t1810 * t1955 - t1864 * t1959;
t1870 = -t1999 * t1915 - t1986;
t1872 = t1999 * t1913 + t1965;
t1840 = t1958 * t1870 - t1954 * t1872;
t1880 = -t2014 - t2015;
t1973 = t1840 * t1955 - t1880 * t1959;
t1888 = -t2013 - t2015;
t1850 = -t1954 * t1882 + t1958 * t1888;
t1869 = (qJD(3) - t1935) * t1915 + t1986;
t1972 = t1850 * t1955 - t1869 * t1959;
t1881 = t1992 - t2004;
t1898 = -t2013 - t2014;
t1858 = t1958 * t1881 - t1954 * t1898;
t1871 = t1887 + t2005;
t1971 = t1858 * t1955 - t1871 * t1959;
t1889 = -g(3) * t2000 - t1987;
t1890 = -g(3) * t2001 + t1997;
t1970 = t1889 * t1959 + t1890 * t1955;
t1902 = t1984 - t1923;
t1926 = t1945 * t1994;
t1903 = t1926 + t1996;
t1969 = t1902 * t1959 + t1903 * t1955;
t1947 = t1955 ^ 2;
t1910 = -t1947 * t2002 - t2012;
t1934 = t1959 * t1955 * t2002;
t1921 = t1934 - t1985;
t1968 = t1910 * t1959 + t1921 * t1955;
t1920 = t1934 + t1985;
t1948 = t1959 ^ 2;
t1924 = -t1948 * t2002 - t2012;
t1967 = t1920 * t1959 + t1924 * t1955;
t1933 = -t1956 * qJDD(1) - t1960 * t1961;
t1932 = t1960 * qJDD(1) - t1956 * t1961;
t1925 = (-t1947 - t1948) * t2002;
t1906 = -t1950 * t1964 - t2010;
t1904 = t1926 - t1996;
t1901 = t1984 + t1923;
t1899 = -t1955 * t1920 + t1959 * t1924;
t1897 = -t1955 * t1910 + t1959 * t1921;
t1876 = -t1955 * t1902 + t1959 * t1903;
t1874 = -t1950 * t1904 + t1967 * t1952;
t1873 = t1952 * t1904 + t1967 * t1950;
t1868 = -t1950 * t1901 + t1968 * t1952;
t1867 = t1952 * t1901 + t1968 * t1950;
t1863 = -t1950 * t1925 + t1969 * t1952;
t1862 = t1952 * t1925 + t1969 * t1950;
t1857 = t1954 * t1881 + t1958 * t1898;
t1855 = -t1955 * t1889 + t1959 * t1890;
t1849 = t1958 * t1882 + t1954 * t1888;
t1848 = -t1950 * t1906 + t1970 * t1952;
t1847 = t1952 * t1906 + t1970 * t1950;
t1839 = t1954 * t1870 + t1958 * t1872;
t1832 = t1959 * t1858 + t1955 * t1871;
t1831 = t1959 * t1850 + t1955 * t1869;
t1827 = t1959 * t1840 + t1955 * t1880;
t1823 = -t1950 * t1857 + t1971 * t1952;
t1822 = t1952 * t1857 + t1971 * t1950;
t1816 = -t1950 * t1849 + t1972 * t1952;
t1815 = t1952 * t1849 + t1972 * t1950;
t1809 = t1958 * t1836 + t1954 * t1837;
t1807 = t1958 * t1834 + t1954 * t1835;
t1804 = -t1950 * t1839 + t1973 * t1952;
t1803 = t1952 * t1839 + t1973 * t1950;
t1799 = t1958 * t1828 + t1954 * t1829;
t1798 = t1959 * t1810 + t1955 * t1864;
t1797 = t1959 * t1808 - t1955 * t1990;
t1793 = t1959 * t1800 + t1955 * t1841;
t1789 = t1958 * t1811 + t1954 * t1812;
t1782 = t1959 * t1790 + t1955 * t1846;
t1781 = -t1950 * t1809 + t1974 * t1952;
t1780 = t1952 * t1809 + t1974 * t1950;
t1777 = -t1950 * t1807 + t1975 * t1952;
t1776 = t1952 * t1807 + t1975 * t1950;
t1775 = -t1950 * t1799 + t1976 * t1952;
t1774 = t1952 * t1799 + t1976 * t1950;
t1768 = t1958 * t1787 + t1954 * t1788;
t1766 = t1958 * t1785 + t1954 * t1786;
t1765 = -t1950 * t1789 + t1977 * t1952;
t1764 = t1952 * t1789 + t1977 * t1950;
t1762 = t1958 * t1783 + t1954 * t1784;
t1761 = t1959 * t1769 + t1955 * t1805;
t1760 = t1959 * t1767 + t1955 * t1801;
t1759 = t1959 * t1763 + t1955 * t1794;
t1757 = t1958 * t1772 + t1954 * t1773;
t1754 = t1959 * t1758 + t1955 * t1830;
t1753 = -t1950 * t1768 + t1978 * t1952;
t1752 = t1952 * t1768 + t1978 * t1950;
t1751 = -t1950 * t1766 + t1979 * t1952;
t1750 = t1952 * t1766 + t1979 * t1950;
t1747 = -t1950 * t1762 + t1980 * t1952;
t1746 = t1952 * t1762 + t1980 * t1950;
t1745 = -t1950 * t1757 + t1981 * t1952;
t1744 = t1952 * t1757 + t1981 * t1950;
t1742 = t1958 * t1748 + t1954 * t1749;
t1741 = t1959 * t1743 + t1955 * t1755;
t1740 = -t1950 * t1742 + t1982 * t1952;
t1739 = t1952 * t1742 + t1982 * t1950;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1933, -t1932, 0, -t1956 * t1937 + t1960 * t1938, 0, 0, 0, 0, 0, 0, -t1956 * t1874 + t1960 * t1899, -t1956 * t1868 + t1960 * t1897, -t1956 * t1863 + t1960 * t1876, -t1956 * t1848 + t1960 * t1855, 0, 0, 0, 0, 0, 0, -t1956 * t1816 + t1960 * t1831, -t1956 * t1823 + t1960 * t1832, -t1956 * t1804 + t1960 * t1827, -t1956 * t1781 + t1960 * t1798, 0, 0, 0, 0, 0, 0, -t1956 * t1775 + t1960 * t1793, -t1956 * t1777 + t1960 * t1797, -t1956 * t1765 + t1960 * t1782, -t1956 * t1745 + t1960 * t1754, 0, 0, 0, 0, 0, 0, -t1956 * t1751 + t1960 * t1760, -t1956 * t1753 + t1960 * t1761, -t1956 * t1747 + t1960 * t1759, -t1956 * t1740 + t1960 * t1741; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1932, t1933, 0, t1960 * t1937 + t1956 * t1938, 0, 0, 0, 0, 0, 0, t1960 * t1874 + t1956 * t1899, t1960 * t1868 + t1956 * t1897, t1960 * t1863 + t1956 * t1876, t1960 * t1848 + t1956 * t1855, 0, 0, 0, 0, 0, 0, t1960 * t1816 + t1956 * t1831, t1960 * t1823 + t1956 * t1832, t1960 * t1804 + t1956 * t1827, t1960 * t1781 + t1956 * t1798, 0, 0, 0, 0, 0, 0, t1960 * t1775 + t1956 * t1793, t1960 * t1777 + t1956 * t1797, t1960 * t1765 + t1956 * t1782, t1960 * t1745 + t1956 * t1754, 0, 0, 0, 0, 0, 0, t1960 * t1751 + t1956 * t1760, t1960 * t1753 + t1956 * t1761, t1960 * t1747 + t1956 * t1759, t1960 * t1740 + t1956 * t1741; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1873, t1867, t1862, t1847, 0, 0, 0, 0, 0, 0, t1815, t1822, t1803, t1780, 0, 0, 0, 0, 0, 0, t1774, t1776, t1764, t1744, 0, 0, 0, 0, 0, 0, t1750, t1752, t1746, t1739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1961, -qJDD(1), 0, t1938, 0, 0, 0, 0, 0, 0, t1899, t1897, t1876, t1855, 0, 0, 0, 0, 0, 0, t1831, t1832, t1827, t1798, 0, 0, 0, 0, 0, 0, t1793, t1797, t1782, t1754, 0, 0, 0, 0, 0, 0, t1760, t1761, t1759, t1741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1961, 0, t1937, 0, 0, 0, 0, 0, 0, t1874, t1868, t1863, t1848, 0, 0, 0, 0, 0, 0, t1816, t1823, t1804, t1781, 0, 0, 0, 0, 0, 0, t1775, t1777, t1765, t1745, 0, 0, 0, 0, 0, 0, t1751, t1753, t1747, t1740; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1873, t1867, t1862, t1847, 0, 0, 0, 0, 0, 0, t1815, t1822, t1803, t1780, 0, 0, 0, 0, 0, 0, t1774, t1776, t1764, t1744, 0, 0, 0, 0, 0, 0, t1750, t1752, t1746, t1739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1924, t1921, t1903, t1890, 0, 0, 0, 0, 0, 0, t1850, t1858, t1840, t1810, 0, 0, 0, 0, 0, 0, t1800, t1808, t1790, t1758, 0, 0, 0, 0, 0, 0, t1767, t1769, t1763, t1743; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1920, t1910, t1902, t1889, 0, 0, 0, 0, 0, 0, -t1869, -t1871, -t1880, -t1864, 0, 0, 0, 0, 0, 0, -t1841, t1990, -t1846, -t1830, 0, 0, 0, 0, 0, 0, -t1801, -t1805, -t1794, -t1755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1904, t1901, t1925, t1906, 0, 0, 0, 0, 0, 0, t1849, t1857, t1839, t1809, 0, 0, 0, 0, 0, 0, t1799, t1807, t1789, t1757, 0, 0, 0, 0, 0, 0, t1766, t1768, t1762, t1742; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1888, t1881, t1870, t1837, 0, 0, 0, 0, 0, 0, t1829, t1835, t1812, t1773, 0, 0, 0, 0, 0, 0, t1786, t1788, t1784, t1749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1882, t1898, t1872, t1836, 0, 0, 0, 0, 0, 0, t1828, t1834, t1811, t1772, 0, 0, 0, 0, 0, 0, t1785, t1787, t1783, t1748; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1869, t1871, t1880, t1864, 0, 0, 0, 0, 0, 0, t1841, -t1990, t1846, t1830, 0, 0, 0, 0, 0, 0, t1801, t1805, t1794, t1755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1856, t1851, t1842, t1792, 0, 0, 0, 0, 0, 0, t1802, t1806, t1795, t1756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1852, t1875, t1844, t1791, 0, 0, 0, 0, 0, 0, -t1817, -t1819, -t1833, -t1778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1841, -t1990, t1846, t1830, 0, 0, 0, 0, 0, 0, t1801, t1805, t1794, t1755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1838, t1826, t1818, t1771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1825, t1845, t1820, t1770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1817, t1819, t1833, t1778;];
f_new_reg = t1;