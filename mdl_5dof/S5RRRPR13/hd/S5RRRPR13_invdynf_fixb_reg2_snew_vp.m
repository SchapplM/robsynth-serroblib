% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRPR13_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:48:08
% EndTime: 2019-12-31 21:48:12
% DurationCPUTime: 4.65s
% Computational Cost: add. (17520->286), mult. (38529->383), div. (0->0), fcn. (29755->10), ass. (0->233)
t1988 = sin(pkin(5));
t1992 = sin(qJ(2));
t2058 = t1992 * qJD(1);
t2048 = t1988 * t2058;
t1980 = qJD(2) * t2048;
t1996 = cos(qJ(2));
t2049 = qJDD(1) * t1996;
t2035 = t1988 * t2049 - t1980;
t2016 = -qJDD(3) + t2035;
t1991 = sin(qJ(3));
t1995 = cos(qJ(3));
t1989 = cos(pkin(5));
t2042 = t1989 * qJD(1) + qJD(2);
t1952 = t1991 * t2048 - t1995 * t2042;
t1954 = t1991 * t2042 + t1995 * t2048;
t2061 = t1954 * t1952;
t1919 = t2016 + t2061;
t1951 = t1952 ^ 2;
t2056 = t1996 * qJD(1);
t2047 = t1988 * t2056;
t1976 = -qJD(3) + t2047;
t1972 = t1976 ^ 2;
t1925 = -t1972 - t1951;
t1886 = t1991 * t1919 + t1995 * t1925;
t2086 = t1886 * t1992;
t1885 = t1995 * t1919 - t1991 * t1925;
t2085 = t1988 * t1885;
t2084 = t1989 * t1885;
t2083 = t1996 * t1886;
t2011 = t2016 - t2061;
t2068 = t1954 ^ 2;
t2045 = -t1972 - t2068;
t1893 = t1991 * t2045 - t1995 * t2011;
t2082 = t1893 * t1992;
t1890 = t1991 * t2011 + t1995 * t2045;
t2081 = t1988 * t1890;
t2080 = t1989 * t1890;
t2050 = qJDD(1) * t1988;
t1963 = qJD(2) * t2047 + t1992 * t2050;
t2038 = t1989 * qJDD(1) + qJDD(2);
t2012 = t1995 * t1963 + t1991 * t2038;
t2054 = qJD(3) + t1976;
t1911 = t2054 * t1952 - t2012;
t2079 = t1991 * t1911;
t2078 = t1995 * t1911;
t2077 = t1996 * t1893;
t2008 = -t1952 * qJD(3) + t2012;
t2062 = t1952 * t1976;
t1910 = t2008 + t2062;
t2074 = -t1951 - t2068;
t2076 = t1992 * t2074;
t2075 = t1996 * t2074;
t2037 = t2042 ^ 2;
t2073 = -t2038 * pkin(2) - t2037 * pkin(8);
t1990 = sin(qJ(5));
t1994 = cos(qJ(5));
t1962 = (-t1996 * pkin(2) - pkin(8) * t1992) * t1988 * qJD(1);
t1993 = sin(qJ(1));
t1997 = cos(qJ(1));
t1979 = -t1997 * g(1) - t1993 * g(2);
t1998 = qJD(1) ^ 2;
t1959 = -t1998 * pkin(1) + pkin(7) * t2050 + t1979;
t1978 = t1993 * g(1) - t1997 * g(2);
t2009 = t1998 * t1988 * pkin(7) + qJDD(1) * pkin(1) + t1978;
t2006 = t1989 * t2009;
t2051 = t1996 * t1959 + t1992 * t2006;
t2064 = t1992 * g(3);
t1899 = t2038 * pkin(8) - t2037 * pkin(2) + (t1962 * t2056 - t2064) * t1988 + t2051;
t2034 = qJD(1) * t2042;
t2014 = t1992 * t2034;
t2015 = t1996 * t2034;
t2065 = t1989 * g(3);
t1900 = t1980 * pkin(2) - t1963 * pkin(8) - t2065 + (-pkin(8) * t2015 + (t2014 - t2049) * pkin(2) - t2009) * t1988;
t1870 = -t1991 * t1899 + t1995 * t1900;
t2013 = pkin(3) * t2016 - t1972 * qJ(4) + qJDD(4) - t1870;
t2001 = t2016 * pkin(9) + t2013 + (t2008 - t2062) * pkin(4);
t2039 = t1991 * t1963 - t1995 * t2038;
t1924 = t1954 * qJD(3) + t2039;
t2040 = t1992 * t1959 - t1996 * t2006;
t2066 = g(3) * t1996;
t1926 = -t1988 * t2066 - t2040;
t2072 = t1924 * pkin(3) - t1910 * qJ(4);
t2002 = -t1951 * pkin(4) + t1924 * pkin(9) + t1962 * t2048 - t1926 + t2072 + t2073;
t1936 = t1954 * pkin(4) + t1976 * pkin(9);
t2067 = -2 * qJD(4);
t2044 = -pkin(3) * t1976 + t2067;
t2036 = -t1936 + t2044;
t1929 = t1952 * pkin(3) - t1954 * qJ(4);
t2043 = t1952 * pkin(9) + t1929;
t1999 = -t1990 * t2002 + t1994 * t2001 + (-t1990 * t2036 + t1994 * t2043) * t1954;
t1933 = -t1994 * t1952 - t1990 * t1976;
t2071 = t1933 ^ 2;
t1935 = t1990 * t1952 - t1994 * t1976;
t2070 = t1935 ^ 2;
t1950 = qJD(5) + t1954;
t2069 = t1950 ^ 2;
t2063 = t1935 * t1933;
t2060 = t1988 ^ 2 * t1998;
t2055 = qJD(3) - t1976;
t2053 = qJD(5) - t1950;
t2052 = qJD(5) + t1950;
t2046 = -t2069 - t2070;
t2041 = t1994 * t1924 + t1990 * t2016;
t1827 = t1994 * t2002 + t1990 * t2001 + (t1990 * t2043 + t1994 * t2036) * t1954;
t1813 = t1990 * t1827 + t1999 * t1994;
t1871 = t1995 * t1899 + t1991 * t1900;
t2007 = -t1972 * pkin(3) - qJ(4) * t2016 - t1952 * t1929 + t1871;
t1840 = -t1924 * pkin(4) - t1951 * pkin(9) + (t2067 - t1936) * t1976 + t2007;
t1812 = t1991 * t1813 + t1995 * t1840;
t1814 = t1994 * t1827 - t1990 * t1999;
t2033 = t1812 * t1992 - t1814 * t1996;
t1857 = t1976 * t2067 + t2007;
t1858 = t1954 * t1929 + t2013;
t1832 = t1995 * t1857 + t1991 * t1858;
t1898 = (t1962 * t2058 + t2066) * t1988 + t2040 + t2073;
t1859 = t2044 * t1954 + t1898 + t2072;
t2032 = t1832 * t1992 - t1859 * t1996;
t1874 = -t2053 * t1935 + t2041;
t2021 = -t1990 * t1924 + t1994 * t2016;
t2010 = t2053 * t1933 + t2021;
t1847 = t1990 * t1874 + t1994 * t2010;
t1888 = -t2070 - t2071;
t1835 = t1991 * t1847 + t1995 * t1888;
t1848 = t1994 * t1874 - t1990 * t2010;
t2031 = t1835 * t1992 - t1848 * t1996;
t1894 = -t2069 - t2071;
t2005 = -qJDD(5) - t2008;
t2003 = -t2005 - t2063;
t1862 = t1990 * t1894 + t1994 * t2003;
t1873 = t2052 * t1935 - t2041;
t1842 = t1991 * t1862 + t1995 * t1873;
t1863 = t1994 * t1894 - t1990 * t2003;
t2030 = t1842 * t1992 - t1863 * t1996;
t1883 = t2005 - t2063;
t1864 = t1990 * t1883 + t1994 * t2046;
t1875 = -t2052 * t1933 - t2021;
t1844 = t1991 * t1864 + t1995 * t1875;
t1865 = t1994 * t1883 - t1990 * t2046;
t2029 = t1844 * t1992 - t1865 * t1996;
t1846 = -t1991 * t1870 + t1995 * t1871;
t2028 = t1846 * t1992 - t1898 * t1996;
t1905 = t2054 * t1954 + t2039;
t1878 = -t1995 * t1905 - t2079;
t2027 = t1878 * t1992 - t2075;
t1943 = t1976 * t1954;
t1907 = -t1924 - t1943;
t1879 = t1995 * t1907 - t2079;
t2026 = t1879 * t1992 - t2075;
t1904 = t2055 * t1954 + t2039;
t2025 = -t1904 * t1996 + t2086;
t1906 = t1924 - t1943;
t2024 = t1906 * t1996 - t2086;
t2023 = -t1910 * t1996 - t2082;
t1908 = -t2055 * t1952 + t2012;
t2022 = t1908 * t1996 + t2082;
t1927 = -t1988 * t2064 + t2051;
t2020 = t1926 * t1996 + t1927 * t1992;
t1939 = t1988 * t2015 - t1963;
t1966 = t1988 * t2014;
t1940 = t1966 + t2035;
t2019 = t1939 * t1996 + t1940 * t1992;
t1986 = t1992 ^ 2;
t1949 = -t1986 * t2060 - t2037;
t1975 = t1996 * t1992 * t2060;
t1961 = t1975 - t2038;
t2018 = t1949 * t1996 + t1961 * t1992;
t1960 = t1975 + t2038;
t1987 = t1996 ^ 2;
t1964 = -t1987 * t2060 - t2037;
t2017 = t1960 * t1996 + t1964 * t1992;
t1974 = -t1993 * qJDD(1) - t1997 * t1998;
t1973 = t1997 * qJDD(1) - t1993 * t1998;
t1965 = (-t1986 - t1987) * t2060;
t1944 = -t1988 * t2009 - t2065;
t1941 = t1966 - t2035;
t1938 = t2042 * t2047 + t1963;
t1932 = -t1992 * t1960 + t1996 * t1964;
t1928 = -t1992 * t1949 + t1996 * t1961;
t1914 = -t1992 * t1939 + t1996 * t1940;
t1913 = -t1988 * t1941 + t2017 * t1989;
t1912 = t1989 * t1941 + t2017 * t1988;
t1902 = -t1988 * t1938 + t2018 * t1989;
t1901 = t1989 * t1938 + t2018 * t1988;
t1897 = -t1988 * t1965 + t2019 * t1989;
t1896 = t1989 * t1965 + t2019 * t1988;
t1889 = -t1992 * t1926 + t1996 * t1927;
t1881 = -t1988 * t1944 + t2020 * t1989;
t1880 = t1989 * t1944 + t2020 * t1988;
t1877 = t1991 * t1907 + t2078;
t1876 = -t1991 * t1905 + t2078;
t1869 = -t1992 * t1908 + t2077;
t1868 = t1992 * t1910 - t2077;
t1867 = -t1992 * t1906 - t2083;
t1866 = t1992 * t1904 + t2083;
t1861 = t1996 * t1879 + t2076;
t1860 = t1996 * t1878 + t2076;
t1856 = t2022 * t1989 + t2081;
t1855 = t2023 * t1989 - t2081;
t1854 = t2022 * t1988 - t2080;
t1853 = t2023 * t1988 + t2080;
t1852 = t2024 * t1989 - t2085;
t1851 = t2025 * t1989 + t2085;
t1850 = t2024 * t1988 + t2084;
t1849 = t2025 * t1988 - t2084;
t1845 = t1995 * t1870 + t1991 * t1871;
t1843 = -t1995 * t1864 + t1991 * t1875;
t1841 = -t1995 * t1862 + t1991 * t1873;
t1839 = -t1988 * t1877 + t2026 * t1989;
t1838 = -t1988 * t1876 + t2027 * t1989;
t1837 = t1989 * t1877 + t2026 * t1988;
t1836 = t1989 * t1876 + t2027 * t1988;
t1834 = -t1995 * t1847 + t1991 * t1888;
t1833 = t1996 * t1846 + t1992 * t1898;
t1831 = t1991 * t1857 - t1995 * t1858;
t1830 = t1996 * t1844 + t1992 * t1865;
t1829 = t1996 * t1842 + t1992 * t1863;
t1828 = t1996 * t1835 + t1992 * t1848;
t1825 = -t1988 * t1845 + t2028 * t1989;
t1824 = t1989 * t1845 + t2028 * t1988;
t1823 = t1996 * t1832 + t1992 * t1859;
t1822 = -t1988 * t1843 + t2029 * t1989;
t1821 = t1989 * t1843 + t2029 * t1988;
t1820 = -t1988 * t1841 + t2030 * t1989;
t1819 = t1989 * t1841 + t2030 * t1988;
t1818 = -t1988 * t1834 + t2031 * t1989;
t1817 = t1989 * t1834 + t2031 * t1988;
t1816 = -t1988 * t1831 + t2032 * t1989;
t1815 = t1989 * t1831 + t2032 * t1988;
t1811 = -t1995 * t1813 + t1991 * t1840;
t1810 = t1996 * t1812 + t1992 * t1814;
t1809 = -t1988 * t1811 + t2033 * t1989;
t1808 = t1989 * t1811 + t2033 * t1988;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1974, -t1973, 0, -t1993 * t1978 + t1997 * t1979, 0, 0, 0, 0, 0, 0, -t1993 * t1913 + t1997 * t1932, -t1993 * t1902 + t1997 * t1928, -t1993 * t1897 + t1997 * t1914, -t1993 * t1881 + t1997 * t1889, 0, 0, 0, 0, 0, 0, -t1993 * t1851 + t1997 * t1866, -t1993 * t1855 + t1997 * t1868, -t1993 * t1839 + t1997 * t1861, -t1993 * t1825 + t1997 * t1833, 0, 0, 0, 0, 0, 0, -t1993 * t1838 + t1997 * t1860, -t1993 * t1852 + t1997 * t1867, -t1993 * t1856 + t1997 * t1869, -t1993 * t1816 + t1997 * t1823, 0, 0, 0, 0, 0, 0, -t1993 * t1820 + t1997 * t1829, -t1993 * t1822 + t1997 * t1830, -t1993 * t1818 + t1997 * t1828, -t1993 * t1809 + t1997 * t1810; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1973, t1974, 0, t1997 * t1978 + t1993 * t1979, 0, 0, 0, 0, 0, 0, t1997 * t1913 + t1993 * t1932, t1997 * t1902 + t1993 * t1928, t1997 * t1897 + t1993 * t1914, t1997 * t1881 + t1993 * t1889, 0, 0, 0, 0, 0, 0, t1997 * t1851 + t1993 * t1866, t1997 * t1855 + t1993 * t1868, t1997 * t1839 + t1993 * t1861, t1997 * t1825 + t1993 * t1833, 0, 0, 0, 0, 0, 0, t1997 * t1838 + t1993 * t1860, t1997 * t1852 + t1993 * t1867, t1997 * t1856 + t1993 * t1869, t1997 * t1816 + t1993 * t1823, 0, 0, 0, 0, 0, 0, t1997 * t1820 + t1993 * t1829, t1997 * t1822 + t1993 * t1830, t1997 * t1818 + t1993 * t1828, t1997 * t1809 + t1993 * t1810; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1912, t1901, t1896, t1880, 0, 0, 0, 0, 0, 0, t1849, t1853, t1837, t1824, 0, 0, 0, 0, 0, 0, t1836, t1850, t1854, t1815, 0, 0, 0, 0, 0, 0, t1819, t1821, t1817, t1808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1998, -qJDD(1), 0, t1979, 0, 0, 0, 0, 0, 0, t1932, t1928, t1914, t1889, 0, 0, 0, 0, 0, 0, t1866, t1868, t1861, t1833, 0, 0, 0, 0, 0, 0, t1860, t1867, t1869, t1823, 0, 0, 0, 0, 0, 0, t1829, t1830, t1828, t1810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1998, 0, t1978, 0, 0, 0, 0, 0, 0, t1913, t1902, t1897, t1881, 0, 0, 0, 0, 0, 0, t1851, t1855, t1839, t1825, 0, 0, 0, 0, 0, 0, t1838, t1852, t1856, t1816, 0, 0, 0, 0, 0, 0, t1820, t1822, t1818, t1809; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1912, t1901, t1896, t1880, 0, 0, 0, 0, 0, 0, t1849, t1853, t1837, t1824, 0, 0, 0, 0, 0, 0, t1836, t1850, t1854, t1815, 0, 0, 0, 0, 0, 0, t1819, t1821, t1817, t1808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1964, t1961, t1940, t1927, 0, 0, 0, 0, 0, 0, t1886, -t1893, t1879, t1846, 0, 0, 0, 0, 0, 0, t1878, -t1886, t1893, t1832, 0, 0, 0, 0, 0, 0, t1842, t1844, t1835, t1812; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1960, t1949, t1939, t1926, 0, 0, 0, 0, 0, 0, -t1904, -t1910, -t2074, -t1898, 0, 0, 0, 0, 0, 0, -t2074, t1906, t1908, -t1859, 0, 0, 0, 0, 0, 0, -t1863, -t1865, -t1848, -t1814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1941, t1938, t1965, t1944, 0, 0, 0, 0, 0, 0, -t1885, t1890, t1877, t1845, 0, 0, 0, 0, 0, 0, t1876, t1885, -t1890, t1831, 0, 0, 0, 0, 0, 0, t1841, t1843, t1834, t1811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1925, t2011, t1907, t1871, 0, 0, 0, 0, 0, 0, -t1905, -t1925, -t2011, t1857, 0, 0, 0, 0, 0, 0, t1873, t1875, t1888, t1840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1919, t2045, t1911, t1870, 0, 0, 0, 0, 0, 0, t1911, t1919, -t2045, -t1858, 0, 0, 0, 0, 0, 0, -t1862, -t1864, -t1847, -t1813; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1904, t1910, t2074, t1898, 0, 0, 0, 0, 0, 0, t2074, -t1906, -t1908, t1859, 0, 0, 0, 0, 0, 0, t1863, t1865, t1848, t1814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2074, -t1906, -t1908, t1859, 0, 0, 0, 0, 0, 0, t1863, t1865, t1848, t1814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1905, t1925, t2011, -t1857, 0, 0, 0, 0, 0, 0, -t1873, -t1875, -t1888, -t1840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1911, -t1919, t2045, t1858, 0, 0, 0, 0, 0, 0, t1862, t1864, t1847, t1813; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1894, t1883, t1874, t1827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2003, t2046, t2010, t1999; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1873, t1875, t1888, t1840;];
f_new_reg = t1;
