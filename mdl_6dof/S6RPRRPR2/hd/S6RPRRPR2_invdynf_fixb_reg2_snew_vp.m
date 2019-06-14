% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRRPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:08:15
% EndTime: 2019-05-05 22:08:22
% DurationCPUTime: 7.93s
% Computational Cost: add. (63547->306), mult. (128900->437), div. (0->0), fcn. (89846->12), ass. (0->245)
t2092 = sin(qJ(3));
t2122 = qJD(1) * t2092;
t2074 = qJD(3) * t2122;
t2096 = cos(qJ(3));
t2077 = t2096 * qJDD(1);
t2109 = t2077 - t2074 - qJDD(4);
t2091 = sin(qJ(4));
t2095 = cos(qJ(4));
t2046 = -t2095 * qJD(3) + t2091 * t2122;
t2047 = qJD(3) * t2091 + t2095 * t2122;
t2116 = t2046 * t2047;
t2013 = -t2109 - t2116;
t2084 = sin(pkin(11));
t2086 = cos(pkin(11));
t2024 = t2086 * t2046 + t2047 * t2084;
t2026 = -t2046 * t2084 + t2047 * t2086;
t2119 = t2024 * t2026;
t1982 = -t2109 - t2119;
t2121 = qJD(1) * t2096;
t2069 = -qJD(4) + t2121;
t2063 = -qJD(6) + t2069;
t2134 = qJD(6) - t2063;
t2133 = qJD(3) ^ 2;
t2090 = sin(qJ(6));
t2094 = cos(qJ(6));
t1990 = t2094 * t2024 + t2026 * t2090;
t2132 = t1990 ^ 2;
t1992 = -t2024 * t2090 + t2026 * t2094;
t2131 = t1992 ^ 2;
t2130 = t2024 ^ 2;
t2129 = t2026 ^ 2;
t2128 = t2046 ^ 2;
t2127 = t2047 ^ 2;
t2126 = t2063 ^ 2;
t2125 = t2069 ^ 2;
t2124 = -2 * qJD(5);
t2123 = -g(3) + qJDD(2);
t2120 = t1990 * t1992;
t2118 = t2024 * t2069;
t2117 = t2026 * t2069;
t2115 = t2046 * t2069;
t2114 = qJD(4) + t2069;
t2093 = sin(qJ(1));
t2097 = cos(qJ(1));
t2062 = -g(1) * t2097 - g(2) * t2093;
t2098 = qJD(1) ^ 2;
t2048 = -pkin(1) * t2098 + t2062;
t2085 = sin(pkin(10));
t2087 = cos(pkin(10));
t2061 = t2093 * g(1) - g(2) * t2097;
t2101 = qJDD(1) * pkin(1) + t2061;
t2027 = -t2085 * t2048 + t2087 * t2101;
t2015 = -qJDD(1) * pkin(2) - t2098 * pkin(7) - t2027;
t2110 = qJD(3) * t2121;
t2111 = t2092 * qJDD(1);
t2051 = t2110 + t2111;
t2052 = t2077 - 0.2e1 * t2074;
t1984 = (-t2051 - t2110) * pkin(8) - t2052 * pkin(3) + t2015;
t2028 = t2087 * t2048 + t2085 * t2101;
t2016 = -pkin(2) * t2098 + qJDD(1) * pkin(7) + t2028;
t2003 = t2096 * t2016 + t2092 * t2123;
t2049 = (-pkin(3) * t2096 - pkin(8) * t2092) * qJD(1);
t1996 = -pkin(3) * t2133 + qJDD(3) * pkin(8) + t2049 * t2121 + t2003;
t1956 = t2091 * t1984 + t2095 * t1996;
t2113 = t2063 + qJD(6);
t2080 = t2092 ^ 2;
t2081 = t2096 ^ 2;
t2112 = t2080 + t2081;
t1955 = t2095 * t1984 - t2091 * t1996;
t2102 = -t2091 * qJDD(3) - t2095 * t2051;
t2021 = -qJD(4) * t2046 - t2102;
t1937 = (-t2021 + t2115) * qJ(5) + t2013 * pkin(4) + t1955;
t2036 = -pkin(4) * t2069 - qJ(5) * t2047;
t2105 = -t2095 * qJDD(3) + t2091 * t2051;
t2100 = -qJD(4) * t2047 - t2105;
t1941 = -pkin(4) * t2128 + qJ(5) * t2100 + t2069 * t2036 + t1956;
t1907 = t2084 * t1937 + t2086 * t1941 + t2024 * t2124;
t1985 = t2086 * t2021 + t2084 * t2100;
t2107 = t2021 * t2084 - t2086 * t2100;
t2108 = -t2090 * t1985 - t2094 * t2107;
t1968 = -t1985 + t2118;
t2053 = -qJDD(1) * t2085 - t2087 * t2098;
t2054 = qJDD(1) * t2087 - t2085 * t2098;
t2106 = t2097 * t2053 - t2054 * t2093;
t2104 = -qJDD(6) + t2109;
t1906 = t2086 * t1937 - t2084 * t1941 + t2026 * t2124;
t2103 = t2053 * t2093 + t2054 * t2097;
t2099 = -t2094 * t1985 + t2090 * t2107;
t2073 = t2096 * t2123;
t1995 = -t2073 - qJDD(3) * pkin(3) - t2133 * pkin(8) + (qJD(1) * t2049 + t2016) * t2092;
t1953 = -t2100 * pkin(4) - t2128 * qJ(5) + t2047 * t2036 + qJDD(5) + t1995;
t2068 = t2096 * t2098 * t2092;
t2066 = -t2081 * t2098 - t2133;
t2065 = -t2080 * t2098 - t2133;
t2060 = -qJDD(3) + t2068;
t2059 = qJDD(3) + t2068;
t2058 = t2112 * t2098;
t2057 = -qJDD(1) * t2093 - t2097 * t2098;
t2056 = qJDD(1) * t2097 - t2093 * t2098;
t2055 = t2112 * qJDD(1);
t2050 = 0.2e1 * t2110 + t2111;
t2035 = t2060 * t2096 - t2065 * t2092;
t2034 = -t2059 * t2092 + t2066 * t2096;
t2033 = t2060 * t2092 + t2065 * t2096;
t2032 = t2059 * t2096 + t2066 * t2092;
t2031 = -t2125 - t2127;
t2030 = t2055 * t2087 - t2058 * t2085;
t2029 = t2055 * t2085 + t2058 * t2087;
t2023 = -t2125 - t2128;
t2012 = t2109 - t2116;
t2011 = -t2127 - t2128;
t2008 = -pkin(5) * t2069 - pkin(9) * t2026;
t2007 = t2035 * t2087 + t2050 * t2085;
t2006 = t2034 * t2087 - t2052 * t2085;
t2005 = t2035 * t2085 - t2050 * t2087;
t2004 = t2034 * t2085 + t2052 * t2087;
t2002 = -t2016 * t2092 + t2073;
t2001 = -t2125 - t2129;
t2000 = t2046 * t2114 + t2102;
t1999 = t2021 + t2115;
t1998 = -t2047 * t2114 - t2105;
t1997 = (qJD(4) - t2069) * t2047 + t2105;
t1994 = -t2027 * t2085 + t2028 * t2087;
t1993 = t2027 * t2087 + t2028 * t2085;
t1988 = t2012 * t2095 - t2031 * t2091;
t1987 = t2012 * t2091 + t2031 * t2095;
t1986 = -t2125 - t2130;
t1981 = t2109 - t2119;
t1980 = -t2013 * t2091 + t2023 * t2095;
t1979 = t2013 * t2095 + t2023 * t2091;
t1974 = -t2126 - t2131;
t1973 = -t2129 - t2130;
t1972 = -t2002 * t2092 + t2003 * t2096;
t1971 = t2002 * t2096 + t2003 * t2092;
t1970 = t1998 * t2095 - t2000 * t2091;
t1969 = t1998 * t2091 + t2000 * t2095;
t1967 = t1985 + t2118;
t1966 = -t2107 - t2117;
t1965 = t2107 - t2117;
t1964 = t1988 * t2096 + t1999 * t2092;
t1963 = t1981 * t2086 - t2001 * t2084;
t1962 = t1988 * t2092 - t1999 * t2096;
t1961 = t1981 * t2084 + t2001 * t2086;
t1960 = t1980 * t2096 + t1997 * t2092;
t1959 = t1980 * t2092 - t1997 * t2096;
t1958 = t1972 * t2087 + t2015 * t2085;
t1957 = t1972 * t2085 - t2015 * t2087;
t1954 = -t2126 - t2132;
t1952 = -t2104 - t2120;
t1951 = t2104 - t2120;
t1950 = -t1982 * t2084 + t1986 * t2086;
t1949 = t1982 * t2086 + t1986 * t2084;
t1948 = t1970 * t2096 + t2011 * t2092;
t1947 = t1970 * t2092 - t2011 * t2096;
t1946 = t1964 * t2087 + t1987 * t2085;
t1945 = t1964 * t2085 - t1987 * t2087;
t1944 = -t2131 - t2132;
t1943 = t1960 * t2087 + t1979 * t2085;
t1942 = t1960 * t2085 - t1979 * t2087;
t1939 = t1966 * t2086 - t1968 * t2084;
t1938 = t1966 * t2084 + t1968 * t2086;
t1934 = t1951 * t2094 - t1974 * t2090;
t1933 = t1951 * t2090 + t1974 * t2094;
t1932 = -t1961 * t2091 + t1963 * t2095;
t1931 = t1961 * t2095 + t1963 * t2091;
t1930 = t1948 * t2087 + t1969 * t2085;
t1929 = t1948 * t2085 - t1969 * t2087;
t1928 = -t1955 * t2091 + t1956 * t2095;
t1927 = t1955 * t2095 + t1956 * t2091;
t1926 = t1990 * t2113 + t2099;
t1925 = -t1990 * t2134 - t2099;
t1924 = -t1992 * t2113 + t2108;
t1923 = t1992 * t2134 - t2108;
t1922 = -t1952 * t2090 + t1954 * t2094;
t1921 = t1952 * t2094 + t1954 * t2090;
t1920 = -t1949 * t2091 + t1950 * t2095;
t1919 = t1949 * t2095 + t1950 * t2091;
t1918 = pkin(5) * t2107 - pkin(9) * t2130 + t2026 * t2008 + t1953;
t1917 = t1928 * t2096 + t1995 * t2092;
t1916 = t1928 * t2092 - t1995 * t2096;
t1915 = t1932 * t2096 + t1967 * t2092;
t1914 = t1932 * t2092 - t1967 * t2096;
t1913 = t1920 * t2096 + t1965 * t2092;
t1912 = t1920 * t2092 - t1965 * t2096;
t1911 = -t1938 * t2091 + t1939 * t2095;
t1910 = t1938 * t2095 + t1939 * t2091;
t1909 = -t1933 * t2084 + t1934 * t2086;
t1908 = t1933 * t2086 + t1934 * t2084;
t1905 = t1911 * t2096 + t1973 * t2092;
t1904 = t1911 * t2092 - t1973 * t2096;
t1903 = t1924 * t2094 - t1926 * t2090;
t1902 = t1924 * t2090 + t1926 * t2094;
t1901 = -t1921 * t2084 + t1922 * t2086;
t1900 = t1921 * t2086 + t1922 * t2084;
t1899 = t1915 * t2087 + t1931 * t2085;
t1898 = t1915 * t2085 - t1931 * t2087;
t1897 = -pkin(5) * t2130 - pkin(9) * t2107 + t2069 * t2008 + t1907;
t1896 = t1917 * t2087 + t1927 * t2085;
t1895 = t1917 * t2085 - t1927 * t2087;
t1894 = pkin(5) * t1982 + t1968 * pkin(9) + t1906;
t1893 = t1913 * t2087 + t1919 * t2085;
t1892 = t1913 * t2085 - t1919 * t2087;
t1891 = -t1908 * t2091 + t1909 * t2095;
t1890 = t1908 * t2095 + t1909 * t2091;
t1889 = -t1906 * t2084 + t1907 * t2086;
t1888 = t1906 * t2086 + t1907 * t2084;
t1887 = t1905 * t2087 + t1910 * t2085;
t1886 = t1905 * t2085 - t1910 * t2087;
t1885 = -t1902 * t2084 + t1903 * t2086;
t1884 = t1902 * t2086 + t1903 * t2084;
t1883 = -t1900 * t2091 + t1901 * t2095;
t1882 = t1900 * t2095 + t1901 * t2091;
t1881 = t1891 * t2096 + t1925 * t2092;
t1880 = t1891 * t2092 - t1925 * t2096;
t1879 = t1894 * t2090 + t1897 * t2094;
t1878 = t1894 * t2094 - t1897 * t2090;
t1877 = t1883 * t2096 + t1923 * t2092;
t1876 = t1883 * t2092 - t1923 * t2096;
t1875 = -t1888 * t2091 + t1889 * t2095;
t1874 = t1888 * t2095 + t1889 * t2091;
t1873 = t1875 * t2096 + t1953 * t2092;
t1872 = t1875 * t2092 - t1953 * t2096;
t1871 = t1881 * t2087 + t1890 * t2085;
t1870 = t1881 * t2085 - t1890 * t2087;
t1869 = -t1884 * t2091 + t1885 * t2095;
t1868 = t1884 * t2095 + t1885 * t2091;
t1867 = t1869 * t2096 + t1944 * t2092;
t1866 = t1869 * t2092 - t1944 * t2096;
t1865 = -t1878 * t2090 + t1879 * t2094;
t1864 = t1878 * t2094 + t1879 * t2090;
t1863 = t1877 * t2087 + t1882 * t2085;
t1862 = t1877 * t2085 - t1882 * t2087;
t1861 = t1873 * t2087 + t1874 * t2085;
t1860 = t1873 * t2085 - t1874 * t2087;
t1859 = t1867 * t2087 + t1868 * t2085;
t1858 = t1867 * t2085 - t1868 * t2087;
t1857 = -t1864 * t2084 + t1865 * t2086;
t1856 = t1864 * t2086 + t1865 * t2084;
t1855 = -t1856 * t2091 + t1857 * t2095;
t1854 = t1856 * t2095 + t1857 * t2091;
t1853 = t1855 * t2096 + t1918 * t2092;
t1852 = t1855 * t2092 - t1918 * t2096;
t1851 = t1853 * t2087 + t1854 * t2085;
t1850 = t1853 * t2085 - t1854 * t2087;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2057, -t2056, 0, -t2061 * t2093 + t2062 * t2097, 0, 0, 0, 0, 0, 0, t2106, -t2103, 0, -t1993 * t2093 + t1994 * t2097, 0, 0, 0, 0, 0, 0, -t2004 * t2093 + t2006 * t2097, -t2005 * t2093 + t2007 * t2097, -t2029 * t2093 + t2030 * t2097, -t1957 * t2093 + t1958 * t2097, 0, 0, 0, 0, 0, 0, -t1942 * t2093 + t1943 * t2097, -t1945 * t2093 + t1946 * t2097, -t1929 * t2093 + t1930 * t2097, -t1895 * t2093 + t1896 * t2097, 0, 0, 0, 0, 0, 0, -t1892 * t2093 + t1893 * t2097, -t1898 * t2093 + t1899 * t2097, -t1886 * t2093 + t1887 * t2097, -t1860 * t2093 + t1861 * t2097, 0, 0, 0, 0, 0, 0, -t1862 * t2093 + t1863 * t2097, -t1870 * t2093 + t1871 * t2097, -t1858 * t2093 + t1859 * t2097, -t1850 * t2093 + t1851 * t2097; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2056, t2057, 0, t2061 * t2097 + t2062 * t2093, 0, 0, 0, 0, 0, 0, t2103, t2106, 0, t1993 * t2097 + t1994 * t2093, 0, 0, 0, 0, 0, 0, t2004 * t2097 + t2006 * t2093, t2005 * t2097 + t2007 * t2093, t2029 * t2097 + t2030 * t2093, t1957 * t2097 + t1958 * t2093, 0, 0, 0, 0, 0, 0, t1942 * t2097 + t1943 * t2093, t1945 * t2097 + t1946 * t2093, t1929 * t2097 + t1930 * t2093, t1895 * t2097 + t1896 * t2093, 0, 0, 0, 0, 0, 0, t1892 * t2097 + t1893 * t2093, t1898 * t2097 + t1899 * t2093, t1886 * t2097 + t1887 * t2093, t1860 * t2097 + t1861 * t2093, 0, 0, 0, 0, 0, 0, t1862 * t2097 + t1863 * t2093, t1870 * t2097 + t1871 * t2093, t1858 * t2097 + t1859 * t2093, t1850 * t2097 + t1851 * t2093; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2123, 0, 0, 0, 0, 0, 0, t2032, t2033, 0, t1971, 0, 0, 0, 0, 0, 0, t1959, t1962, t1947, t1916, 0, 0, 0, 0, 0, 0, t1912, t1914, t1904, t1872, 0, 0, 0, 0, 0, 0, t1876, t1880, t1866, t1852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2098, -qJDD(1), 0, t2062, 0, 0, 0, 0, 0, 0, t2053, -t2054, 0, t1994, 0, 0, 0, 0, 0, 0, t2006, t2007, t2030, t1958, 0, 0, 0, 0, 0, 0, t1943, t1946, t1930, t1896, 0, 0, 0, 0, 0, 0, t1893, t1899, t1887, t1861, 0, 0, 0, 0, 0, 0, t1863, t1871, t1859, t1851; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2098, 0, t2061, 0, 0, 0, 0, 0, 0, t2054, t2053, 0, t1993, 0, 0, 0, 0, 0, 0, t2004, t2005, t2029, t1957, 0, 0, 0, 0, 0, 0, t1942, t1945, t1929, t1895, 0, 0, 0, 0, 0, 0, t1892, t1898, t1886, t1860, 0, 0, 0, 0, 0, 0, t1862, t1870, t1858, t1850; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2123, 0, 0, 0, 0, 0, 0, t2032, t2033, 0, t1971, 0, 0, 0, 0, 0, 0, t1959, t1962, t1947, t1916, 0, 0, 0, 0, 0, 0, t1912, t1914, t1904, t1872, 0, 0, 0, 0, 0, 0, t1876, t1880, t1866, t1852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2098, -qJDD(1), 0, t2028, 0, 0, 0, 0, 0, 0, t2034, t2035, t2055, t1972, 0, 0, 0, 0, 0, 0, t1960, t1964, t1948, t1917, 0, 0, 0, 0, 0, 0, t1913, t1915, t1905, t1873, 0, 0, 0, 0, 0, 0, t1877, t1881, t1867, t1853; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2098, 0, t2027, 0, 0, 0, 0, 0, 0, t2052, -t2050, t2058, -t2015, 0, 0, 0, 0, 0, 0, -t1979, -t1987, -t1969, -t1927, 0, 0, 0, 0, 0, 0, -t1919, -t1931, -t1910, -t1874, 0, 0, 0, 0, 0, 0, -t1882, -t1890, -t1868, -t1854; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2123, 0, 0, 0, 0, 0, 0, t2032, t2033, 0, t1971, 0, 0, 0, 0, 0, 0, t1959, t1962, t1947, t1916, 0, 0, 0, 0, 0, 0, t1912, t1914, t1904, t1872, 0, 0, 0, 0, 0, 0, t1876, t1880, t1866, t1852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2066, t2060, t2077, t2003, 0, 0, 0, 0, 0, 0, t1980, t1988, t1970, t1928, 0, 0, 0, 0, 0, 0, t1920, t1932, t1911, t1875, 0, 0, 0, 0, 0, 0, t1883, t1891, t1869, t1855; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2059, t2065, -t2111, t2002, 0, 0, 0, 0, 0, 0, -t1997, -t1999, -t2011, -t1995, 0, 0, 0, 0, 0, 0, -t1965, -t1967, -t1973, -t1953, 0, 0, 0, 0, 0, 0, -t1923, -t1925, -t1944, -t1918; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2052, t2050, -t2058, t2015, 0, 0, 0, 0, 0, 0, t1979, t1987, t1969, t1927, 0, 0, 0, 0, 0, 0, t1919, t1931, t1910, t1874, 0, 0, 0, 0, 0, 0, t1882, t1890, t1868, t1854; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2023, t2012, t1998, t1956, 0, 0, 0, 0, 0, 0, t1950, t1963, t1939, t1889, 0, 0, 0, 0, 0, 0, t1901, t1909, t1885, t1857; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2013, t2031, t2000, t1955, 0, 0, 0, 0, 0, 0, t1949, t1961, t1938, t1888, 0, 0, 0, 0, 0, 0, t1900, t1908, t1884, t1856; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1997, t1999, t2011, t1995, 0, 0, 0, 0, 0, 0, t1965, t1967, t1973, t1953, 0, 0, 0, 0, 0, 0, t1923, t1925, t1944, t1918; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1986, t1981, t1966, t1907, 0, 0, 0, 0, 0, 0, t1922, t1934, t1903, t1865; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1982, t2001, t1968, t1906, 0, 0, 0, 0, 0, 0, t1921, t1933, t1902, t1864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1965, t1967, t1973, t1953, 0, 0, 0, 0, 0, 0, t1923, t1925, t1944, t1918; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1954, t1951, t1924, t1879; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1952, t1974, t1926, t1878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1923, t1925, t1944, t1918;];
f_new_reg  = t1;