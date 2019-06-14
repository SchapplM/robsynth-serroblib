% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPRRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:34:41
% EndTime: 2019-05-06 17:34:47
% DurationCPUTime: 6.50s
% Computational Cost: add. (46164->289), mult. (106215->390), div. (0->0), fcn. (78200->10), ass. (0->225)
t2068 = sin(pkin(10));
t2069 = cos(pkin(10));
t2077 = cos(qJ(2));
t2110 = qJD(1) * t2077;
t2073 = sin(qJ(2));
t2111 = qJD(1) * t2073;
t2034 = t2068 * t2111 - t2069 * t2110;
t2033 = qJD(4) + t2034;
t2028 = qJD(5) + t2033;
t2123 = qJD(5) + t2028;
t2096 = qJD(2) * t2110;
t2098 = t2073 * qJDD(1);
t2044 = t2096 + t2098;
t2063 = t2077 * qJDD(1);
t2097 = qJD(2) * t2111;
t2045 = t2063 - t2097;
t2016 = t2069 * t2044 + t2068 * t2045;
t2032 = qJD(2) * t2034;
t2009 = t2032 - t2016;
t2066 = t2077 ^ 2;
t2079 = qJD(1) ^ 2;
t2087 = qJD(2) * pkin(2) - qJ(3) * t2111;
t2074 = sin(qJ(1));
t2078 = cos(qJ(1));
t2053 = t2074 * g(1) - t2078 * g(2);
t2088 = qJDD(1) * pkin(1) + t2053;
t2004 = t2045 * pkin(2) + (qJ(3) * t2066 + pkin(7)) * t2079 - t2087 * t2111 - qJDD(3) + t2088;
t2122 = qJD(2) ^ 2;
t2036 = t2068 * t2110 + t2069 * t2111;
t2072 = sin(qJ(4));
t2076 = cos(qJ(4));
t2019 = -qJD(2) * t2076 + t2036 * t2072;
t2020 = qJD(2) * t2072 + t2036 * t2076;
t2071 = sin(qJ(5));
t2075 = cos(qJ(5));
t1998 = t2075 * t2019 + t2020 * t2071;
t2121 = t1998 ^ 2;
t2000 = -t2019 * t2071 + t2020 * t2075;
t2120 = t2000 ^ 2;
t2119 = t2019 ^ 2;
t2118 = t2020 ^ 2;
t2117 = t2028 ^ 2;
t2116 = t2033 ^ 2;
t2115 = t2034 ^ 2;
t2114 = t2036 ^ 2;
t2113 = -2 * qJD(3);
t2112 = -2 * qJD(6);
t2109 = qJD(2) * t2036;
t2108 = t1998 * t2000;
t2107 = t2019 * t2020;
t2106 = t2033 * t2019;
t2105 = t2034 * t2036;
t2104 = t2066 * t2079;
t2103 = t2073 * t2079;
t2102 = qJD(4) - t2033;
t2101 = qJD(5) - t2028;
t2054 = -g(1) * t2078 - g(2) * t2074;
t2084 = -pkin(1) * t2079 + qJDD(1) * pkin(7) + t2054;
t2026 = -t2073 * g(3) + t2077 * t2084;
t2001 = -pkin(2) * t2104 + t2045 * qJ(3) - qJD(2) * t2087 + t2026;
t2082 = t2073 * t2084;
t2080 = -t2082 - t2044 * qJ(3) + qJDD(2) * pkin(2) + (qJ(3) * qJD(1) * qJD(2) + pkin(2) * t2103 - g(3)) * t2077;
t1969 = t2069 * t2001 + t2034 * t2113 + t2068 * t2080;
t2012 = pkin(3) * t2034 - pkin(8) * t2036;
t1955 = -pkin(3) * t2122 + qJDD(2) * pkin(8) - t2012 * t2034 + t1969;
t2100 = t2068 * t2044 - t2069 * t2045;
t2006 = t2100 + t2109;
t1963 = t2006 * pkin(3) + pkin(8) * t2009 - t2004;
t1930 = -t2072 * t1955 + t2076 * t1963;
t2094 = -qJDD(4) - t2100;
t1980 = -t2094 - t2107;
t2089 = -t2072 * qJDD(2) - t2076 * t2016;
t1990 = -qJD(4) * t2019 - t2089;
t1922 = (-t1990 - t2106) * pkin(9) + t1980 * pkin(4) + t1930;
t1931 = t2076 * t1955 + t2072 * t1963;
t2010 = pkin(4) * t2033 - pkin(9) * t2020;
t2091 = -t2076 * qJDD(2) + t2072 * t2016;
t2086 = -qJD(4) * t2020 - t2091;
t1926 = -pkin(4) * t2119 + pkin(9) * t2086 - t2033 * t2010 + t1931;
t1900 = t2071 * t1922 + t2075 * t1926;
t2065 = t2073 ^ 2;
t2099 = t2065 + t2066;
t1899 = t2075 * t1922 - t2071 * t1926;
t2093 = t2071 * t1990 - t2075 * t2086;
t2092 = t2068 * t2001 - t2069 * t2080;
t2090 = -qJDD(5) + t2094;
t2085 = -qJD(5) * t2000 - t2093;
t1966 = -t2090 - t2108;
t1954 = -qJDD(2) * pkin(3) - t2122 * pkin(8) + ((2 * qJD(3)) + t2012) * t2036 + t2092;
t2081 = -t2075 * t1990 - t2071 * t2086;
t1929 = -t2086 * pkin(4) - t2119 * pkin(9) + t2020 * t2010 + t1954;
t1943 = t1998 * t2101 + t2081;
t2060 = t2077 * t2103;
t2059 = -t2104 - t2122;
t2058 = -t2065 * t2079 - t2122;
t2052 = -qJDD(2) + t2060;
t2051 = qJDD(2) + t2060;
t2050 = t2099 * t2079;
t2049 = -qJDD(1) * t2074 - t2078 * t2079;
t2048 = qJDD(1) * t2078 - t2074 * t2079;
t2047 = t2099 * qJDD(1);
t2046 = t2063 - 0.2e1 * t2097;
t2043 = 0.2e1 * t2096 + t2098;
t2041 = t2079 * pkin(7) + t2088;
t2027 = -t2114 - t2122;
t2025 = -t2077 * g(3) - t2082;
t2024 = t2052 * t2077 - t2058 * t2073;
t2023 = -t2051 * t2073 + t2059 * t2077;
t2022 = t2052 * t2073 + t2058 * t2077;
t2021 = t2051 * t2077 + t2059 * t2073;
t2015 = -qJDD(2) - t2105;
t2014 = qJDD(2) - t2105;
t2013 = -t2115 - t2122;
t2008 = -t2032 - t2016;
t2007 = -t2100 + t2109;
t2005 = -t2114 - t2115;
t2003 = -t2025 * t2073 + t2026 * t2077;
t2002 = t2025 * t2077 + t2026 * t2073;
t1993 = -t2116 - t2118;
t1992 = t2015 * t2069 - t2027 * t2068;
t1991 = t2015 * t2068 + t2027 * t2069;
t1988 = -t2116 - t2119;
t1985 = -t2118 - t2119;
t1984 = t2013 * t2069 - t2014 * t2068;
t1983 = t2013 * t2068 + t2014 * t2069;
t1982 = pkin(5) * t2028 - qJ(6) * t2000;
t1981 = t2094 - t2107;
t1979 = -t2117 - t2120;
t1978 = t2007 * t2069 - t2008 * t2068;
t1977 = t2007 * t2068 + t2008 * t2069;
t1976 = t2019 * t2102 + t2089;
t1975 = t1990 - t2106;
t1974 = -t2020 * t2102 - t2091;
t1973 = (qJD(4) + t2033) * t2020 + t2091;
t1972 = -t1991 * t2073 + t1992 * t2077;
t1971 = t1991 * t2077 + t1992 * t2073;
t1970 = -t2117 - t2121;
t1968 = t2036 * t2113 - t2092;
t1967 = t2090 - t2108;
t1965 = t1981 * t2076 - t1993 * t2072;
t1964 = t1981 * t2072 + t1993 * t2076;
t1960 = -t1980 * t2072 + t1988 * t2076;
t1959 = t1980 * t2076 + t1988 * t2072;
t1958 = -t1983 * t2073 + t1984 * t2077;
t1957 = t1983 * t2077 + t1984 * t2073;
t1956 = -t2120 - t2121;
t1951 = -t1977 * t2073 + t1978 * t2077;
t1950 = t1977 * t2077 + t1978 * t2073;
t1949 = t1974 * t2076 - t1976 * t2072;
t1948 = t1974 * t2072 + t1976 * t2076;
t1947 = t1967 * t2075 - t1979 * t2071;
t1946 = t1967 * t2071 + t1979 * t2075;
t1945 = t1965 * t2069 + t1975 * t2068;
t1944 = t1965 * t2068 - t1975 * t2069;
t1942 = -t1998 * t2123 - t2081;
t1941 = -t2000 * t2101 - t2093;
t1940 = t2000 * t2123 + t2093;
t1939 = t1960 * t2069 + t1973 * t2068;
t1938 = t1960 * t2068 - t1973 * t2069;
t1937 = -t1968 * t2068 + t1969 * t2069;
t1936 = t1968 * t2069 + t1969 * t2068;
t1935 = -t1966 * t2071 + t1970 * t2075;
t1934 = t1966 * t2075 + t1970 * t2071;
t1933 = t1949 * t2069 + t1985 * t2068;
t1932 = t1949 * t2068 - t1985 * t2069;
t1928 = -t1946 * t2072 + t1947 * t2076;
t1927 = t1946 * t2076 + t1947 * t2072;
t1924 = -t1944 * t2073 + t1945 * t2077;
t1923 = t1944 * t2077 + t1945 * t2073;
t1919 = t1941 * t2075 - t1943 * t2071;
t1918 = t1941 * t2071 + t1943 * t2075;
t1917 = -t1938 * t2073 + t1939 * t2077;
t1916 = t1938 * t2077 + t1939 * t2073;
t1915 = -t1936 * t2073 + t1937 * t2077;
t1914 = t1936 * t2077 + t1937 * t2073;
t1913 = -t1934 * t2072 + t1935 * t2076;
t1912 = t1934 * t2076 + t1935 * t2072;
t1911 = -t1932 * t2073 + t1933 * t2077;
t1910 = t1932 * t2077 + t1933 * t2073;
t1909 = -t1930 * t2072 + t1931 * t2076;
t1908 = t1930 * t2076 + t1931 * t2072;
t1907 = -pkin(5) * t2085 - qJ(6) * t2121 + t2000 * t1982 + qJDD(6) + t1929;
t1906 = t1928 * t2069 + t1942 * t2068;
t1905 = t1928 * t2068 - t1942 * t2069;
t1904 = t1913 * t2069 + t1940 * t2068;
t1903 = t1913 * t2068 - t1940 * t2069;
t1902 = t1909 * t2069 + t1954 * t2068;
t1901 = t1909 * t2068 - t1954 * t2069;
t1898 = -t1918 * t2072 + t1919 * t2076;
t1897 = t1918 * t2076 + t1919 * t2072;
t1896 = t1898 * t2069 + t1956 * t2068;
t1895 = t1898 * t2068 - t1956 * t2069;
t1894 = -t2028 * t1982 + t2085 * qJ(6) + (-pkin(5) * t1998 + t2112) * t1998 + t1900;
t1893 = pkin(5) * t1966 + qJ(6) * t1943 + t2000 * t2112 + t1899;
t1892 = -t1905 * t2073 + t1906 * t2077;
t1891 = t1905 * t2077 + t1906 * t2073;
t1890 = -t1903 * t2073 + t1904 * t2077;
t1889 = t1903 * t2077 + t1904 * t2073;
t1888 = -t1901 * t2073 + t1902 * t2077;
t1887 = t1901 * t2077 + t1902 * t2073;
t1886 = t1892 * t2078 + t1927 * t2074;
t1885 = t1892 * t2074 - t1927 * t2078;
t1884 = -t1899 * t2071 + t1900 * t2075;
t1883 = t1899 * t2075 + t1900 * t2071;
t1882 = t1890 * t2078 + t1912 * t2074;
t1881 = t1890 * t2074 - t1912 * t2078;
t1880 = -t1895 * t2073 + t1896 * t2077;
t1879 = t1895 * t2077 + t1896 * t2073;
t1878 = -t1893 * t2071 + t1894 * t2075;
t1877 = t1893 * t2075 + t1894 * t2071;
t1876 = t1880 * t2078 + t1897 * t2074;
t1875 = t1880 * t2074 - t1897 * t2078;
t1874 = -t1883 * t2072 + t1884 * t2076;
t1873 = t1883 * t2076 + t1884 * t2072;
t1872 = t1874 * t2069 + t1929 * t2068;
t1871 = t1874 * t2068 - t1929 * t2069;
t1870 = -t1877 * t2072 + t1878 * t2076;
t1869 = t1877 * t2076 + t1878 * t2072;
t1868 = t1870 * t2069 + t1907 * t2068;
t1867 = t1870 * t2068 - t1907 * t2069;
t1866 = -t1871 * t2073 + t1872 * t2077;
t1865 = t1871 * t2077 + t1872 * t2073;
t1864 = -t1867 * t2073 + t1868 * t2077;
t1863 = t1867 * t2077 + t1868 * t2073;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2049, -t2048, 0, -t2053 * t2074 + t2054 * t2078, 0, 0, 0, 0, 0, 0, t2023 * t2078 - t2046 * t2074, t2024 * t2078 + t2043 * t2074, t2047 * t2078 - t2050 * t2074, t2003 * t2078 - t2041 * t2074, 0, 0, 0, 0, 0, 0, t1958 * t2078 + t2006 * t2074, t1972 * t2078 - t2009 * t2074, t1951 * t2078 + t2005 * t2074, t1915 * t2078 - t2004 * t2074, 0, 0, 0, 0, 0, 0, t1917 * t2078 + t1959 * t2074, t1924 * t2078 + t1964 * t2074, t1911 * t2078 + t1948 * t2074, t1888 * t2078 + t1908 * t2074, 0, 0, 0, 0, 0, 0, t1882, t1886, t1876, t1866 * t2078 + t1873 * t2074, 0, 0, 0, 0, 0, 0, t1882, t1886, t1876, t1864 * t2078 + t1869 * t2074; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2048, t2049, 0, t2053 * t2078 + t2054 * t2074, 0, 0, 0, 0, 0, 0, t2023 * t2074 + t2046 * t2078, t2024 * t2074 - t2043 * t2078, t2047 * t2074 + t2050 * t2078, t2003 * t2074 + t2041 * t2078, 0, 0, 0, 0, 0, 0, t1958 * t2074 - t2006 * t2078, t1972 * t2074 + t2009 * t2078, t1951 * t2074 - t2005 * t2078, t1915 * t2074 + t2004 * t2078, 0, 0, 0, 0, 0, 0, t1917 * t2074 - t1959 * t2078, t1924 * t2074 - t1964 * t2078, t1911 * t2074 - t1948 * t2078, t1888 * t2074 - t1908 * t2078, 0, 0, 0, 0, 0, 0, t1881, t1885, t1875, t1866 * t2074 - t1873 * t2078, 0, 0, 0, 0, 0, 0, t1881, t1885, t1875, t1864 * t2074 - t1869 * t2078; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2021, t2022, 0, t2002, 0, 0, 0, 0, 0, 0, t1957, t1971, t1950, t1914, 0, 0, 0, 0, 0, 0, t1916, t1923, t1910, t1887, 0, 0, 0, 0, 0, 0, t1889, t1891, t1879, t1865, 0, 0, 0, 0, 0, 0, t1889, t1891, t1879, t1863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2079, -qJDD(1), 0, t2054, 0, 0, 0, 0, 0, 0, t2023, t2024, t2047, t2003, 0, 0, 0, 0, 0, 0, t1958, t1972, t1951, t1915, 0, 0, 0, 0, 0, 0, t1917, t1924, t1911, t1888, 0, 0, 0, 0, 0, 0, t1890, t1892, t1880, t1866, 0, 0, 0, 0, 0, 0, t1890, t1892, t1880, t1864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2079, 0, t2053, 0, 0, 0, 0, 0, 0, t2046, -t2043, t2050, t2041, 0, 0, 0, 0, 0, 0, -t2006, t2009, -t2005, t2004, 0, 0, 0, 0, 0, 0, -t1959, -t1964, -t1948, -t1908, 0, 0, 0, 0, 0, 0, -t1912, -t1927, -t1897, -t1873, 0, 0, 0, 0, 0, 0, -t1912, -t1927, -t1897, -t1869; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2021, t2022, 0, t2002, 0, 0, 0, 0, 0, 0, t1957, t1971, t1950, t1914, 0, 0, 0, 0, 0, 0, t1916, t1923, t1910, t1887, 0, 0, 0, 0, 0, 0, t1889, t1891, t1879, t1865, 0, 0, 0, 0, 0, 0, t1889, t1891, t1879, t1863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2059, t2052, t2063, t2026, 0, 0, 0, 0, 0, 0, t1984, t1992, t1978, t1937, 0, 0, 0, 0, 0, 0, t1939, t1945, t1933, t1902, 0, 0, 0, 0, 0, 0, t1904, t1906, t1896, t1872, 0, 0, 0, 0, 0, 0, t1904, t1906, t1896, t1868; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2051, t2058, -t2098, t2025, 0, 0, 0, 0, 0, 0, t1983, t1991, t1977, t1936, 0, 0, 0, 0, 0, 0, t1938, t1944, t1932, t1901, 0, 0, 0, 0, 0, 0, t1903, t1905, t1895, t1871, 0, 0, 0, 0, 0, 0, t1903, t1905, t1895, t1867; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2046, t2043, -t2050, -t2041, 0, 0, 0, 0, 0, 0, t2006, -t2009, t2005, -t2004, 0, 0, 0, 0, 0, 0, t1959, t1964, t1948, t1908, 0, 0, 0, 0, 0, 0, t1912, t1927, t1897, t1873, 0, 0, 0, 0, 0, 0, t1912, t1927, t1897, t1869; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2013, t2015, t2007, t1969, 0, 0, 0, 0, 0, 0, t1960, t1965, t1949, t1909, 0, 0, 0, 0, 0, 0, t1913, t1928, t1898, t1874, 0, 0, 0, 0, 0, 0, t1913, t1928, t1898, t1870; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2014, t2027, t2008, t1968, 0, 0, 0, 0, 0, 0, -t1973, -t1975, -t1985, -t1954, 0, 0, 0, 0, 0, 0, -t1940, -t1942, -t1956, -t1929, 0, 0, 0, 0, 0, 0, -t1940, -t1942, -t1956, -t1907; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2006, -t2009, t2005, -t2004, 0, 0, 0, 0, 0, 0, t1959, t1964, t1948, t1908, 0, 0, 0, 0, 0, 0, t1912, t1927, t1897, t1873, 0, 0, 0, 0, 0, 0, t1912, t1927, t1897, t1869; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1988, t1981, t1974, t1931, 0, 0, 0, 0, 0, 0, t1935, t1947, t1919, t1884, 0, 0, 0, 0, 0, 0, t1935, t1947, t1919, t1878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1980, t1993, t1976, t1930, 0, 0, 0, 0, 0, 0, t1934, t1946, t1918, t1883, 0, 0, 0, 0, 0, 0, t1934, t1946, t1918, t1877; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1973, t1975, t1985, t1954, 0, 0, 0, 0, 0, 0, t1940, t1942, t1956, t1929, 0, 0, 0, 0, 0, 0, t1940, t1942, t1956, t1907; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1970, t1967, t1941, t1900, 0, 0, 0, 0, 0, 0, t1970, t1967, t1941, t1894; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1966, t1979, t1943, t1899, 0, 0, 0, 0, 0, 0, t1966, t1979, t1943, t1893; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1940, t1942, t1956, t1929, 0, 0, 0, 0, 0, 0, t1940, t1942, t1956, t1907; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1970, t1967, t1941, t1894; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1966, t1979, t1943, t1893; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1940, t1942, t1956, t1907;];
f_new_reg  = t1;