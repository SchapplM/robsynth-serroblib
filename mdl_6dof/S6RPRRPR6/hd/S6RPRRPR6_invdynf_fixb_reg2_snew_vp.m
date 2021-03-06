% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRRPR6
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRRPR6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:52:58
% EndTime: 2019-05-05 22:53:07
% DurationCPUTime: 9.61s
% Computational Cost: add. (95588->332), mult. (231691->465), div. (0->0), fcn. (180992->12), ass. (0->259)
t2207 = qJD(1) ^ 2;
t2202 = sin(qJ(1));
t2206 = cos(qJ(1));
t2180 = t2202 * g(1) - t2206 * g(2);
t2218 = -qJDD(2) + t2180;
t2197 = cos(pkin(10));
t2224 = pkin(2) * t2197 + pkin(1);
t2195 = sin(pkin(10));
t2191 = t2195 ^ 2;
t2192 = t2197 ^ 2;
t2227 = t2191 + t2192;
t2149 = t2224 * qJDD(1) + (pkin(7) * t2227 + qJ(2)) * t2207 + t2218;
t2177 = t2227 * t2207;
t2201 = sin(qJ(3));
t2205 = cos(qJ(3));
t2171 = (t2195 * t2201 - t2197 * t2205) * qJD(1);
t2168 = qJD(4) + t2171;
t2163 = qJD(6) + t2168;
t2257 = qJD(6) + t2163;
t2217 = t2195 * t2205 + t2197 * t2201;
t2256 = t2217 * qJDD(1);
t2253 = qJD(3) ^ 2;
t2173 = t2217 * qJD(1);
t2200 = sin(qJ(4));
t2204 = cos(qJ(4));
t2159 = -t2204 * qJD(3) + t2173 * t2200;
t2160 = qJD(3) * t2200 + t2173 * t2204;
t2194 = sin(pkin(11));
t2196 = cos(pkin(11));
t2134 = t2196 * t2159 + t2160 * t2194;
t2136 = -t2159 * t2194 + t2160 * t2196;
t2199 = sin(qJ(6));
t2203 = cos(qJ(6));
t2102 = t2203 * t2134 + t2136 * t2199;
t2252 = t2102 ^ 2;
t2104 = -t2134 * t2199 + t2136 * t2203;
t2251 = t2104 ^ 2;
t2250 = t2134 ^ 2;
t2249 = t2136 ^ 2;
t2248 = t2159 ^ 2;
t2247 = t2160 ^ 2;
t2246 = t2163 ^ 2;
t2245 = t2168 ^ 2;
t2244 = t2171 ^ 2;
t2243 = t2173 ^ 2;
t2242 = -2 * qJD(5);
t2241 = t2197 * g(3);
t2240 = qJD(2) * qJD(1);
t2239 = t2102 * t2104;
t2238 = t2134 * t2136;
t2237 = t2136 * t2168;
t2236 = t2159 * t2160;
t2235 = t2168 * t2134;
t2234 = t2168 * t2159;
t2233 = t2171 * qJD(3);
t2232 = t2171 * t2173;
t2166 = t2173 * qJD(3);
t2231 = t2192 * t2207;
t2230 = qJD(4) - t2168;
t2229 = qJD(6) - t2163;
t2219 = t2206 * g(1) + t2202 * g(2);
t2211 = -t2207 * pkin(1) + qJDD(1) * qJ(2) - t2219;
t2156 = -g(3) * t2195 + (t2211 + 0.2e1 * t2240) * t2197;
t2187 = t2197 * qJDD(1);
t2143 = -pkin(2) * t2231 + pkin(7) * t2187 + t2156;
t2225 = -0.2e1 * t2240;
t2210 = -t2241 + (t2225 + t2224 * t2207 + (-pkin(7) - qJ(2)) * qJDD(1) + t2219) * t2195;
t2110 = t2205 * t2143 + t2201 * t2210;
t2146 = pkin(3) * t2171 - pkin(8) * t2173;
t2089 = -pkin(3) * t2253 + qJDD(3) * pkin(8) - t2146 * t2171 + t2110;
t2154 = t2256 - t2233;
t2226 = t2195 * qJDD(1);
t2228 = -t2205 * t2187 + t2201 * t2226;
t2223 = -t2166 - t2228;
t2097 = (-t2154 + t2233) * pkin(8) + (-t2223 + t2166) * pkin(3) - t2149;
t2057 = t2204 * t2089 + t2200 * t2097;
t2056 = -t2200 * t2089 + t2204 * t2097;
t2145 = qJDD(4) - t2223;
t2112 = t2145 - t2236;
t2216 = -t2200 * qJDD(3) - t2204 * t2154;
t2123 = -qJD(4) * t2159 - t2216;
t2043 = (-t2123 - t2234) * qJ(5) + t2112 * pkin(4) + t2056;
t2141 = pkin(4) * t2168 - qJ(5) * t2160;
t2220 = -t2204 * qJDD(3) + t2200 * t2154;
t2214 = -qJD(4) * t2160 - t2220;
t2045 = -pkin(4) * t2248 + qJ(5) * t2214 - t2168 * t2141 + t2057;
t2011 = t2194 * t2043 + t2196 * t2045 + t2134 * t2242;
t2094 = t2196 * t2123 + t2194 * t2214;
t2221 = t2123 * t2194 - t2196 * t2214;
t2222 = -t2199 * t2094 - t2203 * t2221;
t2074 = -t2094 - t2235;
t2090 = t2145 - t2238;
t2109 = -t2201 * t2143 + t2205 * t2210;
t2010 = t2196 * t2043 - t2194 * t2045 + t2136 * t2242;
t2215 = -qJDD(6) - t2145;
t2212 = -t2203 * t2094 + t2199 * t2221;
t2088 = -qJDD(3) * pkin(3) - t2253 * pkin(8) + t2173 * t2146 - t2109;
t2053 = -t2214 * pkin(4) - t2248 * qJ(5) + t2160 * t2141 + qJDD(5) + t2088;
t2182 = t2195 * t2207 * t2197;
t2179 = -qJDD(1) * t2202 - t2206 * t2207;
t2178 = qJDD(1) * t2206 - t2202 * t2207;
t2176 = t2227 * qJDD(1);
t2175 = t2197 * t2177;
t2174 = t2195 * t2177;
t2169 = qJDD(1) * pkin(1) + t2207 * qJ(2) + t2218;
t2161 = -t2243 - t2253;
t2155 = -t2241 + (-t2211 + t2225) * t2195;
t2153 = t2256 - 0.2e1 * t2233;
t2152 = 0.2e1 * t2166 + t2228;
t2151 = -qJDD(3) - t2232;
t2150 = qJDD(3) - t2232;
t2147 = -t2244 - t2253;
t2138 = -t2243 - t2244;
t2137 = t2223 + t2166;
t2129 = -t2245 - t2247;
t2128 = t2151 * t2205 - t2161 * t2201;
t2127 = t2151 * t2201 + t2161 * t2205;
t2126 = -t2245 - t2248;
t2125 = -t2155 * t2195 + t2156 * t2197;
t2124 = t2155 * t2197 + t2156 * t2195;
t2119 = -t2247 - t2248;
t2118 = t2137 * t2205 + t2201 * t2256;
t2117 = t2137 * t2201 - t2205 * t2256;
t2116 = t2147 * t2205 - t2150 * t2201;
t2115 = t2147 * t2201 + t2150 * t2205;
t2114 = pkin(5) * t2168 - pkin(9) * t2136;
t2113 = -t2145 - t2236;
t2111 = -t2245 - t2249;
t2108 = t2159 * t2230 + t2216;
t2107 = t2123 - t2234;
t2106 = -t2160 * t2230 - t2220;
t2105 = (qJD(4) + t2168) * t2160 + t2220;
t2101 = -t2127 * t2195 + t2128 * t2197;
t2100 = t2127 * t2197 + t2128 * t2195;
t2099 = -t2245 - t2250;
t2091 = -t2145 - t2238;
t2087 = t2113 * t2204 - t2129 * t2200;
t2086 = t2113 * t2200 + t2129 * t2204;
t2084 = -t2117 * t2195 + t2118 * t2197;
t2083 = t2117 * t2197 + t2118 * t2195;
t2082 = -t2246 - t2251;
t2081 = -t2112 * t2200 + t2126 * t2204;
t2080 = t2112 * t2204 + t2126 * t2200;
t2079 = -t2115 * t2195 + t2116 * t2197;
t2078 = t2115 * t2197 + t2116 * t2195;
t2077 = -t2249 - t2250;
t2076 = -t2109 * t2201 + t2110 * t2205;
t2075 = t2109 * t2205 + t2110 * t2201;
t2073 = t2094 - t2235;
t2072 = -t2221 + t2237;
t2071 = t2221 + t2237;
t2070 = t2106 * t2204 - t2108 * t2200;
t2069 = t2106 * t2200 + t2108 * t2204;
t2068 = t2091 * t2196 - t2111 * t2194;
t2067 = t2091 * t2194 + t2111 * t2196;
t2066 = -t2246 - t2252;
t2065 = t2087 * t2205 + t2107 * t2201;
t2064 = t2087 * t2201 - t2107 * t2205;
t2063 = t2081 * t2205 + t2105 * t2201;
t2062 = t2081 * t2201 - t2105 * t2205;
t2061 = t2215 - t2239;
t2060 = -t2215 - t2239;
t2059 = -t2090 * t2194 + t2099 * t2196;
t2058 = t2090 * t2196 + t2099 * t2194;
t2055 = t2070 * t2205 + t2119 * t2201;
t2054 = t2070 * t2201 - t2119 * t2205;
t2052 = -t2251 - t2252;
t2051 = -t2075 * t2195 + t2076 * t2197;
t2050 = t2075 * t2197 + t2076 * t2195;
t2049 = t2061 * t2203 - t2082 * t2199;
t2048 = t2061 * t2199 + t2082 * t2203;
t2047 = t2072 * t2196 - t2074 * t2194;
t2046 = t2072 * t2194 + t2074 * t2196;
t2040 = -t2067 * t2200 + t2068 * t2204;
t2039 = t2067 * t2204 + t2068 * t2200;
t2038 = -t2060 * t2199 + t2066 * t2203;
t2037 = t2060 * t2203 + t2066 * t2199;
t2036 = -t2064 * t2195 + t2065 * t2197;
t2035 = t2064 * t2197 + t2065 * t2195;
t2034 = t2102 * t2229 + t2212;
t2033 = -t2102 * t2257 - t2212;
t2032 = -t2104 * t2229 + t2222;
t2031 = t2104 * t2257 - t2222;
t2030 = -t2062 * t2195 + t2063 * t2197;
t2029 = t2062 * t2197 + t2063 * t2195;
t2028 = -t2058 * t2200 + t2059 * t2204;
t2027 = t2058 * t2204 + t2059 * t2200;
t2026 = -t2056 * t2200 + t2057 * t2204;
t2025 = t2056 * t2204 + t2057 * t2200;
t2024 = -t2054 * t2195 + t2055 * t2197;
t2023 = t2054 * t2197 + t2055 * t2195;
t2022 = pkin(5) * t2221 - pkin(9) * t2250 + t2136 * t2114 + t2053;
t2021 = t2040 * t2205 + t2073 * t2201;
t2020 = t2040 * t2201 - t2073 * t2205;
t2019 = t2026 * t2205 + t2088 * t2201;
t2018 = t2026 * t2201 - t2088 * t2205;
t2017 = t2028 * t2205 + t2071 * t2201;
t2016 = t2028 * t2201 - t2071 * t2205;
t2015 = -t2048 * t2194 + t2049 * t2196;
t2014 = t2048 * t2196 + t2049 * t2194;
t2013 = -t2046 * t2200 + t2047 * t2204;
t2012 = t2046 * t2204 + t2047 * t2200;
t2009 = t2013 * t2205 + t2077 * t2201;
t2008 = t2013 * t2201 - t2077 * t2205;
t2007 = -t2037 * t2194 + t2038 * t2196;
t2006 = t2037 * t2196 + t2038 * t2194;
t2005 = t2032 * t2203 - t2034 * t2199;
t2004 = t2032 * t2199 + t2034 * t2203;
t2003 = -pkin(5) * t2250 - pkin(9) * t2221 - t2168 * t2114 + t2011;
t2002 = pkin(5) * t2090 + pkin(9) * t2074 + t2010;
t2001 = -t2020 * t2195 + t2021 * t2197;
t2000 = t2020 * t2197 + t2021 * t2195;
t1999 = -t2018 * t2195 + t2019 * t2197;
t1998 = t2018 * t2197 + t2019 * t2195;
t1997 = -t2016 * t2195 + t2017 * t2197;
t1996 = t2016 * t2197 + t2017 * t2195;
t1995 = -t2014 * t2200 + t2015 * t2204;
t1994 = t2014 * t2204 + t2015 * t2200;
t1993 = -t2010 * t2194 + t2011 * t2196;
t1992 = t2010 * t2196 + t2011 * t2194;
t1991 = -t2008 * t2195 + t2009 * t2197;
t1990 = t2008 * t2197 + t2009 * t2195;
t1989 = -t2006 * t2200 + t2007 * t2204;
t1988 = t2006 * t2204 + t2007 * t2200;
t1987 = -t2004 * t2194 + t2005 * t2196;
t1986 = t2004 * t2196 + t2005 * t2194;
t1985 = t1995 * t2205 + t2033 * t2201;
t1984 = t1995 * t2201 - t2033 * t2205;
t1983 = t2002 * t2199 + t2003 * t2203;
t1982 = t2002 * t2203 - t2003 * t2199;
t1981 = t1989 * t2205 + t2031 * t2201;
t1980 = t1989 * t2201 - t2031 * t2205;
t1979 = -t1992 * t2200 + t1993 * t2204;
t1978 = t1992 * t2204 + t1993 * t2200;
t1977 = t1979 * t2205 + t2053 * t2201;
t1976 = t1979 * t2201 - t2053 * t2205;
t1975 = -t1986 * t2200 + t1987 * t2204;
t1974 = t1986 * t2204 + t1987 * t2200;
t1973 = -t1984 * t2195 + t1985 * t2197;
t1972 = t1984 * t2197 + t1985 * t2195;
t1971 = t1975 * t2205 + t2052 * t2201;
t1970 = t1975 * t2201 - t2052 * t2205;
t1969 = -t1982 * t2199 + t1983 * t2203;
t1968 = t1982 * t2203 + t1983 * t2199;
t1967 = -t1980 * t2195 + t1981 * t2197;
t1966 = t1980 * t2197 + t1981 * t2195;
t1965 = -t1976 * t2195 + t1977 * t2197;
t1964 = t1976 * t2197 + t1977 * t2195;
t1963 = -t1970 * t2195 + t1971 * t2197;
t1962 = t1970 * t2197 + t1971 * t2195;
t1961 = -t1968 * t2194 + t1969 * t2196;
t1960 = t1968 * t2196 + t1969 * t2194;
t1959 = -t1960 * t2200 + t1961 * t2204;
t1958 = t1960 * t2204 + t1961 * t2200;
t1957 = t1959 * t2205 + t2022 * t2201;
t1956 = t1959 * t2201 - t2022 * t2205;
t1955 = -t1956 * t2195 + t1957 * t2197;
t1954 = t1956 * t2197 + t1957 * t2195;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2179, -t2178, 0, -t2180 * t2202 - t2206 * t2219, 0, 0, 0, 0, 0, 0, -t2175 * t2206 - t2187 * t2202, t2174 * t2206 + t2202 * t2226, t2176 * t2206 - t2177 * t2202, t2125 * t2206 - t2169 * t2202, 0, 0, 0, 0, 0, 0, t2079 * t2206 + t2152 * t2202, t2101 * t2206 + t2153 * t2202, t2084 * t2206 + t2138 * t2202, t2051 * t2206 - t2149 * t2202, 0, 0, 0, 0, 0, 0, t2030 * t2206 + t2080 * t2202, t2036 * t2206 + t2086 * t2202, t2024 * t2206 + t2069 * t2202, t1999 * t2206 + t2025 * t2202, 0, 0, 0, 0, 0, 0, t1997 * t2206 + t2027 * t2202, t2001 * t2206 + t2039 * t2202, t1991 * t2206 + t2012 * t2202, t1965 * t2206 + t1978 * t2202, 0, 0, 0, 0, 0, 0, t1967 * t2206 + t1988 * t2202, t1973 * t2206 + t1994 * t2202, t1963 * t2206 + t1974 * t2202, t1955 * t2206 + t1958 * t2202; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2178, t2179, 0, t2180 * t2206 - t2202 * t2219, 0, 0, 0, 0, 0, 0, -t2175 * t2202 + t2187 * t2206, t2174 * t2202 - t2206 * t2226, t2176 * t2202 + t2177 * t2206, t2125 * t2202 + t2169 * t2206, 0, 0, 0, 0, 0, 0, t2079 * t2202 - t2152 * t2206, t2101 * t2202 - t2153 * t2206, t2084 * t2202 - t2138 * t2206, t2051 * t2202 + t2149 * t2206, 0, 0, 0, 0, 0, 0, t2030 * t2202 - t2080 * t2206, t2036 * t2202 - t2086 * t2206, t2024 * t2202 - t2069 * t2206, t1999 * t2202 - t2025 * t2206, 0, 0, 0, 0, 0, 0, t1997 * t2202 - t2027 * t2206, t2001 * t2202 - t2039 * t2206, t1991 * t2202 - t2012 * t2206, t1965 * t2202 - t1978 * t2206, 0, 0, 0, 0, 0, 0, t1967 * t2202 - t1988 * t2206, t1973 * t2202 - t1994 * t2206, t1963 * t2202 - t1974 * t2206, t1955 * t2202 - t1958 * t2206; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2124, 0, 0, 0, 0, 0, 0, t2078, t2100, t2083, t2050, 0, 0, 0, 0, 0, 0, t2029, t2035, t2023, t1998, 0, 0, 0, 0, 0, 0, t1996, t2000, t1990, t1964, 0, 0, 0, 0, 0, 0, t1966, t1972, t1962, t1954; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2207, -qJDD(1), 0, -t2219, 0, 0, 0, 0, 0, 0, -t2175, t2174, t2176, t2125, 0, 0, 0, 0, 0, 0, t2079, t2101, t2084, t2051, 0, 0, 0, 0, 0, 0, t2030, t2036, t2024, t1999, 0, 0, 0, 0, 0, 0, t1997, t2001, t1991, t1965, 0, 0, 0, 0, 0, 0, t1967, t1973, t1963, t1955; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2207, 0, t2180, 0, 0, 0, 0, 0, 0, t2187, -t2226, t2177, t2169, 0, 0, 0, 0, 0, 0, -t2152, -t2153, -t2138, t2149, 0, 0, 0, 0, 0, 0, -t2080, -t2086, -t2069, -t2025, 0, 0, 0, 0, 0, 0, -t2027, -t2039, -t2012, -t1978, 0, 0, 0, 0, 0, 0, -t1988, -t1994, -t1974, -t1958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2124, 0, 0, 0, 0, 0, 0, t2078, t2100, t2083, t2050, 0, 0, 0, 0, 0, 0, t2029, t2035, t2023, t1998, 0, 0, 0, 0, 0, 0, t1996, t2000, t1990, t1964, 0, 0, 0, 0, 0, 0, t1966, t1972, t1962, t1954; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2231, t2182, t2187, t2156, 0, 0, 0, 0, 0, 0, t2116, t2128, t2118, t2076, 0, 0, 0, 0, 0, 0, t2063, t2065, t2055, t2019, 0, 0, 0, 0, 0, 0, t2017, t2021, t2009, t1977, 0, 0, 0, 0, 0, 0, t1981, t1985, t1971, t1957; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2182, -t2191 * t2207, -t2226, t2155, 0, 0, 0, 0, 0, 0, t2115, t2127, t2117, t2075, 0, 0, 0, 0, 0, 0, t2062, t2064, t2054, t2018, 0, 0, 0, 0, 0, 0, t2016, t2020, t2008, t1976, 0, 0, 0, 0, 0, 0, t1980, t1984, t1970, t1956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2187, t2226, -t2177, -t2169, 0, 0, 0, 0, 0, 0, t2152, t2153, t2138, -t2149, 0, 0, 0, 0, 0, 0, t2080, t2086, t2069, t2025, 0, 0, 0, 0, 0, 0, t2027, t2039, t2012, t1978, 0, 0, 0, 0, 0, 0, t1988, t1994, t1974, t1958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2147, t2151, t2137, t2110, 0, 0, 0, 0, 0, 0, t2081, t2087, t2070, t2026, 0, 0, 0, 0, 0, 0, t2028, t2040, t2013, t1979, 0, 0, 0, 0, 0, 0, t1989, t1995, t1975, t1959; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2150, t2161, -t2256, t2109, 0, 0, 0, 0, 0, 0, -t2105, -t2107, -t2119, -t2088, 0, 0, 0, 0, 0, 0, -t2071, -t2073, -t2077, -t2053, 0, 0, 0, 0, 0, 0, -t2031, -t2033, -t2052, -t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2152, t2153, t2138, -t2149, 0, 0, 0, 0, 0, 0, t2080, t2086, t2069, t2025, 0, 0, 0, 0, 0, 0, t2027, t2039, t2012, t1978, 0, 0, 0, 0, 0, 0, t1988, t1994, t1974, t1958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2126, t2113, t2106, t2057, 0, 0, 0, 0, 0, 0, t2059, t2068, t2047, t1993, 0, 0, 0, 0, 0, 0, t2007, t2015, t1987, t1961; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2112, t2129, t2108, t2056, 0, 0, 0, 0, 0, 0, t2058, t2067, t2046, t1992, 0, 0, 0, 0, 0, 0, t2006, t2014, t1986, t1960; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2105, t2107, t2119, t2088, 0, 0, 0, 0, 0, 0, t2071, t2073, t2077, t2053, 0, 0, 0, 0, 0, 0, t2031, t2033, t2052, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2099, t2091, t2072, t2011, 0, 0, 0, 0, 0, 0, t2038, t2049, t2005, t1969; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2090, t2111, t2074, t2010, 0, 0, 0, 0, 0, 0, t2037, t2048, t2004, t1968; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2071, t2073, t2077, t2053, 0, 0, 0, 0, 0, 0, t2031, t2033, t2052, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2066, t2061, t2032, t1983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2060, t2082, t2034, t1982; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2031, t2033, t2052, t2022;];
f_new_reg  = t1;
