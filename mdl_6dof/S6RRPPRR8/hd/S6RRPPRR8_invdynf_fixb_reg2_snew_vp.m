% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPPRR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:23:18
% EndTime: 2019-05-06 11:23:25
% DurationCPUTime: 8.17s
% Computational Cost: add. (37880->322), mult. (84218->383), div. (0->0), fcn. (58249->10), ass. (0->222)
t2211 = sin(qJ(2));
t2253 = qJD(1) * t2211;
t2197 = qJD(2) * t2253;
t2215 = cos(qJ(2));
t2234 = t2215 * qJDD(1);
t2179 = -t2197 + t2234;
t2207 = sin(pkin(10));
t2208 = cos(pkin(10));
t2173 = -t2208 * qJD(2) + t2207 * t2253;
t2175 = t2207 * qJD(2) + t2208 * t2253;
t2243 = t2173 * t2175;
t2133 = t2179 - t2243;
t2171 = t2175 ^ 2;
t2205 = t2215 ^ 2;
t2217 = qJD(1) ^ 2;
t2201 = t2205 * t2217;
t2264 = -t2171 - t2201;
t2112 = t2133 * t2208 - t2207 * t2264;
t2252 = qJD(1) * t2215;
t2232 = qJD(2) * t2252;
t2235 = t2211 * qJDD(1);
t2223 = t2232 + t2235;
t2220 = t2207 * qJDD(2) + t2208 * t2223;
t2242 = t2173 * t2215;
t2233 = qJD(1) * t2242;
t2128 = t2220 + t2233;
t2087 = t2112 * t2215 + t2128 * t2211;
t2110 = t2133 * t2207 + t2208 * t2264;
t2212 = sin(qJ(1));
t2216 = cos(qJ(1));
t2288 = t2087 * t2212 - t2110 * t2216;
t2287 = t2087 * t2216 + t2110 * t2212;
t2085 = t2112 * t2211 - t2128 * t2215;
t2132 = t2220 - t2233;
t2163 = t2175 * t2252;
t2221 = -t2208 * qJDD(2) + t2207 * t2223;
t2265 = -t2163 - t2221;
t2269 = -t2132 * t2208 + t2207 * t2265;
t2258 = t2173 ^ 2;
t2125 = t2171 + t2258;
t2267 = t2132 * t2207 + t2208 * t2265;
t2273 = -t2125 * t2211 + t2215 * t2267;
t2284 = t2212 * t2273 - t2216 * t2269;
t2134 = t2179 + t2243;
t2263 = -t2258 - t2201;
t2270 = -t2134 * t2208 + t2207 * t2263;
t2129 = -t2163 + t2221;
t2268 = t2134 * t2207 + t2208 * t2263;
t2275 = t2129 * t2211 + t2215 * t2268;
t2283 = t2212 * t2275 - t2216 * t2270;
t2282 = t2212 * t2269 + t2216 * t2273;
t2281 = t2212 * t2270 + t2216 * t2275;
t2276 = -t2215 * t2129 + t2211 * t2268;
t2274 = t2125 * t2215 + t2211 * t2267;
t2194 = qJD(5) + t2252;
t2189 = qJD(6) + t2194;
t2266 = qJD(6) + t2189;
t2210 = sin(qJ(5));
t2214 = cos(qJ(5));
t2141 = -t2214 * t2173 + t2175 * t2210;
t2097 = -t2141 * qJD(5) + t2210 * t2221 + t2214 * t2220;
t2244 = t2141 * t2194;
t2078 = -t2097 - t2244;
t2262 = qJD(2) ^ 2;
t2143 = t2173 * t2210 + t2175 * t2214;
t2209 = sin(qJ(6));
t2213 = cos(qJ(6));
t2107 = t2213 * t2141 + t2143 * t2209;
t2261 = t2107 ^ 2;
t2109 = -t2141 * t2209 + t2143 * t2213;
t2260 = t2109 ^ 2;
t2140 = t2141 ^ 2;
t2259 = t2143 ^ 2;
t2257 = t2189 ^ 2;
t2256 = t2194 ^ 2;
t2255 = -0.2e1 * t2175;
t2254 = t2215 * g(3);
t2251 = t2107 * t2109;
t2245 = t2141 * t2143;
t2177 = (-pkin(2) * t2215 - qJ(3) * t2211) * qJD(1);
t2241 = t2177 * t2211;
t2240 = t2194 * t2143;
t2204 = t2211 ^ 2;
t2239 = t2204 * t2217;
t2188 = -g(1) * t2216 - g(2) * t2212;
t2169 = -pkin(1) * t2217 + qJDD(1) * pkin(7) + t2188;
t2238 = t2211 * t2169;
t2146 = pkin(3) * t2173 - qJ(4) * t2175;
t2237 = (2 * qJD(3)) + t2146;
t2236 = qJD(6) - t2189;
t2156 = -g(3) * t2211 + t2215 * t2169;
t2121 = -pkin(2) * t2262 + qJDD(2) * qJ(3) + t2177 * t2252 + t2156;
t2187 = t2212 * g(1) - t2216 * g(2);
t2168 = qJDD(1) * pkin(1) + t2217 * pkin(7) + t2187;
t2178 = 0.2e1 * t2232 + t2235;
t2218 = (-t2179 + t2197) * pkin(2) - t2178 * qJ(3) - t2168;
t2230 = t2207 * t2121 - t2208 * t2218;
t2222 = t2179 * pkin(3) - qJ(4) * t2201 + qJDD(4) + t2230;
t2050 = t2179 * pkin(4) - t2132 * pkin(8) + (pkin(4) * t2173 + t2237) * t2175 + t2222;
t2084 = -0.2e1 * qJD(3) * t2173 + t2208 * t2121 + t2207 * t2218;
t2061 = -pkin(3) * t2201 - t2179 * qJ(4) - 0.2e1 * qJD(4) * t2252 - t2173 * t2146 + t2084;
t2226 = pkin(4) * t2252 - pkin(8) * t2175;
t2058 = -pkin(4) * t2258 + pkin(8) * t2221 - t2226 * t2252 + t2061;
t2024 = t2210 * t2050 + t2214 * t2058;
t2023 = t2214 * t2050 - t2210 * t2058;
t2096 = -t2143 * qJD(5) - t2210 * t2220 + t2214 * t2221;
t2231 = t2213 * t2096 - t2209 * t2097;
t2229 = pkin(5) * t2194 - pkin(9) * t2143;
t2228 = qJDD(5) + t2179;
t2227 = -t2209 * t2096 - t2213 * t2097;
t2225 = -qJDD(2) * pkin(2) - qJ(3) * t2262 + qJDD(3) + t2254;
t2224 = -qJDD(6) - t2228;
t2099 = t2228 - t2245;
t2219 = pkin(3) * t2129 - t2220 * qJ(4) + qJD(4) * t2255 + t2225 + t2238;
t2060 = pkin(4) * t2221 + pkin(8) * t2258 - qJ(4) * t2233 + qJD(1) * t2241 - t2175 * t2226 + t2219;
t2193 = t2215 * t2217 * t2211;
t2192 = -t2201 - t2262;
t2191 = -t2239 - t2262;
t2186 = -qJDD(2) + t2193;
t2185 = qJDD(2) + t2193;
t2184 = t2201 + t2239;
t2183 = -qJDD(1) * t2212 - t2216 * t2217;
t2182 = qJDD(1) * t2216 - t2212 * t2217;
t2181 = (t2204 + t2205) * qJDD(1);
t2180 = -0.2e1 * t2197 + t2234;
t2155 = -t2238 - t2254;
t2154 = t2186 * t2215 - t2191 * t2211;
t2153 = -t2185 * t2211 + t2192 * t2215;
t2152 = t2186 * t2211 + t2191 * t2215;
t2151 = t2185 * t2215 + t2192 * t2211;
t2120 = (qJD(1) * t2177 + t2169) * t2211 + t2225;
t2119 = -t2256 - t2259;
t2117 = -t2155 * t2211 + t2156 * t2215;
t2116 = t2155 * t2215 + t2156 * t2211;
t2106 = -t2256 - t2140;
t2100 = -t2228 - t2245;
t2090 = -t2257 - t2260;
t2089 = -t2140 - t2259;
t2083 = qJD(3) * t2255 - t2230;
t2077 = -t2097 + t2244;
t2076 = t2096 + t2240;
t2075 = t2096 - t2240;
t2070 = t2100 * t2214 - t2119 * t2210;
t2069 = t2100 * t2210 + t2119 * t2214;
t2068 = (-qJ(4) * t2242 + t2241) * qJD(1) + t2219;
t2067 = -t2257 - t2261;
t2066 = -t2099 * t2210 + t2106 * t2214;
t2065 = t2099 * t2214 + t2106 * t2210;
t2064 = t2224 - t2251;
t2063 = -t2224 - t2251;
t2062 = t2237 * t2175 + t2222;
t2059 = -t2260 - t2261;
t2057 = -t2083 * t2207 + t2084 * t2208;
t2056 = t2083 * t2208 + t2084 * t2207;
t2055 = t2076 * t2214 - t2078 * t2210;
t2054 = t2076 * t2210 + t2078 * t2214;
t2052 = t2064 * t2213 - t2090 * t2209;
t2051 = t2064 * t2209 + t2090 * t2213;
t2047 = t2069 * t2207 + t2070 * t2208;
t2046 = -t2069 * t2208 + t2070 * t2207;
t2045 = t2057 * t2215 + t2120 * t2211;
t2044 = t2057 * t2211 - t2120 * t2215;
t2043 = -t2063 * t2209 + t2067 * t2213;
t2042 = t2063 * t2213 + t2067 * t2209;
t2041 = t2065 * t2207 + t2066 * t2208;
t2040 = -t2065 * t2208 + t2066 * t2207;
t2039 = t2107 * t2236 + t2227;
t2038 = -t2107 * t2266 - t2227;
t2037 = -t2109 * t2236 + t2231;
t2036 = t2109 * t2266 - t2231;
t2035 = t2061 * t2208 + t2062 * t2207;
t2034 = t2061 * t2207 - t2062 * t2208;
t2033 = t2047 * t2215 + t2077 * t2211;
t2032 = t2047 * t2211 - t2077 * t2215;
t2031 = t2096 * pkin(5) + t2140 * pkin(9) - t2143 * t2229 + t2060;
t2030 = t2041 * t2215 + t2075 * t2211;
t2029 = t2041 * t2211 - t2075 * t2215;
t2028 = t2054 * t2207 + t2055 * t2208;
t2027 = -t2054 * t2208 + t2055 * t2207;
t2026 = -t2051 * t2210 + t2052 * t2214;
t2025 = t2051 * t2214 + t2052 * t2210;
t2022 = t2035 * t2215 + t2068 * t2211;
t2021 = t2035 * t2211 - t2068 * t2215;
t2020 = t2028 * t2215 - t2089 * t2211;
t2019 = t2028 * t2211 + t2089 * t2215;
t2018 = -t2042 * t2210 + t2043 * t2214;
t2017 = t2042 * t2214 + t2043 * t2210;
t2016 = t2037 * t2213 - t2039 * t2209;
t2015 = t2037 * t2209 + t2039 * t2213;
t2014 = -t2140 * pkin(5) + t2096 * pkin(9) - t2194 * t2229 + t2024;
t2013 = t2099 * pkin(5) + pkin(9) * t2078 + t2023;
t2012 = -t2023 * t2210 + t2024 * t2214;
t2011 = t2023 * t2214 + t2024 * t2210;
t2010 = t2025 * t2207 + t2026 * t2208;
t2009 = -t2025 * t2208 + t2026 * t2207;
t2008 = t2017 * t2207 + t2018 * t2208;
t2007 = -t2017 * t2208 + t2018 * t2207;
t2006 = -t2015 * t2210 + t2016 * t2214;
t2005 = t2015 * t2214 + t2016 * t2210;
t2004 = t2010 * t2215 - t2038 * t2211;
t2003 = t2010 * t2211 + t2038 * t2215;
t2002 = t2013 * t2209 + t2014 * t2213;
t2001 = t2013 * t2213 - t2014 * t2209;
t2000 = t2008 * t2215 - t2036 * t2211;
t1999 = t2008 * t2211 + t2036 * t2215;
t1998 = t2011 * t2207 + t2012 * t2208;
t1997 = -t2011 * t2208 + t2012 * t2207;
t1996 = t1998 * t2215 + t2060 * t2211;
t1995 = t1998 * t2211 - t2060 * t2215;
t1994 = t2005 * t2207 + t2006 * t2208;
t1993 = -t2005 * t2208 + t2006 * t2207;
t1992 = -t2001 * t2209 + t2002 * t2213;
t1991 = t2001 * t2213 + t2002 * t2209;
t1990 = t1994 * t2215 - t2059 * t2211;
t1989 = t1994 * t2211 + t2059 * t2215;
t1988 = -t1991 * t2210 + t1992 * t2214;
t1987 = t1991 * t2214 + t1992 * t2210;
t1986 = t1987 * t2207 + t1988 * t2208;
t1985 = -t1987 * t2208 + t1988 * t2207;
t1984 = t1986 * t2215 + t2031 * t2211;
t1983 = t1986 * t2211 - t2031 * t2215;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2183, -t2182, 0, -t2187 * t2212 + t2188 * t2216, 0, 0, 0, 0, 0, 0, t2153 * t2216 - t2180 * t2212, t2154 * t2216 + t2178 * t2212, t2181 * t2216 - t2184 * t2212, t2117 * t2216 - t2168 * t2212, 0, 0, 0, 0, 0, 0, t2281, t2287, t2282, t2045 * t2216 + t2056 * t2212, 0, 0, 0, 0, 0, 0, t2281, t2282, -t2287, t2022 * t2216 + t2034 * t2212, 0, 0, 0, 0, 0, 0, t2030 * t2216 + t2040 * t2212, t2033 * t2216 + t2046 * t2212, t2020 * t2216 + t2027 * t2212, t1996 * t2216 + t1997 * t2212, 0, 0, 0, 0, 0, 0, t2000 * t2216 + t2007 * t2212, t2004 * t2216 + t2009 * t2212, t1990 * t2216 + t1993 * t2212, t1984 * t2216 + t1985 * t2212; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2182, t2183, 0, t2187 * t2216 + t2188 * t2212, 0, 0, 0, 0, 0, 0, t2153 * t2212 + t2180 * t2216, t2154 * t2212 - t2178 * t2216, t2181 * t2212 + t2184 * t2216, t2117 * t2212 + t2168 * t2216, 0, 0, 0, 0, 0, 0, t2283, t2288, t2284, t2045 * t2212 - t2056 * t2216, 0, 0, 0, 0, 0, 0, t2283, t2284, -t2288, t2022 * t2212 - t2034 * t2216, 0, 0, 0, 0, 0, 0, t2030 * t2212 - t2040 * t2216, t2033 * t2212 - t2046 * t2216, t2020 * t2212 - t2027 * t2216, t1996 * t2212 - t1997 * t2216, 0, 0, 0, 0, 0, 0, t2000 * t2212 - t2007 * t2216, t2004 * t2212 - t2009 * t2216, t1990 * t2212 - t1993 * t2216, t1984 * t2212 - t1985 * t2216; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2151, t2152, 0, t2116, 0, 0, 0, 0, 0, 0, t2276, t2085, t2274, t2044, 0, 0, 0, 0, 0, 0, t2276, t2274, -t2085, t2021, 0, 0, 0, 0, 0, 0, t2029, t2032, t2019, t1995, 0, 0, 0, 0, 0, 0, t1999, t2003, t1989, t1983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2217, -qJDD(1), 0, t2188, 0, 0, 0, 0, 0, 0, t2153, t2154, t2181, t2117, 0, 0, 0, 0, 0, 0, t2275, t2087, t2273, t2045, 0, 0, 0, 0, 0, 0, t2275, t2273, -t2087, t2022, 0, 0, 0, 0, 0, 0, t2030, t2033, t2020, t1996, 0, 0, 0, 0, 0, 0, t2000, t2004, t1990, t1984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2217, 0, t2187, 0, 0, 0, 0, 0, 0, t2180, -t2178, t2184, t2168, 0, 0, 0, 0, 0, 0, -t2270, -t2110, -t2269, -t2056, 0, 0, 0, 0, 0, 0, -t2270, -t2269, t2110, -t2034, 0, 0, 0, 0, 0, 0, -t2040, -t2046, -t2027, -t1997, 0, 0, 0, 0, 0, 0, -t2007, -t2009, -t1993, -t1985; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2151, t2152, 0, t2116, 0, 0, 0, 0, 0, 0, t2276, t2085, t2274, t2044, 0, 0, 0, 0, 0, 0, t2276, t2274, -t2085, t2021, 0, 0, 0, 0, 0, 0, t2029, t2032, t2019, t1995, 0, 0, 0, 0, 0, 0, t1999, t2003, t1989, t1983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2192, t2186, t2234, t2156, 0, 0, 0, 0, 0, 0, t2268, t2112, t2267, t2057, 0, 0, 0, 0, 0, 0, t2268, t2267, -t2112, t2035, 0, 0, 0, 0, 0, 0, t2041, t2047, t2028, t1998, 0, 0, 0, 0, 0, 0, t2008, t2010, t1994, t1986; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2185, t2191, -t2235, t2155, 0, 0, 0, 0, 0, 0, -t2129, -t2128, t2125, -t2120, 0, 0, 0, 0, 0, 0, -t2129, t2125, t2128, -t2068, 0, 0, 0, 0, 0, 0, -t2075, -t2077, t2089, -t2060, 0, 0, 0, 0, 0, 0, t2036, t2038, t2059, -t2031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2180, t2178, -t2184, -t2168, 0, 0, 0, 0, 0, 0, t2270, t2110, t2269, t2056, 0, 0, 0, 0, 0, 0, t2270, t2269, -t2110, t2034, 0, 0, 0, 0, 0, 0, t2040, t2046, t2027, t1997, 0, 0, 0, 0, 0, 0, t2007, t2009, t1993, t1985; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2263, t2133, t2265, t2084, 0, 0, 0, 0, 0, 0, t2263, t2265, -t2133, t2061, 0, 0, 0, 0, 0, 0, t2066, t2070, t2055, t2012, 0, 0, 0, 0, 0, 0, t2018, t2026, t2006, t1988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2134, t2264, -t2132, t2083, 0, 0, 0, 0, 0, 0, -t2134, -t2132, -t2264, -t2062, 0, 0, 0, 0, 0, 0, -t2065, -t2069, -t2054, -t2011, 0, 0, 0, 0, 0, 0, -t2017, -t2025, -t2005, -t1987; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2129, t2128, -t2125, t2120, 0, 0, 0, 0, 0, 0, t2129, -t2125, -t2128, t2068, 0, 0, 0, 0, 0, 0, t2075, t2077, -t2089, t2060, 0, 0, 0, 0, 0, 0, -t2036, -t2038, -t2059, t2031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2263, t2265, -t2133, t2061, 0, 0, 0, 0, 0, 0, t2066, t2070, t2055, t2012, 0, 0, 0, 0, 0, 0, t2018, t2026, t2006, t1988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2129, -t2125, -t2128, t2068, 0, 0, 0, 0, 0, 0, t2075, t2077, -t2089, t2060, 0, 0, 0, 0, 0, 0, -t2036, -t2038, -t2059, t2031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2134, t2132, t2264, t2062, 0, 0, 0, 0, 0, 0, t2065, t2069, t2054, t2011, 0, 0, 0, 0, 0, 0, t2017, t2025, t2005, t1987; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2106, t2100, t2076, t2024, 0, 0, 0, 0, 0, 0, t2043, t2052, t2016, t1992; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2099, t2119, t2078, t2023, 0, 0, 0, 0, 0, 0, t2042, t2051, t2015, t1991; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2075, -t2077, t2089, -t2060, 0, 0, 0, 0, 0, 0, t2036, t2038, t2059, -t2031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2067, t2064, t2037, t2002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2063, t2090, t2039, t2001; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2036, t2038, t2059, -t2031;];
f_new_reg  = t1;
