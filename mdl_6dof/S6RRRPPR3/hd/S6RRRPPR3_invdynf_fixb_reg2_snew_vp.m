% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRPPR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:34:14
% EndTime: 2019-05-07 04:34:20
% DurationCPUTime: 7.24s
% Computational Cost: add. (14040->317), mult. (29757->319), div. (0->0), fcn. (20332->8), ass. (0->196)
t2213 = cos(qJ(1));
t2207 = sin(qJ(3));
t2208 = sin(qJ(2));
t2211 = cos(qJ(3));
t2212 = cos(qJ(2));
t2167 = (t2207 * t2212 + t2208 * t2211) * qJD(1);
t2258 = qJD(1) * t2212;
t2239 = qJD(2) * t2258;
t2241 = t2208 * qJDD(1);
t2174 = t2239 + t2241;
t2197 = t2212 * qJDD(1);
t2259 = qJD(1) * t2208;
t2240 = qJD(2) * t2259;
t2175 = t2197 - t2240;
t2234 = t2207 * t2174 - t2211 * t2175;
t2226 = -qJD(3) * t2167 - t2234;
t2202 = qJD(2) + qJD(3);
t2251 = t2202 * t2167;
t2276 = -t2226 + t2251;
t2263 = t2202 ^ 2;
t2165 = t2207 * t2259 - t2211 * t2258;
t2264 = t2165 ^ 2;
t2268 = -t2264 - t2263;
t2201 = qJDD(2) + qJDD(3);
t2252 = t2167 * t2165;
t2277 = t2201 - t2252;
t2108 = t2207 * t2277 - t2268 * t2211;
t2279 = t2268 * t2207 + t2211 * t2277;
t2074 = t2108 * t2212 + t2208 * t2279;
t2209 = sin(qJ(1));
t2289 = t2074 * t2209;
t2294 = t2213 * t2276 + t2289;
t2288 = t2074 * t2213;
t2293 = t2209 * t2276 - t2288;
t2151 = t2202 * t2165;
t2230 = -t2211 * t2174 - t2207 * t2175;
t2224 = qJD(3) * t2165 + t2230;
t2102 = t2224 + t2151;
t2164 = t2167 ^ 2;
t2269 = -t2164 - t2263;
t2278 = -t2201 - t2252;
t2103 = t2207 * t2278 + t2211 * t2269;
t2106 = -t2207 * t2269 + t2211 * t2278;
t2290 = t2103 * t2208 - t2106 * t2212;
t2292 = -t2102 * t2213 + t2209 * t2290;
t2291 = t2102 * t2209 + t2213 * t2290;
t2069 = t2103 * t2212 + t2106 * t2208;
t2071 = t2108 * t2208 - t2212 * t2279;
t2247 = qJD(3) - t2202;
t2095 = t2167 * t2247 + t2234;
t2283 = t2095 * t2207;
t2282 = t2095 * t2211;
t2122 = t2164 + t2264;
t2281 = t2122 * t2209;
t2280 = t2122 * t2213;
t2270 = t2224 - t2151;
t2262 = 2 * qJD(4);
t2206 = sin(qJ(6));
t2210 = cos(qJ(6));
t2137 = t2165 * t2206 + t2210 * t2202;
t2267 = t2137 ^ 2;
t2139 = t2165 * t2210 - t2202 * t2206;
t2266 = t2139 ^ 2;
t2162 = qJD(6) + t2167;
t2265 = t2162 ^ 2;
t2261 = t2226 * pkin(3);
t2260 = t2201 * qJ(4);
t2101 = t2165 * t2247 + t2230;
t2257 = t2101 * t2207;
t2255 = t2137 * t2139;
t2205 = t2212 ^ 2;
t2215 = qJD(1) ^ 2;
t2250 = t2205 * t2215;
t2184 = -g(1) * t2213 - g(2) * t2209;
t2170 = -pkin(1) * t2215 + qJDD(1) * pkin(7) + t2184;
t2249 = t2208 * t2170;
t2248 = t2208 * t2215;
t2246 = -pkin(4) * t2202 - qJ(5) * t2167 + t2262;
t2133 = pkin(3) * t2165 - qJ(4) * t2167;
t2245 = (2 * qJD(5)) - t2133;
t2244 = qJD(6) - t2162;
t2243 = qJD(6) + t2162;
t2147 = -t2208 * g(3) + t2212 * t2170;
t2228 = qJD(2) * pkin(2) - pkin(8) * t2259;
t2115 = -pkin(2) * t2250 + t2175 * pkin(8) - qJD(2) * t2228 + t2147;
t2218 = qJDD(2) * pkin(2) - t2174 * pkin(8) - t2249 + (pkin(8) * qJD(1) * qJD(2) + pkin(2) * t2248 - g(3)) * t2212;
t2080 = t2211 * t2115 + t2207 * t2218;
t2204 = t2208 ^ 2;
t2242 = t2204 + t2205;
t2238 = -pkin(3) * t2263 + t2080;
t2237 = -t2265 - t2266;
t2183 = t2209 * g(1) - t2213 * g(2);
t2079 = -t2207 * t2115 + t2211 * t2218;
t2236 = -t2210 * t2201 + t2206 * t2226;
t2233 = pkin(4) * t2165 - t2245;
t2227 = qJDD(1) * pkin(1) + t2183;
t2123 = t2175 * pkin(2) - t2228 * t2259 + (pkin(8) * t2205 + pkin(7)) * t2215 + t2227;
t2217 = -pkin(3) * t2251 - qJ(4) * t2102 + t2123;
t2216 = pkin(4) * t2226 - qJ(5) * t2264 + qJDD(5) + t2217;
t2026 = -t2102 * pkin(5) - (-pkin(3) - pkin(9)) * t2226 + (-pkin(9) * t2202 + t2246) * t2167 + t2216;
t2134 = pkin(5) * t2167 - pkin(9) * t2165;
t2225 = -t2201 * pkin(3) - qJ(4) * t2263 + qJDD(4) - t2079;
t2221 = -t2201 * pkin(4) + qJ(5) * t2270 + t2225;
t2029 = -t2263 * pkin(5) - t2201 * pkin(9) + (-t2134 + t2233) * t2167 + t2221;
t2232 = t2026 * t2210 - t2029 * t2206;
t2231 = t2206 * t2201 + t2210 * t2226;
t2223 = qJDD(6) - t2224;
t2222 = t2137 * t2244 + t2231;
t2220 = -pkin(4) * t2264 - qJ(5) * t2226 + t2202 * t2246 + t2238;
t2219 = t2223 - t2255;
t2214 = qJD(2) ^ 2;
t2188 = t2212 * t2248;
t2186 = -t2214 - t2250;
t2185 = -t2204 * t2215 - t2214;
t2182 = -qJDD(2) + t2188;
t2181 = qJDD(2) + t2188;
t2180 = t2242 * t2215;
t2179 = -qJDD(1) * t2209 - t2213 * t2215;
t2178 = qJDD(1) * t2213 - t2209 * t2215;
t2177 = t2242 * qJDD(1);
t2176 = t2197 - 0.2e1 * t2240;
t2173 = 0.2e1 * t2239 + t2241;
t2169 = t2215 * pkin(7) + t2227;
t2146 = -t2212 * g(3) - t2249;
t2143 = t2182 * t2212 - t2185 * t2208;
t2142 = -t2181 * t2208 + t2186 * t2212;
t2141 = t2182 * t2208 + t2185 * t2212;
t2140 = t2181 * t2212 + t2186 * t2208;
t2114 = -t2146 * t2208 + t2147 * t2212;
t2113 = t2146 * t2212 + t2147 * t2208;
t2097 = t2251 + t2226;
t2093 = (qJD(3) + t2202) * t2167 + t2234;
t2091 = t2211 * t2101;
t2090 = -t2265 - t2267;
t2083 = -t2266 - t2267;
t2082 = -t2223 - t2255;
t2078 = -t2137 * t2243 - t2231;
t2077 = -t2139 * t2244 + t2236;
t2076 = t2139 * t2243 - t2236;
t2068 = t2097 * t2211 - t2257;
t2067 = -t2257 - t2282;
t2066 = t2207 * t2270 + t2282;
t2065 = t2097 * t2207 + t2091;
t2064 = t2091 - t2283;
t2063 = -t2211 * t2270 + t2283;
t2056 = t2082 * t2210 - t2206 * t2237;
t2055 = t2206 * t2082 + t2210 * t2237;
t2054 = t2133 * t2167 + t2225;
t2053 = t2090 * t2210 - t2206 * t2219;
t2052 = t2206 * t2090 + t2210 * t2219;
t2051 = -t2165 * t2133 + t2202 * t2262 + t2238 + t2260;
t2050 = t2167 * t2262 + t2217 + t2261;
t2049 = -t2079 * t2207 + t2080 * t2211;
t2048 = t2079 * t2211 + t2080 * t2207;
t2047 = t2077 * t2210 - t2206 * t2222;
t2046 = t2206 * t2077 + t2210 * t2222;
t2045 = -t2065 * t2208 + t2068 * t2212;
t2044 = -t2064 * t2208 + t2067 * t2212;
t2043 = -t2063 * t2208 + t2066 * t2212;
t2042 = t2065 * t2212 + t2068 * t2208;
t2041 = t2064 * t2212 + t2067 * t2208;
t2040 = t2063 * t2212 + t2066 * t2208;
t2039 = t2165 * t2245 + t2220 + t2260;
t2038 = t2167 * t2233 + t2221;
t2037 = t2056 * t2207 + t2078 * t2211;
t2036 = -t2056 * t2211 + t2078 * t2207;
t2035 = t2246 * t2167 + t2216 + t2261;
t2034 = t2053 * t2207 + t2076 * t2211;
t2033 = -t2053 * t2211 + t2076 * t2207;
t2032 = -t2263 * pkin(9) + (pkin(5) + qJ(4)) * t2201 + (t2134 + t2245) * t2165 + t2220;
t2031 = t2047 * t2207 + t2083 * t2211;
t2030 = -t2047 * t2211 + t2083 * t2207;
t2028 = t2051 * t2211 + t2054 * t2207;
t2027 = t2051 * t2207 - t2054 * t2211;
t2025 = -t2048 * t2208 + t2049 * t2212;
t2024 = t2048 * t2212 + t2049 * t2208;
t2023 = t2038 * t2207 + t2039 * t2211;
t2022 = -t2038 * t2211 + t2039 * t2207;
t2021 = -t2036 * t2208 + t2037 * t2212;
t2020 = t2036 * t2212 + t2037 * t2208;
t2019 = -t2033 * t2208 + t2034 * t2212;
t2018 = t2033 * t2212 + t2034 * t2208;
t2017 = -t2030 * t2208 + t2031 * t2212;
t2016 = t2030 * t2212 + t2031 * t2208;
t2015 = -t2027 * t2208 + t2028 * t2212;
t2014 = t2027 * t2212 + t2028 * t2208;
t2013 = t2026 * t2206 + t2029 * t2210;
t2011 = -t2022 * t2208 + t2023 * t2212;
t2010 = t2022 * t2212 + t2023 * t2208;
t2009 = t2013 * t2210 - t2206 * t2232;
t2008 = t2206 * t2013 + t2210 * t2232;
t2007 = t2009 * t2207 + t2032 * t2211;
t2006 = -t2009 * t2211 + t2032 * t2207;
t2005 = -t2006 * t2208 + t2007 * t2212;
t2004 = t2006 * t2212 + t2007 * t2208;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2179, -t2178, 0, -t2183 * t2209 + t2184 * t2213, 0, 0, 0, 0, 0, 0, t2142 * t2213 - t2176 * t2209, t2143 * t2213 + t2173 * t2209, t2177 * t2213 - t2180 * t2209, t2114 * t2213 - t2169 * t2209, 0, 0, 0, 0, 0, 0, t2093 * t2209 - t2288, -t2291, t2045 * t2213 - t2281, t2025 * t2213 - t2123 * t2209, 0, 0, 0, 0, 0, 0, t2293, t2044 * t2213 - t2281, t2291, t2015 * t2213 - t2050 * t2209, 0, 0, 0, 0, 0, 0, t2291, -t2293, t2043 * t2213 + t2281, t2011 * t2213 - t2035 * t2209, 0, 0, 0, 0, 0, 0, t2019 * t2213 - t2052 * t2209, t2021 * t2213 - t2055 * t2209, t2017 * t2213 - t2046 * t2209, t2005 * t2213 - t2008 * t2209; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2178, t2179, 0, t2183 * t2213 + t2184 * t2209, 0, 0, 0, 0, 0, 0, t2142 * t2209 + t2176 * t2213, t2143 * t2209 - t2173 * t2213, t2177 * t2209 + t2180 * t2213, t2114 * t2209 + t2169 * t2213, 0, 0, 0, 0, 0, 0, -t2093 * t2213 - t2289, -t2292, t2045 * t2209 + t2280, t2025 * t2209 + t2123 * t2213, 0, 0, 0, 0, 0, 0, -t2294, t2044 * t2209 + t2280, t2292, t2015 * t2209 + t2050 * t2213, 0, 0, 0, 0, 0, 0, t2292, t2294, t2043 * t2209 - t2280, t2011 * t2209 + t2035 * t2213, 0, 0, 0, 0, 0, 0, t2019 * t2209 + t2052 * t2213, t2021 * t2209 + t2055 * t2213, t2017 * t2209 + t2046 * t2213, t2005 * t2209 + t2008 * t2213; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2140, t2141, 0, t2113, 0, 0, 0, 0, 0, 0, -t2071, t2069, t2042, t2024, 0, 0, 0, 0, 0, 0, -t2071, t2041, -t2069, t2014, 0, 0, 0, 0, 0, 0, -t2069, t2071, t2040, t2010, 0, 0, 0, 0, 0, 0, t2018, t2020, t2016, t2004; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2215, -qJDD(1), 0, t2184, 0, 0, 0, 0, 0, 0, t2142, t2143, t2177, t2114, 0, 0, 0, 0, 0, 0, -t2074, -t2290, t2045, t2025, 0, 0, 0, 0, 0, 0, -t2074, t2044, t2290, t2015, 0, 0, 0, 0, 0, 0, t2290, t2074, t2043, t2011, 0, 0, 0, 0, 0, 0, t2019, t2021, t2017, t2005; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2215, 0, t2183, 0, 0, 0, 0, 0, 0, t2176, -t2173, t2180, t2169, 0, 0, 0, 0, 0, 0, -t2093, t2102, t2122, t2123, 0, 0, 0, 0, 0, 0, -t2276, t2122, -t2102, t2050, 0, 0, 0, 0, 0, 0, -t2102, t2276, -t2122, t2035, 0, 0, 0, 0, 0, 0, t2052, t2055, t2046, t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2140, t2141, 0, t2113, 0, 0, 0, 0, 0, 0, -t2071, t2069, t2042, t2024, 0, 0, 0, 0, 0, 0, -t2071, t2041, -t2069, t2014, 0, 0, 0, 0, 0, 0, -t2069, t2071, t2040, t2010, 0, 0, 0, 0, 0, 0, t2018, t2020, t2016, t2004; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2186, t2182, t2197, t2147, 0, 0, 0, 0, 0, 0, -t2108, t2106, t2068, t2049, 0, 0, 0, 0, 0, 0, -t2108, t2067, -t2106, t2028, 0, 0, 0, 0, 0, 0, -t2106, t2108, t2066, t2023, 0, 0, 0, 0, 0, 0, t2034, t2037, t2031, t2007; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2181, t2185, -t2241, t2146, 0, 0, 0, 0, 0, 0, t2279, t2103, t2065, t2048, 0, 0, 0, 0, 0, 0, t2279, t2064, -t2103, t2027, 0, 0, 0, 0, 0, 0, -t2103, -t2279, t2063, t2022, 0, 0, 0, 0, 0, 0, t2033, t2036, t2030, t2006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2176, t2173, -t2180, -t2169, 0, 0, 0, 0, 0, 0, t2093, -t2102, -t2122, -t2123, 0, 0, 0, 0, 0, 0, t2276, -t2122, t2102, -t2050, 0, 0, 0, 0, 0, 0, t2102, -t2276, t2122, -t2035, 0, 0, 0, 0, 0, 0, -t2052, -t2055, -t2046, -t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2268, t2278, t2097, t2080, 0, 0, 0, 0, 0, 0, t2268, -t2095, -t2278, t2051, 0, 0, 0, 0, 0, 0, -t2278, -t2268, t2095, t2039, 0, 0, 0, 0, 0, 0, t2076, t2078, t2083, t2032; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2277, t2269, t2101, t2079, 0, 0, 0, 0, 0, 0, t2277, t2101, -t2269, -t2054, 0, 0, 0, 0, 0, 0, -t2269, -t2277, -t2270, -t2038, 0, 0, 0, 0, 0, 0, -t2053, -t2056, -t2047, -t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2093, -t2102, -t2122, -t2123, 0, 0, 0, 0, 0, 0, t2276, -t2122, t2102, -t2050, 0, 0, 0, 0, 0, 0, t2102, -t2276, t2122, -t2035, 0, 0, 0, 0, 0, 0, -t2052, -t2055, -t2046, -t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2268, -t2095, -t2278, t2051, 0, 0, 0, 0, 0, 0, -t2278, -t2268, t2095, t2039, 0, 0, 0, 0, 0, 0, t2076, t2078, t2083, t2032; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2276, -t2122, t2102, -t2050, 0, 0, 0, 0, 0, 0, t2102, -t2276, t2122, -t2035, 0, 0, 0, 0, 0, 0, -t2052, -t2055, -t2046, -t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2277, -t2101, t2269, t2054, 0, 0, 0, 0, 0, 0, t2269, t2277, t2270, t2038, 0, 0, 0, 0, 0, 0, t2053, t2056, t2047, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2269, t2277, t2270, t2038, 0, 0, 0, 0, 0, 0, t2053, t2056, t2047, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2278, t2268, -t2095, -t2039, 0, 0, 0, 0, 0, 0, -t2076, -t2078, -t2083, -t2032; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2102, t2276, -t2122, t2035, 0, 0, 0, 0, 0, 0, t2052, t2055, t2046, t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2090, t2082, t2077, t2013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2219, t2237, t2222, t2232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2076, t2078, t2083, t2032;];
f_new_reg  = t1;