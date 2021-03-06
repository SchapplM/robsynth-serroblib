% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6PRRRRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:28:15
% EndTime: 2019-05-05 10:28:24
% DurationCPUTime: 10.19s
% Computational Cost: add. (80294->329), mult. (168101->510), div. (0->0), fcn. (129834->14), ass. (0->283)
t2258 = sin(pkin(6));
t2260 = cos(pkin(6));
t2257 = sin(pkin(12));
t2259 = cos(pkin(12));
t2300 = g(1) * t2257 - g(2) * t2259;
t2322 = -g(3) + qJDD(1);
t2334 = t2258 * t2322 + t2260 * t2300;
t2264 = sin(qJ(4));
t2265 = sin(qJ(3));
t2269 = cos(qJ(4));
t2270 = cos(qJ(3));
t2225 = (-t2264 * t2265 + t2269 * t2270) * qJD(2);
t2227 = (t2264 * t2270 + t2265 * t2269) * qJD(2);
t2263 = sin(qJ(5));
t2268 = cos(qJ(5));
t2198 = -t2268 * t2225 + t2227 * t2263;
t2195 = qJD(6) + t2198;
t2333 = qJD(6) + t2195;
t2200 = t2225 * t2263 + t2227 * t2268;
t2254 = qJD(3) + qJD(4);
t2253 = qJD(5) + t2254;
t2262 = sin(qJ(6));
t2267 = cos(qJ(6));
t2181 = t2200 * t2262 - t2267 * t2253;
t2332 = t2181 ^ 2;
t2183 = t2200 * t2267 + t2253 * t2262;
t2331 = t2183 ^ 2;
t2330 = t2195 ^ 2;
t2329 = t2198 ^ 2;
t2328 = t2200 ^ 2;
t2327 = t2225 ^ 2;
t2326 = t2227 ^ 2;
t2325 = t2253 ^ 2;
t2324 = t2254 ^ 2;
t2323 = t2270 ^ 2;
t2321 = qJD(2) * t2265;
t2320 = t2181 * t2183;
t2319 = t2198 * t2200;
t2318 = t2225 * t2227;
t2317 = t2254 * t2225;
t2273 = qJD(2) ^ 2;
t2316 = t2323 * t2273;
t2315 = qJD(4) - t2254;
t2314 = qJD(5) - t2253;
t2313 = qJD(6) - t2195;
t2238 = -g(1) * t2259 - g(2) * t2257;
t2266 = sin(qJ(2));
t2271 = cos(qJ(2));
t2197 = t2271 * t2238 + t2266 * t2334;
t2191 = -pkin(2) * t2273 + qJDD(2) * pkin(8) + t2197;
t2215 = -t2258 * t2300 + t2260 * t2322;
t2172 = -t2265 * t2191 + t2270 * t2215;
t2307 = t2270 * qJD(2) * qJD(3);
t2312 = t2265 * qJDD(2);
t2232 = t2307 + t2312;
t2248 = t2265 * t2273 * t2270;
t2239 = qJDD(3) + t2248;
t2157 = (-t2232 + t2307) * pkin(9) + t2239 * pkin(3) + t2172;
t2173 = t2270 * t2191 + t2265 * t2215;
t2243 = qJD(3) * pkin(3) - pkin(9) * t2321;
t2308 = qJD(3) * t2321;
t2311 = t2270 * qJDD(2);
t2278 = -t2308 + t2311;
t2158 = -pkin(3) * t2316 + pkin(9) * t2278 - qJD(3) * t2243 + t2173;
t2128 = t2264 * t2157 + t2269 * t2158;
t2301 = t2264 * t2232 - t2269 * t2278;
t2187 = -qJD(4) * t2227 - t2301;
t2214 = pkin(4) * t2254 - pkin(10) * t2227;
t2115 = -pkin(4) * t2327 + pkin(10) * t2187 - t2214 * t2254 + t2128;
t2127 = t2269 * t2157 - t2264 * t2158;
t2275 = -t2269 * t2232 - t2264 * t2278;
t2188 = t2225 * qJD(4) - t2275;
t2310 = qJDD(3) + qJDD(4);
t2202 = t2310 + t2318;
t2274 = (-t2188 + t2317) * pkin(10) + t2202 * pkin(4) + t2127;
t2089 = t2268 * t2115 + t2263 * t2274;
t2255 = t2265 ^ 2;
t2309 = t2255 + t2323;
t2305 = qJDD(5) + t2310;
t2088 = -t2115 * t2263 + t2268 * t2274;
t2285 = -t2263 * t2187 - t2268 * t2188;
t2143 = -qJD(5) * t2198 - t2285;
t2304 = t2253 * t2198 - t2143;
t2303 = -t2262 * t2143 + t2267 * t2305;
t2302 = -t2268 * t2187 + t2263 * t2188;
t2299 = t2266 * t2238 - t2271 * t2334;
t2164 = pkin(5) * t2198 - pkin(11) * t2200;
t2075 = -pkin(5) * t2325 + pkin(11) * t2305 - t2198 * t2164 + t2089;
t2131 = (qJD(5) + t2253) * t2200 + t2302;
t2190 = -qJDD(2) * pkin(2) - t2273 * pkin(8) + t2299;
t2167 = -t2278 * pkin(3) - pkin(9) * t2316 + t2243 * t2321 + t2190;
t2135 = -t2187 * pkin(4) - t2327 * pkin(10) + t2227 * t2214 + t2167;
t2094 = pkin(5) * t2131 + pkin(11) * t2304 + t2135;
t2065 = -t2075 * t2262 + t2094 * t2267;
t2066 = t2075 * t2267 + t2094 * t2262;
t2047 = -t2065 * t2262 + t2066 * t2267;
t2074 = -pkin(5) * t2305 - pkin(11) * t2325 + t2164 * t2200 - t2088;
t2034 = t2047 * t2263 - t2074 * t2268;
t2035 = t2047 * t2268 + t2074 * t2263;
t2022 = t2034 * t2269 + t2035 * t2264;
t2023 = -t2034 * t2264 + t2035 * t2269;
t2019 = -t2022 * t2265 + t2023 * t2270;
t2046 = t2065 * t2267 + t2066 * t2262;
t2298 = t2019 * t2266 - t2046 * t2271;
t2067 = t2088 * t2268 + t2089 * t2263;
t2068 = -t2088 * t2263 + t2089 * t2268;
t2048 = t2067 * t2269 + t2068 * t2264;
t2049 = -t2067 * t2264 + t2068 * t2269;
t2032 = -t2048 * t2265 + t2049 * t2270;
t2297 = t2032 * t2266 - t2135 * t2271;
t2117 = -t2183 * t2313 + t2303;
t2276 = -t2267 * t2143 - t2262 * t2305;
t2119 = t2181 * t2313 + t2276;
t2091 = t2117 * t2267 - t2119 * t2262;
t2142 = -t2331 - t2332;
t2078 = t2091 * t2263 - t2142 * t2268;
t2079 = t2091 * t2268 + t2142 * t2263;
t2057 = t2078 * t2269 + t2079 * t2264;
t2058 = -t2078 * t2264 + t2079 * t2269;
t2041 = -t2057 * t2265 + t2058 * t2270;
t2090 = t2117 * t2262 + t2119 * t2267;
t2296 = t2041 * t2266 - t2090 * t2271;
t2277 = -qJD(5) * t2200 - qJDD(6) - t2302;
t2124 = -t2277 - t2320;
t2146 = -t2330 - t2332;
t2100 = -t2124 * t2262 + t2146 * t2267;
t2116 = t2183 * t2333 - t2303;
t2082 = t2100 * t2263 - t2116 * t2268;
t2083 = t2100 * t2268 + t2116 * t2263;
t2059 = t2082 * t2269 + t2083 * t2264;
t2060 = -t2082 * t2264 + t2083 * t2269;
t2043 = -t2059 * t2265 + t2060 * t2270;
t2099 = t2124 * t2267 + t2146 * t2262;
t2295 = t2043 * t2266 - t2099 * t2271;
t2125 = t2277 - t2320;
t2151 = -t2330 - t2331;
t2102 = t2125 * t2267 - t2151 * t2262;
t2118 = -t2181 * t2333 - t2276;
t2084 = t2102 * t2263 - t2118 * t2268;
t2085 = t2102 * t2268 + t2118 * t2263;
t2061 = t2084 * t2269 + t2085 * t2264;
t2062 = -t2084 * t2264 + t2085 * t2269;
t2045 = -t2061 * t2265 + t2062 * t2270;
t2101 = t2125 * t2262 + t2151 * t2267;
t2294 = t2045 * t2266 - t2101 * t2271;
t2132 = -t2200 * t2314 - t2302;
t2134 = t2198 * t2314 + t2285;
t2103 = t2132 * t2263 + t2134 * t2268;
t2104 = t2132 * t2268 - t2134 * t2263;
t2076 = t2103 * t2269 + t2104 * t2264;
t2077 = -t2103 * t2264 + t2104 * t2269;
t2056 = -t2076 * t2265 + t2077 * t2270;
t2154 = -t2328 - t2329;
t2293 = t2056 * t2266 - t2154 * t2271;
t2095 = t2127 * t2269 + t2128 * t2264;
t2096 = -t2127 * t2264 + t2128 * t2269;
t2072 = -t2095 * t2265 + t2096 * t2270;
t2292 = t2072 * t2266 - t2167 * t2271;
t2161 = -t2325 - t2329;
t2162 = t2305 - t2319;
t2136 = t2161 * t2263 + t2162 * t2268;
t2137 = t2161 * t2268 - t2162 * t2263;
t2109 = t2136 * t2269 + t2137 * t2264;
t2110 = -t2136 * t2264 + t2137 * t2269;
t2081 = -t2109 * t2265 + t2110 * t2270;
t2291 = t2081 * t2266 - t2131 * t2271;
t2163 = -t2305 - t2319;
t2184 = -t2325 - t2328;
t2149 = t2163 * t2263 + t2184 * t2268;
t2150 = t2163 * t2268 - t2184 * t2263;
t2122 = t2149 * t2269 + t2150 * t2264;
t2123 = -t2149 * t2264 + t2150 * t2269;
t2093 = -t2122 * t2265 + t2123 * t2270;
t2290 = t2093 * t2266 + t2271 * t2304;
t2175 = -t2227 * t2315 - t2301;
t2177 = -t2225 * t2315 + t2275;
t2144 = t2175 * t2264 + t2177 * t2269;
t2145 = t2175 * t2269 - t2177 * t2264;
t2121 = -t2144 * t2265 + t2145 * t2270;
t2189 = -t2326 - t2327;
t2289 = t2121 * t2266 - t2189 * t2271;
t2201 = -t2324 - t2327;
t2165 = t2201 * t2264 + t2202 * t2269;
t2166 = t2201 * t2269 - t2202 * t2264;
t2139 = -t2165 * t2265 + t2166 * t2270;
t2174 = (qJD(4) + t2254) * t2227 + t2301;
t2288 = t2139 * t2266 - t2174 * t2271;
t2141 = -t2172 * t2265 + t2173 * t2270;
t2287 = t2141 * t2266 - t2190 * t2271;
t2203 = -t2310 + t2318;
t2211 = -t2324 - t2326;
t2178 = t2203 * t2264 + t2211 * t2269;
t2179 = t2203 * t2269 - t2211 * t2264;
t2148 = -t2178 * t2265 + t2179 * t2270;
t2176 = t2188 + t2317;
t2286 = t2148 * t2266 - t2176 * t2271;
t2284 = t2197 * t2266 - t2271 * t2299;
t2272 = qJD(3) ^ 2;
t2247 = -t2272 - t2316;
t2209 = -t2239 * t2265 + t2247 * t2270;
t2233 = -0.2e1 * t2308 + t2311;
t2283 = t2209 * t2266 + t2233 * t2271;
t2240 = -qJDD(3) + t2248;
t2246 = -t2255 * t2273 - t2272;
t2210 = t2240 * t2270 - t2246 * t2265;
t2231 = 0.2e1 * t2307 + t2312;
t2282 = t2210 * t2266 - t2231 * t2271;
t2234 = t2309 * qJDD(2);
t2237 = t2309 * t2273;
t2281 = t2234 * t2266 + t2237 * t2271;
t2280 = qJDD(2) * t2271 - t2266 * t2273;
t2236 = -qJDD(2) * t2266 - t2271 * t2273;
t2222 = t2280 * t2260;
t2221 = t2236 * t2260;
t2220 = t2280 * t2258;
t2219 = t2236 * t2258;
t2208 = t2240 * t2265 + t2246 * t2270;
t2207 = t2239 * t2270 + t2247 * t2265;
t2206 = t2234 * t2271 - t2237 * t2266;
t2205 = t2281 * t2260;
t2204 = t2281 * t2258;
t2194 = t2210 * t2271 + t2231 * t2266;
t2193 = t2209 * t2271 - t2233 * t2266;
t2171 = -t2258 * t2208 + t2260 * t2282;
t2170 = -t2258 * t2207 + t2260 * t2283;
t2169 = t2260 * t2208 + t2258 * t2282;
t2168 = t2260 * t2207 + t2258 * t2283;
t2160 = t2197 * t2271 + t2266 * t2299;
t2153 = -t2258 * t2215 + t2260 * t2284;
t2152 = t2260 * t2215 + t2258 * t2284;
t2147 = t2178 * t2270 + t2179 * t2265;
t2140 = t2172 * t2270 + t2173 * t2265;
t2138 = t2165 * t2270 + t2166 * t2265;
t2130 = t2141 * t2271 + t2190 * t2266;
t2129 = t2148 * t2271 + t2176 * t2266;
t2126 = t2139 * t2271 + t2174 * t2266;
t2120 = t2144 * t2270 + t2145 * t2265;
t2111 = t2121 * t2271 + t2189 * t2266;
t2108 = -t2258 * t2147 + t2260 * t2286;
t2107 = t2260 * t2147 + t2258 * t2286;
t2106 = -t2258 * t2140 + t2260 * t2287;
t2105 = t2260 * t2140 + t2258 * t2287;
t2098 = -t2258 * t2138 + t2260 * t2288;
t2097 = t2260 * t2138 + t2258 * t2288;
t2092 = t2122 * t2270 + t2123 * t2265;
t2087 = -t2258 * t2120 + t2260 * t2289;
t2086 = t2260 * t2120 + t2258 * t2289;
t2080 = t2109 * t2270 + t2110 * t2265;
t2073 = t2093 * t2271 - t2266 * t2304;
t2071 = t2095 * t2270 + t2096 * t2265;
t2070 = t2081 * t2271 + t2131 * t2266;
t2069 = t2072 * t2271 + t2167 * t2266;
t2064 = -t2258 * t2092 + t2260 * t2290;
t2063 = t2260 * t2092 + t2258 * t2290;
t2055 = t2076 * t2270 + t2077 * t2265;
t2054 = -t2258 * t2080 + t2260 * t2291;
t2053 = t2260 * t2080 + t2258 * t2291;
t2052 = t2056 * t2271 + t2154 * t2266;
t2051 = -t2258 * t2071 + t2260 * t2292;
t2050 = t2260 * t2071 + t2258 * t2292;
t2044 = t2061 * t2270 + t2062 * t2265;
t2042 = t2059 * t2270 + t2060 * t2265;
t2040 = t2057 * t2270 + t2058 * t2265;
t2039 = t2045 * t2271 + t2101 * t2266;
t2038 = -t2258 * t2055 + t2260 * t2293;
t2037 = t2260 * t2055 + t2258 * t2293;
t2036 = t2043 * t2271 + t2099 * t2266;
t2033 = t2041 * t2271 + t2090 * t2266;
t2031 = t2048 * t2270 + t2049 * t2265;
t2030 = t2032 * t2271 + t2135 * t2266;
t2029 = -t2258 * t2044 + t2260 * t2294;
t2028 = t2260 * t2044 + t2258 * t2294;
t2027 = -t2258 * t2042 + t2260 * t2295;
t2026 = t2260 * t2042 + t2258 * t2295;
t2025 = -t2258 * t2040 + t2260 * t2296;
t2024 = t2260 * t2040 + t2258 * t2296;
t2021 = -t2258 * t2031 + t2260 * t2297;
t2020 = t2260 * t2031 + t2258 * t2297;
t2018 = t2022 * t2270 + t2023 * t2265;
t2017 = t2019 * t2271 + t2046 * t2266;
t2016 = -t2258 * t2018 + t2260 * t2298;
t2015 = t2260 * t2018 + t2258 * t2298;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2259 * t2238 - t2257 * t2300, 0, 0, 0, 0, 0, 0, -t2222 * t2257 + t2236 * t2259, -t2221 * t2257 - t2259 * t2280, 0, -t2153 * t2257 + t2160 * t2259, 0, 0, 0, 0, 0, 0, -t2170 * t2257 + t2193 * t2259, -t2171 * t2257 + t2194 * t2259, -t2205 * t2257 + t2206 * t2259, -t2106 * t2257 + t2130 * t2259, 0, 0, 0, 0, 0, 0, -t2098 * t2257 + t2126 * t2259, -t2108 * t2257 + t2129 * t2259, -t2087 * t2257 + t2111 * t2259, -t2051 * t2257 + t2069 * t2259, 0, 0, 0, 0, 0, 0, -t2054 * t2257 + t2070 * t2259, -t2064 * t2257 + t2073 * t2259, -t2038 * t2257 + t2052 * t2259, -t2021 * t2257 + t2030 * t2259, 0, 0, 0, 0, 0, 0, -t2027 * t2257 + t2036 * t2259, -t2029 * t2257 + t2039 * t2259, -t2025 * t2257 + t2033 * t2259, -t2016 * t2257 + t2017 * t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2257 * t2238 + t2259 * t2300, 0, 0, 0, 0, 0, 0, t2222 * t2259 + t2236 * t2257, t2221 * t2259 - t2257 * t2280, 0, t2153 * t2259 + t2160 * t2257, 0, 0, 0, 0, 0, 0, t2170 * t2259 + t2193 * t2257, t2171 * t2259 + t2194 * t2257, t2205 * t2259 + t2206 * t2257, t2106 * t2259 + t2130 * t2257, 0, 0, 0, 0, 0, 0, t2098 * t2259 + t2126 * t2257, t2108 * t2259 + t2129 * t2257, t2087 * t2259 + t2111 * t2257, t2051 * t2259 + t2069 * t2257, 0, 0, 0, 0, 0, 0, t2054 * t2259 + t2070 * t2257, t2064 * t2259 + t2073 * t2257, t2038 * t2259 + t2052 * t2257, t2021 * t2259 + t2030 * t2257, 0, 0, 0, 0, 0, 0, t2027 * t2259 + t2036 * t2257, t2029 * t2259 + t2039 * t2257, t2025 * t2259 + t2033 * t2257, t2016 * t2259 + t2017 * t2257; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2322, 0, 0, 0, 0, 0, 0, t2220, t2219, 0, t2152, 0, 0, 0, 0, 0, 0, t2168, t2169, t2204, t2105, 0, 0, 0, 0, 0, 0, t2097, t2107, t2086, t2050, 0, 0, 0, 0, 0, 0, t2053, t2063, t2037, t2020, 0, 0, 0, 0, 0, 0, t2026, t2028, t2024, t2015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2238, 0, 0, 0, 0, 0, 0, t2236, -t2280, 0, t2160, 0, 0, 0, 0, 0, 0, t2193, t2194, t2206, t2130, 0, 0, 0, 0, 0, 0, t2126, t2129, t2111, t2069, 0, 0, 0, 0, 0, 0, t2070, t2073, t2052, t2030, 0, 0, 0, 0, 0, 0, t2036, t2039, t2033, t2017; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2300, 0, 0, 0, 0, 0, 0, t2222, t2221, 0, t2153, 0, 0, 0, 0, 0, 0, t2170, t2171, t2205, t2106, 0, 0, 0, 0, 0, 0, t2098, t2108, t2087, t2051, 0, 0, 0, 0, 0, 0, t2054, t2064, t2038, t2021, 0, 0, 0, 0, 0, 0, t2027, t2029, t2025, t2016; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2322, 0, 0, 0, 0, 0, 0, t2220, t2219, 0, t2152, 0, 0, 0, 0, 0, 0, t2168, t2169, t2204, t2105, 0, 0, 0, 0, 0, 0, t2097, t2107, t2086, t2050, 0, 0, 0, 0, 0, 0, t2053, t2063, t2037, t2020, 0, 0, 0, 0, 0, 0, t2026, t2028, t2024, t2015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2273, -qJDD(2), 0, t2197, 0, 0, 0, 0, 0, 0, t2209, t2210, t2234, t2141, 0, 0, 0, 0, 0, 0, t2139, t2148, t2121, t2072, 0, 0, 0, 0, 0, 0, t2081, t2093, t2056, t2032, 0, 0, 0, 0, 0, 0, t2043, t2045, t2041, t2019; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t2273, 0, -t2299, 0, 0, 0, 0, 0, 0, t2233, -t2231, t2237, -t2190, 0, 0, 0, 0, 0, 0, -t2174, -t2176, -t2189, -t2167, 0, 0, 0, 0, 0, 0, -t2131, t2304, -t2154, -t2135, 0, 0, 0, 0, 0, 0, -t2099, -t2101, -t2090, -t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2215, 0, 0, 0, 0, 0, 0, t2207, t2208, 0, t2140, 0, 0, 0, 0, 0, 0, t2138, t2147, t2120, t2071, 0, 0, 0, 0, 0, 0, t2080, t2092, t2055, t2031, 0, 0, 0, 0, 0, 0, t2042, t2044, t2040, t2018; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2247, t2240, t2311, t2173, 0, 0, 0, 0, 0, 0, t2166, t2179, t2145, t2096, 0, 0, 0, 0, 0, 0, t2110, t2123, t2077, t2049, 0, 0, 0, 0, 0, 0, t2060, t2062, t2058, t2023; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2239, t2246, -t2312, t2172, 0, 0, 0, 0, 0, 0, t2165, t2178, t2144, t2095, 0, 0, 0, 0, 0, 0, t2109, t2122, t2076, t2048, 0, 0, 0, 0, 0, 0, t2059, t2061, t2057, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2233, t2231, -t2237, t2190, 0, 0, 0, 0, 0, 0, t2174, t2176, t2189, t2167, 0, 0, 0, 0, 0, 0, t2131, -t2304, t2154, t2135, 0, 0, 0, 0, 0, 0, t2099, t2101, t2090, t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2201, t2203, t2175, t2128, 0, 0, 0, 0, 0, 0, t2137, t2150, t2104, t2068, 0, 0, 0, 0, 0, 0, t2083, t2085, t2079, t2035; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2202, t2211, t2177, t2127, 0, 0, 0, 0, 0, 0, t2136, t2149, t2103, t2067, 0, 0, 0, 0, 0, 0, t2082, t2084, t2078, t2034; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2174, t2176, t2189, t2167, 0, 0, 0, 0, 0, 0, t2131, -t2304, t2154, t2135, 0, 0, 0, 0, 0, 0, t2099, t2101, t2090, t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2161, t2163, t2132, t2089, 0, 0, 0, 0, 0, 0, t2100, t2102, t2091, t2047; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2162, t2184, t2134, t2088, 0, 0, 0, 0, 0, 0, -t2116, -t2118, -t2142, -t2074; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2131, -t2304, t2154, t2135, 0, 0, 0, 0, 0, 0, t2099, t2101, t2090, t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2146, t2125, t2117, t2066; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2124, t2151, t2119, t2065; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2116, t2118, t2142, t2074;];
f_new_reg  = t1;
