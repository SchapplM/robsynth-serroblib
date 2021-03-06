% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRRRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:31:50
% EndTime: 2019-05-08 08:32:01
% DurationCPUTime: 11.61s
% Computational Cost: add. (117654->349), mult. (254126->489), div. (0->0), fcn. (196811->12), ass. (0->280)
t2331 = cos(qJ(2));
t2316 = t2331 * qJDD(1);
t2325 = sin(qJ(2));
t2379 = qJD(1) * t2325;
t2359 = qJD(2) * t2379;
t2296 = t2316 - t2359;
t2320 = t2331 ^ 2;
t2334 = qJD(1) ^ 2;
t2326 = sin(qJ(1));
t2332 = cos(qJ(1));
t2304 = t2326 * g(1) - t2332 * g(2);
t2346 = qJDD(1) * pkin(1) + t2304;
t2347 = qJD(2) * pkin(2) - pkin(8) * t2379;
t2267 = t2296 * pkin(2) - t2347 * t2379 + t2346 + (pkin(8) * t2320 + pkin(7)) * t2334;
t2324 = sin(qJ(3));
t2330 = cos(qJ(3));
t2289 = (t2324 * t2331 + t2325 * t2330) * qJD(1);
t2378 = qJD(1) * t2331;
t2358 = qJD(2) * t2378;
t2362 = t2325 * qJDD(1);
t2295 = t2358 + t2362;
t2352 = t2324 * t2295 - t2330 * t2296;
t2264 = -qJD(3) * t2289 - t2352;
t2288 = -t2324 * t2379 + t2330 * t2378;
t2287 = t2288 ^ 2;
t2318 = qJD(2) + qJD(3);
t2350 = pkin(3) * t2318 - pkin(9) * t2289;
t2221 = t2264 * pkin(3) + t2287 * pkin(9) - t2289 * t2350 + t2267;
t2323 = sin(qJ(4));
t2329 = cos(qJ(4));
t2272 = -t2329 * t2288 + t2289 * t2323;
t2270 = qJD(5) + t2272;
t2269 = qJD(6) + t2270;
t2392 = qJD(6) + t2269;
t2274 = t2288 * t2323 + t2289 * t2329;
t2315 = qJD(4) + t2318;
t2322 = sin(qJ(5));
t2328 = cos(qJ(5));
t2260 = t2274 * t2322 - t2328 * t2315;
t2262 = t2274 * t2328 + t2315 * t2322;
t2321 = sin(qJ(6));
t2327 = cos(qJ(6));
t2236 = t2327 * t2260 + t2262 * t2321;
t2391 = t2236 ^ 2;
t2238 = -t2260 * t2321 + t2262 * t2327;
t2390 = t2238 ^ 2;
t2389 = t2260 ^ 2;
t2388 = t2262 ^ 2;
t2387 = t2269 ^ 2;
t2386 = t2270 ^ 2;
t2385 = t2272 ^ 2;
t2384 = t2274 ^ 2;
t2383 = t2289 ^ 2;
t2382 = t2315 ^ 2;
t2381 = t2318 ^ 2;
t2380 = t2325 * g(3);
t2377 = t2236 * t2238;
t2376 = t2260 * t2262;
t2375 = t2270 * t2260;
t2374 = t2272 * t2274;
t2373 = t2288 * t2289;
t2372 = t2318 * t2288;
t2371 = t2320 * t2334;
t2305 = -g(1) * t2332 - g(2) * t2326;
t2291 = -pkin(1) * t2334 + qJDD(1) * pkin(7) + t2305;
t2370 = t2325 * t2291;
t2369 = t2325 * t2334;
t2368 = t2331 * t2291;
t2367 = -qJD(3) + t2318;
t2366 = qJD(4) - t2315;
t2365 = qJD(5) - t2270;
t2364 = qJD(6) - t2269;
t2336 = -pkin(2) * t2371 + t2296 * pkin(8) - qJD(2) * t2347 - t2380;
t2343 = pkin(8) * qJD(1) * qJD(2) + pkin(2) * t2369 - g(3);
t2344 = qJDD(2) * pkin(2) - t2295 * pkin(8) - t2370;
t2235 = t2330 * (t2336 + t2368) + t2324 * (t2331 * t2343 + t2344);
t2210 = -t2287 * pkin(3) + t2264 * pkin(9) - t2318 * t2350 + t2235;
t2234 = -t2324 * t2336 + t2330 * t2344 + (-t2324 * t2291 + t2330 * t2343) * t2331;
t2348 = -t2330 * t2295 - t2324 * t2296;
t2265 = qJD(3) * t2288 - t2348;
t2361 = qJDD(2) + qJDD(3);
t2276 = t2361 + t2373;
t2335 = (-t2265 + t2372) * pkin(9) + t2276 * pkin(3) + t2234;
t2180 = t2329 * t2210 + t2323 * t2335;
t2244 = pkin(4) * t2272 - pkin(10) * t2274;
t2357 = qJDD(4) + t2361;
t2171 = -pkin(4) * t2382 + pkin(10) * t2357 - t2272 * t2244 + t2180;
t2353 = -t2329 * t2264 + t2323 * t2265;
t2213 = (qJD(4) + t2315) * t2274 + t2353;
t2349 = -t2323 * t2264 - t2329 * t2265;
t2224 = -qJD(4) * t2272 - t2349;
t2355 = t2315 * t2272 - t2224;
t2176 = pkin(4) * t2213 + pkin(10) * t2355 - t2221;
t2145 = t2328 * t2171 + t2322 * t2176;
t2319 = t2325 ^ 2;
t2363 = t2319 + t2320;
t2144 = -t2322 * t2171 + t2328 * t2176;
t2179 = -t2210 * t2323 + t2329 * t2335;
t2341 = -t2328 * t2224 - t2322 * t2357;
t2211 = -t2260 * qJD(5) - t2341;
t2354 = t2322 * t2224 - t2328 * t2357;
t2345 = qJD(5) * t2262 + t2354;
t2356 = -t2321 * t2211 - t2327 * t2345;
t2342 = qJD(4) * t2274 + qJDD(5) + t2353;
t2170 = -t2357 * pkin(4) - t2382 * pkin(10) + t2244 * t2274 - t2179;
t2340 = -qJDD(6) - t2342;
t2201 = t2342 - t2376;
t2338 = -t2327 * t2211 + t2321 * t2345;
t2333 = qJD(2) ^ 2;
t2311 = t2331 * t2369;
t2309 = -t2333 - t2371;
t2308 = -t2319 * t2334 - t2333;
t2303 = -qJDD(2) + t2311;
t2302 = qJDD(2) + t2311;
t2301 = t2363 * t2334;
t2300 = -qJDD(1) * t2326 - t2332 * t2334;
t2299 = qJDD(1) * t2332 - t2326 * t2334;
t2298 = t2363 * qJDD(1);
t2297 = t2316 - 0.2e1 * t2359;
t2294 = 0.2e1 * t2358 + t2362;
t2290 = t2334 * pkin(7) + t2346;
t2284 = t2368 - t2380;
t2283 = -t2331 * g(3) - t2370;
t2282 = -t2381 - t2383;
t2281 = t2303 * t2331 - t2308 * t2325;
t2280 = -t2302 * t2325 + t2309 * t2331;
t2279 = t2303 * t2325 + t2308 * t2331;
t2278 = t2302 * t2331 + t2309 * t2325;
t2277 = -t2361 + t2373;
t2275 = -t2381 - t2287;
t2266 = -t2287 - t2383;
t2263 = -t2382 - t2384;
t2256 = -t2283 * t2325 + t2284 * t2331;
t2255 = t2283 * t2331 + t2284 * t2325;
t2252 = t2277 * t2330 - t2282 * t2324;
t2251 = t2277 * t2324 + t2282 * t2330;
t2250 = t2288 * t2367 + t2348;
t2249 = t2265 + t2372;
t2248 = t2289 * t2367 - t2352;
t2247 = (qJD(3) + t2318) * t2289 + t2352;
t2246 = t2275 * t2330 - t2276 * t2324;
t2245 = t2275 * t2324 + t2276 * t2330;
t2243 = -t2357 - t2374;
t2242 = t2357 - t2374;
t2241 = -t2382 - t2385;
t2239 = pkin(5) * t2270 - pkin(11) * t2262;
t2233 = -t2384 - t2385;
t2232 = -t2386 - t2388;
t2231 = t2243 * t2329 - t2263 * t2323;
t2230 = t2243 * t2323 + t2263 * t2329;
t2229 = -t2251 * t2325 + t2252 * t2331;
t2228 = t2251 * t2331 + t2252 * t2325;
t2227 = -t2386 - t2389;
t2226 = t2248 * t2330 - t2250 * t2324;
t2225 = t2248 * t2324 + t2250 * t2330;
t2222 = -t2388 - t2389;
t2220 = -t2245 * t2325 + t2246 * t2331;
t2219 = t2245 * t2331 + t2246 * t2325;
t2218 = t2241 * t2329 - t2242 * t2323;
t2217 = t2241 * t2323 + t2242 * t2329;
t2216 = t2272 * t2366 + t2349;
t2214 = -t2274 * t2366 - t2353;
t2212 = -t2387 - t2390;
t2204 = -t2234 * t2324 + t2235 * t2330;
t2203 = t2234 * t2330 + t2235 * t2324;
t2202 = -t2342 - t2376;
t2200 = -t2387 - t2391;
t2199 = -t2230 * t2324 + t2231 * t2330;
t2198 = t2230 * t2330 + t2231 * t2324;
t2197 = -t2225 * t2325 + t2226 * t2331;
t2196 = t2225 * t2331 + t2226 * t2325;
t2195 = t2260 * t2365 + t2341;
t2194 = t2211 - t2375;
t2193 = -t2262 * t2365 - t2354;
t2192 = (qJD(5) + t2270) * t2262 + t2354;
t2191 = -t2390 - t2391;
t2190 = -t2217 * t2324 + t2218 * t2330;
t2189 = t2217 * t2330 + t2218 * t2324;
t2188 = t2340 - t2377;
t2187 = -t2340 - t2377;
t2186 = t2214 * t2329 - t2216 * t2323;
t2185 = t2214 * t2323 + t2216 * t2329;
t2184 = t2202 * t2328 - t2232 * t2322;
t2183 = t2202 * t2322 + t2232 * t2328;
t2182 = -t2201 * t2322 + t2227 * t2328;
t2181 = t2201 * t2328 + t2227 * t2322;
t2178 = -t2203 * t2325 + t2204 * t2331;
t2177 = t2203 * t2331 + t2204 * t2325;
t2173 = t2188 * t2327 - t2212 * t2321;
t2172 = t2188 * t2321 + t2212 * t2327;
t2168 = -t2198 * t2325 + t2199 * t2331;
t2167 = t2198 * t2331 + t2199 * t2325;
t2166 = t2193 * t2328 - t2195 * t2322;
t2165 = t2193 * t2322 + t2195 * t2328;
t2164 = -t2187 * t2321 + t2200 * t2327;
t2163 = t2187 * t2327 + t2200 * t2321;
t2162 = t2236 * t2364 + t2338;
t2161 = -t2236 * t2392 - t2338;
t2160 = -t2238 * t2364 + t2356;
t2159 = t2238 * t2392 - t2356;
t2158 = t2184 * t2329 + t2194 * t2323;
t2157 = t2184 * t2323 - t2194 * t2329;
t2156 = t2182 * t2329 + t2192 * t2323;
t2155 = t2182 * t2323 - t2192 * t2329;
t2154 = -t2189 * t2325 + t2190 * t2331;
t2153 = t2189 * t2331 + t2190 * t2325;
t2152 = t2166 * t2329 + t2222 * t2323;
t2151 = t2166 * t2323 - t2222 * t2329;
t2150 = -t2185 * t2324 + t2186 * t2330;
t2149 = t2185 * t2330 + t2186 * t2324;
t2148 = -t2179 * t2323 + t2180 * t2329;
t2147 = t2179 * t2329 + t2180 * t2323;
t2146 = pkin(5) * t2345 - pkin(11) * t2389 + t2239 * t2262 + t2170;
t2143 = -t2172 * t2322 + t2173 * t2328;
t2142 = t2172 * t2328 + t2173 * t2322;
t2141 = -t2163 * t2322 + t2164 * t2328;
t2140 = t2163 * t2328 + t2164 * t2322;
t2139 = t2160 * t2327 - t2162 * t2321;
t2138 = t2160 * t2321 + t2162 * t2327;
t2137 = -t2157 * t2324 + t2158 * t2330;
t2136 = t2157 * t2330 + t2158 * t2324;
t2135 = -t2155 * t2324 + t2156 * t2330;
t2134 = t2155 * t2330 + t2156 * t2324;
t2133 = -pkin(5) * t2389 - pkin(11) * t2345 - t2270 * t2239 + t2145;
t2132 = -t2151 * t2324 + t2152 * t2330;
t2131 = t2151 * t2330 + t2152 * t2324;
t2130 = -t2149 * t2325 + t2150 * t2331;
t2129 = t2149 * t2331 + t2150 * t2325;
t2128 = (-t2211 - t2375) * pkin(11) + t2201 * pkin(5) + t2144;
t2127 = -t2147 * t2324 + t2148 * t2330;
t2126 = t2147 * t2330 + t2148 * t2324;
t2125 = t2143 * t2329 + t2161 * t2323;
t2124 = t2143 * t2323 - t2161 * t2329;
t2123 = t2141 * t2329 + t2159 * t2323;
t2122 = t2141 * t2323 - t2159 * t2329;
t2121 = -t2144 * t2322 + t2145 * t2328;
t2120 = t2144 * t2328 + t2145 * t2322;
t2119 = -t2138 * t2322 + t2139 * t2328;
t2118 = t2138 * t2328 + t2139 * t2322;
t2117 = -t2136 * t2325 + t2137 * t2331;
t2116 = t2136 * t2331 + t2137 * t2325;
t2115 = t2121 * t2329 + t2170 * t2323;
t2114 = t2121 * t2323 - t2170 * t2329;
t2113 = -t2134 * t2325 + t2135 * t2331;
t2112 = t2134 * t2331 + t2135 * t2325;
t2111 = -t2131 * t2325 + t2132 * t2331;
t2110 = t2131 * t2331 + t2132 * t2325;
t2109 = t2119 * t2329 + t2191 * t2323;
t2108 = t2119 * t2323 - t2191 * t2329;
t2107 = t2128 * t2321 + t2133 * t2327;
t2106 = t2128 * t2327 - t2133 * t2321;
t2105 = -t2126 * t2325 + t2127 * t2331;
t2104 = t2126 * t2331 + t2127 * t2325;
t2103 = -t2124 * t2324 + t2125 * t2330;
t2102 = t2124 * t2330 + t2125 * t2324;
t2101 = -t2122 * t2324 + t2123 * t2330;
t2100 = t2122 * t2330 + t2123 * t2324;
t2099 = -t2114 * t2324 + t2115 * t2330;
t2098 = t2114 * t2330 + t2115 * t2324;
t2097 = -t2108 * t2324 + t2109 * t2330;
t2096 = t2108 * t2330 + t2109 * t2324;
t2095 = -t2106 * t2321 + t2107 * t2327;
t2094 = t2106 * t2327 + t2107 * t2321;
t2093 = -t2102 * t2325 + t2103 * t2331;
t2092 = t2102 * t2331 + t2103 * t2325;
t2091 = -t2100 * t2325 + t2101 * t2331;
t2090 = t2100 * t2331 + t2101 * t2325;
t2089 = -t2098 * t2325 + t2099 * t2331;
t2088 = t2098 * t2331 + t2099 * t2325;
t2087 = -t2096 * t2325 + t2097 * t2331;
t2086 = t2096 * t2331 + t2097 * t2325;
t2085 = -t2094 * t2322 + t2095 * t2328;
t2084 = t2094 * t2328 + t2095 * t2322;
t2083 = t2085 * t2329 + t2146 * t2323;
t2082 = t2085 * t2323 - t2146 * t2329;
t2081 = -t2082 * t2324 + t2083 * t2330;
t2080 = t2082 * t2330 + t2083 * t2324;
t2079 = -t2080 * t2325 + t2081 * t2331;
t2078 = t2080 * t2331 + t2081 * t2325;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2300, -t2299, 0, -t2304 * t2326 + t2305 * t2332, 0, 0, 0, 0, 0, 0, t2280 * t2332 - t2297 * t2326, t2281 * t2332 + t2294 * t2326, t2298 * t2332 - t2301 * t2326, t2256 * t2332 - t2290 * t2326, 0, 0, 0, 0, 0, 0, t2220 * t2332 + t2247 * t2326, t2229 * t2332 + t2249 * t2326, t2197 * t2332 + t2266 * t2326, t2178 * t2332 - t2267 * t2326, 0, 0, 0, 0, 0, 0, t2154 * t2332 + t2213 * t2326, t2168 * t2332 - t2326 * t2355, t2130 * t2332 + t2233 * t2326, t2105 * t2332 - t2221 * t2326, 0, 0, 0, 0, 0, 0, t2113 * t2332 + t2181 * t2326, t2117 * t2332 + t2183 * t2326, t2111 * t2332 + t2165 * t2326, t2089 * t2332 + t2120 * t2326, 0, 0, 0, 0, 0, 0, t2091 * t2332 + t2140 * t2326, t2093 * t2332 + t2142 * t2326, t2087 * t2332 + t2118 * t2326, t2079 * t2332 + t2084 * t2326; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2299, t2300, 0, t2304 * t2332 + t2305 * t2326, 0, 0, 0, 0, 0, 0, t2280 * t2326 + t2297 * t2332, t2281 * t2326 - t2294 * t2332, t2298 * t2326 + t2301 * t2332, t2256 * t2326 + t2290 * t2332, 0, 0, 0, 0, 0, 0, t2220 * t2326 - t2247 * t2332, t2229 * t2326 - t2249 * t2332, t2197 * t2326 - t2266 * t2332, t2178 * t2326 + t2267 * t2332, 0, 0, 0, 0, 0, 0, t2154 * t2326 - t2213 * t2332, t2168 * t2326 + t2332 * t2355, t2130 * t2326 - t2233 * t2332, t2105 * t2326 + t2221 * t2332, 0, 0, 0, 0, 0, 0, t2113 * t2326 - t2181 * t2332, t2117 * t2326 - t2183 * t2332, t2111 * t2326 - t2165 * t2332, t2089 * t2326 - t2120 * t2332, 0, 0, 0, 0, 0, 0, t2091 * t2326 - t2140 * t2332, t2093 * t2326 - t2142 * t2332, t2087 * t2326 - t2118 * t2332, t2079 * t2326 - t2084 * t2332; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2278, t2279, 0, t2255, 0, 0, 0, 0, 0, 0, t2219, t2228, t2196, t2177, 0, 0, 0, 0, 0, 0, t2153, t2167, t2129, t2104, 0, 0, 0, 0, 0, 0, t2112, t2116, t2110, t2088, 0, 0, 0, 0, 0, 0, t2090, t2092, t2086, t2078; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2334, -qJDD(1), 0, t2305, 0, 0, 0, 0, 0, 0, t2280, t2281, t2298, t2256, 0, 0, 0, 0, 0, 0, t2220, t2229, t2197, t2178, 0, 0, 0, 0, 0, 0, t2154, t2168, t2130, t2105, 0, 0, 0, 0, 0, 0, t2113, t2117, t2111, t2089, 0, 0, 0, 0, 0, 0, t2091, t2093, t2087, t2079; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2334, 0, t2304, 0, 0, 0, 0, 0, 0, t2297, -t2294, t2301, t2290, 0, 0, 0, 0, 0, 0, -t2247, -t2249, -t2266, t2267, 0, 0, 0, 0, 0, 0, -t2213, t2355, -t2233, t2221, 0, 0, 0, 0, 0, 0, -t2181, -t2183, -t2165, -t2120, 0, 0, 0, 0, 0, 0, -t2140, -t2142, -t2118, -t2084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2278, t2279, 0, t2255, 0, 0, 0, 0, 0, 0, t2219, t2228, t2196, t2177, 0, 0, 0, 0, 0, 0, t2153, t2167, t2129, t2104, 0, 0, 0, 0, 0, 0, t2112, t2116, t2110, t2088, 0, 0, 0, 0, 0, 0, t2090, t2092, t2086, t2078; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2309, t2303, t2316, t2284, 0, 0, 0, 0, 0, 0, t2246, t2252, t2226, t2204, 0, 0, 0, 0, 0, 0, t2190, t2199, t2150, t2127, 0, 0, 0, 0, 0, 0, t2135, t2137, t2132, t2099, 0, 0, 0, 0, 0, 0, t2101, t2103, t2097, t2081; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2302, t2308, -t2362, t2283, 0, 0, 0, 0, 0, 0, t2245, t2251, t2225, t2203, 0, 0, 0, 0, 0, 0, t2189, t2198, t2149, t2126, 0, 0, 0, 0, 0, 0, t2134, t2136, t2131, t2098, 0, 0, 0, 0, 0, 0, t2100, t2102, t2096, t2080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2297, t2294, -t2301, -t2290, 0, 0, 0, 0, 0, 0, t2247, t2249, t2266, -t2267, 0, 0, 0, 0, 0, 0, t2213, -t2355, t2233, -t2221, 0, 0, 0, 0, 0, 0, t2181, t2183, t2165, t2120, 0, 0, 0, 0, 0, 0, t2140, t2142, t2118, t2084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2275, t2277, t2248, t2235, 0, 0, 0, 0, 0, 0, t2218, t2231, t2186, t2148, 0, 0, 0, 0, 0, 0, t2156, t2158, t2152, t2115, 0, 0, 0, 0, 0, 0, t2123, t2125, t2109, t2083; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2276, t2282, t2250, t2234, 0, 0, 0, 0, 0, 0, t2217, t2230, t2185, t2147, 0, 0, 0, 0, 0, 0, t2155, t2157, t2151, t2114, 0, 0, 0, 0, 0, 0, t2122, t2124, t2108, t2082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2247, t2249, t2266, -t2267, 0, 0, 0, 0, 0, 0, t2213, -t2355, t2233, -t2221, 0, 0, 0, 0, 0, 0, t2181, t2183, t2165, t2120, 0, 0, 0, 0, 0, 0, t2140, t2142, t2118, t2084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2241, t2243, t2214, t2180, 0, 0, 0, 0, 0, 0, t2182, t2184, t2166, t2121, 0, 0, 0, 0, 0, 0, t2141, t2143, t2119, t2085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2242, t2263, t2216, t2179, 0, 0, 0, 0, 0, 0, -t2192, -t2194, -t2222, -t2170, 0, 0, 0, 0, 0, 0, -t2159, -t2161, -t2191, -t2146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2213, -t2355, t2233, -t2221, 0, 0, 0, 0, 0, 0, t2181, t2183, t2165, t2120, 0, 0, 0, 0, 0, 0, t2140, t2142, t2118, t2084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2227, t2202, t2193, t2145, 0, 0, 0, 0, 0, 0, t2164, t2173, t2139, t2095; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2201, t2232, t2195, t2144, 0, 0, 0, 0, 0, 0, t2163, t2172, t2138, t2094; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2192, t2194, t2222, t2170, 0, 0, 0, 0, 0, 0, t2159, t2161, t2191, t2146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2200, t2188, t2160, t2107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2187, t2212, t2162, t2106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2159, t2161, t2191, t2146;];
f_new_reg  = t1;
