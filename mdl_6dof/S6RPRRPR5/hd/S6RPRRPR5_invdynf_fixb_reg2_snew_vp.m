% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRRPR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:41:54
% EndTime: 2019-05-05 22:42:02
% DurationCPUTime: 9.19s
% Computational Cost: add. (34446->322), mult. (84552->391), div. (0->0), fcn. (66458->10), ass. (0->223)
t2212 = qJDD(3) + qJDD(4);
t2220 = sin(qJ(3));
t2224 = cos(qJ(3));
t2217 = cos(pkin(10));
t2266 = qJD(1) * t2217;
t2216 = sin(pkin(10));
t2267 = qJD(1) * t2216;
t2189 = -t2220 * t2267 + t2224 * t2266;
t2235 = t2216 * t2224 + t2217 * t2220;
t2190 = t2235 * qJD(1);
t2219 = sin(qJ(4));
t2223 = cos(qJ(4));
t2165 = -t2223 * t2189 + t2190 * t2219;
t2167 = t2189 * t2219 + t2190 * t2223;
t2260 = t2167 * t2165;
t2131 = t2212 + t2260;
t2215 = qJD(3) + qJD(4);
t2211 = t2215 ^ 2;
t2272 = t2167 ^ 2;
t2246 = -t2211 - t2272;
t2115 = t2131 * t2219 - t2223 * t2246;
t2117 = t2131 * t2223 + t2219 * t2246;
t2082 = t2115 * t2224 + t2117 * t2220;
t2085 = t2115 * t2220 - t2117 * t2224;
t2060 = t2082 * t2216 + t2085 * t2217;
t2221 = sin(qJ(1));
t2292 = t2060 * t2221;
t2225 = cos(qJ(1));
t2291 = t2060 * t2225;
t2164 = t2165 ^ 2;
t2130 = -t2211 - t2164;
t2239 = t2212 - t2260;
t2096 = t2130 * t2219 + t2223 * t2239;
t2099 = -t2130 * t2223 + t2219 * t2239;
t2074 = t2096 * t2224 - t2099 * t2220;
t2077 = t2096 * t2220 + t2099 * t2224;
t2045 = t2074 * t2216 + t2077 * t2217;
t2290 = t2045 * t2221;
t2289 = t2045 * t2225;
t2057 = t2082 * t2217 - t2085 * t2216;
t2042 = t2074 * t2217 - t2077 * t2216;
t2208 = t2217 * qJDD(1);
t2249 = t2216 * qJDD(1);
t2145 = t2224 * t2208 - t2220 * t2249;
t2257 = t2190 * qJD(3);
t2176 = t2145 - t2257;
t2259 = t2189 * qJD(3);
t2280 = t2235 * qJDD(1);
t2178 = t2280 + t2259;
t2236 = t2219 * t2176 + t2223 * t2178;
t2254 = qJD(4) - t2215;
t2110 = t2165 * t2254 - t2236;
t2288 = t2110 * t2219;
t2287 = t2110 * t2223;
t2227 = qJD(1) ^ 2;
t2198 = t2221 * g(1) - t2225 * g(2);
t2240 = -qJDD(2) + t2198;
t2213 = t2216 ^ 2;
t2214 = t2217 ^ 2;
t2250 = t2213 + t2214;
t2172 = (pkin(7) * t2250 + qJ(2)) * t2227 + (pkin(2) * t2217 + pkin(1)) * qJDD(1) + t2240;
t2187 = t2189 ^ 2;
t2241 = qJD(3) * pkin(3) - pkin(8) * t2190;
t2122 = t2176 * pkin(3) + t2187 * pkin(8) - t2190 * t2241 + t2172;
t2159 = t2215 * t2167;
t2127 = -qJD(4) * t2165 + t2236;
t2261 = t2165 * t2215;
t2243 = -t2127 + t2261;
t2286 = pkin(4) * t2159 + t2243 * qJ(5) - 0.2e1 * qJD(5) * t2167 - t2122;
t2281 = -t2164 - t2272;
t2285 = t2221 * t2281;
t2284 = t2225 * t2281;
t2195 = t2250 * t2227;
t2218 = sin(qJ(6));
t2222 = cos(qJ(6));
t2146 = -t2222 * t2165 + t2215 * t2218;
t2275 = t2146 ^ 2;
t2148 = t2165 * t2218 + t2215 * t2222;
t2274 = t2148 ^ 2;
t2163 = qJD(6) + t2167;
t2273 = t2163 ^ 2;
t2271 = t2190 ^ 2;
t2270 = 0.2e1 * qJD(5);
t2262 = t2146 * t2148;
t2258 = t2189 * t2190;
t2256 = t2214 * t2227;
t2255 = t2217 * t2227;
t2253 = qJD(4) + t2215;
t2252 = qJD(6) - t2163;
t2251 = qJD(6) + t2163;
t2199 = -g(1) * t2225 - g(2) * t2221;
t2191 = -pkin(1) * t2227 + qJDD(1) * qJ(2) + t2199;
t2245 = -t2217 * g(3) - 0.2e1 * qJD(2) * t2267;
t2156 = (pkin(2) * t2255 - pkin(7) * qJDD(1) - t2191) * t2216 + t2245;
t2180 = -g(3) * t2216 + 0.2e1 * qJD(2) * t2266 + t2217 * t2191;
t2160 = -pkin(2) * t2256 + pkin(7) * t2208 + t2180;
t2129 = t2220 * t2156 + t2224 * t2160;
t2247 = -t2273 - t2274;
t2128 = t2224 * t2156 - t2220 * t2160;
t2173 = qJDD(3) + t2258;
t2095 = (-t2178 + t2259) * pkin(8) + t2173 * pkin(3) + t2128;
t2100 = -t2187 * pkin(3) + t2176 * pkin(8) - qJD(3) * t2241 + t2129;
t2070 = t2223 * t2095 - t2219 * t2100;
t2242 = -t2223 * t2176 + t2219 * t2178;
t2126 = qJD(4) * t2167 + t2242;
t2244 = t2222 * t2126 - t2218 * t2212;
t2135 = pkin(4) * t2165 - qJ(5) * t2167;
t2064 = -t2212 * pkin(4) - t2211 * qJ(5) + t2167 * t2135 + qJDD(5) - t2070;
t2038 = -t2239 * pkin(9) + (t2127 + t2261) * pkin(5) + t2064;
t2154 = pkin(5) * t2167 - pkin(9) * t2215;
t2041 = -t2164 * pkin(5) - t2167 * t2154 + (pkin(4) + pkin(9)) * t2126 + t2286;
t2238 = t2038 * t2222 - t2041 * t2218;
t2071 = t2219 * t2095 + t2223 * t2100;
t2237 = -t2218 * t2126 - t2222 * t2212;
t2233 = qJDD(6) + t2127;
t2232 = t2146 * t2252 + t2237;
t2231 = -t2211 * pkin(4) + t2212 * qJ(5) - t2165 * t2135 + t2071;
t2230 = t2233 - t2262;
t2226 = qJD(3) ^ 2;
t2200 = t2216 * t2255;
t2197 = -qJDD(1) * t2221 - t2225 * t2227;
t2196 = qJDD(1) * t2225 - t2221 * t2227;
t2194 = t2250 * qJDD(1);
t2193 = t2217 * t2195;
t2192 = t2216 * t2195;
t2186 = qJDD(1) * pkin(1) + t2227 * qJ(2) + t2240;
t2181 = -t2226 - t2271;
t2179 = -t2216 * t2191 + t2245;
t2177 = t2280 + 0.2e1 * t2259;
t2175 = -t2145 + 0.2e1 * t2257;
t2174 = -qJDD(3) + t2258;
t2171 = -t2226 - t2187;
t2149 = -t2187 - t2271;
t2144 = t2174 * t2224 - t2181 * t2220;
t2143 = t2174 * t2220 + t2181 * t2224;
t2142 = -t2179 * t2216 + t2180 * t2217;
t2141 = t2179 * t2217 + t2180 * t2216;
t2140 = t2145 * t2224 + t2220 * t2280;
t2139 = t2145 * t2220 - t2224 * t2280;
t2138 = t2171 * t2224 - t2173 * t2220;
t2137 = t2171 * t2220 + t2173 * t2224;
t2121 = -t2143 * t2216 + t2144 * t2217;
t2120 = t2143 * t2217 + t2144 * t2216;
t2114 = -t2273 - t2275;
t2113 = -t2274 - t2275;
t2112 = -t2139 * t2216 + t2140 * t2217;
t2111 = t2139 * t2217 + t2140 * t2216;
t2108 = -t2165 * t2253 + t2236;
t2106 = -t2126 + t2159;
t2105 = -t2167 * t2254 - t2242;
t2104 = t2126 + t2159;
t2103 = t2167 * t2253 + t2242;
t2102 = -t2137 * t2216 + t2138 * t2217;
t2101 = t2137 * t2217 + t2138 * t2216;
t2093 = -t2128 * t2220 + t2129 * t2224;
t2092 = t2128 * t2224 + t2129 * t2220;
t2091 = -t2233 - t2262;
t2089 = -t2146 * t2251 - t2237;
t2088 = -t2148 * t2252 + t2244;
t2087 = t2148 * t2251 - t2244;
t2081 = t2106 * t2223 - t2288;
t2080 = t2105 * t2223 - t2288;
t2079 = t2106 * t2219 + t2287;
t2078 = t2105 * t2219 + t2287;
t2073 = t2222 * t2091 - t2218 * t2247;
t2072 = t2218 * t2091 + t2222 * t2247;
t2069 = t2222 * t2114 - t2218 * t2230;
t2068 = t2218 * t2114 + t2222 * t2230;
t2067 = -t2092 * t2216 + t2093 * t2217;
t2066 = t2092 * t2217 + t2093 * t2216;
t2065 = -t2126 * pkin(4) - t2286;
t2063 = t2215 * t2270 + t2231;
t2062 = t2222 * t2088 - t2218 * t2232;
t2061 = t2218 * t2088 + t2222 * t2232;
t2056 = t2072 * t2219 + t2089 * t2223;
t2055 = -t2072 * t2223 + t2089 * t2219;
t2054 = t2068 * t2219 + t2087 * t2223;
t2053 = -t2068 * t2223 + t2087 * t2219;
t2052 = -t2079 * t2220 + t2081 * t2224;
t2051 = -t2078 * t2220 + t2080 * t2224;
t2050 = t2079 * t2224 + t2081 * t2220;
t2049 = t2078 * t2224 + t2080 * t2220;
t2048 = t2061 * t2219 + t2113 * t2223;
t2047 = -t2061 * t2223 + t2113 * t2219;
t2046 = -t2126 * pkin(5) - t2164 * pkin(9) + (t2270 + t2154) * t2215 + t2231;
t2040 = -t2070 * t2219 + t2071 * t2223;
t2039 = t2070 * t2223 + t2071 * t2219;
t2037 = t2063 * t2223 + t2064 * t2219;
t2036 = t2063 * t2219 - t2064 * t2223;
t2035 = -t2055 * t2220 + t2056 * t2224;
t2034 = t2055 * t2224 + t2056 * t2220;
t2033 = -t2053 * t2220 + t2054 * t2224;
t2032 = t2053 * t2224 + t2054 * t2220;
t2031 = -t2050 * t2216 + t2052 * t2217;
t2030 = -t2049 * t2216 + t2051 * t2217;
t2029 = t2050 * t2217 + t2052 * t2216;
t2028 = t2049 * t2217 + t2051 * t2216;
t2027 = -t2047 * t2220 + t2048 * t2224;
t2026 = t2047 * t2224 + t2048 * t2220;
t2025 = -t2039 * t2220 + t2040 * t2224;
t2024 = t2039 * t2224 + t2040 * t2220;
t2023 = t2038 * t2218 + t2041 * t2222;
t2021 = -t2036 * t2220 + t2037 * t2224;
t2020 = t2036 * t2224 + t2037 * t2220;
t2019 = -t2034 * t2216 + t2035 * t2217;
t2018 = t2034 * t2217 + t2035 * t2216;
t2017 = -t2032 * t2216 + t2033 * t2217;
t2016 = t2032 * t2217 + t2033 * t2216;
t2015 = -t2026 * t2216 + t2027 * t2217;
t2014 = t2026 * t2217 + t2027 * t2216;
t2013 = -t2024 * t2216 + t2025 * t2217;
t2012 = t2024 * t2217 + t2025 * t2216;
t2011 = t2222 * t2023 - t2218 * t2238;
t2010 = t2218 * t2023 + t2222 * t2238;
t2009 = t2010 * t2219 + t2046 * t2223;
t2008 = -t2010 * t2223 + t2046 * t2219;
t2007 = -t2020 * t2216 + t2021 * t2217;
t2006 = t2020 * t2217 + t2021 * t2216;
t2005 = -t2008 * t2220 + t2009 * t2224;
t2004 = t2008 * t2224 + t2009 * t2220;
t2003 = -t2004 * t2216 + t2005 * t2217;
t2002 = t2004 * t2217 + t2005 * t2216;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2197, -t2196, 0, -t2198 * t2221 + t2199 * t2225, 0, 0, 0, 0, 0, 0, -t2193 * t2225 - t2208 * t2221, t2192 * t2225 + t2221 * t2249, t2194 * t2225 - t2195 * t2221, t2142 * t2225 - t2186 * t2221, 0, 0, 0, 0, 0, 0, t2102 * t2225 + t2175 * t2221, t2121 * t2225 + t2177 * t2221, t2112 * t2225 + t2149 * t2221, t2067 * t2225 - t2172 * t2221, 0, 0, 0, 0, 0, 0, t2103 * t2221 - t2289, -t2221 * t2243 + t2291, t2031 * t2225 + t2285, t2013 * t2225 - t2122 * t2221, 0, 0, 0, 0, 0, 0, t2030 * t2225 + t2285, -t2104 * t2221 + t2289, -t2108 * t2221 - t2291, t2007 * t2225 - t2065 * t2221, 0, 0, 0, 0, 0, 0, t2017 * t2225 + t2069 * t2221, t2019 * t2225 + t2073 * t2221, t2015 * t2225 + t2062 * t2221, t2003 * t2225 + t2011 * t2221; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2196, t2197, 0, t2198 * t2225 + t2199 * t2221, 0, 0, 0, 0, 0, 0, -t2193 * t2221 + t2208 * t2225, t2192 * t2221 - t2225 * t2249, t2194 * t2221 + t2195 * t2225, t2142 * t2221 + t2186 * t2225, 0, 0, 0, 0, 0, 0, t2102 * t2221 - t2175 * t2225, t2121 * t2221 - t2177 * t2225, t2112 * t2221 - t2149 * t2225, t2067 * t2221 + t2172 * t2225, 0, 0, 0, 0, 0, 0, -t2103 * t2225 - t2290, t2225 * t2243 + t2292, t2031 * t2221 - t2284, t2013 * t2221 + t2122 * t2225, 0, 0, 0, 0, 0, 0, t2030 * t2221 - t2284, t2104 * t2225 + t2290, t2108 * t2225 - t2292, t2007 * t2221 + t2065 * t2225, 0, 0, 0, 0, 0, 0, t2017 * t2221 - t2069 * t2225, t2019 * t2221 - t2073 * t2225, t2015 * t2221 - t2062 * t2225, t2003 * t2221 - t2011 * t2225; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2141, 0, 0, 0, 0, 0, 0, t2101, t2120, t2111, t2066, 0, 0, 0, 0, 0, 0, t2042, -t2057, t2029, t2012, 0, 0, 0, 0, 0, 0, t2028, -t2042, t2057, t2006, 0, 0, 0, 0, 0, 0, t2016, t2018, t2014, t2002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2227, -qJDD(1), 0, t2199, 0, 0, 0, 0, 0, 0, -t2193, t2192, t2194, t2142, 0, 0, 0, 0, 0, 0, t2102, t2121, t2112, t2067, 0, 0, 0, 0, 0, 0, -t2045, t2060, t2031, t2013, 0, 0, 0, 0, 0, 0, t2030, t2045, -t2060, t2007, 0, 0, 0, 0, 0, 0, t2017, t2019, t2015, t2003; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2227, 0, t2198, 0, 0, 0, 0, 0, 0, t2208, -t2249, t2195, t2186, 0, 0, 0, 0, 0, 0, -t2175, -t2177, -t2149, t2172, 0, 0, 0, 0, 0, 0, -t2103, t2243, -t2281, t2122, 0, 0, 0, 0, 0, 0, -t2281, t2104, t2108, t2065, 0, 0, 0, 0, 0, 0, -t2069, -t2073, -t2062, -t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2141, 0, 0, 0, 0, 0, 0, t2101, t2120, t2111, t2066, 0, 0, 0, 0, 0, 0, t2042, -t2057, t2029, t2012, 0, 0, 0, 0, 0, 0, t2028, -t2042, t2057, t2006, 0, 0, 0, 0, 0, 0, t2016, t2018, t2014, t2002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2256, t2200, t2208, t2180, 0, 0, 0, 0, 0, 0, t2138, t2144, t2140, t2093, 0, 0, 0, 0, 0, 0, -t2077, t2085, t2052, t2025, 0, 0, 0, 0, 0, 0, t2051, t2077, -t2085, t2021, 0, 0, 0, 0, 0, 0, t2033, t2035, t2027, t2005; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2200, -t2213 * t2227, -t2249, t2179, 0, 0, 0, 0, 0, 0, t2137, t2143, t2139, t2092, 0, 0, 0, 0, 0, 0, t2074, -t2082, t2050, t2024, 0, 0, 0, 0, 0, 0, t2049, -t2074, t2082, t2020, 0, 0, 0, 0, 0, 0, t2032, t2034, t2026, t2004; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2208, t2249, -t2195, -t2186, 0, 0, 0, 0, 0, 0, t2175, t2177, t2149, -t2172, 0, 0, 0, 0, 0, 0, t2103, -t2243, t2281, -t2122, 0, 0, 0, 0, 0, 0, t2281, -t2104, -t2108, -t2065, 0, 0, 0, 0, 0, 0, t2069, t2073, t2062, t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2171, t2174, t2145, t2129, 0, 0, 0, 0, 0, 0, -t2099, -t2117, t2081, t2040, 0, 0, 0, 0, 0, 0, t2080, t2099, t2117, t2037, 0, 0, 0, 0, 0, 0, t2054, t2056, t2048, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2173, t2181, -t2280, t2128, 0, 0, 0, 0, 0, 0, t2096, -t2115, t2079, t2039, 0, 0, 0, 0, 0, 0, t2078, -t2096, t2115, t2036, 0, 0, 0, 0, 0, 0, t2053, t2055, t2047, t2008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2175, t2177, t2149, -t2172, 0, 0, 0, 0, 0, 0, t2103, -t2243, t2281, -t2122, 0, 0, 0, 0, 0, 0, t2281, -t2104, -t2108, -t2065, 0, 0, 0, 0, 0, 0, t2069, t2073, t2062, t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2130, -t2131, t2106, t2071, 0, 0, 0, 0, 0, 0, t2105, -t2130, t2131, t2063, 0, 0, 0, 0, 0, 0, t2087, t2089, t2113, t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2239, t2246, t2110, t2070, 0, 0, 0, 0, 0, 0, t2110, -t2239, -t2246, -t2064, 0, 0, 0, 0, 0, 0, -t2068, -t2072, -t2061, -t2010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2103, -t2243, t2281, -t2122, 0, 0, 0, 0, 0, 0, t2281, -t2104, -t2108, -t2065, 0, 0, 0, 0, 0, 0, t2069, t2073, t2062, t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2281, -t2104, -t2108, -t2065, 0, 0, 0, 0, 0, 0, t2069, t2073, t2062, t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2105, t2130, -t2131, -t2063, 0, 0, 0, 0, 0, 0, -t2087, -t2089, -t2113, -t2046; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2110, t2239, t2246, t2064, 0, 0, 0, 0, 0, 0, t2068, t2072, t2061, t2010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2114, t2091, t2088, t2023; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2230, t2247, t2232, t2238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2087, t2089, t2113, t2046;];
f_new_reg  = t1;
