% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRRPP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:24:03
% EndTime: 2019-05-05 21:24:09
% DurationCPUTime: 6.51s
% Computational Cost: add. (9900->276), mult. (19260->272), div. (0->0), fcn. (12328->8), ass. (0->158)
t2207 = cos(qJ(3));
t2234 = t2207 * qJD(1);
t2186 = -qJD(4) + t2234;
t2161 = t2186 ^ 2;
t2203 = sin(qJ(4));
t2206 = cos(qJ(4));
t2204 = sin(qJ(3));
t2244 = qJD(1) * t2204;
t2166 = t2203 * qJD(3) + t2206 * t2244;
t2163 = t2166 ^ 2;
t2138 = t2163 + t2161;
t2164 = -t2206 * qJD(3) + t2203 * t2244;
t2140 = t2166 * t2164;
t2191 = qJD(3) * t2244;
t2230 = t2207 * qJDD(1);
t2225 = t2191 - t2230;
t2224 = -qJDD(4) - t2225;
t2217 = -t2224 + t2140;
t2088 = t2203 * t2138 - t2206 * t2217;
t2228 = qJD(3) * t2234;
t2231 = t2204 * qJDD(1);
t2221 = t2228 + t2231;
t2213 = -t2164 * qJD(4) + t2203 * qJDD(3) + t2206 * t2221;
t2243 = t2186 * t2164;
t2211 = t2213 + t2243;
t2066 = t2207 * t2088 + t2204 * t2211;
t2199 = sin(pkin(9));
t2200 = cos(pkin(9));
t2275 = t2206 * t2138 + t2203 * t2217;
t2040 = t2199 * t2066 + t2200 * t2275;
t2043 = t2200 * t2066 - t2199 * t2275;
t2205 = sin(qJ(1));
t2208 = cos(qJ(1));
t2292 = t2208 * t2040 + t2205 * t2043;
t2291 = t2205 * t2040 - t2208 * t2043;
t2216 = -t2206 * qJDD(3) + t2203 * t2221;
t2099 = (qJD(4) + t2186) * t2166 + t2216;
t2212 = t2213 - t2243;
t2072 = t2206 * t2099 - t2203 * t2212;
t2248 = t2164 ^ 2;
t2117 = t2163 + t2248;
t2050 = t2207 * t2072 + t2204 * t2117;
t2069 = t2203 * t2099 + t2206 * t2212;
t2026 = t2199 * t2050 - t2200 * t2069;
t2029 = t2200 * t2050 + t2199 * t2069;
t2290 = t2208 * t2026 + t2205 * t2029;
t2281 = t2205 * t2026 - t2208 * t2029;
t2063 = t2204 * t2088 - t2207 * t2211;
t2215 = t2166 * qJD(4) + t2216;
t2098 = -t2166 * t2186 + t2215;
t2253 = -t2248 - t2161;
t2255 = t2140 + t2224;
t2267 = t2203 * t2255 + t2206 * t2253;
t2276 = t2207 * t2267;
t2271 = t2204 * t2098 + t2276;
t2268 = t2203 * t2253 - t2206 * t2255;
t2279 = t2199 * t2268;
t2282 = t2200 * t2271 + t2279;
t2278 = t2200 * t2268;
t2283 = t2199 * t2271 - t2278;
t2289 = -t2205 * t2283 + t2208 * t2282;
t2288 = t2205 * t2282 + t2208 * t2283;
t2047 = t2204 * t2072 - t2207 * t2117;
t2277 = t2204 * t2267;
t2274 = -t2207 * t2098 + t2277;
t2256 = qJD(4) - t2186;
t2249 = qJD(3) ^ 2;
t2254 = -qJDD(3) * pkin(3) - t2249 * pkin(8);
t2252 = t2215 * pkin(4) - t2211 * qJ(5);
t2247 = -2 * qJD(5);
t2246 = -g(3) + qJDD(2);
t2245 = t2166 * qJ(6);
t2179 = t2205 * g(1) - t2208 * g(2);
t2219 = qJDD(1) * pkin(1) + t2179;
t2180 = -t2208 * g(1) - t2205 * g(2);
t2209 = qJD(1) ^ 2;
t2220 = -t2209 * pkin(1) + t2180;
t2133 = t2199 * t2219 + t2200 * t2220;
t2124 = -t2209 * pkin(2) + qJDD(1) * pkin(7) + t2133;
t2108 = t2207 * t2124 + t2204 * t2246;
t2167 = (-t2207 * pkin(3) - t2204 * pkin(8)) * qJD(1);
t2095 = -t2249 * pkin(3) + qJDD(3) * pkin(8) + t2167 * t2234 + t2108;
t2132 = -t2199 * t2220 + t2200 * t2219;
t2123 = -qJDD(1) * pkin(2) - t2209 * pkin(7) - t2132;
t2168 = 0.2e1 * t2228 + t2231;
t2210 = (t2225 + t2191) * pkin(3) - t2168 * pkin(8) + t2123;
t2054 = t2206 * t2095 + t2203 * t2210;
t2195 = t2204 ^ 2;
t2196 = t2207 ^ 2;
t2232 = t2195 + t2196;
t2227 = -t2186 * pkin(4) + t2247;
t2053 = -t2203 * t2095 + t2206 * t2210;
t2190 = t2207 * t2246;
t2107 = -t2204 * t2124 + t2190;
t2170 = -t2199 * qJDD(1) - t2200 * t2209;
t2171 = t2200 * qJDD(1) - t2199 * t2209;
t2226 = t2208 * t2170 - t2205 * t2171;
t2134 = t2164 * pkin(4) - t2166 * qJ(5);
t2223 = -t2224 * qJ(5) - t2164 * t2134 + t2186 * t2247 + t2054;
t2222 = t2205 * t2170 + t2208 * t2171;
t2218 = t2224 * pkin(4) - t2161 * qJ(5) + qJDD(5) - t2053;
t2094 = -t2190 + (qJD(1) * t2167 + t2124) * t2204 + t2254;
t2185 = t2207 * t2209 * t2204;
t2183 = -t2196 * t2209 - t2249;
t2182 = -t2195 * t2209 - t2249;
t2177 = -qJDD(3) + t2185;
t2176 = qJDD(3) + t2185;
t2175 = t2232 * t2209;
t2174 = -t2205 * qJDD(1) - t2208 * t2209;
t2173 = t2208 * qJDD(1) - t2205 * t2209;
t2172 = t2232 * qJDD(1);
t2169 = -0.2e1 * t2191 + t2230;
t2144 = t2207 * t2177 - t2204 * t2182;
t2143 = -t2204 * t2176 + t2207 * t2183;
t2142 = t2204 * t2177 + t2207 * t2182;
t2141 = t2207 * t2176 + t2204 * t2183;
t2136 = t2200 * t2172 - t2199 * t2175;
t2135 = t2199 * t2172 + t2200 * t2175;
t2112 = t2200 * t2144 + t2199 * t2168;
t2111 = t2200 * t2143 - t2199 * t2169;
t2110 = t2199 * t2144 - t2200 * t2168;
t2109 = t2199 * t2143 + t2200 * t2169;
t2100 = t2256 * t2166 + t2216;
t2093 = -t2199 * t2132 + t2200 * t2133;
t2092 = t2200 * t2132 + t2199 * t2133;
t2076 = -t2204 * t2107 + t2207 * t2108;
t2075 = t2207 * t2107 + t2204 * t2108;
t2062 = t2204 * t2100 + t2276;
t2059 = -t2207 * t2100 + t2277;
t2056 = t2200 * t2076 + t2199 * t2123;
t2055 = t2199 * t2076 - t2200 * t2123;
t2046 = t2166 * t2227 + t2094 + t2252;
t2039 = t2166 * t2134 + t2218;
t2038 = t2200 * t2062 + t2279;
t2035 = t2199 * t2062 - t2278;
t2032 = -t2161 * pkin(4) + t2223;
t2025 = t2216 * pkin(5) + t2248 * qJ(6) + t2167 * t2244 - qJDD(6) - t2107 + (t2256 * pkin(5) + t2227 + t2245) * t2166 + t2252 + t2254;
t2024 = -t2203 * t2053 + t2206 * t2054;
t2023 = t2206 * t2053 + t2203 * t2054;
t2022 = -t2248 * pkin(5) + t2215 * qJ(6) + 0.2e1 * qJD(6) * t2164 + (t2245 + (-pkin(4) - pkin(5)) * t2186) * t2186 + t2223;
t2021 = t2224 * pkin(5) + (pkin(5) * t2164 - 0.2e1 * qJD(6) + t2134) * t2166 + t2218 - t2212 * qJ(6);
t2020 = t2207 * t2024 + t2204 * t2094;
t2019 = t2204 * t2024 - t2207 * t2094;
t2018 = t2206 * t2032 + t2203 * t2039;
t2017 = t2203 * t2032 - t2206 * t2039;
t2016 = t2207 * t2018 + t2204 * t2046;
t2015 = t2204 * t2018 - t2207 * t2046;
t2014 = t2203 * t2021 + t2206 * t2022;
t2013 = -t2206 * t2021 + t2203 * t2022;
t2012 = t2200 * t2020 + t2199 * t2023;
t2011 = t2199 * t2020 - t2200 * t2023;
t2010 = t2207 * t2014 + t2204 * t2025;
t2009 = t2204 * t2014 - t2207 * t2025;
t2008 = t2200 * t2016 + t2199 * t2017;
t2007 = t2199 * t2016 - t2200 * t2017;
t2006 = t2200 * t2010 + t2199 * t2013;
t2005 = t2199 * t2010 - t2200 * t2013;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2174, -t2173, 0, -t2205 * t2179 + t2208 * t2180, 0, 0, 0, 0, 0, 0, t2226, -t2222, 0, -t2205 * t2092 + t2208 * t2093, 0, 0, 0, 0, 0, 0, -t2205 * t2109 + t2208 * t2111, -t2205 * t2110 + t2208 * t2112, -t2205 * t2135 + t2208 * t2136, -t2205 * t2055 + t2208 * t2056, 0, 0, 0, 0, 0, 0, -t2205 * t2035 + t2208 * t2038, -t2291, t2281, -t2205 * t2011 + t2208 * t2012, 0, 0, 0, 0, 0, 0, t2289, t2281, t2291, -t2205 * t2007 + t2208 * t2008, 0, 0, 0, 0, 0, 0, t2289, t2291, -t2281, -t2205 * t2005 + t2208 * t2006; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2173, t2174, 0, t2208 * t2179 + t2205 * t2180, 0, 0, 0, 0, 0, 0, t2222, t2226, 0, t2208 * t2092 + t2205 * t2093, 0, 0, 0, 0, 0, 0, t2208 * t2109 + t2205 * t2111, t2208 * t2110 + t2205 * t2112, t2208 * t2135 + t2205 * t2136, t2208 * t2055 + t2205 * t2056, 0, 0, 0, 0, 0, 0, t2208 * t2035 + t2205 * t2038, t2292, -t2290, t2208 * t2011 + t2205 * t2012, 0, 0, 0, 0, 0, 0, t2288, -t2290, -t2292, t2208 * t2007 + t2205 * t2008, 0, 0, 0, 0, 0, 0, t2288, -t2292, t2290, t2208 * t2005 + t2205 * t2006; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2246, 0, 0, 0, 0, 0, 0, t2141, t2142, 0, t2075, 0, 0, 0, 0, 0, 0, t2059, t2063, -t2047, t2019, 0, 0, 0, 0, 0, 0, t2274, -t2047, -t2063, t2015, 0, 0, 0, 0, 0, 0, t2274, -t2063, t2047, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2209, -qJDD(1), 0, t2180, 0, 0, 0, 0, 0, 0, t2170, -t2171, 0, t2093, 0, 0, 0, 0, 0, 0, t2111, t2112, t2136, t2056, 0, 0, 0, 0, 0, 0, t2038, t2043, -t2029, t2012, 0, 0, 0, 0, 0, 0, t2282, -t2029, -t2043, t2008, 0, 0, 0, 0, 0, 0, t2282, -t2043, t2029, t2006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2209, 0, t2179, 0, 0, 0, 0, 0, 0, t2171, t2170, 0, t2092, 0, 0, 0, 0, 0, 0, t2109, t2110, t2135, t2055, 0, 0, 0, 0, 0, 0, t2035, t2040, -t2026, t2011, 0, 0, 0, 0, 0, 0, t2283, -t2026, -t2040, t2007, 0, 0, 0, 0, 0, 0, t2283, -t2040, t2026, t2005; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2246, 0, 0, 0, 0, 0, 0, t2141, t2142, 0, t2075, 0, 0, 0, 0, 0, 0, t2059, t2063, -t2047, t2019, 0, 0, 0, 0, 0, 0, t2274, -t2047, -t2063, t2015, 0, 0, 0, 0, 0, 0, t2274, -t2063, t2047, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2209, -qJDD(1), 0, t2133, 0, 0, 0, 0, 0, 0, t2143, t2144, t2172, t2076, 0, 0, 0, 0, 0, 0, t2062, t2066, -t2050, t2020, 0, 0, 0, 0, 0, 0, t2271, -t2050, -t2066, t2016, 0, 0, 0, 0, 0, 0, t2271, -t2066, t2050, t2010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2209, 0, t2132, 0, 0, 0, 0, 0, 0, t2169, -t2168, t2175, -t2123, 0, 0, 0, 0, 0, 0, -t2268, t2275, t2069, -t2023, 0, 0, 0, 0, 0, 0, -t2268, t2069, -t2275, -t2017, 0, 0, 0, 0, 0, 0, -t2268, -t2275, -t2069, -t2013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2246, 0, 0, 0, 0, 0, 0, t2141, t2142, 0, t2075, 0, 0, 0, 0, 0, 0, t2059, t2063, -t2047, t2019, 0, 0, 0, 0, 0, 0, t2274, -t2047, -t2063, t2015, 0, 0, 0, 0, 0, 0, t2274, -t2063, t2047, t2009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2183, t2177, t2230, t2108, 0, 0, 0, 0, 0, 0, t2267, t2088, -t2072, t2024, 0, 0, 0, 0, 0, 0, t2267, -t2072, -t2088, t2018, 0, 0, 0, 0, 0, 0, t2267, -t2088, t2072, t2014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2176, t2182, -t2231, t2107, 0, 0, 0, 0, 0, 0, -t2100, -t2211, t2117, -t2094, 0, 0, 0, 0, 0, 0, -t2098, t2117, t2211, -t2046, 0, 0, 0, 0, 0, 0, -t2098, t2211, -t2117, -t2025; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2169, t2168, -t2175, t2123, 0, 0, 0, 0, 0, 0, t2268, -t2275, -t2069, t2023, 0, 0, 0, 0, 0, 0, t2268, -t2069, t2275, t2017, 0, 0, 0, 0, 0, 0, t2268, t2275, t2069, t2013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2253, -t2217, -t2099, t2054, 0, 0, 0, 0, 0, 0, t2253, -t2099, t2217, t2032, 0, 0, 0, 0, 0, 0, t2253, t2217, t2099, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2255, -t2138, -t2212, t2053, 0, 0, 0, 0, 0, 0, -t2255, -t2212, t2138, -t2039, 0, 0, 0, 0, 0, 0, -t2255, t2138, t2212, -t2021; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2100, t2211, -t2117, t2094, 0, 0, 0, 0, 0, 0, t2098, -t2117, -t2211, t2046, 0, 0, 0, 0, 0, 0, t2098, -t2211, t2117, t2025; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2253, -t2099, t2217, t2032, 0, 0, 0, 0, 0, 0, t2253, t2217, t2099, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2098, -t2117, -t2211, t2046, 0, 0, 0, 0, 0, 0, t2098, -t2211, t2117, t2025; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2255, t2212, -t2138, t2039, 0, 0, 0, 0, 0, 0, t2255, -t2138, -t2212, t2021; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2253, t2217, t2099, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2255, -t2138, -t2212, t2021; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2098, t2211, -t2117, -t2025;];
f_new_reg  = t1;
