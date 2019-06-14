% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 19:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPRRP12_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:04:40
% EndTime: 2019-05-06 19:04:46
% DurationCPUTime: 6.25s
% Computational Cost: add. (18939->297), mult. (39112->302), div. (0->0), fcn. (24788->8), ass. (0->185)
t2168 = cos(qJ(2));
t2218 = qJD(1) * t2168;
t2152 = qJD(2) * t2218;
t2164 = sin(qJ(2));
t2202 = t2164 * qJDD(1);
t2129 = t2152 + t2202;
t2190 = qJDD(4) + t2129;
t2180 = qJDD(5) + t2190;
t2163 = sin(qJ(4));
t2167 = cos(qJ(4));
t2123 = qJD(2) * t2163 + t2167 * t2218;
t2125 = qJD(2) * t2167 - t2163 * t2218;
t2162 = sin(qJ(5));
t2166 = cos(qJ(5));
t2097 = t2166 * t2123 + t2125 * t2162;
t2099 = -t2123 * t2162 + t2125 * t2166;
t2215 = t2097 * t2099;
t2062 = t2180 + t2215;
t2096 = t2099 ^ 2;
t2219 = qJD(1) * t2164;
t2147 = qJD(4) + t2219;
t2142 = qJD(5) + t2147;
t2224 = t2142 ^ 2;
t2231 = -t2096 - t2224;
t2034 = t2062 * t2166 + t2162 * t2231;
t2186 = t2062 * t2162 - t2166 * t2231;
t2011 = t2163 * t2034 + t2167 * t2186;
t2259 = t2011 * t2164;
t2258 = t2011 * t2168;
t2013 = t2167 * t2034 - t2163 * t2186;
t2165 = sin(qJ(1));
t2257 = t2013 * t2165;
t2169 = cos(qJ(1));
t2256 = t2013 * t2169;
t2063 = t2180 - t2215;
t2074 = t2097 ^ 2;
t2230 = -t2224 - t2074;
t2236 = -t2063 * t2162 + t2166 * t2230;
t2241 = t2166 * t2063 + t2162 * t2230;
t2243 = -t2163 * t2241 + t2167 * t2236;
t2151 = qJD(2) * t2219;
t2201 = t2168 * qJDD(1);
t2130 = -t2151 + t2201;
t2192 = t2163 * qJDD(2) + t2167 * t2130;
t2093 = -qJD(4) * t2125 - t2192;
t2183 = t2167 * qJDD(2) - t2163 * t2130;
t2094 = -qJD(4) * t2123 + t2183;
t2193 = t2166 * t2093 - t2162 * t2094;
t2204 = qJD(5) + t2142;
t2177 = t2099 * t2204 - t2193;
t2242 = t2163 * t2236 + t2167 * t2241;
t2251 = t2164 * t2242 + t2168 * t2177;
t2255 = t2165 * t2251 - t2169 * t2243;
t2254 = t2165 * t2243 + t2169 * t2251;
t2184 = t2162 * t2093 + t2166 * t2094;
t2176 = -qJD(5) * t2097 + t2184;
t2214 = t2097 * t2142;
t2041 = t2176 + t2214;
t2178 = (-qJD(5) + t2142) * t2099 + t2193;
t2226 = t2041 * t2162 + t2166 * t2178;
t2228 = -t2166 * t2041 + t2162 * t2178;
t2235 = -t2163 * t2228 + t2167 * t2226;
t2057 = t2096 + t2074;
t2234 = t2163 * t2226 + t2167 * t2228;
t2245 = -t2057 * t2168 + t2164 * t2234;
t2250 = t2165 * t2245 - t2169 * t2235;
t2249 = t2164 * t2177 - t2168 * t2242;
t2248 = t2165 * t2235 + t2169 * t2245;
t2244 = -t2057 * t2164 - t2168 * t2234;
t2171 = qJD(1) ^ 2;
t2207 = t2168 * t2171;
t2200 = t2164 * t2207;
t2137 = -qJDD(2) + t2200;
t2158 = t2164 ^ 2;
t2170 = qJD(2) ^ 2;
t2143 = -t2158 * t2171 - t2170;
t2107 = t2137 * t2168 - t2143 * t2164;
t2128 = 0.2e1 * t2152 + t2202;
t2240 = t2107 * t2165 - t2128 * t2169;
t2239 = t2107 * t2169 + t2128 * t2165;
t2136 = qJDD(2) + t2200;
t2159 = t2168 ^ 2;
t2208 = t2159 * t2171;
t2144 = t2170 + t2208;
t2108 = t2136 * t2164 + t2144 * t2168;
t2127 = -0.2e1 * t2151 + t2201;
t2238 = t2108 * t2165 - t2127 * t2169;
t2237 = t2108 * t2169 + t2127 * t2165;
t2139 = t2165 * g(1) - t2169 * g(2);
t2181 = -qJDD(1) * pkin(1) - t2139;
t2222 = 2 * qJD(3);
t2229 = pkin(2) * t2151 - (t2129 + t2152) * qJ(3) - t2219 * t2222 + t2181;
t2227 = t2097 * t2204 - t2184;
t2121 = t2123 ^ 2;
t2225 = t2125 ^ 2;
t2223 = t2147 ^ 2;
t2221 = t2168 * g(3);
t2220 = t2171 * pkin(7);
t2213 = t2123 * t2125;
t2206 = qJD(4) - t2147;
t2205 = qJD(4) + t2147;
t2138 = pkin(3) * t2219 - qJD(2) * pkin(8);
t2059 = -t2138 * t2219 + (-pkin(3) * t2159 - pkin(7)) * t2171 + (-pkin(2) - pkin(8)) * t2130 + t2229;
t2182 = -qJDD(2) * pkin(2) - t2170 * qJ(3) + qJDD(3) + t2221;
t2140 = -g(1) * t2169 - g(2) * t2165;
t2118 = -pkin(1) * t2171 + qJDD(1) * pkin(7) + t2140;
t2198 = t2171 * (-pkin(2) * t2168 - qJ(3) * t2164) + t2118;
t2068 = -qJDD(2) * pkin(8) + (t2129 - t2152) * pkin(3) + (-pkin(8) * t2207 + t2198) * t2164 + t2182;
t2028 = t2167 * t2059 + t2163 * t2068;
t2113 = pkin(4) * t2147 - pkin(9) * t2125;
t2026 = -pkin(4) * t2121 + pkin(9) * t2093 - t2113 * t2147 + t2028;
t2179 = t2190 - t2213;
t2187 = -t2163 * t2059 + t2167 * t2068;
t2172 = (-t2123 * t2147 - t2094) * pkin(9) + t2179 * pkin(4) + t2187;
t1997 = t2166 * t2026 + t2162 * t2172;
t2203 = t2158 + t2159;
t2199 = -t2223 - t2225;
t1996 = -t2162 * t2026 + t2166 * t2172;
t2073 = pkin(5) * t2097 - qJ(6) * t2099;
t1990 = -pkin(5) * t2224 + qJ(6) * t2180 + 0.2e1 * qJD(6) * t2142 - t2097 * t2073 + t1997;
t1991 = -pkin(5) * t2180 - qJ(6) * t2224 + t2099 * t2073 + qJDD(6) - t1996;
t2189 = t1990 * t2162 - t1991 * t2166;
t2188 = t1996 * t2166 + t1997 * t2162;
t2105 = t2136 * t2168 - t2144 * t2164;
t2103 = t2137 * t2164 + t2143 * t2168;
t2175 = t2123 * t2206 - t2183;
t2155 = t2164 * g(3);
t2173 = -t2170 * pkin(2) + qJDD(2) * qJ(3) + t2168 * t2198 - t2155;
t2067 = t2130 * pkin(3) - pkin(8) * t2208 + (t2222 + t2138) * qJD(2) + t2173;
t2029 = -t2093 * pkin(4) - t2121 * pkin(9) + t2125 * t2113 + t2067;
t2134 = t2203 * t2171;
t2133 = -qJDD(1) * t2165 - t2169 * t2171;
t2132 = qJDD(1) * t2169 - t2165 * t2171;
t2131 = t2203 * qJDD(1);
t2117 = -t2181 + t2220;
t2112 = t2168 * t2118 - t2155;
t2111 = -t2164 * t2118 - t2221;
t2102 = t2131 * t2169 - t2134 * t2165;
t2101 = t2131 * t2165 + t2134 * t2169;
t2095 = -t2223 - t2121;
t2090 = -t2190 - t2213;
t2089 = -t2121 - t2225;
t2085 = t2198 * t2164 + t2182;
t2084 = qJD(2) * t2222 + t2173;
t2081 = -t2111 * t2164 + t2112 * t2168;
t2080 = t2111 * t2168 + t2112 * t2164;
t2079 = -t2123 * t2205 + t2183;
t2078 = -t2125 * t2206 - t2192;
t2077 = t2125 * t2205 + t2192;
t2075 = t2130 * pkin(2) + t2220 - t2229;
t2072 = t2167 * t2090 - t2163 * t2199;
t2071 = t2163 * t2090 + t2167 * t2199;
t2066 = t2167 * t2095 - t2163 * t2179;
t2065 = t2163 * t2095 + t2167 * t2179;
t2055 = t2084 * t2168 + t2085 * t2164;
t2054 = t2084 * t2164 - t2085 * t2168;
t2053 = t2167 * t2078 - t2163 * t2175;
t2052 = t2163 * t2078 + t2167 * t2175;
t2051 = t2071 * t2164 + t2079 * t2168;
t2050 = -t2071 * t2168 + t2079 * t2164;
t2045 = t2065 * t2164 + t2077 * t2168;
t2044 = -t2065 * t2168 + t2077 * t2164;
t2042 = t2176 - t2214;
t2033 = t2052 * t2164 + t2089 * t2168;
t2032 = -t2052 * t2168 + t2089 * t2164;
t2010 = t2167 * t2028 - t2163 * t2187;
t2009 = t2163 * t2028 + t2167 * t2187;
t2008 = -t2168 * t2227 - t2259;
t2006 = -t2164 * t2227 + t2258;
t2004 = t2009 * t2164 + t2067 * t2168;
t2003 = -t2009 * t2168 + t2067 * t2164;
t2002 = -t2042 * t2168 + t2259;
t2000 = -t2042 * t2164 - t2258;
t1998 = t2177 * pkin(5) + qJ(6) * t2227 - 0.2e1 * qJD(6) * t2099 + t2029;
t1985 = -t1996 * t2162 + t1997 * t2166;
t1983 = t1990 * t2166 + t1991 * t2162;
t1981 = t2167 * t1985 - t2163 * t2188;
t1980 = t2163 * t1985 + t2167 * t2188;
t1979 = t1980 * t2164 + t2029 * t2168;
t1978 = -t1980 * t2168 + t2029 * t2164;
t1977 = t2167 * t1983 - t2163 * t2189;
t1976 = t2163 * t1983 + t2167 * t2189;
t1975 = t1976 * t2164 + t1998 * t2168;
t1974 = -t1976 * t2168 + t1998 * t2164;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2133, -t2132, 0, -t2139 * t2165 + t2140 * t2169, 0, 0, 0, 0, 0, 0, -t2237, t2239, t2102, t2081 * t2169 - t2117 * t2165, 0, 0, 0, 0, 0, 0, t2102, t2237, -t2239, t2055 * t2169 - t2075 * t2165, 0, 0, 0, 0, 0, 0, t2045 * t2169 + t2066 * t2165, t2051 * t2169 + t2072 * t2165, t2033 * t2169 + t2053 * t2165, t2004 * t2169 + t2010 * t2165, 0, 0, 0, 0, 0, 0, t2254, t2008 * t2169 - t2257, t2248, t1979 * t2169 + t1981 * t2165, 0, 0, 0, 0, 0, 0, t2254, t2248, t2002 * t2169 + t2257, t1975 * t2169 + t1977 * t2165; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2132, t2133, 0, t2139 * t2169 + t2140 * t2165, 0, 0, 0, 0, 0, 0, -t2238, t2240, t2101, t2081 * t2165 + t2117 * t2169, 0, 0, 0, 0, 0, 0, t2101, t2238, -t2240, t2055 * t2165 + t2075 * t2169, 0, 0, 0, 0, 0, 0, t2045 * t2165 - t2066 * t2169, t2051 * t2165 - t2072 * t2169, t2033 * t2165 - t2053 * t2169, t2004 * t2165 - t2010 * t2169, 0, 0, 0, 0, 0, 0, t2255, t2008 * t2165 + t2256, t2250, t1979 * t2165 - t1981 * t2169, 0, 0, 0, 0, 0, 0, t2255, t2250, t2002 * t2165 - t2256, t1975 * t2165 - t1977 * t2169; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2105, t2103, 0, t2080, 0, 0, 0, 0, 0, 0, 0, -t2105, -t2103, t2054, 0, 0, 0, 0, 0, 0, t2044, t2050, t2032, t2003, 0, 0, 0, 0, 0, 0, t2249, t2006, t2244, t1978, 0, 0, 0, 0, 0, 0, t2249, t2244, t2000, t1974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2171, -qJDD(1), 0, t2140, 0, 0, 0, 0, 0, 0, -t2108, t2107, t2131, t2081, 0, 0, 0, 0, 0, 0, t2131, t2108, -t2107, t2055, 0, 0, 0, 0, 0, 0, t2045, t2051, t2033, t2004, 0, 0, 0, 0, 0, 0, t2251, t2008, t2245, t1979, 0, 0, 0, 0, 0, 0, t2251, t2245, t2002, t1975; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2171, 0, t2139, 0, 0, 0, 0, 0, 0, t2127, -t2128, t2134, t2117, 0, 0, 0, 0, 0, 0, t2134, -t2127, t2128, t2075, 0, 0, 0, 0, 0, 0, -t2066, -t2072, -t2053, -t2010, 0, 0, 0, 0, 0, 0, -t2243, t2013, -t2235, -t1981, 0, 0, 0, 0, 0, 0, -t2243, -t2235, -t2013, -t1977; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2105, t2103, 0, t2080, 0, 0, 0, 0, 0, 0, 0, -t2105, -t2103, t2054, 0, 0, 0, 0, 0, 0, t2044, t2050, t2032, t2003, 0, 0, 0, 0, 0, 0, t2249, t2006, t2244, t1978, 0, 0, 0, 0, 0, 0, t2249, t2244, t2000, t1974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2144, t2137, t2201, t2112, 0, 0, 0, 0, 0, 0, t2201, t2144, -t2137, t2084, 0, 0, 0, 0, 0, 0, t2077, t2079, t2089, t2067, 0, 0, 0, 0, 0, 0, t2177, -t2227, -t2057, t2029, 0, 0, 0, 0, 0, 0, t2177, -t2057, -t2042, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2136, t2143, -t2202, t2111, 0, 0, 0, 0, 0, 0, -t2202, -t2136, -t2143, -t2085, 0, 0, 0, 0, 0, 0, -t2065, -t2071, -t2052, -t2009, 0, 0, 0, 0, 0, 0, -t2242, t2011, -t2234, -t1980, 0, 0, 0, 0, 0, 0, -t2242, -t2234, -t2011, -t1976; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2127, t2128, -t2134, -t2117, 0, 0, 0, 0, 0, 0, -t2134, t2127, -t2128, -t2075, 0, 0, 0, 0, 0, 0, t2066, t2072, t2053, t2010, 0, 0, 0, 0, 0, 0, t2243, -t2013, t2235, t1981, 0, 0, 0, 0, 0, 0, t2243, t2235, t2013, t1977; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2134, t2127, -t2128, -t2075, 0, 0, 0, 0, 0, 0, t2066, t2072, t2053, t2010, 0, 0, 0, 0, 0, 0, t2243, -t2013, t2235, t1981, 0, 0, 0, 0, 0, 0, t2243, t2235, t2013, t1977; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2201, -t2144, t2137, -t2084, 0, 0, 0, 0, 0, 0, -t2077, -t2079, -t2089, -t2067, 0, 0, 0, 0, 0, 0, -t2177, t2227, t2057, -t2029, 0, 0, 0, 0, 0, 0, -t2177, t2057, t2042, -t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2202, t2136, t2143, t2085, 0, 0, 0, 0, 0, 0, t2065, t2071, t2052, t2009, 0, 0, 0, 0, 0, 0, t2242, -t2011, t2234, t1980, 0, 0, 0, 0, 0, 0, t2242, t2234, t2011, t1976; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2095, t2090, t2078, t2028, 0, 0, 0, 0, 0, 0, t2236, -t2034, t2226, t1985, 0, 0, 0, 0, 0, 0, t2236, t2226, t2034, t1983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2179, t2199, t2175, t2187, 0, 0, 0, 0, 0, 0, t2241, -t2186, t2228, t2188, 0, 0, 0, 0, 0, 0, t2241, t2228, t2186, t2189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2077, t2079, t2089, t2067, 0, 0, 0, 0, 0, 0, t2177, -t2227, -t2057, t2029, 0, 0, 0, 0, 0, 0, t2177, -t2057, -t2042, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2230, -t2062, t2178, t1997, 0, 0, 0, 0, 0, 0, t2230, t2178, t2062, t1990; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2063, t2231, -t2041, t1996, 0, 0, 0, 0, 0, 0, t2063, -t2041, -t2231, -t1991; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2177, -t2227, -t2057, t2029, 0, 0, 0, 0, 0, 0, t2177, -t2057, -t2042, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2230, t2178, t2062, t1990; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2177, -t2057, -t2042, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2063, t2041, t2231, t1991;];
f_new_reg  = t1;