% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPPRP5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:29:26
% EndTime: 2019-05-06 09:29:31
% DurationCPUTime: 6.30s
% Computational Cost: add. (17075->293), mult. (37887->299), div. (0->0), fcn. (23662->8), ass. (0->182)
t2162 = cos(qJ(2));
t2210 = qJD(1) * t2162;
t2146 = qJD(2) * t2210;
t2159 = sin(qJ(2));
t2195 = t2159 * qJDD(1);
t2123 = t2146 + t2195;
t2182 = qJDD(5) + t2123;
t2155 = sin(pkin(9));
t2156 = cos(pkin(9));
t2117 = qJD(2) * t2155 + t2156 * t2210;
t2119 = qJD(2) * t2156 - t2155 * t2210;
t2158 = sin(qJ(5));
t2161 = cos(qJ(5));
t2086 = t2161 * t2117 + t2119 * t2158;
t2088 = -t2117 * t2158 + t2119 * t2161;
t2207 = t2086 * t2088;
t2057 = t2182 + t2207;
t2085 = t2088 ^ 2;
t2211 = qJD(1) * t2159;
t2141 = qJD(5) + t2211;
t2216 = t2141 ^ 2;
t2223 = -t2085 - t2216;
t2025 = t2057 * t2161 + t2158 * t2223;
t2179 = t2057 * t2158 - t2161 * t2223;
t2005 = t2155 * t2025 + t2156 * t2179;
t2251 = t2005 * t2159;
t2250 = t2005 * t2162;
t2007 = t2156 * t2025 - t2155 * t2179;
t2160 = sin(qJ(1));
t2249 = t2007 * t2160;
t2163 = cos(qJ(1));
t2248 = t2007 * t2163;
t2058 = t2182 - t2207;
t2066 = t2086 ^ 2;
t2222 = -t2216 - t2066;
t2228 = -t2058 * t2158 + t2161 * t2222;
t2233 = t2161 * t2058 + t2158 * t2222;
t2235 = -t2155 * t2233 + t2156 * t2228;
t2145 = qJD(2) * t2211;
t2194 = t2162 * qJDD(1);
t2124 = -t2145 + t2194;
t2105 = -qJDD(2) * t2155 - t2156 * t2124;
t2106 = qJDD(2) * t2156 - t2124 * t2155;
t2185 = t2161 * t2105 - t2158 * t2106;
t2197 = qJD(5) + t2141;
t2172 = t2088 * t2197 - t2185;
t2234 = t2155 * t2228 + t2156 * t2233;
t2243 = t2159 * t2234 + t2162 * t2172;
t2247 = t2160 * t2243 - t2163 * t2235;
t2246 = t2160 * t2235 + t2163 * t2243;
t2177 = t2158 * t2105 + t2161 * t2106;
t2170 = -qJD(5) * t2086 + t2177;
t2206 = t2086 * t2141;
t2041 = t2170 + t2206;
t2173 = (-qJD(5) + t2141) * t2088 + t2185;
t2218 = t2041 * t2158 + t2161 * t2173;
t2220 = -t2161 * t2041 + t2158 * t2173;
t2227 = -t2155 * t2220 + t2156 * t2218;
t2049 = t2085 + t2066;
t2226 = t2155 * t2218 + t2156 * t2220;
t2237 = -t2049 * t2162 + t2159 * t2226;
t2242 = t2160 * t2237 - t2163 * t2227;
t2241 = t2159 * t2172 - t2162 * t2234;
t2240 = t2160 * t2227 + t2163 * t2237;
t2236 = -t2049 * t2159 - t2162 * t2226;
t2165 = qJD(1) ^ 2;
t2198 = t2162 * t2165;
t2191 = t2159 * t2198;
t2131 = -qJDD(2) + t2191;
t2164 = qJD(2) ^ 2;
t2152 = t2159 ^ 2;
t2200 = t2152 * t2165;
t2137 = -t2164 - t2200;
t2097 = t2131 * t2162 - t2137 * t2159;
t2122 = 0.2e1 * t2146 + t2195;
t2232 = t2097 * t2160 - t2122 * t2163;
t2231 = t2097 * t2163 + t2122 * t2160;
t2130 = qJDD(2) + t2191;
t2153 = t2162 ^ 2;
t2199 = t2153 * t2165;
t2138 = t2164 + t2199;
t2098 = t2130 * t2159 + t2138 * t2162;
t2125 = -0.2e1 * t2145 + t2194;
t2230 = t2098 * t2160 - t2125 * t2163;
t2229 = t2098 * t2163 + t2125 * t2160;
t2133 = t2160 * g(1) - t2163 * g(2);
t2175 = -qJDD(1) * pkin(1) - t2133;
t2215 = 2 * qJD(3);
t2221 = pkin(2) * t2145 - (t2123 + t2146) * qJ(3) - t2211 * t2215 + t2175;
t2219 = t2086 * t2197 - t2177;
t2116 = t2117 ^ 2;
t2217 = t2119 ^ 2;
t2214 = -2 * qJD(4);
t2213 = t2162 * g(3);
t2212 = t2165 * pkin(7);
t2205 = t2117 * t2119;
t2132 = pkin(3) * t2211 - qJD(2) * qJ(4);
t2051 = -t2132 * t2211 + (-pkin(3) * t2153 - pkin(7)) * t2165 + (-pkin(2) - qJ(4)) * t2124 + t2221;
t2176 = -qJDD(2) * pkin(2) - t2164 * qJ(3) + qJDD(3) + t2213;
t2134 = -g(1) * t2163 - g(2) * t2160;
t2112 = -pkin(1) * t2165 + qJDD(1) * pkin(7) + t2134;
t2190 = t2165 * (-pkin(2) * t2162 - qJ(3) * t2159) + t2112;
t2060 = -qJDD(2) * qJ(4) + (t2123 - t2146) * pkin(3) + (-qJ(4) * t2198 + t2190) * t2159 + t2176;
t2022 = t2156 * t2051 + t2155 * t2060 + t2117 * t2214;
t2107 = pkin(4) * t2211 - pkin(8) * t2119;
t2018 = -pkin(4) * t2116 + pkin(8) * t2105 - t2107 * t2211 + t2022;
t2193 = t2117 * t2211;
t2169 = -t2106 - t2193;
t2171 = -t2155 * t2051 + t2156 * t2060 + t2119 * t2214;
t2174 = t2123 - t2205;
t2166 = pkin(4) * t2174 + pkin(8) * t2169 + t2171;
t1987 = t2161 * t2018 + t2158 * t2166;
t2196 = t2152 + t2153;
t2192 = t2119 * t2211;
t1986 = -t2158 * t2018 + t2161 * t2166;
t2184 = -t2200 - t2217;
t2065 = pkin(5) * t2086 - qJ(6) * t2088;
t1980 = -pkin(5) * t2216 + qJ(6) * t2182 + 0.2e1 * qJD(6) * t2141 - t2086 * t2065 + t1987;
t1981 = -pkin(5) * t2182 - qJ(6) * t2216 + t2088 * t2065 + qJDD(6) - t1986;
t2181 = t1980 * t2158 - t1981 * t2161;
t2180 = t1986 * t2161 + t1987 * t2158;
t2095 = t2130 * t2162 - t2138 * t2159;
t2093 = t2131 * t2159 + t2137 * t2162;
t2149 = t2159 * g(3);
t2167 = -t2164 * pkin(2) + qJDD(2) * qJ(3) + t2162 * t2190 - t2149;
t2056 = qJDD(4) + t2124 * pkin(3) - qJ(4) * t2199 + (t2215 + t2132) * qJD(2) + t2167;
t2027 = -t2105 * pkin(4) - t2116 * pkin(8) + t2119 * t2107 + t2056;
t2129 = t2196 * t2165;
t2128 = -qJDD(1) * t2160 - t2163 * t2165;
t2127 = qJDD(1) * t2163 - t2160 * t2165;
t2126 = t2196 * qJDD(1);
t2111 = -t2175 + t2212;
t2102 = t2162 * t2112 - t2149;
t2101 = -t2159 * t2112 - t2213;
t2091 = t2126 * t2163 - t2129 * t2160;
t2090 = t2126 * t2160 + t2129 * t2163;
t2089 = -t2200 - t2116;
t2084 = -t2123 - t2205;
t2083 = t2106 - t2193;
t2082 = t2105 + t2192;
t2081 = -t2105 + t2192;
t2078 = -t2116 - t2217;
t2075 = t2190 * t2159 + t2176;
t2072 = qJD(2) * t2215 + t2167;
t2071 = -t2101 * t2159 + t2102 * t2162;
t2070 = t2101 * t2162 + t2102 * t2159;
t2069 = t2124 * pkin(2) + t2212 - t2221;
t2068 = t2156 * t2084 - t2155 * t2184;
t2067 = t2155 * t2084 + t2156 * t2184;
t2062 = t2156 * t2089 - t2155 * t2174;
t2061 = t2155 * t2089 + t2156 * t2174;
t2053 = t2156 * t2082 - t2155 * t2169;
t2052 = t2155 * t2082 + t2156 * t2169;
t2047 = t2072 * t2162 + t2075 * t2159;
t2046 = t2072 * t2159 - t2075 * t2162;
t2045 = t2067 * t2159 + t2083 * t2162;
t2044 = -t2067 * t2162 + t2083 * t2159;
t2043 = t2061 * t2159 + t2081 * t2162;
t2042 = -t2061 * t2162 + t2081 * t2159;
t2039 = t2170 - t2206;
t2033 = t2052 * t2159 + t2078 * t2162;
t2032 = -t2052 * t2162 + t2078 * t2159;
t2004 = t2156 * t2022 - t2155 * t2171;
t2003 = t2155 * t2022 + t2156 * t2171;
t2002 = -t2162 * t2219 - t2251;
t2000 = -t2159 * t2219 + t2250;
t1998 = t2172 * pkin(5) + qJ(6) * t2219 - 0.2e1 * qJD(6) * t2088 + t2027;
t1997 = -t2039 * t2162 + t2251;
t1995 = -t2039 * t2159 - t2250;
t1993 = t2003 * t2159 + t2056 * t2162;
t1992 = -t2003 * t2162 + t2056 * t2159;
t1979 = -t1986 * t2158 + t1987 * t2161;
t1977 = t1980 * t2161 + t1981 * t2158;
t1975 = t2156 * t1979 - t2155 * t2180;
t1974 = t2155 * t1979 + t2156 * t2180;
t1973 = t1974 * t2159 + t2027 * t2162;
t1972 = -t1974 * t2162 + t2027 * t2159;
t1971 = t2156 * t1977 - t2155 * t2181;
t1970 = t2155 * t1977 + t2156 * t2181;
t1969 = t1970 * t2159 + t1998 * t2162;
t1968 = -t1970 * t2162 + t1998 * t2159;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2128, -t2127, 0, -t2133 * t2160 + t2134 * t2163, 0, 0, 0, 0, 0, 0, -t2229, t2231, t2091, t2071 * t2163 - t2111 * t2160, 0, 0, 0, 0, 0, 0, t2091, t2229, -t2231, t2047 * t2163 - t2069 * t2160, 0, 0, 0, 0, 0, 0, t2043 * t2163 + t2062 * t2160, t2045 * t2163 + t2068 * t2160, t2033 * t2163 + t2053 * t2160, t1993 * t2163 + t2004 * t2160, 0, 0, 0, 0, 0, 0, t2246, t2002 * t2163 - t2249, t2240, t1973 * t2163 + t1975 * t2160, 0, 0, 0, 0, 0, 0, t2246, t2240, t1997 * t2163 + t2249, t1969 * t2163 + t1971 * t2160; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2127, t2128, 0, t2133 * t2163 + t2134 * t2160, 0, 0, 0, 0, 0, 0, -t2230, t2232, t2090, t2071 * t2160 + t2111 * t2163, 0, 0, 0, 0, 0, 0, t2090, t2230, -t2232, t2047 * t2160 + t2069 * t2163, 0, 0, 0, 0, 0, 0, t2043 * t2160 - t2062 * t2163, t2045 * t2160 - t2068 * t2163, t2033 * t2160 - t2053 * t2163, t1993 * t2160 - t2004 * t2163, 0, 0, 0, 0, 0, 0, t2247, t2002 * t2160 + t2248, t2242, t1973 * t2160 - t1975 * t2163, 0, 0, 0, 0, 0, 0, t2247, t2242, t1997 * t2160 - t2248, t1969 * t2160 - t1971 * t2163; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2095, t2093, 0, t2070, 0, 0, 0, 0, 0, 0, 0, -t2095, -t2093, t2046, 0, 0, 0, 0, 0, 0, t2042, t2044, t2032, t1992, 0, 0, 0, 0, 0, 0, t2241, t2000, t2236, t1972, 0, 0, 0, 0, 0, 0, t2241, t2236, t1995, t1968; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2165, -qJDD(1), 0, t2134, 0, 0, 0, 0, 0, 0, -t2098, t2097, t2126, t2071, 0, 0, 0, 0, 0, 0, t2126, t2098, -t2097, t2047, 0, 0, 0, 0, 0, 0, t2043, t2045, t2033, t1993, 0, 0, 0, 0, 0, 0, t2243, t2002, t2237, t1973, 0, 0, 0, 0, 0, 0, t2243, t2237, t1997, t1969; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2165, 0, t2133, 0, 0, 0, 0, 0, 0, t2125, -t2122, t2129, t2111, 0, 0, 0, 0, 0, 0, t2129, -t2125, t2122, t2069, 0, 0, 0, 0, 0, 0, -t2062, -t2068, -t2053, -t2004, 0, 0, 0, 0, 0, 0, -t2235, t2007, -t2227, -t1975, 0, 0, 0, 0, 0, 0, -t2235, -t2227, -t2007, -t1971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2095, t2093, 0, t2070, 0, 0, 0, 0, 0, 0, 0, -t2095, -t2093, t2046, 0, 0, 0, 0, 0, 0, t2042, t2044, t2032, t1992, 0, 0, 0, 0, 0, 0, t2241, t2000, t2236, t1972, 0, 0, 0, 0, 0, 0, t2241, t2236, t1995, t1968; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2138, t2131, t2194, t2102, 0, 0, 0, 0, 0, 0, t2194, t2138, -t2131, t2072, 0, 0, 0, 0, 0, 0, t2081, t2083, t2078, t2056, 0, 0, 0, 0, 0, 0, t2172, -t2219, -t2049, t2027, 0, 0, 0, 0, 0, 0, t2172, -t2049, -t2039, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2130, t2137, -t2195, t2101, 0, 0, 0, 0, 0, 0, -t2195, -t2130, -t2137, -t2075, 0, 0, 0, 0, 0, 0, -t2061, -t2067, -t2052, -t2003, 0, 0, 0, 0, 0, 0, -t2234, t2005, -t2226, -t1974, 0, 0, 0, 0, 0, 0, -t2234, -t2226, -t2005, -t1970; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2125, t2122, -t2129, -t2111, 0, 0, 0, 0, 0, 0, -t2129, t2125, -t2122, -t2069, 0, 0, 0, 0, 0, 0, t2062, t2068, t2053, t2004, 0, 0, 0, 0, 0, 0, t2235, -t2007, t2227, t1975, 0, 0, 0, 0, 0, 0, t2235, t2227, t2007, t1971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2129, t2125, -t2122, -t2069, 0, 0, 0, 0, 0, 0, t2062, t2068, t2053, t2004, 0, 0, 0, 0, 0, 0, t2235, -t2007, t2227, t1975, 0, 0, 0, 0, 0, 0, t2235, t2227, t2007, t1971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2194, -t2138, t2131, -t2072, 0, 0, 0, 0, 0, 0, -t2081, -t2083, -t2078, -t2056, 0, 0, 0, 0, 0, 0, -t2172, t2219, t2049, -t2027, 0, 0, 0, 0, 0, 0, -t2172, t2049, t2039, -t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2195, t2130, t2137, t2075, 0, 0, 0, 0, 0, 0, t2061, t2067, t2052, t2003, 0, 0, 0, 0, 0, 0, t2234, -t2005, t2226, t1974, 0, 0, 0, 0, 0, 0, t2234, t2226, t2005, t1970; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2089, t2084, t2082, t2022, 0, 0, 0, 0, 0, 0, t2228, -t2025, t2218, t1979, 0, 0, 0, 0, 0, 0, t2228, t2218, t2025, t1977; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2174, t2184, t2169, t2171, 0, 0, 0, 0, 0, 0, t2233, -t2179, t2220, t2180, 0, 0, 0, 0, 0, 0, t2233, t2220, t2179, t2181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2081, t2083, t2078, t2056, 0, 0, 0, 0, 0, 0, t2172, -t2219, -t2049, t2027, 0, 0, 0, 0, 0, 0, t2172, -t2049, -t2039, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2222, -t2057, t2173, t1987, 0, 0, 0, 0, 0, 0, t2222, t2173, t2057, t1980; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2058, t2223, -t2041, t1986, 0, 0, 0, 0, 0, 0, t2058, -t2041, -t2223, -t1981; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2172, -t2219, -t2049, t2027, 0, 0, 0, 0, 0, 0, t2172, -t2049, -t2039, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2222, t2173, t2057, t1980; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2172, -t2049, -t2039, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2058, t2041, t2223, t1981;];
f_new_reg  = t1;