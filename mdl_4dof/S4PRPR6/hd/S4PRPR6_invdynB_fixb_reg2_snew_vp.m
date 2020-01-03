% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRPR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:43
% EndTime: 2019-12-31 16:24:47
% DurationCPUTime: 2.89s
% Computational Cost: add. (6541->306), mult. (14725->479), div. (0->0), fcn. (10288->8), ass. (0->219)
t306 = sin(pkin(6));
t308 = cos(pkin(6));
t287 = g(1) * t308 + g(2) * t306;
t310 = sin(qJ(2));
t312 = cos(qJ(2));
t340 = g(3) - qJDD(1);
t257 = -t287 * t312 - t310 * t340;
t357 = qJD(2) ^ 2;
t363 = -pkin(2) * t357 + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t257;
t307 = cos(pkin(7));
t305 = sin(pkin(7));
t301 = t305 ^ 2;
t302 = t307 ^ 2;
t358 = t357 * (t301 + t302);
t273 = t307 * t358;
t330 = t310 * qJDD(2);
t251 = t273 * t312 + t307 * t330;
t367 = t306 * t251;
t366 = t308 * t251;
t309 = sin(qJ(4));
t311 = cos(qJ(4));
t354 = t305 * t309;
t267 = (-t307 * t311 + t354) * qJD(2);
t317 = t305 * t311 + t307 * t309;
t269 = t317 * qJD(2);
t233 = t269 * t267;
t359 = qJDD(4) - t233;
t365 = t309 * t359;
t364 = t311 * t359;
t362 = t306 * t340;
t361 = t308 * t340;
t286 = g(1) * t306 - g(2) * t308;
t276 = t308 * t286;
t243 = -t306 * t287 + t276;
t298 = t302 * t357;
t337 = t301 * t357;
t281 = t298 + t337;
t360 = t317 * qJDD(2);
t263 = t267 ^ 2;
t264 = t269 ^ 2;
t304 = qJDD(2) * pkin(2);
t275 = t307 * t286;
t200 = -t275 + (pkin(3) * t307 * t357 - pkin(5) * qJDD(2) - t363) * t305;
t216 = -t305 * t286 + t307 * t363;
t334 = qJDD(2) * t307;
t203 = -pkin(3) * t298 + pkin(5) * t334 + t216;
t167 = -t200 * t311 + t203 * t309;
t168 = t200 * t309 + t203 * t311;
t138 = -t167 * t311 + t168 * t309;
t356 = t305 * t138;
t355 = t305 * t307;
t329 = t312 * qJDD(2);
t324 = t307 * t329;
t336 = t310 * t357;
t255 = t305 * t324 - t336 * t355;
t353 = t306 * t255;
t285 = t329 - t336;
t352 = t306 * t285;
t351 = t306 * t286;
t349 = t307 * t138;
t348 = t308 * t255;
t347 = t308 * t285;
t256 = -t287 * t310 + t312 * t340;
t242 = -qJ(3) * t357 + qJDD(3) + t256 - t304;
t221 = -pkin(3) * t334 - pkin(5) * t281 + t242;
t346 = t309 * t221;
t226 = qJDD(4) + t233;
t345 = t309 * t226;
t344 = t310 * t242;
t343 = t311 * t221;
t342 = t311 * t226;
t341 = t312 * t242;
t339 = t267 * qJD(4);
t338 = t269 * qJD(4);
t333 = t301 * qJDD(2);
t332 = t306 * qJDD(2);
t331 = t308 * qJDD(2);
t327 = t310 * t233;
t326 = t312 * t233;
t325 = t305 * t334;
t323 = -t242 + t304;
t215 = t305 * t363 + t275;
t185 = t215 * t305 + t216 * t307;
t297 = t302 * qJDD(2);
t279 = t297 + t333;
t237 = t279 * t310 + t281 * t312;
t238 = t279 * t312 - t281 * t310;
t322 = -pkin(1) * t237 - pkin(2) * t281 + qJ(1) * t238 - qJ(3) * t279 - t185;
t272 = t305 * t358;
t246 = t272 * t310 - t305 * t329;
t249 = t272 * t312 + t305 * t330;
t321 = -pkin(1) * t246 + qJ(1) * t249 - qJ(3) * t272 + t305 * t323;
t248 = -t273 * t310 + t324;
t320 = -pkin(1) * t248 - qJ(1) * t251 + qJ(3) * t273 - t307 * t323;
t284 = t312 * t357 + t330;
t319 = pkin(1) * t285 + qJ(1) * t284 - t256;
t318 = -pkin(1) * t284 + qJ(1) * t285 - t257;
t139 = t167 * t309 + t168 * t311;
t220 = t256 * t310 + t257 * t312;
t244 = -t287 * t308 - t351;
t250 = pkin(4) * t284 - t286 * t312;
t247 = -pkin(4) * t285 - t286 * t310;
t265 = qJDD(2) * t354 - t311 * t334;
t184 = t215 * t307 - t216 * t305;
t219 = t256 * t312 - t257 * t310;
t313 = qJD(4) ^ 2;
t282 = t298 - t337;
t280 = t297 - t333;
t271 = t308 * t284;
t270 = t306 * t284;
t260 = -t264 - t313;
t259 = -t264 + t313;
t258 = t263 - t313;
t254 = t284 * t355;
t241 = t308 * t249;
t240 = t306 * t249;
t239 = t280 * t312 - t282 * t310;
t232 = -t264 + t263;
t231 = t360 - t339;
t230 = t360 - 0.2e1 * t339;
t229 = -t265 - t338;
t228 = t265 + 0.2e1 * t338;
t224 = -t313 - t263;
t223 = (-t267 * t311 + t269 * t309) * qJD(4);
t222 = (-t267 * t309 - t269 * t311) * qJD(4);
t217 = -t263 - t264;
t213 = t231 * t311 - t309 * t338;
t212 = t231 * t309 + t311 * t338;
t211 = -t229 * t309 + t311 * t339;
t210 = t229 * t311 + t309 * t339;
t209 = -t260 * t309 - t342;
t208 = -t259 * t309 + t364;
t207 = t258 * t311 - t345;
t206 = t260 * t311 - t345;
t205 = t259 * t311 + t365;
t204 = t258 * t309 + t342;
t202 = t220 * t308 - t351;
t201 = t220 * t306 + t276;
t196 = -t228 * t311 - t230 * t309;
t195 = -t265 * t311 + t309 * t360;
t194 = -t228 * t309 + t230 * t311;
t193 = -t265 * t309 - t311 * t360;
t192 = t224 * t311 - t365;
t191 = t224 * t309 + t364;
t188 = -t222 * t305 + t223 * t307;
t187 = -t222 * t307 - t223 * t305;
t186 = qJDD(4) * t310 + t188 * t312;
t182 = -pkin(5) * t206 + t343;
t181 = -pkin(4) * t246 - t216 * t310 + t307 * t341;
t180 = -pkin(4) * t248 - t215 * t310 + t305 * t341;
t179 = -t212 * t305 + t213 * t307;
t178 = -t210 * t305 + t211 * t307;
t177 = -t212 * t307 - t213 * t305;
t176 = -t210 * t307 - t211 * t305;
t175 = -pkin(5) * t191 + t346;
t174 = -t206 * t305 + t209 * t307;
t173 = -t205 * t305 + t208 * t307;
t172 = -t204 * t305 + t207 * t307;
t171 = t206 * t307 + t209 * t305;
t170 = -t205 * t307 - t208 * t305;
t169 = -t204 * t307 - t207 * t305;
t165 = -pkin(4) * t237 + t184 * t312;
t164 = t185 * t312 + t344;
t163 = t185 * t310 - t341;
t162 = -pkin(3) * t230 + pkin(5) * t209 + t346;
t161 = t173 * t312 + t310 * t360;
t160 = t172 * t312 - t265 * t310;
t159 = -t194 * t305 + t196 * t307;
t158 = -t193 * t305 + t195 * t307;
t157 = -t194 * t307 - t196 * t305;
t156 = t193 * t307 + t195 * t305;
t155 = -t191 * t305 + t192 * t307;
t154 = t191 * t307 + t192 * t305;
t153 = -pkin(3) * t228 + pkin(5) * t192 - t343;
t152 = t179 * t312 + t327;
t151 = t178 * t312 - t327;
t149 = t174 * t312 + t230 * t310;
t148 = t174 * t310 - t230 * t312;
t147 = t159 * t312 - t232 * t310;
t146 = t155 * t312 + t228 * t310;
t145 = t155 * t310 - t228 * t312;
t144 = t158 * t312 + t217 * t310;
t143 = t158 * t310 - t217 * t312;
t142 = -pkin(2) * t156 - pkin(3) * t193;
t141 = t164 * t308 - t184 * t306;
t140 = t164 * t306 + t184 * t308;
t137 = -pkin(1) * t163 + pkin(2) * t242 - qJ(3) * t185;
t136 = t149 * t308 + t171 * t306;
t135 = t149 * t306 - t171 * t308;
t134 = -pkin(3) * t221 + pkin(5) * t139;
t133 = -pkin(2) * t171 - pkin(3) * t206 + t168;
t132 = -pkin(5) * t193 - t138;
t131 = -pkin(4) * t163 - (pkin(2) * t310 - qJ(3) * t312) * t184;
t130 = t146 * t308 + t154 * t306;
t129 = t146 * t306 - t154 * t308;
t128 = -pkin(2) * t154 - pkin(3) * t191 + t167;
t127 = t144 * t308 + t156 * t306;
t126 = t144 * t306 - t156 * t308;
t125 = -qJ(3) * t171 - t162 * t305 + t182 * t307;
t124 = -pkin(3) * t217 + pkin(5) * t195 + t139;
t123 = -qJ(3) * t154 - t153 * t305 + t175 * t307;
t122 = t139 * t307 - t356;
t121 = t139 * t305 + t349;
t120 = -pkin(1) * t148 + pkin(2) * t230 - qJ(3) * t174 - t162 * t307 - t182 * t305;
t119 = t122 * t312 + t221 * t310;
t118 = t122 * t310 - t221 * t312;
t117 = -pkin(1) * t145 + pkin(2) * t228 - qJ(3) * t155 - t153 * t307 - t175 * t305;
t116 = -pkin(2) * t121 - pkin(3) * t138;
t115 = -pkin(4) * t148 + t125 * t312 - t133 * t310;
t114 = -qJ(3) * t156 - t124 * t305 + t132 * t307;
t113 = -pkin(4) * t145 + t123 * t312 - t128 * t310;
t112 = -pkin(5) * t349 - qJ(3) * t121 - t134 * t305;
t111 = t119 * t308 + t121 * t306;
t110 = t119 * t306 - t121 * t308;
t109 = -pkin(1) * t143 + pkin(2) * t217 - qJ(3) * t158 - t124 * t307 - t132 * t305;
t108 = -pkin(4) * t143 + t114 * t312 - t142 * t310;
t107 = -pkin(1) * t118 + pkin(2) * t221 + pkin(5) * t356 - qJ(3) * t122 - t134 * t307;
t106 = -pkin(4) * t118 + t112 * t312 - t116 * t310;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, 0, 0, 0, 0, 0, -t271, -t347, 0, t202, 0, 0, 0, 0, 0, 0, -t366, t241, t308 * t238, t141, 0, 0, 0, 0, 0, 0, t130, t136, t127, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, 0, 0, 0, 0, 0, 0, -t270, -t352, 0, t201, 0, 0, 0, 0, 0, 0, -t367, t240, t306 * t238, t140, 0, 0, 0, 0, 0, 0, t129, t135, t126, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, 0, 0, 0, 0, 0, 0, t285, -t284, 0, -t219, 0, 0, 0, 0, 0, 0, t248, t246, t237, t163, 0, 0, 0, 0, 0, 0, t145, t148, t143, t118; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t362, -t361, -t243, -qJ(1) * t243, 0, 0, t347, 0, -t271, t332, t308 * t247 + t306 * t319, t308 * t250 + t306 * t318, t308 * t219, -qJ(1) * t201 - (pkin(1) * t306 - pkin(4) * t308) * t219, t301 * t332 + t348, t239 * t308 + 0.2e1 * t306 * t325, t241, t302 * t332 - t348, t366, 0, t308 * t180 - t306 * t320, t308 * t181 - t306 * t321, t308 * t165 - t306 * t322, -qJ(1) * t140 + t131 * t308 - t137 * t306, t152 * t308 - t177 * t306, t147 * t308 - t157 * t306, t161 * t308 - t170 * t306, t151 * t308 - t176 * t306, t160 * t308 - t169 * t306, t186 * t308 - t187 * t306, -qJ(1) * t129 + t113 * t308 - t117 * t306, -qJ(1) * t135 + t115 * t308 - t120 * t306, -qJ(1) * t126 + t108 * t308 - t109 * t306, -qJ(1) * t110 + t106 * t308 - t107 * t306; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t361, -t362, t244, qJ(1) * t244, 0, 0, t352, 0, -t270, -t331, t306 * t247 - t308 * t319, t306 * t250 - t308 * t318, t306 * t219, qJ(1) * t202 - (-pkin(1) * t308 - pkin(4) * t306) * t219, -t301 * t331 + t353, t239 * t306 - 0.2e1 * t308 * t325, t240, -t302 * t331 - t353, t367, 0, t306 * t180 + t308 * t320, t306 * t181 + t308 * t321, t306 * t165 + t308 * t322, qJ(1) * t141 + t131 * t306 + t137 * t308, t152 * t306 + t177 * t308, t147 * t306 + t157 * t308, t161 * t306 + t170 * t308, t151 * t306 + t176 * t308, t160 * t306 + t169 * t308, t186 * t306 + t187 * t308, qJ(1) * t130 + t113 * t306 + t117 * t308, qJ(1) * t136 + t115 * t306 + t120 * t308, qJ(1) * t127 + t108 * t306 + t109 * t308, qJ(1) * t111 + t106 * t306 + t107 * t308; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t286, t287, 0, 0, 0, 0, t284, 0, t285, 0, -t250, t247, t220, pkin(1) * t286 + pkin(4) * t220, t254, t280 * t310 + t282 * t312, t246, -t254, -t248, 0, -pkin(4) * t251 + t215 * t312 + t305 * t344, pkin(4) * t249 + t216 * t312 + t307 * t344, pkin(4) * t238 + t184 * t310, pkin(4) * t164 - (-pkin(2) * t312 - qJ(3) * t310 - pkin(1)) * t184, t179 * t310 - t326, t159 * t310 + t232 * t312, t173 * t310 - t312 * t360, t178 * t310 + t326, t172 * t310 + t265 * t312, -qJDD(4) * t312 + t188 * t310, -pkin(1) * t154 + pkin(4) * t146 + t123 * t310 + t128 * t312, -pkin(1) * t171 + pkin(4) * t149 + t125 * t310 + t133 * t312, -pkin(1) * t156 + pkin(4) * t144 + t114 * t310 + t142 * t312, -pkin(1) * t121 + pkin(4) * t119 + t112 * t310 + t116 * t312;];
tauB_reg = t1;
