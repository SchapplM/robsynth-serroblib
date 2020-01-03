% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:16
% EndTime: 2019-12-31 18:01:22
% DurationCPUTime: 3.97s
% Computational Cost: add. (12547->313), mult. (18911->434), div. (0->0), fcn. (8076->8), ass. (0->200)
t305 = -qJD(1) + qJD(4);
t303 = t305 ^ 2;
t316 = cos(qJ(4));
t304 = qJDD(1) - qJDD(4);
t313 = sin(qJ(4));
t348 = t313 * t304;
t271 = t316 * t303 - t348;
t341 = t316 * t304;
t273 = t313 * t303 + t341;
t310 = sin(pkin(8));
t311 = cos(pkin(8));
t229 = t310 * t271 + t311 * t273;
t309 = g(3) + qJDD(3);
t253 = pkin(6) * t271 + t316 * t309;
t367 = pkin(6) * t273 + t313 * t309;
t188 = qJ(3) * t229 + t310 * t253 + t311 * t367;
t314 = sin(qJ(1));
t317 = cos(qJ(1));
t363 = -t311 * t271 + t310 * t273;
t373 = -qJ(3) * t363 + t311 * t253 - t310 * t367;
t374 = t317 * t229 + t314 * t363;
t379 = -pkin(5) * t374 + t317 * t188 - t314 * t373;
t194 = t314 * t229 - t317 * t363;
t378 = -pkin(5) * t194 + t314 * t188 + t317 * t373;
t319 = qJD(1) ^ 2;
t290 = t314 * g(1) - t317 * g(2);
t330 = -qJDD(2) + t290;
t324 = -t319 * qJ(2) - t330;
t358 = pkin(1) + pkin(2);
t322 = -t358 * qJDD(1) + t324;
t306 = qJDD(1) * qJ(2);
t291 = t317 * g(1) + t314 * g(2);
t326 = (2 * qJD(2) * qJD(1)) - t291;
t325 = t306 + t326;
t323 = -t358 * t319 + t325;
t213 = t310 * t322 + t311 * t323;
t207 = -t319 * pkin(3) + t213;
t321 = t310 * t323;
t320 = -t321 + t311 * t324 + (-t311 * t358 - pkin(3)) * qJDD(1);
t172 = t313 * t207 - t316 * t320;
t173 = t316 * t207 + t313 * t320;
t334 = t313 * t172 + t316 * t173;
t148 = t316 * t172 - t313 * t173;
t353 = t311 * t148;
t132 = -t310 * t334 + t353;
t354 = t310 * t148;
t133 = t311 * t334 + t354;
t122 = t317 * t132 + t314 * t133;
t375 = t314 * t132 - t317 * t133;
t281 = -t310 * qJDD(1) + t311 * t319;
t282 = t311 * qJDD(1) + t310 * t319;
t233 = t317 * t281 + t314 * t282;
t257 = qJ(3) * t282 + t310 * t309;
t329 = qJ(3) * t281 + t311 * t309;
t369 = -pkin(5) * t233 + t314 * t257 + t317 * t329;
t333 = -t314 * t281 + t317 * t282;
t365 = -pkin(5) * t333 + t317 * t257 - t314 * t329;
t212 = -t311 * t322 + t321;
t180 = t311 * t212 - t310 * t213;
t181 = t310 * t212 + t311 * t213;
t154 = t317 * t180 + t314 * t181;
t364 = t314 * t180 - t317 * t181;
t357 = qJ(2) * t309;
t356 = qJDD(1) * pkin(1);
t312 = sin(qJ(5));
t307 = t312 ^ 2;
t355 = t307 * t303;
t167 = t304 * pkin(4) - t303 * pkin(7) + t172;
t352 = t312 * t167;
t315 = cos(qJ(5));
t289 = t315 * t303 * t312;
t279 = qJDD(5) + t289;
t351 = t312 * t279;
t280 = qJDD(5) - t289;
t350 = t312 * t280;
t349 = t312 * t304;
t344 = t315 * t167;
t343 = t315 * t279;
t342 = t315 * t280;
t293 = t315 * t304;
t168 = -t303 * pkin(4) - t304 * pkin(7) + t173;
t163 = t315 * t168 + t312 * t309;
t308 = t315 ^ 2;
t339 = t307 + t308;
t338 = qJD(5) * t305;
t337 = t358 * t309;
t336 = t312 * t338;
t335 = t315 * t338;
t162 = t312 * t168 - t315 * t309;
t260 = -t319 * pkin(1) + t325;
t261 = -t324 + t356;
t219 = t317 * t260 - t314 * t261;
t249 = -t314 * t290 - t317 * t291;
t332 = t313 * t289;
t331 = t316 * t289;
t283 = t314 * qJDD(1) + t317 * t319;
t264 = -pkin(5) * t283 + t317 * g(3);
t284 = t317 * qJDD(1) - t314 * t319;
t263 = pkin(5) * t284 + t314 * g(3);
t143 = t315 * t162 - t312 * t163;
t144 = t312 * t162 + t315 * t163;
t218 = t314 * t260 + t317 * t261;
t248 = t317 * t290 - t314 * t291;
t318 = qJD(5) ^ 2;
t294 = t308 * t303;
t288 = -t294 - t318;
t287 = t294 - t318;
t286 = -t318 - t355;
t285 = t318 - t355;
t275 = t294 - t355;
t274 = t294 + t355;
t269 = t339 * t304;
t268 = -t293 - 0.2e1 * t336;
t267 = -t293 - t336;
t266 = t335 - t349;
t265 = 0.2e1 * t335 - t349;
t262 = t339 * t338;
t246 = t313 * qJDD(5) + t316 * t262;
t245 = -t316 * qJDD(5) + t313 * t262;
t244 = t315 * t266 - t307 * t338;
t243 = -t312 * t267 - t308 * t338;
t242 = -t312 * t286 - t342;
t241 = -t312 * t285 + t343;
t240 = t315 * t288 - t351;
t239 = t315 * t287 - t350;
t238 = t315 * t286 - t350;
t237 = t312 * t288 + t343;
t232 = -t316 * t269 - t313 * t274;
t231 = -t313 * t269 + t316 * t274;
t224 = -t312 * t265 + t315 * t268;
t223 = t316 * t241 - t312 * t348;
t222 = t316 * t239 - t313 * t293;
t221 = t313 * t241 + t312 * t341;
t220 = t313 * t239 + t315 * t341;
t217 = t316 * t244 - t332;
t216 = t316 * t243 + t332;
t215 = t313 * t244 + t331;
t214 = t313 * t243 - t331;
t211 = t316 * t242 + t313 * t265;
t210 = t316 * t240 - t313 * t268;
t209 = t313 * t242 - t316 * t265;
t208 = t313 * t240 + t316 * t268;
t203 = -t310 * t245 + t311 * t246;
t202 = -t311 * t245 - t310 * t246;
t201 = t316 * t224 - t313 * t275;
t200 = t313 * t224 + t316 * t275;
t199 = -t310 * t231 + t311 * t232;
t198 = t311 * t231 + t310 * t232;
t193 = -t310 * t221 + t311 * t223;
t192 = -t310 * t220 + t311 * t222;
t191 = -t311 * t221 - t310 * t223;
t190 = -t311 * t220 - t310 * t222;
t185 = -t310 * t215 + t311 * t217;
t184 = -t310 * t214 + t311 * t216;
t183 = -t311 * t215 - t310 * t217;
t182 = -t311 * t214 - t310 * t216;
t177 = -t310 * t209 + t311 * t211;
t176 = -t310 * t208 + t311 * t210;
t175 = t311 * t209 + t310 * t211;
t174 = t311 * t208 + t310 * t210;
t171 = qJ(3) * t180 + t357;
t169 = -qJ(3) * t181 + t337;
t165 = -t310 * t200 + t311 * t201;
t164 = -t311 * t200 - t310 * t201;
t161 = -pkin(7) * t238 + t344;
t160 = -pkin(7) * t237 + t352;
t159 = -pkin(4) * t238 + t163;
t158 = -pkin(4) * t237 + t162;
t157 = t314 * t198 + t317 * t199;
t156 = -t317 * t198 + t314 * t199;
t153 = t314 * t175 + t317 * t177;
t152 = t314 * t174 + t317 * t176;
t151 = -t317 * t175 + t314 * t177;
t150 = -t317 * t174 + t314 * t176;
t145 = -pkin(3) * t309 + pkin(6) * t334;
t141 = -pkin(6) * t231 + t316 * t143;
t140 = pkin(6) * t232 + t313 * t143;
t139 = -pkin(6) * t209 - t313 * t159 + t316 * t161;
t138 = -pkin(6) * t208 - t313 * t158 + t316 * t160;
t137 = t316 * t144 + t313 * t167;
t136 = t313 * t144 - t316 * t167;
t135 = -pkin(3) * t238 + pkin(6) * t211 + t316 * t159 + t313 * t161;
t134 = -pkin(3) * t237 + pkin(6) * t210 + t316 * t158 + t313 * t160;
t129 = -qJ(3) * t198 - t310 * t140 + t311 * t141;
t128 = -qJ(3) * t199 - t311 * t140 - t310 * t141;
t127 = -t310 * t136 + t311 * t137;
t126 = t311 * t136 + t310 * t137;
t125 = qJ(2) * t238 - qJ(3) * t175 - t310 * t135 + t311 * t139;
t124 = qJ(2) * t237 - qJ(3) * t174 - t310 * t134 + t311 * t138;
t121 = pkin(6) * t353 + qJ(3) * t132 - t310 * t145 + t357;
t120 = -pkin(6) * t136 - (pkin(4) * t313 - pkin(7) * t316) * t143;
t119 = -pkin(6) * t354 - qJ(3) * t133 - t311 * t145 + t337;
t118 = -qJ(3) * t177 - t311 * t135 - t310 * t139 + t358 * t238;
t117 = -qJ(3) * t176 - t311 * t134 - t310 * t138 + t358 * t237;
t116 = pkin(6) * t137 - (-pkin(4) * t316 - pkin(7) * t313 - pkin(3)) * t143;
t115 = t314 * t126 + t317 * t127;
t114 = -t317 * t126 + t314 * t127;
t113 = -qJ(2) * t143 - qJ(3) * t126 - t310 * t116 + t311 * t120;
t112 = -qJ(3) * t127 - t311 * t116 - t310 * t120 - t143 * t358;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t283, -t284, 0, t249, 0, 0, 0, 0, 0, 0, -t283, 0, t284, t219, 0, 0, 0, 0, 0, 0, -t233, t333, 0, -t364, 0, 0, 0, 0, 0, 0, -t194, t374, 0, -t375, 0, 0, 0, 0, 0, 0, t152, t153, t157, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t284, -t283, 0, t248, 0, 0, 0, 0, 0, 0, t284, 0, t283, t218, 0, 0, 0, 0, 0, 0, t333, t233, 0, t154, 0, 0, 0, 0, 0, 0, t374, t194, 0, t122, 0, 0, 0, 0, 0, 0, t150, t151, t156, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, 0, 0, 0, 0, 0, 0, -t237, -t238, 0, t143; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t284, 0, -t283, 0, -t263, -t264, -t248, -pkin(5) * t248, 0, t284, 0, 0, t283, 0, -t263, -t218, t264, -pkin(5) * t218 + (-pkin(1) * t314 + qJ(2) * t317) * g(3), 0, 0, -t333, 0, -t233, 0, t365, t369, t154, -pkin(5) * t154 - t314 * t169 + t317 * t171, 0, 0, -t374, 0, -t194, 0, t379, t378, t122, -pkin(5) * t122 - t314 * t119 + t317 * t121, -t314 * t183 + t317 * t185, -t314 * t164 + t317 * t165, -t314 * t191 + t317 * t193, -t314 * t182 + t317 * t184, -t314 * t190 + t317 * t192, -t314 * t202 + t317 * t203, -pkin(5) * t150 - t314 * t117 + t317 * t124, -pkin(5) * t151 - t314 * t118 + t317 * t125, -pkin(5) * t156 - t314 * t128 + t317 * t129, -pkin(5) * t114 - t314 * t112 + t317 * t113; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t283, 0, t284, 0, t264, -t263, t249, pkin(5) * t249, 0, t283, 0, 0, -t284, 0, t264, t219, t263, pkin(5) * t219 + (pkin(1) * t317 + qJ(2) * t314) * g(3), 0, 0, -t233, 0, t333, 0, t369, -t365, t364, -pkin(5) * t364 + t317 * t169 + t314 * t171, 0, 0, -t194, 0, t374, 0, t378, -t379, t375, -pkin(5) * t375 + t317 * t119 + t314 * t121, t317 * t183 + t314 * t185, t317 * t164 + t314 * t165, t317 * t191 + t314 * t193, t317 * t182 + t314 * t184, t317 * t190 + t314 * t192, t317 * t202 + t314 * t203, pkin(5) * t152 + t317 * t117 + t314 * t124, pkin(5) * t153 + t317 * t118 + t314 * t125, pkin(5) * t157 + t317 * t128 + t314 * t129, pkin(5) * t115 + t317 * t112 + t314 * t113; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t290, t291, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t330 + 0.2e1 * t356, 0, 0.2e1 * t306 + t326, pkin(1) * t261 + qJ(2) * t260, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t281 + t358 * t282 + t212, qJ(2) * t282 + t358 * t281 + t213, 0, qJ(2) * t181 + t180 * t358, 0, 0, 0, 0, 0, t304, pkin(3) * t273 + qJ(2) * t363 + t229 * t358 + t172, pkin(3) * t271 + qJ(2) * t229 - t358 * t363 + t173, 0, pkin(3) * t148 + qJ(2) * t133 + t132 * t358, (-t266 - t335) * t312, -t315 * t265 - t312 * t268, -t315 * t285 - t351, (-t267 + t336) * t315, -t312 * t287 - t342, 0, -pkin(3) * t208 - pkin(4) * t268 - pkin(7) * t240 + qJ(2) * t176 - t358 * t174 + t344, -pkin(3) * t209 + pkin(4) * t265 - pkin(7) * t242 + qJ(2) * t177 - t358 * t175 - t352, -pkin(3) * t231 - pkin(4) * t274 + pkin(7) * t269 + qJ(2) * t199 - t358 * t198 - t144, -pkin(3) * t136 + pkin(4) * t167 - pkin(7) * t144 + qJ(2) * t127 - t358 * t126;];
tauB_reg = t1;
