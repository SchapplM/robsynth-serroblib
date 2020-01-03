% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRR11_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:56
% EndTime: 2019-12-31 18:06:01
% DurationCPUTime: 2.27s
% Computational Cost: add. (6139->319), mult. (11645->432), div. (0->0), fcn. (6358->6), ass. (0->215)
t308 = sin(qJ(5));
t311 = cos(qJ(5));
t312 = cos(qJ(4));
t348 = qJD(1) * t312;
t269 = -t311 * qJD(4) + t308 * t348;
t271 = qJD(4) * t308 + t311 * t348;
t238 = t271 * t269;
t346 = qJD(1) * qJD(4);
t334 = t312 * t346;
t309 = sin(qJ(4));
t344 = qJDD(1) * t309;
t274 = -t334 - t344;
t317 = qJDD(5) - t274;
t374 = -t238 + t317;
t376 = t308 * t374;
t375 = t311 * t374;
t292 = qJD(1) * t309 + qJD(5);
t255 = t292 * t269;
t335 = t309 * t346;
t342 = qJDD(1) * t312;
t275 = -t335 + t342;
t336 = t269 * qJD(5) - t308 * qJDD(4) - t311 * t275;
t204 = t336 + t255;
t331 = -t311 * qJDD(4) + t308 * t275;
t200 = (qJD(5) - t292) * t271 + t331;
t315 = qJD(1) ^ 2;
t310 = sin(qJ(1));
t313 = cos(qJ(1));
t285 = g(1) * t313 + g(2) * t310;
t304 = qJDD(1) * qJ(2);
t323 = t285 - t304;
t366 = pkin(1) + qJ(3);
t373 = t366 * t315 - qJDD(3) + t323;
t265 = t269 ^ 2;
t266 = t271 ^ 2;
t291 = t292 ^ 2;
t372 = 2 * qJD(3);
t371 = pkin(2) + pkin(3);
t370 = pkin(4) * t309;
t369 = pkin(4) * t312;
t341 = qJDD(1) * t313;
t278 = -t310 * t315 + t341;
t368 = pkin(5) * t278;
t343 = qJDD(1) * t310;
t279 = t313 * t315 + t343;
t367 = pkin(5) * t279;
t365 = qJ(2) - pkin(6);
t364 = qJDD(1) * pkin(1);
t345 = qJD(2) * qJD(1);
t302 = 0.2e1 * t345;
t236 = -qJDD(1) * pkin(6) + t302 - t373;
t223 = g(3) * t309 + t236 * t312;
t314 = qJD(4) ^ 2;
t327 = -pkin(7) * t312 + t370;
t322 = t315 * t327;
t198 = qJDD(4) * pkin(4) + pkin(7) * t314 - t312 * t322 + t223;
t363 = t198 * t308;
t362 = t198 * t311;
t218 = t238 + t317;
t361 = t218 * t308;
t360 = t218 * t311;
t337 = t309 * t312 * t315;
t282 = qJDD(4) + t337;
t359 = t282 * t309;
t358 = t282 * t312;
t283 = qJDD(4) - t337;
t357 = t283 * t309;
t356 = t283 * t312;
t355 = t292 * t308;
t354 = t292 * t311;
t305 = t309 ^ 2;
t353 = t305 * t315;
t306 = t312 ^ 2;
t352 = t306 * t315;
t284 = t310 * g(1) - t313 * g(2);
t324 = -qJDD(2) + t284;
t318 = qJ(2) * t315 + t324;
t332 = t366 * qJDD(1);
t316 = t332 + t318;
t340 = qJD(1) * t372;
t248 = t316 + t340;
t307 = t315 * pkin(6);
t235 = -t307 + t248;
t351 = t309 * t235;
t350 = t312 * t235;
t328 = pkin(7) * t309 + t369;
t191 = -pkin(4) * t274 - pkin(7) * t275 - t307 + (t328 * qJD(4) + t372) * qJD(1) + t316;
t333 = t312 * g(3) - t309 * t236;
t199 = -t314 * pkin(4) + qJDD(4) * pkin(7) - t309 * t322 - t333;
t165 = t308 * t191 + t311 * t199;
t349 = t305 + t306;
t339 = t309 * t238;
t338 = t312 * t238;
t164 = -t311 * t191 + t199 * t308;
t144 = t164 * t308 + t311 * t165;
t182 = -t309 * t223 - t312 * t333;
t249 = -0.2e1 * t345 + t373;
t207 = -t248 * t310 - t313 * t249;
t256 = -pkin(1) * t315 + t302 - t323;
t257 = t318 + t364;
t213 = t313 * t256 - t257 * t310;
t240 = -t284 * t310 - t313 * t285;
t330 = t310 * t337;
t329 = t313 * t337;
t326 = g(3) * t310 + t368;
t325 = g(3) * t313 - t367;
t143 = -t164 * t311 + t165 * t308;
t181 = t223 * t312 - t309 * t333;
t206 = t248 * t313 - t249 * t310;
t210 = t256 * t310 + t257 * t313;
t239 = t284 * t313 - t285 * t310;
t321 = -t285 + 0.2e1 * t304 + t302;
t299 = -pkin(2) * t315 + g(3);
t320 = -pkin(2) * t343 + t313 * t299 - t367;
t290 = -t314 - t352;
t289 = t314 - t352;
t288 = -t314 - t353;
t287 = -t314 + t353;
t281 = (-t305 + t306) * t315;
t280 = t349 * t315;
t277 = t349 * qJDD(1);
t276 = -0.2e1 * t335 + t342;
t273 = 0.2e1 * t334 + t344;
t268 = t349 * t346;
t254 = -t266 + t291;
t253 = t265 - t291;
t252 = t275 * t309 + t306 * t346;
t251 = t274 * t312 + t305 * t346;
t247 = -t290 * t309 - t358;
t246 = t288 * t312 - t357;
t245 = t290 * t312 - t359;
t244 = t289 * t312 + t357;
t243 = t288 * t309 + t356;
t242 = t287 * t309 + t358;
t241 = t371 * t277;
t237 = -pkin(2) * t341 - t299 * t310 - t368;
t233 = -t266 + t265;
t232 = -t277 * t313 + t280 * t310;
t231 = -t277 * t310 - t280 * t313;
t230 = -t266 - t291;
t229 = -pkin(2) * t248 + qJ(2) * g(3);
t227 = -t273 * t309 + t276 * t312;
t226 = -pkin(2) * t249 + t366 * g(3);
t225 = -t291 - t265;
t221 = -qJD(5) * t271 - t331;
t216 = t265 + t266;
t215 = t245 * t313 - t276 * t310;
t214 = t243 * t313 - t273 * t310;
t212 = t245 * t310 + t276 * t313;
t211 = t243 * t310 + t273 * t313;
t209 = (-t269 * t311 + t271 * t308) * t292;
t208 = (-t269 * t308 - t271 * t311) * t292;
t205 = -t255 + t336;
t201 = (-qJD(5) - t292) * t271 - t331;
t196 = -t271 * t355 - t311 * t336;
t195 = t271 * t354 - t308 * t336;
t194 = -t221 * t308 + t269 * t354;
t193 = t221 * t311 + t269 * t355;
t192 = t209 * t309 - t312 * t317;
t190 = t253 * t311 - t361;
t189 = -t254 * t308 + t375;
t188 = t253 * t308 + t360;
t187 = t254 * t311 + t376;
t184 = -t230 * t308 - t360;
t183 = t230 * t311 - t361;
t180 = t225 * t311 - t376;
t179 = t225 * t308 + t375;
t178 = t196 * t309 - t338;
t177 = t194 * t309 + t338;
t176 = t371 * t280 + t182;
t175 = -t365 * t247 - t371 * t276 + t351;
t174 = -t365 * t246 - t371 * t273 - t350;
t173 = t181 * t313 - t235 * t310;
t172 = t181 * t310 + t235 * t313;
t171 = -t200 * t311 - t205 * t308;
t170 = t201 * t311 + t204 * t308;
t169 = -t200 * t308 + t205 * t311;
t168 = t201 * t308 - t204 * t311;
t167 = t190 * t309 + t200 * t312;
t166 = t189 * t309 + t205 * t312;
t162 = t184 * t312 - t204 * t309;
t161 = t184 * t309 + t204 * t312;
t160 = -pkin(7) * t183 - t362;
t159 = t371 * t245 - t366 * t247 + t333;
t158 = t371 * t243 - t366 * t246 + t223;
t157 = t180 * t312 - t201 * t309;
t156 = t180 * t309 + t201 * t312;
t155 = -pkin(7) * t179 - t363;
t154 = t170 * t309 + t233 * t312;
t153 = t171 * t312 - t216 * t309;
t152 = t171 * t309 + t216 * t312;
t151 = -pkin(4) * t183 + t165;
t150 = -pkin(4) * t179 + t164;
t149 = -t365 * t182 - t371 * t235;
t148 = t161 * t313 - t183 * t310;
t147 = t161 * t310 + t183 * t313;
t146 = t156 * t313 - t179 * t310;
t145 = t156 * t310 + t179 * t313;
t142 = t152 * t313 - t169 * t310;
t141 = t152 * t310 + t169 * t313;
t140 = t371 * t181 - t366 * t182;
t139 = t144 * t312 - t198 * t309;
t138 = t144 * t309 + t198 * t312;
t137 = -pkin(7) * t169 - t143;
t136 = t138 * t313 - t143 * t310;
t135 = t138 * t310 + t143 * t313;
t134 = pkin(4) * t204 + pkin(7) * t184 + t371 * t161 - t366 * t162 - t363;
t133 = pkin(4) * t201 + pkin(7) * t180 + t371 * t156 - t366 * t157 + t362;
t132 = t151 * t312 + t160 * t309 - t365 * t162 - t371 * t183;
t131 = t150 * t312 + t155 * t309 - t365 * t157 - t371 * t179;
t130 = t137 * t309 - t365 * t153 + (-t369 - t371) * t169;
t129 = pkin(4) * t216 + pkin(7) * t171 + t371 * t152 - t366 * t153 + t144;
t128 = pkin(4) * t198 + pkin(7) * t144 + t371 * t138 - t366 * t139;
t127 = -t365 * t139 + (-t328 - t371) * t143;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t279, -t278, 0, t240, 0, 0, 0, 0, 0, 0, 0, t279, t278, t213, 0, 0, 0, 0, 0, 0, 0, t278, -t279, t207, 0, 0, 0, 0, 0, 0, t214, t215, t232, t173, 0, 0, 0, 0, 0, 0, t146, t148, t142, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t278, -t279, 0, t239, 0, 0, 0, 0, 0, 0, 0, -t278, t279, t210, 0, 0, 0, 0, 0, 0, 0, t279, t278, t206, 0, 0, 0, 0, 0, 0, t211, t212, t231, t172, 0, 0, 0, 0, 0, 0, t145, t147, t141, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t246, t247, 0, t182, 0, 0, 0, 0, 0, 0, t157, t162, t153, t139; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t278, 0, -t279, 0, -t326, -t325, -t239, -pkin(5) * t239, 0, -t278, t279, 0, 0, 0, -t210, t326, t325, -pkin(5) * t210 + (-pkin(1) * t310 + qJ(2) * t313) * g(3), 0, t279, t278, 0, 0, 0, -t206, t320, t237, -pkin(5) * t206 - t226 * t310 + t229 * t313, t252 * t313 - t330, t227 * t313 - t281 * t310, t244 * t313 - t310 * t342, t251 * t313 + t330, t242 * t313 + t309 * t343, -qJDD(4) * t310 - t268 * t313, -pkin(5) * t211 - t158 * t310 + t174 * t313, -pkin(5) * t212 - t159 * t310 + t175 * t313, -pkin(5) * t231 + t176 * t313 + t241 * t310, -pkin(5) * t172 - t140 * t310 + t149 * t313, t178 * t313 - t195 * t310, t154 * t313 - t168 * t310, t166 * t313 - t187 * t310, t177 * t313 - t193 * t310, t167 * t313 - t188 * t310, t192 * t313 - t208 * t310, -pkin(5) * t145 + t131 * t313 - t133 * t310, -pkin(5) * t147 + t132 * t313 - t134 * t310, -pkin(5) * t141 - t129 * t310 + t130 * t313, -pkin(5) * t135 + t127 * t313 - t128 * t310; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t279, 0, t278, 0, t325, -t326, t240, pkin(5) * t240, 0, -t279, -t278, 0, 0, 0, t213, -t325, t326, pkin(5) * t213 + (pkin(1) * t313 + qJ(2) * t310) * g(3), 0, -t278, t279, 0, 0, 0, t207, -t237, t320, pkin(5) * t207 + t226 * t313 + t229 * t310, t252 * t310 + t329, t227 * t310 + t281 * t313, t244 * t310 + t312 * t341, t251 * t310 - t329, t242 * t310 - t309 * t341, qJDD(4) * t313 - t268 * t310, pkin(5) * t214 + t158 * t313 + t174 * t310, pkin(5) * t215 + t159 * t313 + t175 * t310, pkin(5) * t232 + t176 * t310 - t241 * t313, pkin(5) * t173 + t140 * t313 + t149 * t310, t178 * t310 + t195 * t313, t154 * t310 + t168 * t313, t166 * t310 + t187 * t313, t177 * t310 + t193 * t313, t167 * t310 + t188 * t313, t192 * t310 + t208 * t313, pkin(5) * t146 + t131 * t310 + t133 * t313, pkin(5) * t148 + t132 * t310 + t134 * t313, pkin(5) * t142 + t129 * t313 + t130 * t310, pkin(5) * t136 + t127 * t310 + t128 * t313; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t284, t285, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t324 - 0.2e1 * t364, t321, pkin(1) * t257 + qJ(2) * t256, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t321, t324 + 0.2e1 * t332 + t340, -qJ(2) * t249 + t366 * t248, (t275 - t335) * t312, -t273 * t312 - t276 * t309, -t289 * t309 + t356, (-t274 + t334) * t309, t287 * t312 - t359, 0, t365 * t243 + t366 * t273 + t351, t365 * t245 + t366 * t276 + t350, -t365 * t277 - t366 * t280 - t181, t365 * t181 + t366 * t235, t196 * t312 + t339, t170 * t312 - t233 * t309, t189 * t312 - t205 * t309, t194 * t312 - t339, t190 * t312 - t200 * t309, t209 * t312 + t309 * t317, -t150 * t309 + t155 * t312 + t365 * t156 + t366 * t179, -t151 * t309 + t160 * t312 + t365 * t161 + t366 * t183, t137 * t312 + t365 * t152 + (t366 + t370) * t169, t365 * t138 + (t327 + t366) * t143;];
tauB_reg = t1;
