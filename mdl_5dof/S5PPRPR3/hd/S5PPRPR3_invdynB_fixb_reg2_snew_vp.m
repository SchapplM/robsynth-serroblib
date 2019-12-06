% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:24
% EndTime: 2019-12-05 15:05:32
% DurationCPUTime: 4.84s
% Computational Cost: add. (12284->367), mult. (20086->576), div. (0->0), fcn. (15353->10), ass. (0->251)
t321 = sin(qJ(3));
t315 = sin(pkin(7));
t318 = cos(pkin(7));
t293 = g(1) * t318 + g(2) * t315;
t314 = sin(pkin(8));
t317 = cos(pkin(8));
t347 = g(3) - qJDD(1);
t275 = -t293 * t317 - t314 * t347;
t292 = g(1) * t315 - g(2) * t318;
t284 = -qJDD(2) + t292;
t323 = cos(qJ(3));
t226 = t275 * t323 - t284 * t321;
t325 = qJD(3) ^ 2;
t220 = -pkin(3) * t325 + t226;
t313 = sin(pkin(9));
t316 = cos(pkin(9));
t225 = t275 * t321 + t284 * t323;
t326 = qJDD(3) * pkin(3) - t225;
t180 = t220 * t313 - t316 * t326;
t181 = t220 * t316 + t313 * t326;
t333 = t180 * t313 + t181 * t316;
t141 = t180 * t316 - t181 * t313;
t349 = t323 * t141;
t112 = -t321 * t333 + t349;
t354 = t321 * t141;
t113 = t323 * t333 + t354;
t285 = qJDD(3) * t313 + t316 * t325;
t286 = qJDD(3) * t316 - t313 * t325;
t248 = t285 * t321 - t286 * t323;
t379 = t314 * t248;
t274 = -t293 * t314 + t317 * t347;
t269 = qJDD(4) + t274;
t223 = qJ(4) * t285 + t269 * t316;
t330 = -qJ(4) * t286 + t269 * t313;
t378 = t321 * t223 + t323 * t330;
t377 = t223 * t323 - t321 * t330;
t329 = -t285 * t323 - t286 * t321;
t361 = t317 * t318;
t196 = t248 * t361 + t315 * t329;
t363 = t315 * t317;
t194 = t248 * t363 - t318 * t329;
t375 = t315 * t347;
t373 = t318 * t347;
t256 = t317 * t274;
t209 = -t275 * t314 + t256;
t320 = sin(qJ(5));
t310 = t320 ^ 2;
t368 = t310 * t325;
t367 = t314 * t274;
t289 = qJDD(3) * t323 - t321 * t325;
t365 = t314 * t289;
t364 = t315 * t284;
t362 = t317 * t284;
t360 = t318 * t284;
t288 = qJDD(3) * t321 + t323 * t325;
t359 = t318 * t288;
t358 = t318 * t289;
t173 = -qJDD(3) * pkin(4) - pkin(6) * t325 + t180;
t357 = t320 * t173;
t322 = cos(qJ(5));
t302 = t320 * t325 * t322;
t294 = qJDD(5) + t302;
t356 = t320 * t294;
t295 = qJDD(5) - t302;
t355 = t320 * t295;
t353 = t321 * t274;
t352 = t322 * t173;
t351 = t322 * t294;
t350 = t322 * t295;
t348 = t323 * t274;
t174 = -pkin(4) * t325 + qJDD(3) * pkin(6) + t181;
t155 = t174 * t322 + t269 * t320;
t311 = t322 ^ 2;
t346 = t310 + t311;
t345 = qJD(3) * qJD(5);
t344 = qJDD(3) * t314;
t343 = t317 * qJDD(3);
t342 = t320 * qJDD(3);
t306 = t322 * qJDD(3);
t341 = pkin(1) * t314 + pkin(5);
t340 = t320 * t345;
t339 = t322 * t345;
t154 = t174 * t320 - t269 * t322;
t130 = t154 * t320 + t155 * t322;
t287 = t346 * qJDD(3);
t308 = t311 * t325;
t290 = t308 + t368;
t250 = t287 * t313 + t290 * t316;
t251 = t287 * t316 - t290 * t313;
t200 = t250 * t323 + t251 * t321;
t201 = -t250 * t321 + t251 * t323;
t338 = -pkin(2) * t200 - pkin(3) * t250 - pkin(4) * t290 - pkin(6) * t287 + qJ(2) * t201 - t130;
t337 = -pkin(2) * t329 + pkin(3) * t285 + qJ(2) * t248 + t181;
t336 = pkin(2) * t248 - pkin(3) * t286 + qJ(2) * t329 + t180;
t335 = pkin(2) * t289 + qJ(2) * t288 - t225;
t334 = -pkin(2) * t288 + qJ(2) * t289 - t226;
t210 = t275 * t317 + t367;
t253 = -t292 * t315 - t293 * t318;
t332 = t313 * t302;
t331 = t316 * t302;
t129 = t154 * t322 - t155 * t320;
t184 = t225 * t323 - t226 * t321;
t185 = t225 * t321 + t226 * t323;
t252 = t292 * t318 - t293 * t315;
t328 = t288 * t315 + t317 * t358;
t327 = t289 * t363 - t359;
t324 = qJD(5) ^ 2;
t301 = -t308 - t324;
t300 = t308 - t324;
t299 = -t324 - t368;
t298 = t324 - t368;
t297 = t318 * t344;
t296 = t315 * t344;
t291 = t308 - t368;
t283 = t306 - 0.2e1 * t340;
t282 = t306 - t340;
t281 = t339 + t342;
t280 = 0.2e1 * t339 + t342;
t279 = t346 * t345;
t277 = t314 * t288;
t273 = t281 * t322 - t310 * t345;
t272 = -t282 * t320 - t311 * t345;
t271 = qJDD(5) * t313 + t279 * t316;
t270 = -qJDD(5) * t316 + t279 * t313;
t267 = -t299 * t320 - t350;
t266 = -t298 * t320 + t351;
t265 = t301 * t322 - t356;
t264 = t300 * t322 - t355;
t263 = t299 * t322 - t355;
t262 = -t298 * t322 - t356;
t261 = t301 * t320 + t351;
t260 = -t300 * t320 - t350;
t259 = (-t281 - t339) * t320;
t258 = (-t282 + t340) * t322;
t243 = -t280 * t320 + t283 * t322;
t242 = -t280 * t322 - t283 * t320;
t241 = t289 * t315 - t317 * t359;
t239 = -t288 * t363 - t358;
t237 = t314 * t329;
t236 = t273 * t316 - t332;
t235 = t272 * t316 + t332;
t234 = t273 * t313 + t331;
t233 = t272 * t313 - t331;
t232 = t266 * t316 + t313 * t342;
t231 = t264 * t316 + t306 * t313;
t230 = t266 * t313 - t316 * t342;
t229 = t264 * t313 - t306 * t316;
t228 = pkin(5) * t288 + t348;
t227 = -pkin(5) * t289 + t353;
t219 = t267 * t316 + t280 * t313;
t218 = t265 * t316 - t283 * t313;
t217 = t267 * t313 - t280 * t316;
t216 = t265 * t313 + t283 * t316;
t212 = t243 * t316 - t291 * t313;
t211 = t243 * t313 + t291 * t316;
t205 = t288 * t341 + t348;
t204 = t289 * t341 - t353;
t203 = -t270 * t321 + t271 * t323;
t202 = -t270 * t323 - t271 * t321;
t199 = t210 * t318 - t364;
t198 = t210 * t315 + t360;
t197 = -t248 * t315 + t329 * t361;
t195 = t248 * t318 + t329 * t363;
t193 = -t234 * t321 + t236 * t323;
t192 = -t233 * t321 + t235 * t323;
t191 = -t234 * t323 - t236 * t321;
t190 = -t233 * t323 - t235 * t321;
t189 = -t230 * t321 + t232 * t323;
t188 = -t229 * t321 + t231 * t323;
t187 = -t230 * t323 - t232 * t321;
t186 = -t229 * t323 - t231 * t321;
t179 = -t217 * t321 + t219 * t323;
t178 = -t216 * t321 + t218 * t323;
t177 = t217 * t323 + t219 * t321;
t176 = t216 * t323 + t218 * t321;
t171 = -t211 * t321 + t212 * t323;
t170 = -t211 * t323 - t212 * t321;
t169 = t317 * t228 + t314 * t334;
t168 = t317 * t227 + t314 * t335;
t167 = t193 * t317 - t259 * t314;
t166 = t192 * t317 - t258 * t314;
t165 = t189 * t317 - t262 * t314;
t164 = t188 * t317 - t260 * t314;
t163 = t185 * t317 + t367;
t162 = t185 * t314 - t256;
t161 = -pkin(5) * t329 + t377;
t160 = pkin(5) * t248 + t378;
t159 = t179 * t317 + t263 * t314;
t158 = t178 * t317 + t261 * t314;
t157 = t179 * t314 - t263 * t317;
t156 = t178 * t314 - t261 * t317;
t153 = -pkin(6) * t263 + t352;
t152 = -pkin(6) * t261 + t357;
t149 = t200 * t315 + t201 * t361;
t148 = -t200 * t318 + t201 * t363;
t147 = t171 * t317 - t242 * t314;
t146 = -t248 * t341 - t378;
t145 = -t329 * t341 + t377;
t144 = -pkin(4) * t263 + t155;
t143 = -pkin(4) * t261 + t154;
t138 = t163 * t318 - t184 * t315;
t137 = t163 * t315 + t184 * t318;
t136 = -pkin(3) * t269 + qJ(4) * t333;
t135 = -pkin(1) * t162 + pkin(2) * t274 - pkin(5) * t185;
t134 = t159 * t318 + t177 * t315;
t133 = t158 * t318 + t176 * t315;
t132 = t159 * t315 - t177 * t318;
t131 = t158 * t315 - t176 * t318;
t127 = -pkin(2) * t177 - pkin(3) * t217 + pkin(4) * t280 - pkin(6) * t267 - t357;
t126 = -pkin(2) * t176 - pkin(3) * t216 - pkin(4) * t283 - pkin(6) * t265 + t352;
t125 = t317 * t161 - t314 * t337;
t124 = t317 * t160 - t314 * t336;
t123 = -qJ(2) * t162 - (pkin(2) * t314 - pkin(5) * t317) * t184;
t122 = -qJ(4) * t250 + t129 * t316;
t121 = qJ(4) * t251 + t129 * t313;
t120 = -qJ(4) * t217 - t144 * t313 + t153 * t316;
t119 = -qJ(4) * t216 - t143 * t313 + t152 * t316;
t118 = -pkin(3) * t263 + qJ(4) * t219 + t144 * t316 + t153 * t313;
t117 = -pkin(3) * t261 + qJ(4) * t218 + t143 * t316 + t152 * t313;
t116 = t130 * t316 + t173 * t313;
t115 = t130 * t313 - t173 * t316;
t109 = t113 * t317 + t269 * t314;
t108 = t113 * t314 - t269 * t317;
t107 = pkin(2) * t112 + pkin(3) * t141;
t106 = -pkin(5) * t200 - t121 * t321 + t122 * t323;
t105 = -t323 * t121 - t321 * t122 - t201 * t341;
t104 = -t115 * t321 + t116 * t323;
t103 = t115 * t323 + t116 * t321;
t102 = -pkin(5) * t177 - t118 * t321 + t120 * t323;
t101 = -pkin(5) * t176 - t117 * t321 + t119 * t323;
t100 = pkin(5) * t112 + qJ(4) * t349 - t136 * t321;
t99 = t109 * t318 - t112 * t315;
t98 = t109 * t315 + t112 * t318;
t97 = -qJ(4) * t115 - (pkin(4) * t313 - pkin(6) * t316) * t129;
t96 = -pkin(1) * t157 + pkin(2) * t263 - pkin(5) * t179 - t118 * t323 - t120 * t321;
t95 = -pkin(1) * t156 + pkin(2) * t261 - pkin(5) * t178 - t117 * t323 - t119 * t321;
t94 = t104 * t317 - t129 * t314;
t93 = t104 * t314 + t129 * t317;
t92 = t317 * t106 - t314 * t338;
t91 = qJ(4) * t116 - (-pkin(4) * t316 - pkin(6) * t313 - pkin(3)) * t129;
t90 = -qJ(2) * t157 + t102 * t317 - t127 * t314;
t89 = -qJ(2) * t156 + t101 * t317 - t126 * t314;
t88 = -pkin(1) * t108 + pkin(2) * t269 - pkin(5) * t113 - qJ(4) * t354 - t136 * t323;
t87 = -pkin(2) * t103 - pkin(3) * t115 + pkin(4) * t173 - pkin(6) * t130;
t86 = t103 * t315 + t318 * t94;
t85 = -t103 * t318 + t315 * t94;
t84 = -qJ(2) * t108 + t100 * t317 - t107 * t314;
t83 = -pkin(5) * t103 - t321 * t91 + t323 * t97;
t82 = -pkin(1) * t93 - pkin(2) * t129 - pkin(5) * t104 - t321 * t97 - t323 * t91;
t81 = -qJ(2) * t93 - t314 * t87 + t317 * t83;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, 0, 0, 0, 0, 0, 0, t241, -t328, 0, t138, 0, 0, 0, 0, 0, 0, t197, t196, 0, t99, 0, 0, 0, 0, 0, 0, t133, t134, t149, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, 0, 0, 0, 0, 0, 0, t239, -t327, 0, t137, 0, 0, 0, 0, 0, 0, t195, t194, 0, t98, 0, 0, 0, 0, 0, 0, t131, t132, t148, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, 0, 0, 0, 0, 0, 0, -t277, -t365, 0, t162, 0, 0, 0, 0, 0, 0, t237, t379, 0, t108, 0, 0, 0, 0, 0, 0, t156, t157, t314 * t201, t93; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t375, -t373, -t252, -qJ(1) * t252, 0, 0, 0, 0, 0, 0, -t274 * t315 - t314 * t360, -t275 * t315 - t317 * t360, t318 * t209, -qJ(1) * t198 - (pkin(1) * t315 - qJ(2) * t318) * t209, 0, 0, t328, 0, t241, t297, -qJ(1) * t239 + t168 * t318 - t205 * t315, qJ(1) * t327 + t169 * t318 - t204 * t315, t184 * t361 + t185 * t315, -qJ(1) * t137 + t123 * t318 - t135 * t315, 0, 0, -t196, 0, t197, t297, -qJ(1) * t195 + t124 * t318 - t145 * t315, -qJ(1) * t194 + t125 * t318 - t146 * t315, t112 * t361 + t113 * t315, -qJ(1) * t98 - t315 * t88 + t318 * t84, t167 * t318 - t191 * t315, t147 * t318 - t170 * t315, t165 * t318 - t187 * t315, t166 * t318 - t190 * t315, t164 * t318 - t186 * t315, -t202 * t315 + t203 * t361, -qJ(1) * t131 - t315 * t95 + t318 * t89, -qJ(1) * t132 - t315 * t96 + t318 * t90, -qJ(1) * t148 - t105 * t315 + t318 * t92, -qJ(1) * t85 - t315 * t82 + t318 * t81; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t373, -t375, t253, qJ(1) * t253, 0, 0, 0, 0, 0, 0, t274 * t318 - t314 * t364, t275 * t318 - t315 * t362, t315 * t209, qJ(1) * t199 - (-pkin(1) * t318 - qJ(2) * t315) * t209, 0, 0, t327, 0, t239, t296, qJ(1) * t241 + t168 * t315 + t205 * t318, -qJ(1) * t328 + t169 * t315 + t204 * t318, t184 * t363 - t185 * t318, qJ(1) * t138 + t123 * t315 + t135 * t318, 0, 0, -t194, 0, t195, t296, qJ(1) * t197 + t124 * t315 + t145 * t318, qJ(1) * t196 + t125 * t315 + t146 * t318, t112 * t363 - t113 * t318, qJ(1) * t99 + t315 * t84 + t318 * t88, t167 * t315 + t191 * t318, t147 * t315 + t170 * t318, t165 * t315 + t187 * t318, t166 * t315 + t190 * t318, t164 * t315 + t186 * t318, t202 * t318 + t203 * t363, qJ(1) * t133 + t315 * t89 + t318 * t95, qJ(1) * t134 + t315 * t90 + t318 * t96, qJ(1) * t149 + t105 * t318 + t315 * t92, qJ(1) * t86 + t315 * t81 + t318 * t82; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t292, t293, 0, 0, 0, 0, 0, 0, 0, 0, t362, -t314 * t284, t210, pkin(1) * t284 + qJ(2) * t210, 0, 0, t365, 0, -t277, -t343, -pkin(1) * t289 + t314 * t227 - t317 * t335, pkin(1) * t288 + t314 * t228 - t317 * t334, t314 * t184, qJ(2) * t163 - (-pkin(2) * t317 - pkin(5) * t314 - pkin(1)) * t184, 0, 0, -t379, 0, t237, -t343, pkin(1) * t248 + t314 * t160 + t317 * t336, -pkin(1) * t329 + t314 * t161 + t317 * t337, t314 * t112, pkin(1) * t112 + qJ(2) * t109 + t100 * t314 + t107 * t317, t193 * t314 + t259 * t317, t171 * t314 + t242 * t317, t189 * t314 + t262 * t317, t192 * t314 + t258 * t317, t188 * t314 + t260 * t317, t314 * t203, -pkin(1) * t176 + qJ(2) * t158 + t101 * t314 + t126 * t317, -pkin(1) * t177 + qJ(2) * t159 + t102 * t314 + t127 * t317, -pkin(1) * t200 + t314 * t106 + t317 * t338, -pkin(1) * t103 + qJ(2) * t94 + t314 * t83 + t317 * t87;];
tauB_reg = t1;
