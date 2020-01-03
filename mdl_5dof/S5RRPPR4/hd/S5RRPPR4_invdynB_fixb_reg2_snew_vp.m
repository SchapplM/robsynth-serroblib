% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:51
% EndTime: 2019-12-31 19:27:56
% DurationCPUTime: 3.60s
% Computational Cost: add. (13072->322), mult. (16428->440), div. (0->0), fcn. (8290->8), ass. (0->201)
t349 = g(3) + qJDD(4);
t350 = sin(pkin(8));
t344 = qJDD(1) + qJDD(2);
t351 = cos(pkin(8));
t345 = (qJD(1) + qJD(2));
t396 = t345 ^ 2;
t363 = t351 * t344 + t350 * t396;
t281 = qJ(4) * t363 + t350 * t349;
t353 = sin(qJ(2));
t356 = cos(qJ(2));
t369 = t350 * t344 - t351 * t396;
t398 = t353 * t369 + t356 * t363;
t400 = -qJ(4) * t369 + t351 * t349;
t208 = -pkin(6) * t398 + t356 * t281 - t353 * t400;
t354 = sin(qJ(1));
t357 = cos(qJ(1));
t252 = t353 * t363 - t356 * t369;
t413 = -pkin(6) * t252 + t353 * t281 + t356 * t400;
t414 = -t354 * t252 + t357 * t398;
t424 = -pkin(5) * t414 + t357 * t208 - t354 * t413;
t216 = t357 * t252 + t354 * t398;
t423 = -pkin(5) * t216 + t354 * t208 + t357 * t413;
t334 = t344 * qJ(3);
t325 = t357 * g(1) + t354 * g(2);
t360 = qJD(1) ^ 2;
t314 = -t360 * pkin(1) - t325;
t324 = t354 * g(1) - t357 * g(2);
t362 = qJDD(1) * pkin(1) + t324;
t267 = t356 * t314 + t353 * t362;
t377 = (2 * qJD(3) * t345) + t267;
t366 = t334 + t377;
t395 = pkin(2) + pkin(3);
t232 = -t395 * t396 + t366;
t335 = t344 * pkin(2);
t266 = t353 * t314 - t356 * t362;
t364 = -qJDD(3) - t266;
t245 = -qJ(3) * t396 - t335 - t364;
t361 = -t344 * pkin(3) + t245;
t196 = t350 * t232 - t351 * t361;
t197 = t351 * t232 + t350 * t361;
t171 = t351 * t196 - t350 * t197;
t172 = t350 * t196 + t351 * t197;
t148 = t353 * t171 - t356 * t172;
t407 = t356 * t171 + t353 * t172;
t142 = -t354 * t148 + t357 * t407;
t421 = t357 * t148 + t354 * t407;
t307 = t353 * t344 + t356 * t396;
t310 = -t356 * t344 + t353 * t396;
t259 = t357 * t307 - t354 * t310;
t285 = pkin(6) * t310 - t353 * g(3);
t288 = pkin(6) * t307 - t356 * g(3);
t419 = pkin(5) * t259 - t354 * t285 + t357 * t288;
t263 = t354 * t307 + t357 * t310;
t418 = pkin(5) * t263 + t357 * t285 + t354 * t288;
t372 = t353 * t266 + t356 * t267;
t225 = t356 * t266 - t353 * t267;
t381 = t357 * t225;
t415 = -t354 * t372 + t381;
t386 = t354 * t225;
t182 = t357 * t372 + t386;
t412 = pkin(1) * t307;
t240 = -pkin(2) * t396 + t366;
t199 = t353 * t240 - t356 * t245;
t373 = t356 * t240 + t353 * t245;
t174 = -t354 * t199 + t357 * t373;
t173 = t357 * t199 + t354 * t373;
t352 = sin(qJ(5));
t347 = t352 ^ 2;
t393 = t347 * t396;
t188 = t344 * pkin(4) - pkin(7) * t396 + t196;
t390 = t352 * t188;
t355 = cos(qJ(5));
t323 = t355 * t396 * t352;
t315 = qJDD(5) + t323;
t389 = t352 * t315;
t316 = qJDD(5) - t323;
t388 = t352 * t316;
t387 = t352 * t344;
t384 = t355 * t188;
t383 = t355 * t315;
t382 = t355 * t316;
t329 = t355 * t344;
t189 = -pkin(4) * t396 - t344 * pkin(7) + t197;
t184 = t355 * t189 + t352 * t349;
t348 = t355 ^ 2;
t379 = -t347 - t348;
t378 = qJD(5) * t345;
t376 = t352 * t378;
t375 = t355 * t378;
t183 = t352 * t189 - t355 * t349;
t279 = -t354 * t324 - t357 * t325;
t368 = t350 * t323;
t367 = t351 * t323;
t318 = t357 * qJDD(1) - t354 * t360;
t365 = -pkin(5) * t318 - t354 * g(3);
t161 = t355 * t183 - t352 * t184;
t162 = t352 * t183 + t355 * t184;
t278 = t357 * t324 - t354 * t325;
t359 = qJD(5) ^ 2;
t358 = pkin(1) * g(3);
t330 = t348 * t396;
t322 = -t330 - t359;
t321 = t330 - t359;
t320 = -t359 - t393;
t319 = t359 - t393;
t317 = t354 * qJDD(1) + t357 * t360;
t312 = t330 - t393;
t311 = t330 + t393;
t306 = t379 * t344;
t298 = -t329 + 0.2e1 * t376;
t297 = -t329 + t376;
t296 = -t375 - t387;
t295 = -0.2e1 * t375 - t387;
t294 = pkin(1) * t310;
t293 = -pkin(5) * t317 + t357 * g(3);
t292 = t379 * t378;
t277 = t350 * qJDD(5) + t351 * t292;
t276 = t351 * qJDD(5) - t350 * t292;
t275 = t355 * t296 + t347 * t378;
t274 = -t352 * t297 + t348 * t378;
t273 = -t352 * t320 - t382;
t272 = -t352 * t319 + t383;
t271 = t355 * t322 - t389;
t270 = t355 * t321 - t388;
t269 = t355 * t320 - t388;
t268 = t352 * t322 + t383;
t258 = t351 * t306 - t350 * t311;
t257 = t350 * t306 + t351 * t311;
t250 = -t352 * t295 + t355 * t298;
t249 = t351 * t272 - t350 * t387;
t248 = t351 * t270 - t350 * t329;
t247 = -t350 * t272 - t351 * t387;
t246 = -t350 * t270 - t351 * t329;
t244 = t351 * t275 - t368;
t243 = t351 * t274 + t368;
t242 = -t350 * t275 - t367;
t241 = -t350 * t274 + t367;
t238 = t351 * t273 + t350 * t295;
t237 = t351 * t271 - t350 * t298;
t236 = t350 * t273 - t351 * t295;
t235 = t350 * t271 + t351 * t298;
t230 = -t353 * t276 + t356 * t277;
t229 = t356 * t276 + t353 * t277;
t228 = t351 * t250 - t350 * t312;
t227 = -t350 * t250 - t351 * t312;
t222 = pkin(6) * t372 + t358;
t221 = t353 * t257 + t356 * t258;
t220 = -t356 * t257 + t353 * t258;
t213 = -t353 * t247 + t356 * t249;
t212 = -t353 * t246 + t356 * t248;
t211 = t356 * t247 + t353 * t249;
t210 = t356 * t246 + t353 * t248;
t205 = -t353 * t242 + t356 * t244;
t204 = -t353 * t241 + t356 * t243;
t203 = t356 * t242 + t353 * t244;
t202 = t356 * t241 + t353 * t243;
t195 = -pkin(6) * t199 + (-pkin(2) * t353 + qJ(3) * t356) * g(3);
t194 = t353 * t236 + t356 * t238;
t193 = t353 * t235 + t356 * t237;
t192 = -t356 * t236 + t353 * t238;
t191 = -t356 * t235 + t353 * t237;
t190 = pkin(6) * t373 + t358 + (pkin(2) * t356 + qJ(3) * t353) * g(3);
t186 = -t353 * t227 + t356 * t228;
t185 = t356 * t227 + t353 * t228;
t180 = -pkin(7) * t269 + t384;
t179 = -pkin(7) * t268 + t390;
t178 = -pkin(4) * t269 + t184;
t177 = -pkin(4) * t268 + t183;
t176 = -t354 * t220 + t357 * t221;
t175 = t357 * t220 + t354 * t221;
t168 = -t354 * t192 + t357 * t194;
t167 = -t354 * t191 + t357 * t193;
t166 = t357 * t192 + t354 * t194;
t165 = t357 * t191 + t354 * t193;
t164 = qJ(3) * t349 + qJ(4) * t171;
t163 = -qJ(4) * t172 + t395 * t349;
t159 = -qJ(4) * t257 + t351 * t161;
t158 = -qJ(4) * t258 - t350 * t161;
t157 = t351 * t162 + t350 * t188;
t156 = t350 * t162 - t351 * t188;
t155 = qJ(3) * t269 - qJ(4) * t236 - t350 * t178 + t351 * t180;
t154 = qJ(3) * t268 - qJ(4) * t235 - t350 * t177 + t351 * t179;
t153 = -qJ(4) * t238 - t351 * t178 - t350 * t180 + t395 * t269;
t152 = -qJ(4) * t237 - t351 * t177 - t350 * t179 + t395 * t268;
t147 = -pkin(6) * t220 - t353 * t158 + t356 * t159;
t146 = pkin(6) * t221 + t356 * t158 + t353 * t159;
t145 = t353 * t156 + t356 * t157;
t144 = -t356 * t156 + t353 * t157;
t141 = -pkin(6) * t407 - t353 * t163 + t356 * t164;
t140 = pkin(1) * t349 - pkin(6) * t148 + t356 * t163 + t353 * t164;
t139 = -pkin(6) * t192 - t353 * t153 + t356 * t155;
t138 = -pkin(6) * t191 - t353 * t152 + t356 * t154;
t137 = pkin(1) * t269 + pkin(6) * t194 + t356 * t153 + t353 * t155;
t136 = pkin(1) * t268 + pkin(6) * t193 + t356 * t152 + t353 * t154;
t135 = -qJ(4) * t156 - (pkin(4) * t350 - pkin(7) * t351 + qJ(3)) * t161;
t134 = -qJ(4) * t157 - (pkin(4) * t351 + pkin(7) * t350 + t395) * t161;
t133 = -t354 * t144 + t357 * t145;
t132 = t357 * t144 + t354 * t145;
t131 = -pkin(6) * t144 - t353 * t134 + t356 * t135;
t130 = -pkin(1) * t161 + pkin(6) * t145 + t356 * t134 + t353 * t135;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t317, -t318, 0, t279, 0, 0, 0, 0, 0, 0, -t259, t263, 0, t182, 0, 0, 0, 0, 0, 0, -t259, 0, -t263, t174, 0, 0, 0, 0, 0, 0, -t216, t414, 0, -t421, 0, 0, 0, 0, 0, 0, t167, t168, t176, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t318, -t317, 0, t278, 0, 0, 0, 0, 0, 0, -t263, -t259, 0, -t415, 0, 0, 0, 0, 0, 0, -t263, 0, t259, t173, 0, 0, 0, 0, 0, 0, t414, t216, 0, t142, 0, 0, 0, 0, 0, 0, t165, t166, t175, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349, 0, 0, 0, 0, 0, 0, -t268, -t269, 0, t161; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t318, 0, -t317, 0, t365, -t293, -t278, -pkin(5) * t278, 0, 0, -t263, 0, -t259, 0, t418, t419, t415, pkin(5) * t415 + pkin(6) * t381 - t354 * t222, 0, -t263, 0, 0, t259, 0, t418, -t173, -t419, -pkin(5) * t173 - t354 * t190 + t357 * t195, 0, 0, -t414, 0, -t216, 0, t424, t423, t142, -pkin(5) * t142 - t354 * t140 + t357 * t141, -t354 * t203 + t357 * t205, -t354 * t185 + t357 * t186, -t354 * t211 + t357 * t213, -t354 * t202 + t357 * t204, -t354 * t210 + t357 * t212, -t354 * t229 + t357 * t230, -pkin(5) * t165 - t354 * t136 + t357 * t138, -pkin(5) * t166 - t354 * t137 + t357 * t139, -pkin(5) * t175 - t354 * t146 + t357 * t147, -pkin(5) * t132 - t354 * t130 + t357 * t131; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t317, 0, t318, 0, t293, t365, t279, pkin(5) * t279, 0, 0, t259, 0, -t263, 0, -t419, t418, t182, pkin(5) * t182 + pkin(6) * t386 + t357 * t222, 0, t259, 0, 0, t263, 0, -t419, t174, -t418, pkin(5) * t174 + t357 * t190 + t354 * t195, 0, 0, -t216, 0, t414, 0, t423, -t424, t421, -pkin(5) * t421 + t357 * t140 + t354 * t141, t357 * t203 + t354 * t205, t357 * t185 + t354 * t186, t357 * t211 + t354 * t213, t357 * t202 + t354 * t204, t357 * t210 + t354 * t212, t357 * t229 + t354 * t230, pkin(5) * t167 + t357 * t136 + t354 * t138, pkin(5) * t168 + t357 * t137 + t354 * t139, pkin(5) * t176 + t357 * t146 + t354 * t147, pkin(5) * t133 + t357 * t130 + t354 * t131; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t324, t325, 0, 0, 0, 0, 0, 0, 0, t344, -t266 - t294, -t267 - t412, 0, -pkin(1) * t225, 0, 0, 0, t344, 0, 0, -t294 + 0.2e1 * t335 + t364, 0, 0.2e1 * t334 + t377 + t412, pkin(1) * t199 - pkin(2) * t245 + qJ(3) * t240, 0, 0, 0, 0, 0, t344, pkin(1) * t398 + qJ(3) * t369 + t363 * t395 + t196, pkin(1) * t252 + qJ(3) * t363 - t369 * t395 + t197, 0, pkin(1) * t407 + qJ(3) * t172 + t171 * t395, (-t296 + t375) * t352, -t355 * t295 - t352 * t298, -t355 * t319 - t389, (-t297 - t376) * t355, -t352 * t321 - t382, 0, pkin(1) * t191 - pkin(4) * t298 - pkin(7) * t271 + qJ(3) * t237 - t395 * t235 + t384, pkin(1) * t192 + pkin(4) * t295 - pkin(7) * t273 + qJ(3) * t238 - t395 * t236 - t390, pkin(1) * t220 - pkin(4) * t311 - pkin(7) * t306 + qJ(3) * t258 - t395 * t257 - t162, pkin(1) * t144 + pkin(4) * t188 - pkin(7) * t162 + qJ(3) * t157 - t395 * t156;];
tauB_reg = t1;
