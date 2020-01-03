% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:37
% EndTime: 2019-12-31 16:50:42
% DurationCPUTime: 3.97s
% Computational Cost: add. (9256->370), mult. (18102->551), div. (0->0), fcn. (11410->8), ass. (0->258)
t383 = sin(qJ(1));
t386 = cos(qJ(1));
t351 = t386 * g(1) + t383 * g(2);
t388 = qJD(1) ^ 2;
t335 = -t388 * pkin(1) - t351;
t377 = sin(pkin(7));
t378 = cos(pkin(7));
t350 = g(1) * t383 - t386 * g(2);
t390 = qJDD(1) * pkin(1) + t350;
t288 = t377 * t335 - t378 * t390;
t289 = t378 * t335 + t377 * t390;
t400 = t288 * t377 + t378 * t289;
t235 = t288 * t378 - t289 * t377;
t412 = t386 * t235;
t440 = -t383 * t400 + t412;
t415 = t383 * t235;
t195 = t386 * t400 + t415;
t341 = qJDD(1) * t377 + t388 * t378;
t342 = qJDD(1) * t378 - t377 * t388;
t291 = -t383 * t341 + t386 * t342;
t374 = g(3) - qJDD(2);
t316 = qJ(2) * t341 - t374 * t378;
t391 = -qJ(2) * t342 - t374 * t377;
t439 = -pkin(4) * t291 + t383 * t316 + t386 * t391;
t431 = t386 * t341 + t383 * t342;
t437 = pkin(4) * t431 + t386 * t316 - t383 * t391;
t381 = sin(qJ(4));
t382 = sin(qJ(3));
t406 = qJD(1) * qJD(3);
t363 = t382 * t406;
t385 = cos(qJ(3));
t367 = t385 * qJDD(1);
t339 = t367 - t363;
t328 = -qJDD(4) + t339;
t384 = cos(qJ(4));
t409 = qJD(1) * t382;
t331 = -t384 * qJD(3) + t381 * t409;
t333 = qJD(3) * t381 + t384 * t409;
t421 = t331 * t333;
t389 = -t328 - t421;
t436 = t381 * t389;
t434 = t384 * t389;
t408 = qJD(1) * t385;
t358 = -qJD(4) + t408;
t317 = t331 * t358;
t401 = t385 * t406;
t405 = t382 * qJDD(1);
t338 = t401 + t405;
t402 = t331 * qJD(4) - t381 * qJDD(3) - t384 * t338;
t256 = t402 - t317;
t398 = -t384 * qJDD(3) + t338 * t381;
t252 = (qJD(4) + t358) * t333 + t398;
t326 = t331 ^ 2;
t327 = t333 ^ 2;
t356 = t358 ^ 2;
t429 = pkin(3) * t382;
t428 = pkin(3) * t385;
t275 = -t388 * pkin(2) + qJDD(1) * pkin(5) + t289;
t395 = -pkin(6) * t382 - t428;
t336 = t395 * qJD(1);
t362 = t385 * t374;
t387 = qJD(3) ^ 2;
t237 = t362 - qJDD(3) * pkin(3) - t387 * pkin(6) + (qJD(1) * t336 + t275) * t382;
t427 = t237 * t381;
t426 = t237 * t384;
t271 = t328 - t421;
t425 = t271 * t381;
t424 = t271 * t384;
t274 = -qJDD(1) * pkin(2) - t388 * pkin(5) + t288;
t423 = t274 * t382;
t422 = t274 * t385;
t357 = t385 * t388 * t382;
t349 = qJDD(3) - t357;
t420 = t349 * t385;
t419 = t358 * t381;
t418 = t358 * t384;
t348 = qJDD(3) + t357;
t417 = t382 * t348;
t416 = t382 * t349;
t372 = t382 ^ 2;
t411 = t388 * t372;
t392 = -t339 + t363;
t393 = t338 + t401;
t227 = t392 * pkin(3) - t393 * pkin(6) + t274;
t262 = t385 * t275 - t382 * t374;
t238 = -t387 * pkin(3) + qJDD(3) * pkin(6) + t336 * t408 + t262;
t192 = t381 * t227 + t384 * t238;
t373 = t385 ^ 2;
t410 = t372 + t373;
t404 = t385 * t421;
t403 = t382 * t421;
t191 = -t384 * t227 + t238 * t381;
t260 = t275 * t382 + t362;
t215 = t260 * t382 + t385 * t262;
t300 = -t350 * t383 - t386 * t351;
t397 = t377 * t357;
t396 = t378 * t357;
t345 = t386 * qJDD(1) - t383 * t388;
t394 = -pkin(4) * t345 - g(3) * t383;
t162 = -t191 * t384 + t192 * t381;
t163 = t191 * t381 + t192 * t384;
t214 = t260 * t385 - t382 * t262;
t299 = t386 * t350 - t383 * t351;
t370 = t373 * t388;
t355 = -t370 - t387;
t354 = t370 - t387;
t353 = -t387 - t411;
t352 = t387 - t411;
t347 = t370 - t411;
t346 = t370 + t411;
t344 = t383 * qJDD(1) + t386 * t388;
t343 = t410 * qJDD(1);
t340 = t367 - 0.2e1 * t363;
t337 = 0.2e1 * t401 + t405;
t330 = t385 * t348;
t329 = t410 * t406;
t318 = -pkin(4) * t344 + t386 * g(3);
t312 = -t327 + t356;
t311 = t326 - t356;
t310 = t338 * t385 - t372 * t406;
t309 = -t339 * t382 - t373 * t406;
t308 = qJDD(3) * t377 + t329 * t378;
t307 = -qJDD(3) * t378 + t329 * t377;
t306 = -t382 * t353 - t420;
t305 = -t352 * t382 + t330;
t304 = t355 * t385 - t417;
t303 = t354 * t385 - t416;
t302 = t353 * t385 - t416;
t301 = t355 * t382 + t330;
t298 = -t327 + t326;
t297 = -t327 - t356;
t296 = t343 * t378 - t346 * t377;
t295 = t343 * t377 + t346 * t378;
t290 = -t382 * t337 + t340 * t385;
t287 = -t356 - t326;
t285 = -qJD(4) * t333 - t398;
t283 = t378 * t310 - t397;
t282 = t378 * t309 + t397;
t281 = t377 * t310 + t396;
t280 = t377 * t309 - t396;
t279 = t305 * t378 + t377 * t405;
t278 = t378 * t303 + t377 * t367;
t277 = t305 * t377 - t378 * t405;
t276 = t377 * t303 - t378 * t367;
t270 = t326 + t327;
t268 = t306 * t378 + t337 * t377;
t267 = t304 * t378 - t340 * t377;
t266 = t306 * t377 - t337 * t378;
t265 = t304 * t377 + t340 * t378;
t264 = (t331 * t384 - t333 * t381) * t358;
t263 = (-t331 * t381 - t333 * t384) * t358;
t261 = t290 * t378 - t347 * t377;
t259 = t290 * t377 + t347 * t378;
t257 = t317 + t402;
t253 = (-qJD(4) + t358) * t333 - t398;
t251 = t333 * t419 - t384 * t402;
t250 = t333 * t418 + t381 * t402;
t249 = -t285 * t381 - t331 * t418;
t248 = -t285 * t384 + t331 * t419;
t247 = -t383 * t295 + t386 * t296;
t246 = t386 * t295 + t383 * t296;
t245 = t264 * t385 - t382 * t328;
t244 = -pkin(5) * t302 + t422;
t243 = t311 * t384 + t425;
t242 = -t312 * t381 + t434;
t241 = -pkin(5) * t301 + t423;
t240 = -t311 * t381 + t424;
t239 = -t312 * t384 - t436;
t231 = -pkin(2) * t302 + t262;
t230 = -pkin(2) * t301 + t260;
t229 = -t297 * t381 + t424;
t228 = t297 * t384 + t425;
t226 = pkin(1) * t374 + qJ(2) * t400;
t225 = t287 * t384 - t436;
t224 = t287 * t381 + t434;
t221 = t251 * t385 + t403;
t220 = t249 * t385 - t403;
t219 = -t383 * t266 + t386 * t268;
t218 = -t383 * t265 + t386 * t267;
t217 = t386 * t266 + t383 * t268;
t216 = t386 * t265 + t383 * t267;
t212 = -t252 * t384 - t257 * t381;
t211 = t253 * t384 + t256 * t381;
t210 = -t252 * t381 + t257 * t384;
t209 = -t253 * t381 + t256 * t384;
t208 = t245 * t378 - t263 * t377;
t207 = t245 * t377 + t263 * t378;
t206 = t243 * t385 - t382 * t252;
t205 = t242 * t385 - t382 * t257;
t204 = t229 * t385 - t382 * t256;
t203 = t382 * t229 + t256 * t385;
t202 = -qJ(2) * t295 + t214 * t378;
t201 = qJ(2) * t296 + t214 * t377;
t200 = t225 * t385 - t382 * t253;
t199 = t382 * t225 + t253 * t385;
t198 = t215 * t378 + t274 * t377;
t197 = t215 * t377 - t274 * t378;
t196 = -pkin(6) * t228 + t426;
t193 = t211 * t385 - t382 * t298;
t190 = -pkin(6) * t224 + t427;
t189 = t221 * t378 - t250 * t377;
t188 = t220 * t378 - t248 * t377;
t187 = t221 * t377 + t250 * t378;
t186 = t220 * t377 + t248 * t378;
t185 = t212 * t385 - t382 * t270;
t184 = t382 * t212 + t270 * t385;
t183 = -qJ(2) * t266 - t231 * t377 + t244 * t378;
t182 = -qJ(2) * t265 - t230 * t377 + t241 * t378;
t181 = t206 * t378 - t240 * t377;
t180 = t205 * t378 - t239 * t377;
t179 = t206 * t377 + t240 * t378;
t178 = t205 * t377 + t239 * t378;
t177 = t204 * t378 + t228 * t377;
t176 = t204 * t377 - t228 * t378;
t175 = -pkin(1) * t302 + qJ(2) * t268 + t231 * t378 + t244 * t377;
t174 = -pkin(1) * t301 + qJ(2) * t267 + t230 * t378 + t241 * t377;
t173 = t200 * t378 + t224 * t377;
t172 = t200 * t377 - t224 * t378;
t171 = -pkin(3) * t228 + t192;
t170 = -pkin(3) * t224 + t191;
t169 = t193 * t378 - t209 * t377;
t168 = t193 * t377 + t209 * t378;
t167 = t185 * t378 + t210 * t377;
t166 = t185 * t377 - t210 * t378;
t165 = -t383 * t197 + t386 * t198;
t164 = t386 * t197 + t383 * t198;
t161 = -pkin(2) * t203 - pkin(3) * t256 - pkin(6) * t229 - t427;
t160 = -pkin(2) * t199 - pkin(3) * t253 - pkin(6) * t225 + t426;
t159 = -qJ(2) * t197 - (pkin(2) * t377 - pkin(5) * t378) * t214;
t158 = t163 * t385 + t382 * t237;
t157 = t382 * t163 - t237 * t385;
t156 = -pkin(6) * t210 - t162;
t155 = -t383 * t176 + t386 * t177;
t154 = t386 * t176 + t383 * t177;
t153 = -t383 * t172 + t386 * t173;
t152 = t386 * t172 + t383 * t173;
t151 = qJ(2) * t198 - (-pkin(2) * t378 - pkin(5) * t377 - pkin(1)) * t214;
t150 = -pkin(5) * t203 - t382 * t171 + t196 * t385;
t149 = -pkin(5) * t199 - t382 * t170 + t190 * t385;
t148 = -t383 * t166 + t386 * t167;
t147 = t386 * t166 + t383 * t167;
t146 = -pkin(2) * t184 - pkin(3) * t270 - pkin(6) * t212 - t163;
t145 = t158 * t378 + t162 * t377;
t144 = t158 * t377 - t162 * t378;
t143 = -pkin(5) * t184 + t156 * t385 + t210 * t429;
t142 = -pkin(2) * t157 + pkin(3) * t237 - pkin(6) * t163;
t141 = -pkin(5) * t157 + (-pkin(6) * t385 + t429) * t162;
t140 = -qJ(2) * t176 + t150 * t378 - t161 * t377;
t139 = -qJ(2) * t172 + t149 * t378 - t160 * t377;
t138 = -pkin(1) * t203 + qJ(2) * t177 + t150 * t377 + t161 * t378;
t137 = -pkin(1) * t199 + qJ(2) * t173 + t149 * t377 + t160 * t378;
t136 = -t383 * t144 + t386 * t145;
t135 = t386 * t144 + t383 * t145;
t134 = -qJ(2) * t166 + t143 * t378 - t146 * t377;
t133 = -pkin(1) * t184 + qJ(2) * t167 + t143 * t377 + t146 * t378;
t132 = -qJ(2) * t144 + t141 * t378 - t142 * t377;
t131 = -pkin(1) * t157 + qJ(2) * t145 + t141 * t377 + t142 * t378;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t344, -t345, 0, t300, 0, 0, 0, 0, 0, 0, -t431, -t291, 0, t195, 0, 0, 0, 0, 0, 0, t218, t219, t247, t165, 0, 0, 0, 0, 0, 0, t153, t155, t148, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t345, -t344, 0, t299, 0, 0, 0, 0, 0, 0, t291, -t431, 0, -t440, 0, 0, 0, 0, 0, 0, t216, t217, t246, t164, 0, 0, 0, 0, 0, 0, t152, t154, t147, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t374, 0, 0, 0, 0, 0, 0, t301, t302, 0, -t214, 0, 0, 0, 0, 0, 0, t199, t203, t184, t157; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t345, 0, -t344, 0, t394, -t318, -t299, -pkin(4) * t299, 0, 0, t291, 0, -t431, 0, t439, t437, t440, pkin(4) * t440 + qJ(2) * t412 - t383 * t226, -t383 * t281 + t386 * t283, -t383 * t259 + t386 * t261, -t383 * t277 + t386 * t279, -t383 * t280 + t386 * t282, -t383 * t276 + t386 * t278, -t383 * t307 + t386 * t308, -pkin(4) * t216 - t383 * t174 + t386 * t182, -pkin(4) * t217 - t383 * t175 + t386 * t183, -pkin(4) * t246 - t383 * t201 + t386 * t202, -pkin(4) * t164 - t383 * t151 + t386 * t159, -t383 * t187 + t386 * t189, -t383 * t168 + t386 * t169, -t383 * t178 + t386 * t180, -t383 * t186 + t386 * t188, -t383 * t179 + t386 * t181, -t383 * t207 + t386 * t208, -pkin(4) * t152 - t383 * t137 + t386 * t139, -pkin(4) * t154 - t383 * t138 + t386 * t140, -pkin(4) * t147 - t383 * t133 + t386 * t134, -pkin(4) * t135 - t383 * t131 + t386 * t132; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t344, 0, t345, 0, t318, t394, t300, pkin(4) * t300, 0, 0, t431, 0, t291, 0, -t437, t439, t195, pkin(4) * t195 + qJ(2) * t415 + t386 * t226, t386 * t281 + t383 * t283, t386 * t259 + t383 * t261, t386 * t277 + t383 * t279, t386 * t280 + t383 * t282, t386 * t276 + t383 * t278, t386 * t307 + t383 * t308, pkin(4) * t218 + t386 * t174 + t383 * t182, pkin(4) * t219 + t386 * t175 + t383 * t183, pkin(4) * t247 + t386 * t201 + t383 * t202, pkin(4) * t165 + t386 * t151 + t383 * t159, t386 * t187 + t383 * t189, t386 * t168 + t383 * t169, t386 * t178 + t383 * t180, t386 * t186 + t383 * t188, t386 * t179 + t383 * t181, t386 * t207 + t383 * t208, pkin(4) * t153 + t386 * t137 + t383 * t139, pkin(4) * t155 + t386 * t138 + t383 * t140, pkin(4) * t148 + t386 * t133 + t383 * t134, pkin(4) * t136 + t386 * t131 + t383 * t132; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t350, t351, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t342 - t288, -pkin(1) * t341 - t289, 0, -pkin(1) * t235, t393 * t382, t337 * t385 + t382 * t340, t352 * t385 + t417, -t392 * t385, t382 * t354 + t420, 0, pkin(1) * t265 + pkin(2) * t340 + pkin(5) * t304 - t422, pkin(1) * t266 - pkin(2) * t337 + pkin(5) * t306 + t423, pkin(1) * t295 + pkin(2) * t346 + pkin(5) * t343 + t215, pkin(1) * t197 - pkin(2) * t274 + pkin(5) * t215, t382 * t251 - t404, t382 * t211 + t298 * t385, t382 * t242 + t257 * t385, t382 * t249 + t404, t382 * t243 + t252 * t385, t382 * t264 + t328 * t385, pkin(1) * t172 - pkin(2) * t224 + pkin(5) * t200 + t170 * t385 + t382 * t190, pkin(1) * t176 - pkin(2) * t228 + pkin(5) * t204 + t171 * t385 + t382 * t196, pkin(1) * t166 + pkin(5) * t185 + t382 * t156 + (-pkin(2) - t428) * t210, pkin(1) * t144 + pkin(5) * t158 + (-pkin(2) + t395) * t162;];
tauB_reg = t1;
