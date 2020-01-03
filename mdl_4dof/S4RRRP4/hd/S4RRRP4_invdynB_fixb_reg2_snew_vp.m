% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:45
% EndTime: 2019-12-31 17:15:48
% DurationCPUTime: 3.16s
% Computational Cost: add. (9318->344), mult. (20311->474), div. (0->0), fcn. (13266->6), ass. (0->263)
t398 = sin(qJ(3));
t401 = cos(qJ(3));
t402 = cos(qJ(2));
t399 = sin(qJ(2));
t433 = qJD(1) * t399;
t359 = -t401 * t402 * qJD(1) + t398 * t433;
t361 = (t398 * t402 + t399 * t401) * qJD(1);
t330 = t361 * t359;
t394 = qJDD(2) + qJDD(3);
t417 = -t394 + t330;
t476 = t417 * pkin(3);
t452 = t398 * t417;
t442 = t401 * t417;
t430 = qJD(1) * qJD(2);
t421 = t402 * t430;
t429 = t399 * qJDD(1);
t370 = t421 + t429;
t389 = t402 * qJDD(1);
t422 = t399 * t430;
t371 = t389 - t422;
t416 = t398 * t370 - t401 * t371;
t312 = -t361 * qJD(3) - t416;
t397 = t402 ^ 2;
t405 = qJD(1) ^ 2;
t400 = sin(qJ(1));
t403 = cos(qJ(1));
t380 = t400 * g(1) - t403 * g(2);
t411 = qJDD(1) * pkin(1) + t380;
t412 = qJD(2) * pkin(2) - pkin(6) * t433;
t315 = t371 * pkin(2) + (pkin(6) * t397 + pkin(5)) * t405 - t412 * t433 + t411;
t395 = qJD(2) + qJD(3);
t346 = t395 * pkin(3) - t361 * qJ(4);
t406 = -t361 * t346 - qJDD(4) + t315;
t477 = t312 * pkin(3) + t406;
t352 = t395 * t359;
t408 = t359 * qJD(3) - t401 * t370 - t398 * t371;
t475 = -t352 - t408;
t358 = t361 ^ 2;
t393 = t395 ^ 2;
t343 = -t358 - t393;
t327 = t330 + t394;
t453 = t398 * t327;
t292 = t401 * t343 - t453;
t443 = t401 * t327;
t293 = -t398 * t343 - t443;
t258 = t402 * t292 + t399 * t293;
t474 = -pkin(1) * t258 - pkin(2) * t292;
t357 = t359 ^ 2;
t325 = -t393 - t357;
t275 = t398 * t325 - t442;
t276 = t401 * t325 + t452;
t245 = t402 * t275 + t399 * t276;
t473 = -pkin(1) * t245 - pkin(2) * t275;
t287 = (qJD(3) - t395) * t361 + t416;
t291 = -t352 + t408;
t255 = -t287 * t398 + t401 * t291;
t257 = -t287 * t401 - t398 * t291;
t220 = -t399 * t255 + t402 * t257;
t314 = -t357 - t358;
t209 = t400 * t220 - t403 * t314;
t468 = pkin(4) * t209;
t246 = -t399 * t275 + t402 * t276;
t286 = (qJD(3) + t395) * t361 + t416;
t225 = t400 * t246 - t403 * t286;
t467 = pkin(4) * t225;
t259 = -t399 * t292 + t402 * t293;
t231 = t400 * t259 - t403 * t475;
t466 = pkin(4) * t231;
t218 = t402 * t255 + t399 * t257;
t465 = pkin(5) * t218;
t464 = pkin(5) * t245;
t463 = pkin(5) * t258;
t462 = pkin(6) * t255;
t461 = pkin(6) * t275;
t460 = pkin(6) * t292;
t458 = t395 * t398;
t457 = t395 * t401;
t396 = t399 ^ 2;
t456 = t396 * t405;
t391 = t397 * t405;
t446 = t399 * t405;
t381 = t403 * g(1) + t400 * g(2);
t363 = -t405 * pkin(1) + qJDD(1) * pkin(5) - t381;
t449 = t399 * t363;
t302 = qJDD(2) * pkin(2) - t370 * pkin(6) - t449 + (pkin(2) * t446 + pkin(6) * t430 - g(3)) * t402;
t345 = -t399 * g(3) + t402 * t363;
t306 = -pkin(2) * t391 + t371 * pkin(6) - qJD(2) * t412 + t345;
t267 = -t401 * t302 + t398 * t306;
t424 = -qJ(4) * t408 + t267;
t410 = -qJ(4) * t352 - t424;
t432 = qJD(4) * t361;
t223 = t410 - 0.2e1 * t432 - t476;
t455 = t398 * t223;
t454 = t398 * t315;
t268 = t398 * t302 + t401 * t306;
t227 = -t401 * t267 + t398 * t268;
t451 = t399 * t227;
t362 = t405 * pkin(5) + t411;
t450 = t399 * t362;
t387 = t402 * t446;
t378 = qJDD(2) + t387;
t448 = t399 * t378;
t379 = qJDD(2) - t387;
t447 = t399 * t379;
t445 = t401 * t223;
t444 = t401 * t315;
t441 = t402 * t227;
t440 = t402 * t362;
t439 = t402 * t379;
t438 = -pkin(1) * t314 + pkin(5) * t220;
t437 = -pkin(1) * t286 + pkin(5) * t246;
t436 = -pkin(1) * t475 + pkin(5) * t259;
t435 = -t343 - t357;
t434 = t396 + t397;
t428 = t400 * qJDD(1);
t427 = t403 * qJDD(1);
t426 = t400 * t330;
t425 = t403 * t330;
t420 = -pkin(2) * t314 + pkin(6) * t257;
t419 = -pkin(2) * t286 + pkin(6) * t276;
t418 = -pkin(2) * t475 + pkin(6) * t293;
t228 = t398 * t267 + t401 * t268;
t344 = t402 * g(3) + t449;
t305 = t399 * t344 + t402 * t345;
t336 = -t400 * t380 - t403 * t381;
t415 = t400 * t387;
t414 = t403 * t387;
t200 = -pkin(1) * t218 - pkin(2) * t255;
t375 = -t400 * t405 + t427;
t413 = -pkin(4) * t375 - t400 * g(3);
t304 = t402 * t344 - t399 * t345;
t335 = t403 * t380 - t400 * t381;
t409 = t312 * qJ(4) - 0.2e1 * qJD(4) * t359 - t395 * t346 + t268;
t404 = qJD(2) ^ 2;
t385 = -t391 - t404;
t384 = t391 - t404;
t383 = -t404 - t456;
t382 = t404 - t456;
t377 = t391 - t456;
t376 = t391 + t456;
t374 = t403 * t405 + t428;
t373 = t434 * qJDD(1);
t372 = t389 - 0.2e1 * t422;
t369 = 0.2e1 * t421 + t429;
t367 = t402 * t378;
t366 = t434 * t430;
t356 = -pkin(4) * t374 + t403 * g(3);
t355 = 0.2e1 * t432;
t350 = -t358 + t393;
t349 = t357 - t393;
t348 = t402 * t370 - t396 * t430;
t347 = -t399 * t371 - t397 * t430;
t342 = -t399 * t383 - t439;
t341 = -t399 * t382 + t367;
t340 = t402 * t385 - t448;
t339 = t402 * t384 - t447;
t338 = t402 * t383 - t447;
t337 = t399 * t385 + t367;
t333 = t403 * t373 - t400 * t376;
t332 = t400 * t373 + t403 * t376;
t331 = -t399 * t369 + t402 * t372;
t329 = t357 - t358;
t323 = t403 * t342 + t400 * t369;
t322 = t403 * t340 - t400 * t372;
t321 = t400 * t342 - t403 * t369;
t320 = t400 * t340 + t403 * t372;
t319 = (-t359 * t401 + t361 * t398) * t395;
t318 = (-t359 * t398 - t361 * t401) * t395;
t317 = -pkin(5) * t338 - t440;
t316 = -pkin(5) * t337 - t450;
t311 = -pkin(1) * t338 + t345;
t310 = -pkin(1) * t337 + t344;
t297 = t401 * t349 - t453;
t296 = -t398 * t350 - t442;
t295 = t398 * t349 + t443;
t294 = t401 * t350 - t452;
t282 = -t361 * t458 - t401 * t408;
t281 = t361 * t457 - t398 * t408;
t280 = -t398 * t312 + t359 * t457;
t279 = t401 * t312 + t359 * t458;
t278 = t403 * t305 - t400 * t362;
t277 = t400 * t305 + t403 * t362;
t273 = -t399 * t318 + t402 * t319;
t272 = t402 * t318 + t399 * t319;
t271 = -pkin(3) * t475 - qJ(4) * t327;
t270 = t403 * t273 + t400 * t394;
t269 = t400 * t273 - t403 * t394;
t265 = -t444 - t460;
t264 = -t454 - t461;
t263 = -t399 * t295 + t402 * t297;
t262 = -t399 * t294 + t402 * t296;
t261 = t402 * t295 + t399 * t297;
t260 = t402 * t294 + t399 * t296;
t256 = -t401 * t286 - t398 * t475;
t254 = -t398 * t286 + t401 * t475;
t251 = t357 * qJ(4) + t477;
t250 = -t399 * t281 + t402 * t282;
t249 = -t399 * t279 + t402 * t280;
t248 = t402 * t281 + t399 * t282;
t247 = t402 * t279 + t399 * t280;
t243 = t435 * qJ(4) - t477;
t242 = t418 - t454;
t241 = t403 * t250 + t426;
t240 = t403 * t249 - t426;
t239 = t400 * t250 - t425;
t238 = t400 * t249 + t425;
t237 = t419 + t444;
t236 = t403 * t263 - t400 * t287;
t235 = t403 * t262 - t400 * t291;
t234 = t400 * t263 + t403 * t287;
t233 = t400 * t262 + t403 * t291;
t232 = t403 * t259 + t400 * t475;
t230 = pkin(4) * t232;
t229 = -t357 * pkin(3) + t409;
t226 = t403 * t246 + t400 * t286;
t224 = pkin(4) * t226;
t222 = (t325 + t357) * qJ(4) + (-t286 + t312) * pkin(3) + t406;
t221 = pkin(2) * t315 + pkin(6) * t228;
t219 = -t399 * t254 + t402 * t256;
t217 = t402 * t254 + t399 * t256;
t215 = t355 + (-t291 + t352) * qJ(4) + t476 + t424;
t214 = t268 + t474;
t213 = t403 * t219 - t400 * t329;
t212 = t400 * t219 + t403 * t329;
t211 = -qJ(4) * t287 + (-t314 - t357) * pkin(3) + t409;
t210 = t403 * t220 + t400 * t314;
t208 = pkin(4) * t210;
t207 = t267 + t473;
t206 = t401 * t243 - t398 * t271 - t460;
t205 = -t227 - t462;
t204 = pkin(3) * t251 + qJ(4) * t229;
t203 = qJ(4) * t442 - t398 * t222 - t461;
t202 = t228 + t420;
t201 = t398 * t243 + t401 * t271 + t418;
t199 = qJ(4) * t452 + t401 * t222 + t419;
t198 = -t399 * t242 + t402 * t265 - t463;
t197 = t435 * pkin(3) + t409 + t474;
t196 = t401 * t229 - t455;
t195 = t398 * t229 + t445;
t194 = t402 * t228 - t451;
t193 = t399 * t228 + t441;
t192 = -t399 * t237 + t402 * t264 - t464;
t191 = -pkin(3) * t291 + t200;
t190 = t355 - t410 + t473 + 0.2e1 * t476;
t189 = t403 * t194 - t400 * t315;
t188 = t400 * t194 + t403 * t315;
t187 = -pkin(1) * t193 - pkin(2) * t227;
t186 = -t398 * t211 + t401 * t215 - t462;
t185 = t401 * t211 + t398 * t215 + t420;
t184 = -t399 * t201 + t402 * t206 - t463;
t183 = -t399 * t195 + t402 * t196;
t182 = t402 * t195 + t399 * t196;
t181 = -t399 * t199 + t402 * t203 - t464;
t180 = -pkin(5) * t193 - pkin(6) * t441 - t399 * t221;
t179 = -pkin(6) * t195 - qJ(4) * t445 - t398 * t204;
t178 = -t399 * t202 + t402 * t205 - t465;
t177 = t403 * t183 - t400 * t251;
t176 = t400 * t183 + t403 * t251;
t175 = pkin(2) * t251 + pkin(6) * t196 - qJ(4) * t455 + t401 * t204;
t174 = -pkin(1) * t182 - pkin(2) * t195 - pkin(3) * t223;
t173 = -t399 * t185 + t402 * t186 - t465;
t172 = -pkin(5) * t182 - t399 * t175 + t402 * t179;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t374, -t375, 0, t336, 0, 0, 0, 0, 0, 0, t322, t323, t333, t278, 0, 0, 0, 0, 0, 0, t226, t232, t210, t189, 0, 0, 0, 0, 0, 0, t226, t232, t210, t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t375, -t374, 0, t335, 0, 0, 0, 0, 0, 0, t320, t321, t332, t277, 0, 0, 0, 0, 0, 0, t225, t231, t209, t188, 0, 0, 0, 0, 0, 0, t225, t231, t209, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t337, t338, 0, -t304, 0, 0, 0, 0, 0, 0, t245, t258, t218, t193, 0, 0, 0, 0, 0, 0, t245, t258, t218, t182; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t375, 0, -t374, 0, t413, -t356, -t335, -pkin(4) * t335, t403 * t348 - t415, t403 * t331 - t400 * t377, t403 * t341 + t399 * t428, t403 * t347 + t415, t403 * t339 + t400 * t389, t400 * qJDD(2) + t403 * t366, -pkin(4) * t320 - t400 * t310 + t403 * t316, -pkin(4) * t321 - t400 * t311 + t403 * t317, -pkin(4) * t332 + t403 * t304, -pkin(4) * t277 - (pkin(1) * t400 - pkin(5) * t403) * t304, t241, t213, t235, t240, t236, t270, t403 * t192 - t400 * t207 - t467, t403 * t198 - t400 * t214 - t466, t403 * t178 - t400 * t200 - t468, -pkin(4) * t188 + t403 * t180 - t400 * t187, t241, t213, t235, t240, t236, t270, t403 * t181 - t400 * t190 - t467, t403 * t184 - t400 * t197 - t466, t403 * t173 - t400 * t191 - t468, -pkin(4) * t176 + t403 * t172 - t400 * t174; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t374, 0, t375, 0, t356, t413, t336, pkin(4) * t336, t400 * t348 + t414, t400 * t331 + t403 * t377, t400 * t341 - t399 * t427, t400 * t347 - t414, t400 * t339 - t402 * t427, -t403 * qJDD(2) + t400 * t366, pkin(4) * t322 + t403 * t310 + t400 * t316, pkin(4) * t323 + t403 * t311 + t400 * t317, pkin(4) * t333 + t400 * t304, pkin(4) * t278 - (-pkin(1) * t403 - pkin(5) * t400) * t304, t239, t212, t233, t238, t234, t269, t400 * t192 + t403 * t207 + t224, t400 * t198 + t403 * t214 + t230, t400 * t178 + t403 * t200 + t208, pkin(4) * t189 + t400 * t180 + t403 * t187, t239, t212, t233, t238, t234, t269, t400 * t181 + t403 * t190 + t224, t400 * t184 + t403 * t197 + t230, t400 * t173 + t403 * t191 + t208, pkin(4) * t177 + t400 * t172 + t403 * t174; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t380, t381, 0, 0, (t370 + t421) * t399, t402 * t369 + t399 * t372, t402 * t382 + t448, (t371 - t422) * t402, t399 * t384 + t439, 0, pkin(1) * t372 + pkin(5) * t340 + t440, -pkin(1) * t369 + pkin(5) * t342 - t450, pkin(1) * t376 + pkin(5) * t373 + t305, pkin(1) * t362 + pkin(5) * t305, t248, t217, t260, t247, t261, t272, t402 * t237 + t399 * t264 + t437, t402 * t242 + t399 * t265 + t436, t402 * t202 + t399 * t205 + t438, pkin(1) * t315 + pkin(5) * t194 - pkin(6) * t451 + t402 * t221, t248, t217, t260, t247, t261, t272, t402 * t199 + t399 * t203 + t437, t402 * t201 + t399 * t206 + t436, t402 * t185 + t399 * t186 + t438, pkin(1) * t251 + pkin(5) * t183 + t402 * t175 + t399 * t179;];
tauB_reg = t1;
