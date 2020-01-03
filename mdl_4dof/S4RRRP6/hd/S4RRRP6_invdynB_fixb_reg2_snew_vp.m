% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:14
% EndTime: 2019-12-31 17:19:18
% DurationCPUTime: 3.13s
% Computational Cost: add. (8187->339), mult. (16714->468), div. (0->0), fcn. (10688->6), ass. (0->260)
t410 = sin(qJ(1));
t413 = cos(qJ(1));
t388 = g(1) * t413 + g(2) * t410;
t414 = qJD(1) ^ 2;
t364 = -pkin(1) * t414 + qJDD(1) * pkin(5) - t388;
t409 = sin(qJ(2));
t412 = cos(qJ(2));
t478 = pkin(2) * t412;
t427 = -pkin(6) * t409 - t478;
t375 = t427 * qJD(1);
t468 = g(3) * t412;
t482 = qJD(2) ^ 2;
t318 = (qJD(1) * t375 + t364) * t409 - qJDD(2) * pkin(2) - t482 * pkin(6) + t468;
t442 = qJD(1) * qJD(2);
t398 = t409 * t442;
t439 = qJDD(1) * t412;
t378 = -t398 + t439;
t368 = -qJDD(3) + t378;
t408 = sin(qJ(3));
t411 = cos(qJ(3));
t446 = qJD(1) * t409;
t371 = -t411 * qJD(2) + t408 * t446;
t373 = qJD(2) * t408 + t411 * t446;
t457 = t371 * t373;
t433 = t368 + t457;
t484 = t433 * pkin(3);
t434 = t412 * t442;
t441 = qJDD(1) * t409;
t377 = t434 + t441;
t432 = -t411 * qJDD(2) + t408 * t377;
t333 = -qJD(3) * t373 - t432;
t445 = qJD(1) * t412;
t395 = -qJD(3) + t445;
t350 = -pkin(3) * t395 - qJ(4) * t373;
t366 = t371 ^ 2;
t261 = -t333 * pkin(3) - t366 * qJ(4) + t350 * t373 + qJDD(4) + t318;
t461 = t433 * t408;
t460 = t433 * t411;
t358 = t371 * t395;
t416 = qJD(3) * t371 - qJDD(2) * t408 - t377 * t411;
t483 = t358 - t416;
t304 = (qJD(3) + t395) * t373 + t432;
t367 = t373 ^ 2;
t393 = t395 ^ 2;
t335 = -t393 - t366;
t283 = t335 * t408 - t460;
t481 = pkin(2) * t283;
t337 = -t367 - t393;
t327 = t368 - t457;
t463 = t327 * t408;
t288 = t337 * t411 + t463;
t480 = pkin(2) * t288;
t479 = pkin(2) * t409;
t308 = t358 + t416;
t274 = -t304 * t411 - t308 * t408;
t326 = -t366 - t367;
t248 = t274 * t412 + t326 * t409;
t272 = -t304 * t408 + t308 * t411;
t220 = t248 * t410 - t272 * t413;
t477 = pkin(4) * t220;
t284 = t335 * t411 + t461;
t303 = (qJD(3) - t395) * t373 + t432;
t257 = t284 * t412 + t303 * t409;
t228 = t257 * t410 - t283 * t413;
t476 = pkin(4) * t228;
t462 = t327 * t411;
t289 = -t337 * t408 + t462;
t260 = t289 * t412 + t409 * t483;
t235 = t260 * t410 - t288 * t413;
t475 = pkin(4) * t235;
t247 = t274 * t409 - t326 * t412;
t474 = pkin(5) * t247;
t256 = t284 * t409 - t303 * t412;
t473 = pkin(5) * t256;
t259 = t289 * t409 - t412 * t483;
t472 = pkin(5) * t259;
t471 = pkin(6) * t272;
t470 = pkin(6) * t283;
t469 = pkin(6) * t288;
t387 = g(1) * t410 - t413 * g(2);
t363 = qJDD(1) * pkin(1) + pkin(5) * t414 + t387;
t424 = -t378 + t398;
t425 = t377 + t434;
t302 = pkin(2) * t424 - pkin(6) * t425 - t363;
t352 = -g(3) * t409 + t412 * t364;
t319 = -t482 * pkin(2) + qJDD(2) * pkin(6) + t375 * t445 + t352;
t275 = -t411 * t302 + t408 * t319;
t435 = -qJ(4) * t416 + t275;
t422 = qJ(4) * t358 - t435;
t444 = qJD(4) * t373;
t233 = t422 - 0.2e1 * t444 - t484;
t467 = t233 * t408;
t466 = t233 * t411;
t465 = t318 * t408;
t464 = t318 * t411;
t459 = t363 * t409;
t458 = t363 * t412;
t394 = t412 * t414 * t409;
t385 = qJDD(2) + t394;
t456 = t385 * t409;
t386 = qJDD(2) - t394;
t455 = t386 * t409;
t454 = t386 * t412;
t453 = t395 * t408;
t452 = t395 * t411;
t404 = t409 ^ 2;
t451 = t404 * t414;
t450 = -pkin(1) * t272 + pkin(5) * t248;
t449 = -pkin(1) * t283 + pkin(5) * t257;
t448 = -pkin(1) * t288 + pkin(5) * t260;
t276 = t408 * t302 + t411 * t319;
t405 = t412 ^ 2;
t447 = t404 + t405;
t440 = qJDD(1) * t410;
t438 = qJDD(1) * t413;
t437 = t409 * t457;
t436 = t412 * t457;
t351 = t364 * t409 + t468;
t311 = t351 * t409 + t412 * t352;
t343 = -t387 * t410 - t413 * t388;
t430 = t410 * t394;
t429 = t413 * t394;
t382 = -t410 * t414 + t438;
t426 = -pkin(4) * t382 - g(3) * t410;
t231 = -t275 * t411 + t276 * t408;
t232 = t275 * t408 + t276 * t411;
t310 = t351 * t412 - t352 * t409;
t342 = t387 * t413 - t388 * t410;
t421 = t333 * qJ(4) - 0.2e1 * qJD(4) * t371 + t350 * t395 + t276;
t420 = -pkin(1) * t247 + pkin(2) * t326 - pkin(6) * t274;
t419 = -pkin(1) * t256 + pkin(2) * t303 - pkin(6) * t284;
t418 = -pkin(1) * t259 + pkin(2) * t483 - pkin(6) * t289;
t402 = t405 * t414;
t392 = -t402 - t482;
t391 = t402 - t482;
t390 = -t451 - t482;
t389 = -t451 + t482;
t384 = t402 - t451;
t383 = t402 + t451;
t381 = t413 * t414 + t440;
t380 = t447 * qJDD(1);
t379 = -0.2e1 * t398 + t439;
t376 = 0.2e1 * t434 + t441;
t370 = t412 * t385;
t369 = t447 * t442;
t362 = 0.2e1 * t444;
t360 = -pkin(4) * t381 + g(3) * t413;
t356 = -t367 + t393;
t355 = t366 - t393;
t354 = t377 * t412 - t404 * t442;
t353 = -t378 * t409 - t405 * t442;
t349 = -t390 * t409 - t454;
t348 = -t389 * t409 + t370;
t347 = t392 * t412 - t456;
t346 = t391 * t412 - t455;
t345 = t390 * t412 - t455;
t344 = t392 * t409 + t370;
t340 = -t367 + t366;
t339 = t380 * t413 - t383 * t410;
t338 = t380 * t410 + t383 * t413;
t336 = -t376 * t409 + t379 * t412;
t325 = t349 * t413 + t376 * t410;
t324 = t347 * t413 - t379 * t410;
t323 = t349 * t410 - t376 * t413;
t322 = t347 * t410 + t379 * t413;
t321 = -pkin(5) * t345 - t458;
t320 = -pkin(5) * t344 - t459;
t317 = (t371 * t411 - t373 * t408) * t395;
t316 = (-t371 * t408 - t373 * t411) * t395;
t313 = -pkin(1) * t345 + t352;
t312 = -pkin(1) * t344 + t351;
t299 = t373 * t453 - t411 * t416;
t298 = t373 * t452 + t408 * t416;
t297 = -t333 * t408 - t371 * t452;
t296 = -t333 * t411 + t371 * t453;
t295 = t317 * t412 - t368 * t409;
t294 = t317 * t409 + t368 * t412;
t293 = t355 * t411 + t463;
t292 = -t356 * t408 - t460;
t291 = -t355 * t408 + t462;
t290 = -t356 * t411 + t461;
t287 = t311 * t413 - t363 * t410;
t286 = t311 * t410 + t363 * t413;
t281 = t299 * t412 + t437;
t280 = t297 * t412 - t437;
t279 = t299 * t409 - t436;
t278 = t297 * t409 + t436;
t277 = -pkin(3) * t483 + qJ(4) * t327;
t273 = -t303 * t411 - t408 * t483;
t271 = t303 * t408 - t411 * t483;
t269 = t295 * t413 - t316 * t410;
t268 = t295 * t410 + t316 * t413;
t267 = t464 - t469;
t266 = t293 * t412 - t304 * t409;
t265 = t292 * t412 - t308 * t409;
t264 = t293 * t409 + t304 * t412;
t263 = t292 * t409 + t308 * t412;
t262 = t465 - t470;
t254 = t273 * t412 - t340 * t409;
t253 = t273 * t409 + t340 * t412;
t252 = t281 * t413 - t298 * t410;
t251 = t280 * t413 - t296 * t410;
t250 = t281 * t410 + t298 * t413;
t249 = t280 * t410 + t296 * t413;
t246 = -qJ(4) * t337 + t261;
t244 = -pkin(2) * t272 - pkin(3) * t308;
t243 = t276 - t480;
t242 = t275 - t481;
t241 = -pkin(3) * t366 + t421;
t240 = t266 * t413 - t291 * t410;
t239 = t265 * t413 - t290 * t410;
t238 = t266 * t410 + t291 * t413;
t237 = t265 * t410 + t290 * t413;
t236 = t260 * t413 + t288 * t410;
t234 = pkin(4) * t236;
t230 = -pkin(3) * t303 + qJ(4) * t335 - t261;
t229 = t257 * t413 + t283 * t410;
t227 = pkin(4) * t229;
t226 = t232 * t412 + t318 * t409;
t225 = t232 * t409 - t318 * t412;
t224 = t254 * t413 - t271 * t410;
t223 = t254 * t410 + t271 * t413;
t222 = t362 + (-t308 - t358) * qJ(4) + t484 + t435;
t221 = t248 * t413 + t272 * t410;
t219 = pkin(4) * t221;
t218 = -qJ(4) * t304 + (-t326 - t366) * pkin(3) + t421;
t217 = -t480 + (-t337 - t366) * pkin(3) + t421;
t216 = t418 - t465;
t215 = t419 + t464;
t214 = -t231 - t471;
t213 = -pkin(3) * t261 + qJ(4) * t241;
t212 = t362 - t422 - t481 + 0.2e1 * t484;
t211 = t246 * t411 - t277 * t408 - t469;
t210 = qJ(4) * t460 - t230 * t408 - t470;
t209 = t241 * t411 - t467;
t208 = t241 * t408 + t466;
t207 = -t243 * t409 + t267 * t412 - t472;
t206 = -t242 * t409 + t262 * t412 - t473;
t205 = t226 * t413 + t231 * t410;
t204 = t226 * t410 - t231 * t413;
t203 = -t246 * t408 - t277 * t411 + t418;
t202 = -t232 + t420;
t201 = -pkin(1) * t225 + pkin(2) * t318 - pkin(6) * t232;
t200 = -qJ(4) * t461 - t230 * t411 + t419;
t199 = t209 * t412 + t261 * t409;
t198 = t209 * t409 - t261 * t412;
t197 = t214 * t412 + t272 * t479 - t474;
t196 = -pkin(2) * t208 - pkin(3) * t233;
t195 = -t218 * t408 + t222 * t411 - t471;
t194 = -pkin(5) * t225 + (-pkin(6) * t412 + t479) * t231;
t193 = t211 * t412 - t217 * t409 - t472;
t192 = t210 * t412 - t212 * t409 - t473;
t191 = -t218 * t411 - t222 * t408 + t420;
t190 = -pkin(6) * t208 - qJ(4) * t466 - t213 * t408;
t189 = t199 * t413 + t208 * t410;
t188 = t199 * t410 - t208 * t413;
t187 = t195 * t412 - t244 * t409 - t474;
t186 = -pkin(1) * t198 + pkin(2) * t261 - pkin(6) * t209 + qJ(4) * t467 - t213 * t411;
t185 = -pkin(5) * t198 + t190 * t412 - t196 * t409;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t381, -t382, 0, t343, 0, 0, 0, 0, 0, 0, t324, t325, t339, t287, 0, 0, 0, 0, 0, 0, t229, t236, t221, t205, 0, 0, 0, 0, 0, 0, t229, t236, t221, t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t382, -t381, 0, t342, 0, 0, 0, 0, 0, 0, t322, t323, t338, t286, 0, 0, 0, 0, 0, 0, t228, t235, t220, t204, 0, 0, 0, 0, 0, 0, t228, t235, t220, t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t344, t345, 0, -t310, 0, 0, 0, 0, 0, 0, t256, t259, t247, t225, 0, 0, 0, 0, 0, 0, t256, t259, t247, t198; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t382, 0, -t381, 0, t426, -t360, -t342, -pkin(4) * t342, t354 * t413 - t430, t336 * t413 - t384 * t410, t348 * t413 + t409 * t440, t353 * t413 + t430, t346 * t413 + t410 * t439, qJDD(2) * t410 + t369 * t413, -pkin(4) * t322 - t312 * t410 + t320 * t413, -pkin(4) * t323 - t313 * t410 + t321 * t413, -pkin(4) * t338 + t310 * t413, -pkin(4) * t286 - (pkin(1) * t410 - pkin(5) * t413) * t310, t252, t224, t239, t251, t240, t269, t206 * t413 - t215 * t410 - t476, t207 * t413 - t216 * t410 - t475, t197 * t413 - t202 * t410 - t477, -pkin(4) * t204 + t194 * t413 - t201 * t410, t252, t224, t239, t251, t240, t269, t192 * t413 - t200 * t410 - t476, t193 * t413 - t203 * t410 - t475, t187 * t413 - t191 * t410 - t477, -pkin(4) * t188 + t185 * t413 - t186 * t410; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t381, 0, t382, 0, t360, t426, t343, pkin(4) * t343, t354 * t410 + t429, t336 * t410 + t384 * t413, t348 * t410 - t409 * t438, t353 * t410 - t429, t346 * t410 - t412 * t438, -qJDD(2) * t413 + t369 * t410, pkin(4) * t324 + t312 * t413 + t320 * t410, pkin(4) * t325 + t313 * t413 + t321 * t410, pkin(4) * t339 + t310 * t410, pkin(4) * t287 - (-pkin(1) * t413 - pkin(5) * t410) * t310, t250, t223, t237, t249, t238, t268, t206 * t410 + t215 * t413 + t227, t207 * t410 + t216 * t413 + t234, t197 * t410 + t202 * t413 + t219, pkin(4) * t205 + t194 * t410 + t201 * t413, t250, t223, t237, t249, t238, t268, t192 * t410 + t200 * t413 + t227, t193 * t410 + t203 * t413 + t234, t187 * t410 + t191 * t413 + t219, pkin(4) * t189 + t185 * t410 + t186 * t413; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t387, t388, 0, 0, t425 * t409, t376 * t412 + t379 * t409, t389 * t412 + t456, -t424 * t412, t391 * t409 + t454, 0, pkin(1) * t379 + pkin(5) * t347 + t458, -pkin(1) * t376 + pkin(5) * t349 - t459, pkin(1) * t383 + pkin(5) * t380 + t311, pkin(1) * t363 + pkin(5) * t311, t279, t253, t263, t278, t264, t294, t242 * t412 + t262 * t409 + t449, t243 * t412 + t267 * t409 + t448, t214 * t409 - t272 * t478 + t450, pkin(5) * t226 + (-pkin(1) + t427) * t231, t279, t253, t263, t278, t264, t294, t210 * t409 + t212 * t412 + t449, t211 * t409 + t217 * t412 + t448, t195 * t409 + t244 * t412 + t450, -pkin(1) * t208 + pkin(5) * t199 + t190 * t409 + t196 * t412;];
tauB_reg = t1;
